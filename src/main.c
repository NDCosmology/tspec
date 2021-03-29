/************************************************
Title: tspec
Author: Jared Coughlin
Date: 2/25/16
Purpose: Calculates particle temperatures
Notes:   * Uses amiga halo files
         * 6/16/16: Added general T0(z) instead
           of having to just use Bertone's fixed
           value. Did this using Bolton's table
           that Bertone sent me
         * 5/26/17: Started adding mpi
         * 6/19/17: Finished parallelizing the 
           read-in and changing the timing.
           Starting to parallelize the temp calc.
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



int main(int argc, char **argv)
{
   int i;
   int ngas;
   int n_this_task;
   double start;
   double end;
   double tot_time_local;
   double tot_time_global;
   PARTICLE_DATA *P;
   PARTICLE_DATA *All_P;

   // Set up MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &thistask);

   // Get start time
   start = MPI_Wtime();

   // Initialize
   if(thistask == 0)
   {
      printf("Initializing...\n");
      fflush(stdout);
   }
   init(argc, argv[1]);

   // Load the halo information
   if(thistask == 0)
   {
      printf("Loading halos...\n");
      fflush(stdout);
   }
   load_halos();

   // Read the snapshot file
   if(thistask == 0)
   {
      printf("Loading snapshot...\n");
      fflush(stdout);
      header = load_header();
      All_P = load_snapshot(&ngas);
      
      // Sort by id in ascending order. I'm not sure how to do this
      // in parallel, which is why it's in serial
      qsort(All_P, ngas, sizeof(PARTICLE_DATA), pid_cmp);
   }

   // Bcast the header
   MPI_Bcast(&header, 1, mpi_header_type, 0, MPI_COMM_WORLD);  
   
   // Set cosmological parameters from header
   LITTLE_H = header.HubbleParam;
   OMEGA_R0 = 0.0;
   OMEGA_M0 = header.Omega0;
   OMEGA_DE0 = header.OmegaLambda;
   OMEGA_K0 = 1.0 - header.Omega0 - header.OmegaLambda;

   // Flag halo particles
   if(thistask == 0)
   {
      printf("Flagging halo particles...\n");
      fflush(stdout);
   }
   flag_halo_parts(All_P);

   // Free halo resources
   for(i = 0; i < nhalos_max; i++)
   {   
      // If there's only one AHF file set then H is only allocated on root     
      if(n_halo_tasks == 1)
      {
         if(thistask == 0)
         {
            free(H[i].plist);
            free(H[i].sublist);
         }
      }

      else
      {
         free(H[i].plist);
         free(H[i].sublist);
      }
   }

   // If there's only one AHF file set then H is only allocated on root
   if(n_halo_tasks == 1)
   {
      if(thistask == 0)
      {
         free(H);
      }
   }

   else
   {
      free(H);
   }

   #ifdef DEBUGGING
     if(thistask == 0)
     {
        write_flagged_particles(All_P, ngas);
     }
   #endif

   // Calculate temperatures
   if(thistask == 0)
   {
      printf("Calculating temperatures...\n");
      fflush(stdout);
   }
   P = get_temperatures(All_P, ngas, &n_this_task);

   // Write data
   if(thistask == 0)
   {
      printf("Writing data...\n");
      fflush(stdout);
   }
   write_particle_data(P, n_this_task);

   // Get end time
   end = MPI_Wtime();

   // Get max time across all processors (probably root)
   tot_time_local = end - start;
   MPI_Allreduce(&tot_time_local, &tot_time_global, 1, MPI_DOUBLE, MPI_MAX, 
                 MPI_COMM_WORLD); 

   if(thistask == 0)
   {
      printf("Total number of halos: %d\n", nhalos_max);
      printf("Total time to run: %lf\n", tot_time_global);
      printf("Done.\n");
   }

   // Clean up mpi
   MPI_Finalize();

   return 0;
}
