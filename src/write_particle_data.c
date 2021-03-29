/************************************************
Title: write.c
Purpose: Contains functions used for, well, writing
         the data, including the temp
Notes:   This function was ripped from the original dspec
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/***********************
         write
***********************/
void write_particle_data(PARTICLE_DATA *P, int n_this_task)
{
   // Does as the name says, really. With the data.

   FILE *fd;
   char tspec_file[256];
   int blksize;
   int k;
   int n_with_mass = 0; // This is very bad. I should recalc this, but
                        // I know it's zero cos they're dm particles.
   PARTICLE_DATA *SP;
   int *rcnts;
   int *displs;

   // First regather all the particles onto root, as I just want a single
   // output file, and this is the easiest way to get that. I want a single
   // output file because it means no changes have to be made to my other 
   // codes.

   // Allocate memory for recvcounts, displs, and SP
   if(thistask == 0)
   {
      if(!(rcnts = calloc(ntasks, sizeof(int))))
      {
         printf("Error, count not allocate memory for rcnts when gathering to root!\n");
         exit(EXIT_FAILURE);
      }

      if(!(displs = calloc(ntasks, sizeof(int))))
      {
         printf("Error, count not allocate memory for displs when gathering to root!\n");
         exit(EXIT_FAILURE);
      }

      if(!(SP = calloc(header.npartTotal[1], sizeof(PARTICLE_DATA))))
      {
         printf("Error, could not allocate memory for SP when gathering to root!\n");
         exit(EXIT_FAILURE);
      }
   }

   // Get n_this_task from each processor
   MPI_Gather(&n_this_task, 1, MPI_INT, rcnts, 1, MPI_INT, 0, MPI_COMM_WORLD);

   // Set up displs
   if(thistask == 0)
   {
      displs[0] = 0;
      for(k = 1; k < ntasks; k++)
      {
         displs[k] = displs[k - 1] + rcnts[k - 1];
      }
   }

   MPI_Gatherv(P, n_this_task, mpi_particle_type, SP, rcnts, displs, mpi_particle_type, 0, 
               MPI_COMM_WORLD);

   // Free
   if(thistask == 0)
   {
      free(displs);
      free(rcnts);
   }

   free(P);

   // Now continue on with writing the file on root
   if(thistask == 0)
   {
      // Open file for writing
      sprintf(tspec_file, "%s-tspec", snapfile);

      if(!(fd = fopen(tspec_file, "wb")))
      {
         printf("Error, could not open file for writing dm data!\n");
         exit(EXIT_FAILURE);
      }

      // Write the header
      blksize = 256;

      // Because I changed this code to read all the data at once and am
      // now writing it all to just one file, I need to change header.npart
      // to be the same as header.npartTotal. Only need to do for dm, since 
      // all others should be 0
      header.npart[1] = header.npartTotal[1];

      // Also have to have num_files for the same reason
      header.num_files = 1;

      my_fwrite(&blksize, sizeof(int), 1, fd);
      my_fwrite(&header, sizeof(IO_HEADER), 1, fd);
      my_fwrite(&blksize, sizeof(int), 1, fd);

      // Positions
      blksize = header.npartTotal[1] * 3 * sizeof(float);
      my_fwrite(&blksize, sizeof(int), 1, fd);
      for(k = 0; k < header.npartTotal[1]; k++)
      {
         my_fwrite(&SP[k].pos[0], sizeof(float), 3, fd);
      }
      my_fwrite(&blksize, sizeof(int), 1, fd);
   
      // Velocities
      blksize = header.npartTotal[1] * 3 * sizeof(float);
      my_fwrite(&blksize, sizeof(int), 1, fd);
      for(k = 0; k < header.npartTotal[1]; k++)
      {
         my_fwrite(&SP[k].vel[0], sizeof(float), 3, fd);
      }
      my_fwrite(&blksize, sizeof(int), 1, fd);

      // IDs
      blksize = header.npartTotal[1] * sizeof(int);
      my_fwrite(&blksize, sizeof(int), 1, fd);
      for(k = 0; k < header.npartTotal[1]; k++)
      {
         my_fwrite(&SP[k].id, sizeof(int), 1, fd);
      }
      my_fwrite(&blksize, sizeof(int), 1, fd);
   
      // Masses
      blksize = n_with_mass * sizeof(float);

      if(n_with_mass > 0)
      {
         my_fwrite(&blksize, sizeof(int), 1, fd);

         for(k = 0; k < header.npartTotal[1]; k++)
         {
            my_fwrite(&SP[k].mass, sizeof(float), 1, fd);
         }

         my_fwrite(&blksize, sizeof(int), 1, fd);
      }

      // Sph properties
      // Temp
      blksize = header.npartTotal[1] * sizeof(float);
      my_fwrite(&blksize, sizeof(int), 1, fd);
      for(k = 0; k < header.npartTotal[1]; k++)
      {
         my_fwrite(&SP[k].temp, sizeof(float), 1, fd);
      }
      my_fwrite(&blksize, sizeof(int), 1, fd);

      // Density
      blksize = header.npartTotal[1] * sizeof(float);
      my_fwrite(&blksize, sizeof(int), 1, fd);
      for(k = 0; k < header.npartTotal[1]; k++)
      {
         my_fwrite(&SP[k].density, sizeof(float), 1, fd);
      }
      my_fwrite(&blksize, sizeof(int), 1, fd);

      // Hsml
      blksize = header.npartTotal[1] * sizeof(float);
      my_fwrite(&blksize, sizeof(int), 1, fd);
      for(k = 0; k < header.npartTotal[1]; k++)
      {
         my_fwrite(&SP[k].hsml, sizeof(float), 1, fd);
      }
      my_fwrite(&blksize, sizeof(int), 1, fd);

      // Close the file
      fclose(fd);

      // Free
      free(SP);
   }
}



/***********************
       my_fwrite
***********************/
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
   // Taken from gadget

   size_t nwritten;

   if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
   {
      printf("Error writing to output file!\n");
      exit(EXIT_FAILURE);
   }

   return nwritten;
}
