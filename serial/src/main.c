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
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"



int main(int argc, char **argv)
{
   int i;
   int ngas;
   char buffer[256];
   clock_t start;
   clock_t end;
   FILE *fb;
   PARTICLE_DATA *P;

   char buffer2[256];

   // Check args
   printf("Initializing...\n");
   fflush(stdout);

   if(argc != 2)
   {
      printf("Usage: ./tspec tspec_param.param\n");
      exit(EXIT_FAILURE);
   }

   // Get start time
   start = clock();

   // Set args by reading parameter file
   if(!(fb = fopen(argv[1], "r")))
   {
      printf("Error, could not open parameter file for reading!\n");
      exit(EXIT_FAILURE);
   }

   fgets(buffer, sizeof(buffer), fb);
   sscanf(buffer, "Snapshot%s", snapfile);

   fgets(buffer, sizeof(buffer), fb);
   sscanf(buffer, "HaloPartsFile%s", part_file);

   fgets(buffer, sizeof(buffer), fb);
   sscanf(buffer, "HaloFile%s", halo_file);

   fgets(buffer, sizeof(buffer), fb);
   sscanf(buffer, "HaloSubFile%s", subfile);

   fgets(buffer,sizeof(buffer),fb);
   sscanf(buffer, "GUL_IN_CM%s", buffer2);
   GUL_IN_CM = atof(buffer2);
   
   fgets(buffer,sizeof(buffer),fb);
   sscanf(buffer, "GUV_IN_CM_PER_S%s", buffer2);
   GUV_IN_CM_PER_S = atof(buffer2);

   fgets(buffer,sizeof(buffer),fb);
   sscanf(buffer, "GUM_IN_G%s", buffer2);
   GUM_IN_G = atof(buffer2);

   fgets(buffer,sizeof(buffer),fb);
   sscanf(buffer, "DE%s", buffer2);
   de_flag = atoi(buffer2);

   fgets(buffer, sizeof(buffer), fb);
   sscanf(buffer, "DE_W0%s", buffer2);
   DE_W0 = atof(buffer2);
   
   fgets(buffer, sizeof(buffer), fb);
   sscanf(buffer, "DE_WA%s", buffer2);
   DE_WA = atof(buffer2);

   // Close file
   fclose(fb);

   // Do a DE error check
   de_error_check();

   // Load the snapshot header
   header = load_header();

   // Load the halo information
   printf("Loading halos...\n");
   fflush(stdout);
   load_halos();

   // Read the snapshot file
   printf("Loading snapshot...\n");
   fflush(stdout);
   P = load_snapshot(&ngas);

   LITTLE_H = header.HubbleParam;
   OMEGA_R0 = 0.0;
   OMEGA_M0 = header.Omega0;
   OMEGA_DE0 = header.OmegaLambda;
   OMEGA_K0 = 1.0 - header.Omega0 - header.OmegaLambda;

   // Sort by id in ascending order
   qsort(P, ngas, sizeof(PARTICLE_DATA), pid_cmp);

   // Flag halo particles
   printf("Flagging halo particles...\n");
   fflush(stdout);
   flag_halo_parts(P, ngas);

   /**********************
          debugging
   **********************/
   #ifdef DEBUGGING
      // Write all of the particles flagged as being in halos to a file, since gdb
      // is very slow at this. Also write all of the halo particles to a file (sans-duplicates)
      // because I believe remove_duplicates works
      FILE *fd;
      int j;
      if(!(fd = fopen("flaggedparts_code.txt", "w")))
      {
         printf("Error, cold not open file for writing flagged particles!\n");
         exit(EXIT_FAILURE);
      }
      for(i = 0; i < ngas; i++)
      {
         if(P[i].in_halo == 1)
         {
            fprintf(fd, "%d\n", P[i].id);
         }
      }
      fclose(fd);

      if(!(fd = fopen("amiga_halo_parts.txt", "w")))
      {
         printf("Error, could not open file for writing amiga halo parts!\n");
         exit(EXIT_FAILURE);
      }
      for(i = 0; i < nhalos; i++)
      {
         for(j = 0; j < H[i].npart; j++)
         {
            if(H[i].plist[j] != -1)
            {
               fprintf(fd, "%d\n", H[i].plist[j]);
            }
         }
      }
      fclose(fd);
   #endif
   /**88********************
        end_debugging
   ************************/

   // Calculate temperatures
   printf("Calculating temperatures...\n");
   fflush(stdout);
   get_temperatures(P, ngas);

   // Write data
   printf("Writing data...\n");
   fflush(stdout);
   write(P, ngas);

   // Free particle data
   free(P);

   // Free halo resources
   for(i = 0; i < nhalos; i++)
   {
      free(H[i].plist);
      free(H[i].sublist);
   }

   free(H);
   
   // Get end time
   end = clock();

   printf("Total time to run: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
   printf("Done.\n");

   return 0;
}
