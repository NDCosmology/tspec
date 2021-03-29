/************************************************
Title: write.c
Purpose: Contains functions used for, well, writing
         the data, including the temp
Notes:   This function was ripped from the original dspec
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"



/***********************
         write
***********************/
void write(PARTICLE_DATA *SP, int ngas)
{
   // Does as the name says, really. With the data.

   FILE *fd;
   char tspec_file[256];
   int blksize;
   int k;
   int n_with_mass = 0; // This is very bad. I should recacl this, but
                        // I know it's zero cos they're dm particles.

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
   blksize = ngas * 3 * sizeof(float);
   my_fwrite(&blksize, sizeof(int), 1, fd);
   for(k = 0; k < ngas; k++)
   {
      my_fwrite(&SP[k].pos[0], sizeof(float), 3, fd);
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);
   
   // Velocities
   blksize = ngas * 3 * sizeof(float);
   my_fwrite(&blksize, sizeof(int), 1, fd);
   for(k = 0; k < ngas; k++)
   {
      my_fwrite(&SP[k].vel[0], sizeof(float), 3, fd);
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);

   // IDs
   blksize = ngas * sizeof(int);
   my_fwrite(&blksize, sizeof(int), 1, fd);
   for(k = 0; k < ngas; k++)
   {
      my_fwrite(&SP[k].id, sizeof(int), 1, fd);
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);
   
   // Masses
   blksize = n_with_mass * sizeof(float);

   if(n_with_mass > 0)
   {
      my_fwrite(&blksize, sizeof(int), 1, fd);

      for(k = 0; k < ngas; k++)
      {
         my_fwrite(&SP[k].mass, sizeof(float), 1, fd);
      }

      my_fwrite(&blksize, sizeof(int), 1, fd);
   }

   // Sph properties
   // Temp
   blksize = ngas * sizeof(float);
   my_fwrite(&blksize, sizeof(int), 1, fd);
   for(k = 0; k < ngas; k++)
   {
      my_fwrite(&SP[k].temp, sizeof(float), 1, fd);
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);

   // Density
   blksize = ngas * sizeof(float);
   my_fwrite(&blksize, sizeof(int), 1, fd);
   for(k = 0; k < ngas; k++)
   {
      my_fwrite(&SP[k].density, sizeof(float), 1, fd);
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);

   // Hsml
   blksize = ngas * sizeof(float);
   my_fwrite(&blksize, sizeof(int), 1, fd);
   for(k = 0; k < ngas; k++)
   {
      my_fwrite(&SP[k].hsml, sizeof(float), 1, fd);
   }
   my_fwrite(&blksize, sizeof(int), 1, fd);

   // Close the file
   fclose(fd);

   // If desired, write the hsml and rho for use with splash
   #ifdef GSPLASH_DATA
      write_hsml(SP, ngas);
      write_dens(SP, ngas);
   #endif
}



/***********************
      write_hsml
***********************/
void write_hsml(PARTICLE_DATA *D, int ngas)
{
   // Does as the name says. For use with gsplash

   #ifdef GSPLASH_DATA
      int i;
      char fname[100];
      FILE *fd;

      sprintf(fname, "dm_hsml.txt");

      // Open the file
      if(!(fd = fopen(fname, "w")))
      {
         printf("Error, could not open file for writing hsml!\n");
         exit(EXIT_FAILURE);
      }

      for(i = 0; i < ngas; i++)
      {
         fprintf(fd, "%f\n", D[i].hsml); 
      }

      // Close the file
      fclose(fd);
   #endif
}



/***********************
      write_dens
***********************/
void write_dens(PARTICLE_DATA *D, int ngas)
{
   // Splash takes a list of densities as well for dm particles,
   // so this function writes them out

   #ifdef GSPLASH_DATA
      int i;
      char fname[100];
      FILE *fd;

      sprintf(fname, "dm_rho.txt");

      // Open the file
      if(!(fd = fopen(fname, "w")))
      {
         printf("Error, could not open file for writing dens!\n");
         exit(EXIT_FAILURE);
      }

      for(i = 0; i < ngas; i++)
      {
         fprintf(fd, "%f\n", D[i].density);
      }

      // Close the file
      fclose(fd);
   #endif
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
