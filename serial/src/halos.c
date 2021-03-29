/************************************************
Title: halos.c
Purpose: Contains functions related to reading in
         information from amiga's halo files
Notes:
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "allvars.h"
#include "proto.h"



/***********************
      load_halos
***********************/
void load_halos(void)
{
   // Driver function for reading in halos

   // Read in halo particles
   read_halo_particles();

   // Read in substructure info
   read_halo_substruct();

   // Get M_vir for each halo
   read_virial_mass();
}



/***********************
  read_halo_particles
***********************/
void read_halo_particles(void)
{
   // Reads in the particles for each halo

   int i;
   int j;
   int type;
   FILE *fd;

   // Open particles file
   if(!(fd = fopen(part_file, "r")))
   {
      printf("Error, could not open file for reading number of halos!\n");
      exit(EXIT_FAILURE);
   }

   // Get number of halos
   fscanf(fd, " %d\n", &nhalos);

   // Allocate memory for halos struct
   if(!(H = calloc(nhalos, sizeof(HALO_DATA))))
   {
      printf("Error, could not allocate memory for halos!\n");
      exit(EXIT_FAILURE);
   }

   // Initialize halos struct
   for(i = 0; i < nhalos; i++)
   {
      H[i].npart = 0;
      H[i].hid = 0;
      H[i].plist = NULL;
      H[i].nsub = 0;
      H[i].sublist = NULL;
      H[i].m_vir = 0.0;
      H[i].host_id = 0;
   }

   // Read in halo info
   for(i = 0; i < nhalos; i++)
   {
      // Get number of particles and halo id
      fscanf(fd, "%d %ld\n", &H[i].npart, &H[i].hid);

      // Allocate memory for particle list
      if(!(H[i].plist = calloc(H[i].npart, sizeof(int))))
      {
         printf("Error, could not allocate memory for plist!\n");
         exit(EXIT_FAILURE);
      }

      // Read in halo particles
      for(j = 0; j < H[i].npart; j++)
      {
         fscanf(fd, "%d %d\n", &H[i].plist[j], &type);
      }

      // Sort the particle list
      qsort(H[i].plist, H[i].npart, sizeof(int), cmpfunc);
   }

   // Close particles file
   fclose(fd);
}



/***********************
  read_halo_substruct
***********************/
void read_halo_substruct(void)
{
   // Gets list of subhalos and nsub_halos

   FILE *fd;
   int i;
   int j;
   int k;
   int nsub;
   long int hid;
   int nlines = 0;

   // Open file for reading
   if(!(fd = fopen(subfile, "r")))
   {
      printf("Error, could not open file for reading substructure!\n");
      exit(EXIT_FAILURE);
   } 

   // Get the number of lines in the file
   while((i = getc(fd)) != EOF)
   {
      if(i == '\n')
      {
         nlines++;
      }
   }

   // There has to be an even number of lines in the file b/c it writes
   // hid nsub_halos\n sub1 sub2 ...subN\n (that is, there are two lines
   // for every halo: the halo id with nsub_halos, and then the list of
   // subhalos. Only those halos that have sub structure are written)
   nlines = nlines / 2;

   // Go back to beginning
   rewind(fd);

   // Loop over sub structure file
   for(i = 0; i < nlines; i++)
   {
      j = 0;

      // Line up current halo with entry in H
      fscanf(fd, "%ld %d\n", &hid, &nsub);

      while(H[j].hid != hid)
      {
         j++;

         if(j >= nhalos)
         {
            printf("Error, could not find current halo in substruct in H!\n");
            exit(EXIT_FAILURE);
         }
      }

      H[j].nsub = nsub;

      // Allocate memory for sublist
      if(!(H[j].sublist = calloc(H[j].nsub, sizeof(long int))))
      {
         printf("Error, could not allocate memory for sublist!\n");
         exit(EXIT_FAILURE);
      }

      // Read in the sublist
      for(k = 0; k < H[j].nsub; k++)
      {
         if(k != (H[j].nsub - 1))
         {
            fscanf(fd, "%ld ", &H[j].sublist[k]);
         }

         else
         {
            fscanf(fd, "%ld\n", &H[j].sublist[k]);
         }
      }
   }

   // Close file
   fclose(fd);
}



/***********************
   read_virial_mass
***********************/
void read_virial_mass(void)
{
   // Reads the virial mass for each halo from the
   // halos file

   FILE *fd;
   int i;
   int j;
   int nsub;
   long int hid;

   // Open the file for reading
   if(!(fd = fopen(halo_file, "r")))
   {
      printf("Error, could not open halo file for reading!\n");
      exit(EXIT_FAILURE);
   }

   // Skip the header
   while((i = getc(fd)) != '\n')
   {
      continue;
   }

   // Loop over every halo
   for(i = 0; i < nhalos; i++)
   {
      // Read M_vir. It's the 4th column
      fscanf(fd, "%ld %ld %d %f", &hid, &H[i].host_id, &nsub, &H[i].m_vir);

      // Error check
      if((hid != H[i].hid) || (nsub != H[i].nsub))
      {  
         printf("Error reading nsub or halos have different order in halos and particles files!\n");
         exit(EXIT_FAILURE);
      }

      // Read the rest of the line
      while((j = getc(fd)) != '\n')
      {
         continue;
      }        
   }

   // Close file
   fclose(fd);
}



/***********************
        cmpfunc
***********************/
int cmpfunc(const void *p1, const void *p2)
{
   // Used by qsort to sort pids in ascending order
   
   const int *elem1 = p1;
   const int *elem2 = p2;

   return (int)(*elem1 - *elem2);
}
