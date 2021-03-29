/************************************************
Title: halos.c
Purpose: Contains functions related to reading in
         information from amiga's halo files
Notes:   * This file SHOULD be all set to go in
           parallel. See notes about some subhalos
           missing from the halos list, though.
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/***********************
      load_halos
***********************/
void load_halos(void)
{
   // Driver function for reading in halos

   // Read in halo properties for multiple AHF file sets
   if(n_halo_tasks > 1)
   {
      read_halo_particles();

      // Get total number of halos across all file sets
      MPI_Reduce(&nhalos_local, &nhalos_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);  

      // Read in substructure info
      read_halo_substruct();

      // Get M_vir for each halo
      read_virial_mass();
   }

   // Read in halo properties for a single AHF file set. It's the same thing as above,
   // but now only root reads the AHF files
   else
   {
      if(thistask == 0)
      {
         read_halo_particles();
         read_halo_substruct();
         read_virial_mass();
      }
   }
}



/***********************
  read_halo_particles
***********************/
void read_halo_particles(void)
{
   // Reads in the particles for each halo
   // Every halo that's listed in the particles file is also listed in the
   // halos file and vice versa, so I'm good getting the hids from the 
   // particles file

   int i;
   int j;
   int type;
   char *fname;
   FILE *fd = NULL;

   // Get the particle file name from the prefix passed
   // in the parameter file
   fname = get_halo_fname("p");

   // Open particles file
   if(!(fd = fopen(fname, "r")))
   {
      printf("Error, could not open file for reading number of halos!\n");
      exit(EXIT_FAILURE);
   }

   // Get number of halos
   fscanf(fd, " %d\n", &nhalos_local);

   // Get the max number of halos held by a single processor. This is to make communicating
   // between processors easier because every processor will always have something to send.
   // Only do this if we're working with multiple AHF file sets
   if(n_halo_tasks != 1)
   {
      MPI_Allreduce(&nhalos_local, &nhalos_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   }

   // When working with only one set, the max will be the local value
   else
   {
      nhalos_max = nhalos_local;
   }

   // Allocate memory for halos struct
   if(!(H = calloc(nhalos_max, sizeof(HALO_DATA))))
   {
      printf("Error, could not allocate memory for halos!\n");
      exit(EXIT_FAILURE);
   }

   // Initialize halos struct
   for(i = 0; i < nhalos_max; i++)
   {
      H[i].npart = 0;
      H[i].hid = 0;
      H[i].plist = NULL;
      H[i].nsub = 0;
      H[i].sublist = NULL;
      H[i].m_vir = 0.0;
      H[i].host_id = 0;
      H[i].new_id = 0;
   }

   // Read in halo info. We loop over nhalos_max instead of nhalos_local b/c
   // it is much easier to do the communication to root if every task has
   // the same number of halos. See flag.c. If we're outside the range of
   // nhalos_local, I make a 'ghost halo' (ghostlo). This ghostlo has 1
   // particle with id -1 and m_vir = -1 to distinuish it as as ghostlo
   // intead of as a halo. These ghostlos are also padded at the end of H,
   // which means they should not interfere with reading the substruct and
   // viral mass from the other halo files.
   for(i = 0; i < nhalos_max; i++)
   {
      if(i < nhalos_local)
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

      else
      {
         if(!(H[i].plist = calloc(1, sizeof(int))))
         {
            printf("Error, could not allocate memory for ghostlo!\n");
            exit(EXIT_FAILURE);
         }

         H[i].plist[0] = -1;
         H[i].npart = 1;
         H[i].m_vir = -1;
         H[i].hid = -1;
         H[i].nsub = 0;
      }
   }

   // Close particles file
   fclose(fd);

   // Clean
   free(fname);
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
   char *fname;

   // Get file name
   fname = get_halo_fname("s");

   // Open file for reading
   if(!(fd = fopen(fname, "r")))
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
   
   // When mpi is turned on, ahf produces one set of files (halos, subs,
   // parts, profs) for each processor it was run with. From my testing,
   // it seems that every halo with subs that's listed in the halos file
   // for that processor is also in that processor's subs file. Additionally,
   // each host halo that's listed in the subs file is also in that 
   // processor's halos file (that is, there appear to be none that are left out
   // or in another processor's file). However, it now appears that some subhalos
   // that are listed as subs in the subs file are not listed in the halos file for
   // that processor.  

   // Because all this does is read the sublist from the subhalos file, this function
   // will work in parallel
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

         if(j >= nhalos_local)
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

   // Clean
   free(fname);
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
   char *fname;

   // Get filename
   fname = get_halo_fname("h");

   // Open the file for reading
   if(!(fd = fopen(fname, "r")))
   {
      printf("Error, could not open halo file for reading!\n");
      exit(EXIT_FAILURE);
   }

   // Skip the header (only applicable on the root set of files)
   if(thistask == 0)
   {
       while((i = getc(fd)) != '\n')
       {
          continue;
       }
   }

   // Loop over every halo
   for(i = 0; i < nhalos_local; i++)
   {
      // Read M_vir. It's the 4th column
      fscanf(fd, "%ld %ld %d %f", &hid, &H[i].host_id, &nsub, &H[i].m_vir);

      // Error check. I might need to do this after each halo has been read in, as I'm
      // not sure about the ordering when ahf is run in parallel. That is, I might need
      // to read each halo and then sort them by hid in both H and whatever new array I use
      // here, if necessary
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

   // Clean
   free(fname);
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



/***********************
    get_halo_fname
***********************/
char *get_halo_fname(char *ftype)
{
   // This function gets the name of the file specified by ftype
   // based on the prefix passed in the parameter file. This is
   // needed because AHF appends strange numbers to each file that
   // I don't know how to anticipate and can't find in the docs. 
   // In order to have serial compatibility, we need to build the file
   // names based on whether or not we're running with multiple processors.
   // This always assumes that the number of ahf file sets is equal to the
   // number of processors tspec is being run with, as I don't care to sort
   // out the case of different numbers of filesets and processors.

   char fname[100];
   char command[200];
   char f_to_open[200];
   char *buf;
   char *split_str;
   FILE *fp = NULL;

   if(!(buf = calloc(150, sizeof(char))))
   {
      printf("Error, could not allocate memory for buf!\n");
      exit(EXIT_FAILURE);   
   }

   // Build the fname with the wildcard character
   if(strcmp(ftype, "p") == 0)
   {
      // Build parallel name. Make sure there's more than one AHF file set
      if((ntasks > 1) && (n_halo_tasks != 1))
      {
         sprintf(fname, "%s.%04d*.AHF_particles", part_file, thistask);
      }

      // Serial name
      else
      {
         sprintf(fname, "%s.*.AHF_particles", part_file);

      }
   }

   else if(strcmp(ftype, "s") == 0)
   {
      // Parallel name
      if((ntasks > 1) && (n_halo_tasks != 1))
      {
         sprintf(fname, "%s.%04d*.AHF_substructure", subfile, thistask);
      }

      // Serial name
      else
      {
         sprintf(fname, "%s.*.AHF_substructure", subfile);
      
      }
   }

   else if(strcmp(ftype, "h") == 0)
   {
      // Parallel name
      if((ntasks > 1) && (n_halo_tasks != 1))
      {
         sprintf(fname, "%s.%04d*.AHF_halos", halo_file, thistask);
      }

      // Serial name
      else
      {
         sprintf(fname, "%s.*.AHF_halos", halo_file);

      }
   }

   sprintf(command, "/bin/ls %s", fname);

   // Run the ls command. The only matching option should be the file I need.
   // This is run to expand the wildcard character
   if((fp = popen(command, "r")) == NULL)
   {
      printf("Error, did not run command!\n");
      exit(EXIT_FAILURE);
   }

   if(fgets(f_to_open, sizeof(f_to_open) - 1, fp))
   {
      // ls appends a newline, so we have to get rid of it
      split_str = strtok(f_to_open, "\n");
      strcpy(buf, split_str);
   }

   pclose(fp);

   return buf;
}
