/************************************************
Title: load.c
Purpose: Contains functions related to reading the 
         snapshot
Notes: These were mostly developed and tested in
       rs.c
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"



/***********************
      read_header
***********************/
IO_HEADER load_header(void)
{
   // Reads the header from the snapshot. The header is
   // always 256 bytes. It, like every other block in the
   // snapshot, it padded by an integer that holds the
   // size, in bytes, of the block

   FILE *fd;
   int blksize1 = 0;
   int blksize2 = 1;
   enum fields blocknr;
   char read_file[256];
   IO_HEADER h;

   // Open the file
   if(!(fd = fopen(snapfile, "rb")))
   {
      // We could be here either because the file failed to load, or else
      // because there is more than one file per snapshot, so let's try
      // and load snapfile.0 first before quitting
      sprintf(read_file, "%s.%d", snapfile, 0);

      if(!(fd = fopen(read_file, "rb")))
      {
         printf("Error, could not open snapshot for reading header!\n");
         exit(EXIT_FAILURE);
      }
   }

   // Read the header
   my_fread(&blksize1, sizeof(int), 1, fd);
   my_fread(&h, sizeof(IO_HEADER), 1, fd);
   my_fread(&blksize2, sizeof(int), 1, fd);

   // Close the file
   fclose(fd);

   // Check to make sure correct size was read and written
   // Check header size
   if(sizeof(IO_HEADER) != 256)
   {
      printf("Error, incorrect header size! Must be 256 bytes!\n");
      exit(EXIT_FAILURE);
   }

   blocknr = HEADER;
   block_check(blocknr, blksize1, blksize2, h);

   // Do a quick error check on the redshift to make sure it's within
   // the bounds of the HM table
   if(h.redshift > 9.479)
   {
      printf("Error, snapshot's z is outside bounds of HM table!\n");
      exit(EXIT_FAILURE);
   }

   return h;
}



/***********************
     load_snapshot
***********************/
PARTICLE_DATA *load_snapshot(int *ngas)
{
   // Does as the name says. Reads in snapshot segment file_num

   int i;
   int master;
   int k;
   int n;
   int pc = 0;
   int pc_new;
   int pc_sph;
   int blksize1;
   int blksize2;
   int numpart;
   int n_with_masses;
   FILE *fd;
   char filename[256];
   enum fields blocknr;
   IO_HEADER theader;
   PARTICLE_DATA *AP;
   PARTICLE_DATA *D;

   // Open the file
   if(header.num_files == 1)
   {
      if(!(fd = fopen(snapfile, "rb")))
      {
         printf("Error, could not open snapshot for reading!\n");
         exit(EXIT_FAILURE);
      }
   }

   else
   {
      // Allocate memory for D
      if(!(D = calloc(header.npartTotal[1], sizeof(PARTICLE_DATA))))
      {
         printf("Error, could not allocate memory for all particles!\n");
         exit(EXIT_FAILURE);
      }

      // Initialize master index
      master = 0;
   }

   // Loop over every snapshot segment
   for(i = 0; i < header.num_files; i++)
   {
      // Progress bar
      printf("Loading snapshot file %d of %d\n", i + 1, header.num_files);

      // Open file if there's more than one per snapshot
      if(header.num_files > 1)
      {
         sprintf(filename, "%s.%d", snapfile, i);

         if(!(fd = fopen(filename, "rb")))
         {
            printf("Error, could not open snapshot segment for reading!\n");
            exit(EXIT_FAILURE);
         }
      }

      // Read the header
      my_fread(&blksize1, sizeof(int), 1, fd);
      my_fread(&theader, sizeof(IO_HEADER), 1, fd);
      my_fread(&blksize2, sizeof(int), 1, fd);

      blocknr = HEADER;
      block_check(blocknr, blksize1, blksize2, theader);  

      // Get the number of particles in the file
      numpart = theader.npart[0] + theader.npart[1] + theader.npart[2] +
                theader.npart[3] + theader.npart[4] + theader.npart[5];

      // Get the number of particles in the mass block
      for(k = 0, n_with_masses = 0; k < 6; k++)
      {
         if(theader.mass[k] == 0)
         {
            n_with_masses += theader.npart[k];
         }
      }

      // Allocate memory for particles
      if(!(AP = calloc(numpart, sizeof(PARTICLE_DATA))))
      {
         printf("Error, could not allocate memory for particles!\n");
         exit(EXIT_FAILURE);
      }

	   // Read positions and type
	   my_fread(&blksize1, sizeof(int), 1, fd);

	   for(k = 0, pc_new = pc; k < 6; k++)
	   {
		   for(n = 0; n < theader.npart[k]; n++)
		   {
			   my_fread(&AP[pc_new].pos[0], sizeof(float), 3, fd);
			   AP[pc_new].type = k;								
			   pc_new++;
		   }
	   }

	   my_fread(&blksize2, sizeof(int), 1, fd);

      blocknr = POS;
      block_check(blocknr, blksize1, blksize2, theader);

	   // Read velocities
	   my_fread(&blksize1, sizeof(int), 1, fd);

	   for(k = 0, pc_new = pc; k < 6; k++)
	   {
		   for(n = 0; n < theader.npart[k]; n++)
		   {
			   my_fread(&AP[pc_new].vel[0], sizeof(float), 3, fd);
			   pc_new++;
		   }
	   }

	   my_fread(&blksize2, sizeof(int), 1, fd);

      blocknr = VEL;
      block_check(blocknr, blksize1, blksize2, theader);

	   // Read Ids
	   my_fread(&blksize1, sizeof(int), 1, fd);

	   for(k = 0, pc_new = pc; k < 6; k++)
	   {
		   for(n = 0; n < theader.npart[k]; n++)
		   {
			   my_fread(&AP[pc_new].id, sizeof(int), 1, fd);														
			   pc_new++;
		   }
	   }

	   my_fread(&blksize2, sizeof(int), 1, fd);

      blocknr = IDS;
      block_check(blocknr, blksize1, blksize2, theader);

	   // Read masses
	   if(n_with_masses > 0)
	   {
		   my_fread(&blksize1, sizeof(int), 1, fd);
	   }

	   for(k = 0, pc_new = pc; k < 6; k++)
	   {
		   for(n = 0; n < theader.npart[k]; n++)
		   {
			   AP[pc_new].type = k;
			   if(theader.mass[k] == 0)
			   {
				   my_fread(&AP[pc_new].mass,sizeof(float), 1, fd);
			   }
			   else
			   {
				   AP[pc_new].mass = theader.mass[k];
			   }
			   pc_new++;
		   }
	   }

	   if(n_with_masses > 0)
	   {
	      my_fread(&blksize2, sizeof(int), 1, fd);

         blocknr = MASS;
         block_check(blocknr, blksize1, blksize2, theader);
	   }

      // Gas only properties
      if(theader.npart[1] > 0)
      {
         // Skip reading temp since dspec doesn't write it

		   // Read Density
		   my_fread(&blksize1, sizeof(int), 1, fd);

		   for(n = 0, pc_sph = pc; n < theader.npart[1]; n++)
		   {
			   my_fread(&AP[pc_sph].density, sizeof(float), 1, fd);
			   pc_sph++;
		   }

		   my_fread(&blksize2, sizeof(int), 1, fd);

         blocknr = RHO;
         block_check(blocknr, blksize1, blksize2, theader);

		   // Read Hsml
		   my_fread(&blksize1, sizeof(int), 1, fd);

		   for(n = 0, pc_sph = pc; n < theader.npart[1]; n++)
		   {
			   my_fread(&AP[pc_sph].hsml, sizeof(float), 1, fd);

            // Set the m_vir and in_halo flags
            AP[pc_sph].m_vir = 0.0;
            AP[pc_sph].in_halo = 0;
			   pc_sph++;
		   }

		   my_fread(&blksize2, sizeof(int), 1, fd);

         blocknr = HSML;
         block_check(blocknr, blksize1, blksize2, theader);
      }

      // If there's more than one file per snapshot, we need to save the
      // results before moving on
      if(header.num_files > 1)
      {
         for(k = 0; k < theader.npart[1]; k++)
         {
            D[master].pos[0]  = AP[k].pos[0];
            D[master].pos[1]  = AP[k].pos[1];
            D[master].pos[2]  = AP[k].pos[2];
            D[master].vel[0]  = AP[k].vel[0];
            D[master].vel[1]  = AP[k].vel[1];
            D[master].vel[2]  = AP[k].vel[2];
            D[master].mass    = AP[k].mass;
            D[master].density = AP[k].density;
            D[master].hsml    = AP[k].hsml;
            D[master].type    = AP[k].type;
            D[master].id      = AP[k].id;

            // Increment master index
            master++;
         }

         // Free AP
         free(AP);
      }

      // Close the file
      fclose(fd);
   }

   // Update ngas
   *ngas = theader.npartTotal[1];

   if(theader.num_files == 1)
   {
      return AP;
   }

   else
   {
      return D;
   }
}



/***********************
      block_check
***********************/
void block_check(enum fields blocknr, int blksize1, int blksize2, IO_HEADER dummy_header)
{
   // This function gets the size of block and makes sure that
   // the padding values that are written around it (the ones 
   // obtained with SKIP) match that value.

   int size = 0;

   // Make sure the two padding values match
   if(blksize1 != blksize2)
   {
      printf("Paddings don't match! Block: %d\n", blocknr);
      exit(EXIT_FAILURE);
   }

   // Make sure the padding matches the actual size of the block
   size = get_block_size(blocknr, dummy_header);

   if(blksize1 != size)
   {
      printf("Paddings do not match actual block size! Block: %d\n", blocknr);
      exit(EXIT_FAILURE);
   }
}



/***********************
    get_block_size
***********************/
int get_block_size(enum fields blocknr, IO_HEADER h)
{
   // This function gets the actual block size to make sure that the correct
   // values are written to the padding.

   int bsize = 0;
   int i;
   int nmass = 0;

   switch(blocknr)
   {
      case HEADER:
         bsize = sizeof(IO_HEADER);
         break;

      case POS:
      case VEL:
         bsize = sizeof(float) * 3 * (h.npart[0] + h.npart[1] + h.npart[2] + 
                 h.npart[3] + h.npart[4] + h.npart[5]); 
         break;

      case IDS:
         bsize = sizeof(int) * (h.npart[0] + h.npart[1] + h.npart[2] + h.npart[3] +
                 h.npart[4] + h.npart[5]);
         break;

      case MASS:
         for(i = 0; i < 6; i++)
         {
            if(h.mass[i] == 0)
            {
               nmass += h.npart[i];
            }
         }

         bsize = sizeof(float) * nmass;
         break;

      case RHO:
      case HSML:
         bsize = sizeof(float) * h.npart[1];
         break;
   }

   return bsize;
}



/***********************
        my_fread
***********************/
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
   // This function is taken from gadget. It checks to make sure fread has 
   // read the correct number of elements.

   size_t nread;

   if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
   {
      printf("I/O error with fread.\n");
      exit(EXIT_FAILURE);
   }

   return nread;
}



/***********************
        pid_cmp
***********************/
int pid_cmp(const void *p1, const void *p2)
{
   // Used with qsort to sort by id

   const PARTICLE_DATA *elem1 = p1;
   const PARTICLE_DATA *elem2 = p2;

   return (int)(elem1->id - elem2->id);
}
