/************************************************
Title: init.c
Purpose: Read and distribute parameter file
Notes:
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/***********************
         init
***********************/
void init(int nargs, char *paramfile)
{
   char buffer[256];
   char buffer2[256];
   FILE *fb;
 
   // Check args
   if(thistask == 0)
   {
       printf("Initializing...\n");
       fflush(stdout);

       if(nargs != 2)
       {
          printf("Usage: ./tspec tspec_param.param\n");
          exit(EXIT_FAILURE);
       }

       // Set args by reading parameter file
       if(!(fb = fopen(paramfile, "r")))
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

       fgets(buffer, sizeof(buffer), fb);
       sscanf(buffer, "N_Halo_Files%s", buffer2);
       n_halo_tasks = atoi(buffer2);

      // Close file
      fclose(fb);
   }

   // Bcast the parameter file data. This is definitely not the best way to do this
   MPI_Bcast(&snapfile, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
   MPI_Bcast(&part_file, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
   MPI_Bcast(&halo_file, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
   MPI_Bcast(&subfile, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
   MPI_Bcast(&GUL_IN_CM, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&GUV_IN_CM_PER_S, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&GUM_IN_G, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&de_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&DE_W0, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&DE_WA, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n_halo_tasks, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
   // Do a DE error check
   de_error_check();

   // Create custom MPI data types
   make_custom_mpi_type();
}



/***********************
  make_custom_mpi_type
***********************/
void make_custom_mpi_type(void)
{
   // See https://www.rc.colorado.edu/sites/default/files/Datatypes.pdf 

   // Header type variables
   int h_blocks[14] = {6,6,1,1,1,1,6,1,1,1,1,1,1,96};
   MPI_Datatype h_types[14] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT,\
                           MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE,\
                           MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR};
   MPI_Aint h_disp[14];
   MPI_Aint lb;
   MPI_Aint charex, intex, floatex, doublex;

   // Particle type variables
   int p_blocks[10] = {3,3,1,1,1,1,1,1,1,1};
   MPI_Datatype p_types[10] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,\
                            MPI_FLOAT, MPI_FLOAT, MPI_INT, MPI_INT, MPI_INT};
   MPI_Aint p_disp[10];

   // Get the extents
   MPI_Type_get_extent(MPI_CHAR, &lb, &charex);
   MPI_Type_get_extent(MPI_INT, &lb, &intex);
   MPI_Type_get_extent(MPI_FLOAT, &lb, &floatex);
   MPI_Type_get_extent(MPI_DOUBLE, &lb, &doublex);

   // Now fill the header displacements
   h_disp[0] = (MPI_Aint)0;
   h_disp[1] = 6 * intex;
   h_disp[2] = h_disp[1] + (6 * doublex);
   h_disp[3] = h_disp[2] + doublex;
   h_disp[4] = h_disp[3] + doublex;
   h_disp[5] = h_disp[4] + intex;
   h_disp[6] = h_disp[5] + intex;
   h_disp[7] = h_disp[6] + (6 * intex);
   h_disp[8] = h_disp[7] + intex;
   h_disp[9] = h_disp[8] + intex;
   h_disp[10] = h_disp[9] + doublex;
   h_disp[11] = h_disp[10] + doublex;
   h_disp[12] = h_disp[11] + doublex;
   h_disp[13] = h_disp[12] + doublex;

   // Make MPI header structure
   MPI_Type_create_struct(14, h_blocks, h_disp, h_types, &mpi_header_type);
   MPI_Type_commit(&mpi_header_type);

   // Fill the particle displacements
   p_disp[0] = (MPI_Aint)0;
   p_disp[1] = 3 * floatex;
   p_disp[2] = p_disp[1] + (3 * floatex);
   p_disp[3] = p_disp[2] + floatex;
   p_disp[4] = p_disp[3] + floatex;
   p_disp[5] = p_disp[4] + floatex;
   p_disp[6] = p_disp[5] + floatex;
   p_disp[7] = p_disp[6] + floatex;
   p_disp[8] = p_disp[7] + intex;
   p_disp[9] = p_disp[8] + intex;

   // Make the MPI particle type
   MPI_Type_create_struct(10, p_blocks, p_disp, p_types, &mpi_particle_type);
   MPI_Type_commit(&mpi_particle_type);
}
