/************************************************
Title: allvars.c
Purpose: Contains all instantizations of globals
Notes:
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "allvars.h"



// MPI
int thistask;
int ntasks;
MPI_Datatype mpi_header_type;
MPI_Datatype mpi_particle_type;

// Files
char snapfile[256];
char part_file[100];
char halo_file[100];
char subfile[100];

// Cosmology
float a_dot;
int nhalos_local;
int nhalos_tot = 0;
int nhalos_max;
int max_subs_global;
float LITTLE_H;
float OMEGA_R0;
float OMEGA_K0;
float OMEGA_M0;
float OMEGA_DE0;
float BARYON_FRAC = 0.155;
float PROTON_MASS = 1.6726e-24; // in g
float MOL_WEIGHT = 0.588;      // Primordial H and He
float BOLTZMANN = 1.3806e-16;  // erg/K

// Units
double GUL_IN_CM;
double GUV_IN_CM_PER_S;
double GUM_IN_G;

int de_flag;

float DE_W0;
float DE_WA;

int n_halo_tasks; 

// Particle Data
IO_HEADER header;

// Halo Data
HALO_DATA *H;
