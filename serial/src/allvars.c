/************************************************
Title: allvars.c
Purpose: Contains all instantizations of globals
Notes:
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "allvars.h"



// Files
char snapfile[256];
char part_file[100];
char halo_file[100];
char subfile[100];

// Cosmology
float a_dot;
int nhalos;
float LITTLE_H;
float OMEGA_R0;
float OMEGA_K0;
float OMEGA_M0;
float OMEGA_DE0;

// Units
double GUL_IN_CM;
double GUV_IN_CM_PER_S;
double GUM_IN_G;

int de_flag;

float DE_W0;
float DE_WA;

// Particle Data
IO_HEADER header;

// Halo Data
HALO_DATA *H;
