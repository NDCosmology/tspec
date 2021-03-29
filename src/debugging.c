/************************************************
Title: debugging.c
Purpose: Contains functions that are used only
         when the debugging option is turned on.
Notes:
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



void write_flagged_particles(PARTICLE_DATA *P, int ngas)
{
    // Write all of the particles flagged as being in halos to a file, since gdb
    // is very slow at this. Also write all of the halo particles to a file (sans-duplicates)
    // because I believe remove_duplicates works
    FILE *fd;
    int i;

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
}
