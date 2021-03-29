/************************************************
Title: flag.c
Purpose: Contains functions related to flagging 
         halo particles.
Notes:   * I misunderstood what the substruct option
            in amiga does. It only controls if the substruct
            file is written. AHF still finds subhalos. If I 
            want to ignore those, I have to explicitly skip
            them, which is why I commented out the else block
            below. Also, I don't want to skip the halos with
            subhalos, I want to skip the subhalos themselves,
            so I changed the if statement to look at host_id.
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "allvars.h"
#include "proto.h"



/***********************
    flag_halo_parts
***********************/
void flag_halo_parts(PARTICLE_DATA *P, int ngas)
{
   // Flags which particles belong to halos and assigns
   // the appropriate m_vir

   int counter = 0;
   int i;
   int index;

   while(counter < nhalos)
   {
      // Progress bar
      //printf("\r%f", ((float)counter / (float)nhalos) * 100.0);
      //fflush(stdout);

      // Halo has no substructures
      if(H[counter].nsub == 0)
      //if(H[counter].host_id == 0) // See NOTES in header.
      {
         // Only keep those halos > 10^9 M_sun
         //if(H[counter].m_vir >= 1.0e9)
         //{
            // Loop over plist
            for(i = 0; i < H[counter].npart; i++)
            {
               index = H[counter].plist[i] - 1;

               // Flag the particle
               P[index].in_halo = 1;
               P[index].m_vir = H[counter].m_vir;
            }
         //}

         // Update counter
         counter++;
      }

      // Halo has substructures
      else
      {
         // I want to test with no sub struct, so if we somehow get here,
         // that's bad. Fail.
         //printf("Error, there was a halo with substructure! %d\n", H[counter]hid);
         //exit(EXIT_FAILURE);

         // Remove duplicate particles from current halo
         remove_duplicates(counter);
      
         // Loop over plist
         for(i = 0; i < H[counter].npart; i++)
         {
            // Only do non-duplicate particles
            if(H[counter].plist[i] != -1)
            {
               index = H[counter].plist[i] - 1;

               // Flag the particle
               P[index].in_halo = 1;
               P[index].m_vir = H[counter].m_vir;
            }
         }

         // Update the counter
         counter++;
      }
   
      //counter++;
   }
}



/***********************
   remove_duplicates
***********************/
void remove_duplicates(int current)
{
   // Removes those particles that are in the first subhalo
   // level down from the current halo's plist

   int i;
   int j;
   int k;
   int sub_ind;

   // Loop over the subhalos
   for(i = 0; i < H[current].nsub; i++)
   {
      sub_ind = 0;

      // Find subhalo
      while(H[sub_ind].hid != H[current].sublist[i])
      {
         sub_ind++;
      }

      // Loop over every particle in current and compare to subhalo
      for(j = 0; j < H[current].npart; j++)
      {
         for(k = 0; k < H[sub_ind].npart; k++)
         {
            if(H[current].plist[j] == H[sub_ind].plist[k])
            {
               // Flag as a duplicate by changing its id to -1
               H[current].plist[j] = -1;
               continue;
            }
         }
      }
   }
}
