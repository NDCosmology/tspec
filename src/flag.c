/************************************************
Title: flag.c
Purpose: Contains functions related to flagging 
         halo particles.
Notes:   * Particle ids in gadget start at 1, but my array is
           indexed from 0, hence the index = H[counter].plist[i] - 1
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"



/***********************
    flag_halo_parts
***********************/
void flag_halo_parts(PARTICLE_DATA *P)
{
   // This is the driver function for flagging particles that belong to halos. Normally,
   // AHF is run in parallel and with the same number of processors used here in tspec.
   // This means that each processor can read its own file set and we go from there. But,
   // there is either a memory leak in AHF or (as I suspect) AHF simply gives each
   // processor a full copy of P. This causes the memory usage to balloon out of control,
   // and meant that, for my larger sims, I had to run AHF in serial, which is very slow
   // and produces just one file set. As such, flagging can be done slightly differently
   // in the two cases, and there's more MPI overhead in the former case, which is why
   // I've split them into two different functions. The general idea is the same in both,
   // it's just easier to read and maintain if they're separate.

   // Case of multiple AHF file sets
   if(n_halo_tasks != 1)
   {
      flag_halo_parts_mult_file_sets(P);
   }

   // Case of one AHF file set
   else
   {
      flag_halo_parts_single_file_set(P);
   }
}



/***********************
flag_halo_parts_mult_file_sets
***********************/
void flag_halo_parts_mult_file_sets(PARTICLE_DATA *P)
{
   // Flags which particles belong to halos and assigns
   // the appropriate m_vir. In parallel this involves sending
   // the plists to root and having root do the flagging since
   // there isn't enough memory to give each processor it's own
   // copy of P.

   int i;
   int j;
   int k;
   int l;
   int max_subs_local;
   int totcounts;
   int *npart_per_plist;
   int *recvcnts;
   int *displs;
   int *plist;
   int *plist_displs;
   float *mass_list;

   // Allocate memory for gatherv arrays on root (recvbuf, recvcnts, and displs, etc)
   if(thistask == 0)
   {
      if(!(npart_per_plist = calloc(ntasks, sizeof(int))))
      {
         printf("Error, could not allocate memory for recvbuf!\n");
         exit(EXIT_FAILURE);
      }

      if(!(recvcnts = calloc(ntasks, sizeof(int))))
      {
         printf("Error, could not allocate memory for recvcnts!\n");
         exit(EXIT_FAILURE);
      }

      if(!(displs = calloc(ntasks, sizeof(int))))
      {
         printf("Error, could not allocate memory for displs!\n");
         exit(EXIT_FAILURE);
      }

       if(!(plist_displs = calloc(ntasks, sizeof(int))))
      {
         printf("Error, could not allocate memory for plist_displs!\n");
         exit(EXIT_FAILURE);
      }

      if(!(mass_list = calloc(ntasks, sizeof(float))))
      {
        printf("Error, could not allocate memory for mass list!\n");
        exit(EXIT_FAILURE);
      }

      // Only sending one value at first (npart). Displs refers to the spacing between
      // elements in index space, so here it's just i, since every processor is only
      // sending one value.
      for(i = 0; i < ntasks; i++)
      {
         displs[i] = i;
         recvcnts[i] = 1;
      }
   }

   // Get max number of subhalos
   max_subs_local = H[0].nsub;
   for(i = 1; i < nhalos_max; i++)
   {
      if(H[i].nsub > max_subs_local)
      {
         max_subs_local = H[i].nsub;        
      }
   }

   MPI_Allreduce(&max_subs_local, &max_subs_global, 1, MPI_INT, MPI_MAX, 
               MPI_COMM_WORLD);

   // Loop over every halo
   for(i = 0; i < nhalos_max; i++)
   {
      // Check to see if halo has substructure. If it does, we need to remove
      // duplicate particles before doing the communication. The reason that we need to
      // remove duplicates in the multiple file sets case and not the single file sets
      // case is because a host and all of its subhalos are not necessarily in the same
      // file set. This means that when the plists are sent to root for flagging, there
      // is no guarantee that we'll have the subhalo plist passed to root AFTER the host,
      // as is the case for the single file set. This means a particle that's really in a
      // subhalo might end up having the m_vir of the host assigned to it because that
      // plist just happened to be passed after the sub plist.
      remove_duplicates(i);

      // For each halo, we need to send it's plist to root. This means that the memory to
      // hold the plist needs to be allocated on root. This means root needs to know how 
      // many particles are being sent from each processor. 

      // Start by sending the number of particles in the plist
      MPI_Gatherv(&H[i].npart, 1, MPI_INT, npart_per_plist, recvcnts, displs, MPI_INT, 0,
                  MPI_COMM_WORLD);

      // Get the halo virial masses
      MPI_Gatherv(&H[i].m_vir, 1, MPI_FLOAT, mass_list, recvcnts, displs, MPI_FLOAT, 0, 
                  MPI_COMM_WORLD);

      // Now that we have the number of particles in the current halo, we can allocate 
      // memory for the plist on root
      if(thistask == 0)
      {
         totcounts = 0;

         for(j = 0; j < ntasks; j++)
         {
            totcounts += npart_per_plist[j];
         }

         // Allocate memory for plist
         if(!(plist = calloc(totcounts, sizeof(int))))
         {
            printf("Error, could not allocate memory for plist!\n");
            exit(EXIT_FAILURE);
         }

         // Set up the displacements for the plist
         plist_displs[0] = 0;

         for(j = 1; j < ntasks; j++)
         {
            plist_displs[j] = npart_per_plist[j - 1] + plist_displs[j - 1];
         }
      }

      // Now send the plists
      MPI_Gatherv(&H[i].plist[0], H[i].npart, MPI_INT, plist, npart_per_plist, 
                 plist_displs, MPI_INT, 0, MPI_COMM_WORLD);

      // Loop over plist
      if(thistask == 0)
      {
         // Flag particles
         for(j = 0; j < totcounts; j++)
         {
            if(plist[j] != -1)
            {
               P[plist[j] - 1].in_halo = 1;

               // We can just get the virial mass here if there's only one file set
               if(n_halo_tasks == 1)
               {
                  P[plist[j] - 1].m_vir = H[i].m_vir;
               }
            }
         }

         // Assign virial mass. We loop over every task. For every task, we loop over the
         // number of particles sent by that task. We then assign the corresponding m to
         // those particles. We just need l as an ever increasing index.
         for(j = 0, l = 0; j < ntasks; j++)
         {
            for(k = 0; k < npart_per_plist[j]; k++)
            {
               // Skip the ghostlos
               if(plist[l] != -1)
               {   
                  P[plist[l] - 1].m_vir = mass_list[j];
               }
               l++;
            } 
         }

         // Free the plist for the next halo. If there's only one file set, plist only
         // exists on root
         if(thistask == 0)
         {
            free(plist);
         }
      }
   }

   // Free memory for gatherv arrays
   if(thistask == 0)
   {
      free(npart_per_plist);
      free(recvcnts);
      free(displs);
      free(plist_displs);
      free(mass_list);
   }
}



/***********************
flag_halo_parts_single_file_set
***********************/
void flag_halo_parts_single_file_set(PARTICLE_DATA *P)
{
   // Flags which particles belong to halos and assigns
   // the appropriate m_vir. 
   // NOTE: See flag() for the original way of doing this

   int i;
   int j;

   #ifdef PROFILING
      clock_t start;
      clock_t end;
      double tot_remove_duplicates = 0.0;
      double tot_flag = 0.0;
   #endif

   if(thistask == 0)
   {
      // Loop over every halo
      for(i = 0; i < nhalos_max; i++)
      {
        printf("Flagging particles for halo %ld\n", H[i].hid);

         // Remove duplicates
         #ifdef PROFILING
            // Get start time
            start = clock();
         #endif
         //remove_duplicates_single_set_optimized(i);
         remove_duplicates_single_set(i);

         #ifdef PROFILING
            // Get end time
            end = clock();

            // Add to total time for removing duplicates
            tot_remove_duplicates += (double)(end - start) / CLOCKS_PER_SEC;
         #endif

         // Loop over every particle in the halo
         #ifdef PROFILING
            // Get start time
            start = clock();
         #endif
         for(j = 0; j < H[i].npart; j++)
         {
            // Flag particles
            if(H[i].plist[j] != -1)
            {
               P[H[i].plist[j] - 1].in_halo = 1;
               P[H[i].plist[j] - 1].m_vir = H[i].m_vir;
            }
         }

         #ifdef PROFILING
            // Get end time
            end = clock();

            // Add to total time
            tot_flag += (double)(end - start) / CLOCKS_PER_SEC;
         #endif
      }

      #ifdef PROFILING
         // Print totals
         printf("Total time spent removing duplicates: %e secs\n", tot_remove_duplicates);
         printf("Total time spent flagging: %e secs\n", tot_flag);
      #endif
   }
}



/***********************
          flag
************************/
void flag(PARTICLE_DATA *P)
{
   int i;
   int j;

   for(i = 0; i < nhalos_max; i++)
   {
      if(thistask == 0)
      {
         for(j = 0; j < H[i].npart; j++)
         {
            P[H[i].plist[j] - 1].in_halo = 1;
            P[H[i].plist[j]-1].m_vir = H[i].m_vir;
         }
      }
   }
}



/***********************
   remove_duplicates
***********************/
void remove_duplicates(int current)
{
   // Driver function for removing particles in a subhalo from the host halo's
   // plist. This is only necessary to do because some subhalos reside on a
   // different processor than their host. This means that it's possible for the
   // subhalo to get flagged on root before the host, which would overwrite the
   // subhalo's m_vir with the host's m_vir, which is bad, as it screws up the
   // temperature calculation.

   int i;
   int j;
   int k;
   int l;
   int host_index;
   long int host_id;
   int tag;
   long int tag_long;
   int n_mia_local = 0;
   int npart_in_mia = 0;
   int n_mia_subids_sent;
   int *mia_plist = NULL;
   long int *mia_subids_local = NULL;
   long int *mia_subids_rbuf = NULL;
   MPI_Status status;

   // Allocate for rbuf. This is memory inefficient, but I don't really care at
   // this point. I'll change it if I have to
   if(!(mia_subids_rbuf = calloc(max_subs_global, sizeof(long int))))
   {
      printf("Error, could not allocate memory for mia_subids_rbuf!\n");
      exit(EXIT_FAILURE);
   }

   // Get list of missing subhalo ids on current processor
   mia_subids_local = get_mia_subs(current, &n_mia_local);

   // Need every processor to send
   for(i = 0; i < ntasks; i++)
   {
      // This processor's turn to send
      if(i == thistask)
      {
         // Current processor needs to send to every other processor
         for(j = 0; j < ntasks; j++)
         {
            // Don't send to itself
            if(thistask != j)
            {
               MPI_Send(mia_subids_local, n_mia_local, MPI_LONG, j, 0, 
                        MPI_COMM_WORLD);
            }
         }

         // Now we need to wait for the other processors to send the plists of the
         // missing halos back to current task. Each mia halo might be on a different
         // processor than the others, but there has to be n_mia_local recvs (one for
         // each mia halo)
         for(j = 0; j < n_mia_local; j++)
         {
            // I'm tagging the message with the mia halo id because many different 
            // processors could be sending back here all at the same time, and if 
            // every message has the same tag, does that mess things up? Or does 
            // MPI queue them and then accept them in order?
            
            // The mpi tag must be an int, but the hids are long ints and generally
            // larger than INT_MAX. As such, to keep the tag (not even sure if I
            // have to) corresponding to the halo id, I can subtract off INT_MAX
            // until we're less than INT_MAX.
            tag_long = mia_subids_local[j];
            while(tag_long > INT_MAX)
            {
               tag_long -= INT_MAX;
            }

            // I don't think this cast is necessary, but I don't really care
            tag = tag_long;

            MPI_Recv(&npart_in_mia, 1, MPI_INT, MPI_ANY_SOURCE,
                     tag, MPI_COMM_WORLD, &status);

            // Allocate memory for plist
            if(!(mia_plist = calloc(npart_in_mia, sizeof(int))))
            {
               printf("Error, could not allocate memory for mia_plist!\n");
               exit(EXIT_FAILURE);
            }

            MPI_Recv(mia_plist, npart_in_mia, MPI_INT, MPI_ANY_SOURCE, 
                     tag, MPI_COMM_WORLD, &status);

            // Receive the host id
            MPI_Recv(&host_id, 1, MPI_LONG, MPI_ANY_SOURCE, tag, 
                     MPI_COMM_WORLD, &status);

            // Now that we have the plist of the mia halo, we can go through its host
            // plist and remove the duplicates. 
            // Find host index
            for(k = 0; k < nhalos_local; k++)
            {
               // This does not work, as the mia halo id is not on this processor,
               // which is why I had to send for it in the first place.
               if(host_id == H[k].hid)
               {
                  host_index = k;
                  break;
               }
            }
        
            for(k = 0; k < npart_in_mia; k++)
            {
               for(l = 0; l < H[host_index].npart; l++)
               {
                  if(mia_plist[k] == H[host_index].plist[l])
                  {
                     // Here if there's a duplicate. Flag it as being one
                     H[host_index].plist[l] = -1;
                     break;
                  }
               }
            }

            // Free
            free(mia_plist);
         }
      }

      else
      {
         // All the other processors need to listen for the current processor's list of 
         // mia subids
         MPI_Recv(mia_subids_rbuf, max_subs_global, MPI_LONG, MPI_ANY_SOURCE, 0, 
                  MPI_COMM_WORLD, &status);

         // Now we get the number of mia halos sent
         MPI_Get_count(&status, MPI_LONG, &n_mia_subids_sent);

         // Now search to see if current processor contains any of the sending processor's
         //  mia halos
         for(j = 0; j < n_mia_subids_sent; j++)
         {
            for(k = 0; k < nhalos_local; k++)
            {
               if(mia_subids_rbuf[j] == H[k].hid)
               {
                  // We found one! Now we have to send the plist back to the sending 
                  // processor. MPI_Send sends a copy of the data (as far as I can 
                  // tell from my testing), which means that I don't need to send 
                  // the plist back from the searching processor to here, as the mia 
                  // halo plist will never get modified. It's only the host plist 
                  // (which resides on the searching processor) that gets changed. 


                  // Get remapped id to use as tag. See above
                  tag_long = H[k].hid;

                  while(tag_long > INT_MAX)
                  {
                     tag_long -= INT_MAX;
                  }

                  tag = tag_long;

                  // Send number of particles
                  MPI_Send(&H[k].npart, 1, MPI_INT, status.MPI_SOURCE, tag, 
                           MPI_COMM_WORLD);

                  // Send the plist
                  MPI_Send(H[k].plist, H[k].npart, MPI_INT, status.MPI_SOURCE, 
                           tag, MPI_COMM_WORLD);

                  // Send the host
                  MPI_Send(&H[k].host_id, 1, MPI_LONG, status.MPI_SOURCE, tag, 
                           MPI_COMM_WORLD);

                  break;
               }
            }
         }
      }
   }

   // Free
   free(mia_subids_local);
   free(mia_subids_rbuf);
}



/***********************
      get_mia_subs
***********************/
long int *get_mia_subs(int current, int *n_mia_local)
{
   // Give every processor a list of all of the subhalos that are missing from their
   // home set for the current halo.

   int i;
   int j;
   int *found;
   long int *mia_subids = NULL;

   // Allocate memory for found list
   if(!(found = calloc(H[current].nsub, sizeof(int))))
   {
      printf("Error, could not allocate memory for found list!\n");
      exit(EXIT_FAILURE);
   }

   // Initialize number of missing halos. The idea here is that we start
   // with all of them and then subtract off the number that we find.
   // Whatever is left is the number of subhalos that are missing.
   *n_mia_local = H[current].nsub;

   // Loop over every subhalo
   for(i = 0; i < H[current].nsub; i++)
   {
      // Now loop over every halo in the current processor's set
      for(j = 0; j < nhalos_local; j++)
      {
         // Check for a match
         if(H[j].hid == H[current].sublist[i])
         {
            // We're here if the subhalo is in the file set, so we flag it as found
            // and then decrement the number of missing.
            found[i] = 1;
            *n_mia_local -= 1;
            break;
         }
      }
   }

   // Now allocate memory based on the number of subhalos that are not in the current
   // file set. Then add the halos missing from current processor to the list.
   if(*n_mia_local > 0)
   {
      if(!(mia_subids = calloc(*n_mia_local, sizeof(long int))))
      {
         printf("Error, could not allocate memory for mia_subids!\n");
         exit(EXIT_FAILURE);
      }

      // Now we do a second loop to save the subhalo ids that are not in the current set
      // Loop over every subhalo
      for(i = 0, j = 0; i < H[current].nsub; i++)
      {
         if(found[i] == 0)
         {
            mia_subids[j] = H[current].sublist[i];
            j++;
         }
      }
   }

   free(found);

   return mia_subids;
}



/***********************
remove_duplicates_single_set_optimized
***********************/
void remove_duplicates_single_set_optimized(int current)
{
   // Removes those particles that are in the first subhalo
   // level down from the current halo's plist. This is taken from the serial version
   // of this code. 
   // The halos are sorted in the halos file by mass, and it's impossible
   // for a subhalo to be more massive than its host, which means that the subhalos
   // will always be further down the list (read: at a higher index of H[]) than
   // H[current], so we can start our search for the subhalo at H[current + 1].

   int i;
   int j;
   int k;
   int sub_ind;
   int max_iters;

   // If current is the last halo, then there's no point in looking because this
   // halo (H[current]) cannot have any subhalos since it's last in the list and they're
   // sorted by mass and a subhalo must be less massive than its host, which means it
   // will be futher down the list, which, if this is the last element, doesn't exist
   if(current == nhalos_max - 1)
   {
      return;
   }

   // Set our starting index to the current halo
   sub_ind = current;

   // Loop over the subhalos
   for(i = 0; i < H[current].nsub; i++)
   {
      // The closest the subhalo could be is the next element in the H array
      sub_ind = sub_ind + 1;

      // Reset the max iters counter
      max_iters = 0;

      // Find subhalo
      while(H[sub_ind].hid != H[current].sublist[i])
      {
         sub_ind++;
         max_iters++;

         // Check max iters to avoid potential infinite loop
         if(max_iters > nhalos_max)
         {
            printf("Error, maximum number of iterations reached when removing \
               duplicates! Halo: %d, sublist entry: %d\n", current, i);
            for(i = 0; i < nhalos_max; i++)
            {
               free(H[i].plist);
               free(H[i].sublist);
            }
            free(H);

            exit(EXIT_FAILURE);
         }
      }

      // Loop over every particle in current and compare to subhalo
      for(j = 0, k = 0; j < H[current].npart; j++)
      {
         // We don't need to loop over the whole sub_ind plist. Since they're sorted,
         // if H[sub_ind].plist[k] > H[current].plist[j], then we're done, plist[j] is
         // not in H[sub_ind].plist.
         if(H[current].plist[j] < H[sub_ind].plist[k])
         {
            continue;
         }

         if(H[current].plist[j] == H[sub_ind].plist[k])
         {
            // Flag as a duplicate by changing its id to -1
            H[current].plist[j] = -1;

            // Now set k_start to k + 1 so we don't re-loop over elements that we know
            // cannot 
            k++;
            continue;
         }

         // The only way that H[sub_ind].plist[k] < H[current].plist[j] is if we've
         // reached the end of plist[k] but not plist[j]. That is, if
         // H[current].plist = [1,2,17,19,20,33,43,51,57,68] 
         // and H[sub_ind].plist = [17,20,57], then when k is 2 the pid is 57, but we
         // still have pid 68 left in the host plist. This is our exit condition
         if(H[current].plist[j] > H[sub_ind].plist[k])
         {
            break;
         }
      }

      // We can check to make sure every particle in the subhalo was found
      if(k != H[sub_ind].npart)
      {
         printf("Error, did not remove all duplicates from halo index: %d, \
            sub halo index: %d\n", current, sub_ind);
         exit(EXIT_FAILURE);
      }
   }
}



/***********************
remove_duplicates_single_set
***********************/
void remove_duplicates_single_set(int current)
{
   // Removes those particles that are in the first subhalo
   // level down from the current halo's plist. This is taken from the serial version
   // of this code

   int i;
   int j;
   int k;
   int sub_ind;
    int max_iters;

   // Loop over the subhalos
   for(i = 0; i < H[current].nsub; i++)
   {
      sub_ind = 0;
    max_iters = 0;

      // Find subhalo
      while(H[sub_ind].hid != H[current].sublist[i])
      {
         sub_ind++;
         max_iters++;

         // Check max iters to avoid potential infinite loop
         if(max_iters > nhalos_max)
         {
            printf("Error, maximum number of iterations reached when removing \
               duplicates! Halo: %d, sublist entry: %d\n", current, i);
            for(i = 0; i < nhalos_max; i++)
            {
               free(H[i].plist);
               free(H[i].sublist);
            }
            free(H);

            exit(EXIT_FAILURE);
         }
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
