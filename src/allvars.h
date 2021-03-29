/************************************************
Title: allvars.h
Purpose: Contains all global variable declarations
         and structure definitions.
Notes:
************************************************/
#include <mpi.h>

#ifndef ALLVARS_H
   #define ALLVARS_H

    // MPI
    extern int thistask;
    extern int ntasks;
    extern MPI_Datatype mpi_header_type;
    extern MPI_Datatype mpi_particle_type;

   // Files
   extern char snapfile[256];  // Holds file name for snapshot
   extern char part_file[100]; // Holds amiga particles data
   extern char halo_file[100]; // Holds amiga halo data
   extern char subfile[100];   // Holds amiga substructure data

   // Cosmology
   extern float a_dot;
   extern int nhalos_local;    // Number of halos in a give file set
   extern int nhalos_tot;      // Total number of halos across all file sets
   extern int nhalos_max;      // Max number of halos held by a single processor
   extern int max_subs_global; // The max number of subhalos held by a single host
   extern float LITTLE_H;
   extern float OMEGA_R0;
   extern float OMEGA_K0;
   extern float OMEGA_M0;
   extern float OMEGA_DE0;
   extern float BARYON_FRAC;
   extern float PROTON_MASS;
   extern float MOL_WEIGHT;
   extern float BOLTZMANN;
   

   // Units
   extern double GUL_IN_CM;
   extern double GUM_IN_G;
   extern double GUV_IN_CM_PER_S;

   extern int de_flag;

   extern float DE_W0;
   extern float DE_WA;

   extern int n_halo_tasks;   // Because of the memory leak in AHF, large sims only run
                              // in serial. If the sim is smaller then AHF and tspec will
                              // be run with the same number of processors and all is
                              // right with the world. However, if there's just one file
                              // set then I need to know so I can read and distribute the
                              // particles appropriately. This variable therefore holds
                              // the number of AHF file sets, which is the number of
                              // cores AHF was run with.

   // Fields for block checking
   enum fields
   {
      HEADER,
      POS,
      VEL,
      IDS,
      MASS,
      RHO,
      HSML
   };   

   // Snapshot header
   typedef struct IO_HEADER
   {
      int npart[6];
      double mass[6];
      double time;
      double redshift;
      int flag_sfr;
      int flag_feedback;
      int npartTotal[6];
      int flag_cooling;
      int num_files;
      double BoxSize;
      double Omega0;
      double OmegaLambda;
      double HubbleParam;
      char fill[96];
   } IO_HEADER;

   // Particle data
   typedef struct PARTICLE_DATA
   {
      float pos[3];
      float vel[3];
      float mass;
      float density;
      float temp;
      float hsml;
      float m_vir;
      int type;
      int id;
      int in_halo;
   } PARTICLE_DATA;

   // Halos struct
   typedef struct HALO_DATA
   {
      int npart;         // Number of particles in halo
      int new_id;        // Ressigned halo id for use as an mpi tag
      long int hid;      // halo id
      int *plist;        // Holds list of particle ids that are in halo
      int nsub;          // Number of sub halos the halo has
      long int *sublist; // List of sub halo ids
      float m_vir;       // halo's virial mass
      long int host_id;  // ID of halo's host. 0 if there is no host 
   } HALO_DATA;

   // Global structures
   extern IO_HEADER header;   // Master header for the snapshot
   extern HALO_DATA *H;       // Structure for holding all the halo info
#endif
