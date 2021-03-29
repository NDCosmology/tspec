/************************************************
Title: proto.h
Purpose: Contains all function prototypes
Notes: Organized by where function is defined
************************************************/
#ifndef ALLVARS_H
   #include "allvars.h"
#endif



/***********************
         de.c
***********************/
void de_error_check(void);
void get_a_dot(void);



/***********************
       debugging
***********************/
void write_flagged_particles(PARTICLE_DATA *, int);



/***********************
         flag.c
***********************/
void flag_halo_parts(PARTICLE_DATA *);
void flag_halo_parts_mult_file_sets(PARTICLE_DATA *);
void flag_halo_parts_single_file_set(PARTICLE_DATA *);
void flag(PARTICLE_DATA *);
void remove_duplicates(int);
long int *get_mia_subs(int, int *);
void remove_duplicates_single_set_optimized(int);
void remove_duplicates_single_set(int);



/***********************
        halos.c
***********************/
void load_halos(void);
void read_halo_particles(void);
void read_halo_substruct(void);
void read_virial_mass(void);
int cmpfunc(const void *, const void *);
char *get_halo_fname(char *);



/***********************
         init
***********************/
void init(int, char *);
void make_custom_mpi_type(void);



/***********************
        load.c
***********************/
IO_HEADER load_header(void);
PARTICLE_DATA *load_snapshot(int *);
void block_check(enum fields, int, int, IO_HEADER);
int get_block_size(enum fields, IO_HEADER);
size_t my_fread(void *, size_t, size_t, FILE *);
int pid_cmp(const void *, const void *);



/***********************
     temperature.c
***********************/
PARTICLE_DATA *get_temperatures(PARTICLE_DATA *, int, int *);
float get_T0(void);
PARTICLE_DATA *split_particles(PARTICLE_DATA *, int, int *);



/***********************
        write.c
***********************/
void write_particle_data(PARTICLE_DATA *, int);
size_t my_fwrite(void *, size_t, size_t, FILE *);
