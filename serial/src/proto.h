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
         flag.c
***********************/
void flag_halo_parts(PARTICLE_DATA *, int);
void remove_duplicates(int);



/***********************
        halos.c
***********************/
void load_halos(void);
void read_halo_particles(void);
void read_halo_substruct(void);
void read_virial_mass(void);
int cmpfunc(const void *, const void *);



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
void get_temperatures(PARTICLE_DATA *, int);
float get_T0(void);



/***********************
        write.c
***********************/
void write(PARTICLE_DATA *, int);
void write_hsml(PARTICLE_DATA *, int);
void write_dens(PARTICLE_DATA *, int);
size_t my_fwrite(void *, size_t, size_t, FILE *);
