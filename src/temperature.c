/************************************************
Title: temperature.c
Purpose: Contains functions related to calculating
         the temperature of each particle
Notes: Based on Bertone's code
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "allvars.h"
#include "proto.h"



/***********************
   get_temperatures
***********************/
PARTICLE_DATA *get_temperatures(PARTICLE_DATA *All_P, int ngas, int *n_this_task)
{
   // Does as the name says

   int i;
   float T0;
   float rho_c;
   float rho_mean;
   float rho_b;
   float r_vir;
   float v_vir;
   float rho_phys;
   float G = 6.67e-8;
   PARTICLE_DATA *P;

   // Get a_dot
   get_a_dot();

   // Get T0
   if(thistask == 0)
   {
      T0 = get_T0();
   }

   MPI_Bcast(&T0, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

   // Get the critical density of the Universe (cgs units)
   rho_c = 3.0 * pow(a_dot / header.time, 2.0) / (8.0 * M_PI * G);

   // This is Omega_m(z) = Omega_m0 * (1+z)^3 / E^2, E = H(z) / H0. See Mo02 and notes from
   // 12/8/16. rho_mean is now in cgs
   rho_mean = ((OMEGA_M0 * pow(1. + header.redshift, 3.0) * pow(3.2407e-18 * LITTLE_H, 2.0)) / 
              pow(a_dot / header.time, 2.0)) * rho_c;

   // Get the average baryon density by multiplying by baryon_frac
   rho_b = BARYON_FRAC * rho_mean;

   // Divide particles amongst the processors
   P = split_particles(All_P, ngas, n_this_task);

   for(i = 0; i < *n_this_task; i++)
   {
      // Check to see if particle is in halo
      if(P[i].in_halo == 1)
      {
         // Get r_vir (Bertone thesis eq 2.10) G is in cgs, so put cm^3 to km^3 and 
         // g to M_sun. This puts r_vir in km. See notes 3/11/16
         r_vir = pow((1.989e18 * G * P[i].m_vir) / (100.0 * LITTLE_H * (a_dot * a_dot 
                 / (header.time * header.time))), 1.0 / 3.0);

         // Get V_vir (Bertone thesis eq 2.11) in km/s (1.989e18 is the conversion)
         v_vir = sqrt(1.989e18 * G * P[i].m_vir / (r_vir * LITTLE_H));

         // Get the temp (Bertone thesis eq 2.12) (Kelvin) (1e10 is the conversion)
         P[i].temp = 1e10 * MOL_WEIGHT * PROTON_MASS * v_vir * v_vir / (2.0 * BOLTZMANN);
      }

      // Not in halo
      else
      {
         // Bertone thesis eq. 1.4 (Kelvin. Need factor for putting P.dens -> g/cm^3)
         rho_phys = P[i].density / (header.time * header.time * header.time);
         P[i].temp = T0 * pow(rho_phys * BARYON_FRAC * LITTLE_H * LITTLE_H * GUM_IN_G / (rho_b * pow(GUL_IN_CM, 3.0)), 1.0 / 1.7);
      }
   }

   return P;
}



/***********************
        get_T0
***********************/
float get_T0(void)
{
   // This function reads Bolton's table and interpolates to
   // the correct value of T0 for the current redshift

   FILE *fd;
   int i;
   int nlines = 0;
   double dummy;
   float T;
   double *z;
   double *temp;
   gsl_interp_accel *accl;
   gsl_spline *spline;

   // Open the file
   if(!(fd = fopen("./temp_S3.dat", "r")))
   {
      printf("Error, could not open Bolton's file for reading T0!\n");
      exit(EXIT_FAILURE);
   }

   // Get the number of lines in the file
   while((i = getc(fd)) != EOF)
   {
      if(i == '\n')
      {
         nlines++;
      }
   }

   // Get rid of heading line
   nlines = nlines - 1;

   // Go back to beginning of file
   rewind(fd);

   // Skip heading line
   while((i = getc(fd)) != '\n')
   {
      continue;
   }

   // Allocate memory for z and T0 arrays
   if(!(z = calloc(nlines, sizeof(double))))
   {
      printf("Error, could not allocate memory for z array in get_T0!\n");
      exit(EXIT_FAILURE);
   }
   
   if(!(temp = calloc(nlines, sizeof(double))))
   {
      printf("Error, could not allocate memory for T0 array in get_T0!\n");
      exit(EXIT_FAILURE);
   }

   // Read z and T0. Order of data is: z nHI/nH nHeI/nH nHeII/nH T0(K)
   for(i = 0; i < nlines; i++)
   {
      fscanf(fd, "%lf %lf %lf %lf %lf\n", &z[i], &dummy, &dummy, &dummy, &temp[i]);
   }

   // The gsl requires the x values (redshifts, in this case) to be in strictly
   // increasing order, but the temp file goes from past (high z) to present (low
   // z). So I'm going to convert them to scale factors (which increase from past
   // to present) so the gsl can do its thing
   for(i = 0; i < nlines; i++)
   {
      z[i] = 1.0 / (1.0 + z[i]);
   }

   // Do interpolation to correct value of T0 using gsl
   accl = gsl_interp_accel_alloc();
   spline = gsl_spline_alloc(gsl_interp_cspline, nlines);
   gsl_spline_init(spline, z, temp, nlines);
   T = gsl_spline_eval(spline, 1.0 / (1.0 + header.redshift), accl);

   gsl_spline_free(spline);
   gsl_interp_accel_free(accl);

   return T;
}



/***********************
    split_particles
***********************/
PARTICLE_DATA *split_particles(PARTICLE_DATA *All_P, int ngas, int *n_this_task)
{
   int i;
   int n_to_send;
   int nleft;
   int *p_displs;
   int *p_sendcnts;
   PARTICLE_DATA *p_rbuf;
   
   // Broadcast ngas
   MPI_Bcast(&ngas, 1, MPI_INT, 0, MPI_COMM_WORLD);

   // Get number of particles per processor
   n_to_send = floor((float)ngas / (float)ntasks);

   // Figure out how many particles are left
   nleft = ngas - (ntasks * n_to_send);

   // Append these particles to the last processor
   if(thistask == (ntasks - 1))
   {
      n_to_send += nleft;
   }

   // Create the displs, rbuf, and sendcounts
   if(thistask == 0)
   {
      if(!(p_displs = calloc(ntasks, sizeof(int))))
      {
         printf("Error, could not allocate memory for p_displs!\n");
         exit(EXIT_FAILURE);
      }

      if(!(p_sendcnts = calloc(ntasks, sizeof(int))))
      {
         printf("Error, could not allocate memory for p_sendcnts!\n");
         exit(EXIT_FAILURE);
      }
   }

   if(!(p_rbuf = calloc(n_to_send, sizeof(PARTICLE_DATA))))
   {
      printf("Error, could not allocate memory for p_rbuf!\n");
      exit(EXIT_FAILURE); 
   }

   // Now tell root how many particles to send to each processor
   MPI_Gather(&n_to_send, 1, MPI_INT, p_sendcnts, 1, MPI_INT, 0, MPI_COMM_WORLD);

   // Now set up the displacement array
   if(thistask == 0)
   {
      p_displs[0] = 0;

      for(i = 1; i < ntasks; i++)
      {
         p_displs[i] = p_displs[i - 1] + p_sendcnts[i - 1];
      }
   }

   // Send the particles to the other processors
   MPI_Scatterv(All_P, p_sendcnts, p_displs, mpi_particle_type, p_rbuf, n_to_send, 
                mpi_particle_type, 0, MPI_COMM_WORLD);

   // Free memory
   if(thistask == 0)
   {
      free(p_displs);
      free(p_sendcnts);
   }

   // Update n_this_task
   *n_this_task = n_to_send;

   // We no longer need All_P on root, and so we free it
   if(thistask == 0)
   {
      free(All_P);
   }

   return p_rbuf;
}
