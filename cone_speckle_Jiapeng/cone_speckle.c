/**
 * This program is to generate the cone speckle use the monte-carlo simulation.
 *
 * @author - Jiapeng Huang - MPL - 2013.3.8
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <malloc.h>
#include <sys/types.h>
#include <float.h>
#include <errno.h>
#include <error.h>
#include <hdf5.h> /* hdf5 output */
#include <ctype.h>
#include <unistd.h> /* command line arguments, get memory */
#include <getopt.h> /* command line arguments */
#include <assert.h>
#include <complex.h>
#include <gsl/gsl_rng.h> /* random number generator using mt19937 */
#include <errno.h>

/* return a random double */
__inline double random_double_range(gsl_rng *r,double min, double max)
{
	return ((max-min)*gsl_rng_uniform(r)) + min;
}

/* return a random int */
__inline int random_int(gsl_rng *r,int max){
  return gsl_rng_uniform_int(r,max+1);
}

/* define the camera structure*/
typedef struct{
  double cam_lx; /* size of camera ccd x*/
  double cam_ly; /* size of camera ccd y*/
  unsigned int cam_sx; /*camera pixels x*/
  unsigned int cam_sy; /*camera pixels y*/
} camera_t;

/* a scatterer */
typedef struct {
  double x,y,z; /* scatterer position */
} scatterer_t;

/*a wave-number vector k*/
typedef struct {
	double kx,ky,kz;
} k_t;

/**
 * This is a method to get the field on a ccd pixel of a scatterer's scattering light. By calculating ray trace,
 * the method will return the index of the ccd pixel.
 */
__inline int target_pixel_index(scatterer_t scatter, k_t k, camera_t cam, double distance){
	int index = 0;
	return index;
}

/**
 * This is a method to calculate the distance between a certain ccd pixel and a certain scatterer
 *
 */
__inline double distance_ccd_scatter(scatterer_t scatter, camera_t cam, int ccd_index){
	double distance = 0.0;
	return distance;
}

/* Calculate the distance between two scatters*/
__inline int distance_two_scatters(scatterer_t scatt_a,scatterer_t scatt_b){
	int distance = 0.0;
	distance = sqrt((scatt_a.x - scatt_b.x)*(scatt_a.x - scatt_b.x) + ((scatt_a.y - scatt_b.y)*(scatt_a.y - scatt_b.y)) + ((scatt_a.z - scatt_b.z)*(scatt_a.z - scatt_b.z)));
	return distance;
}
/**
 *
 * This is the method to create the exponential probability function for the spp multiple scattering.
 *
 */
__inline bool is_radiation_out(double pathlength, double mean_free_pathlength,gsl_rng *r){
	double random_probability = 0.0;
	double out_probability = 0.0;
	bool out = false;
	random_probability = random_double_range(r, 0.0, 1.0);
	out_probability = cexp(-(pathlength/mean_free_pathlength));
	if(random_probability >= out_probability){
		out = true;
	}
	else {
		out = false;
	}
	return out;
}
/**
 *
 * This is the method to create a random wave vector for the SPP.
 * In this case, since the SPP are confined in the metal-dielectric interface. The wave vector are considered as
 * (kx,ky,0), the kz component should be zero
 * // TODO: is this true, kz = 0?
 */

__inline k_t random_spp_wavevector(double lambda_spp, gsl_rng *r){
	k_t k;
	double k_absolute = 2.0 * M_PI /lambda_spp;
	double theta = random_double_range(r, 0.0, 2.0*M_PI);
	k.kz = 0.0;
	k.kx = k_absolute * sin(theta);
	k.ky = k_absolute * cos(theta);
	return k;
}

/**
 *
 * This is a method to create the angular distribution of the spp radiation. According to the distribution function at
 * the paper, "directional surface plasmon scattering from silver films" , H.J. SIMON and J.K.GUHA.
 *
 * //TODO: In this paper, the equation about the angular performance and the correlation function about the surface roughness
 *  spectrum.
 */

__inline double spp_radiation_co(double lambda_spp, double lambda_light, double radiation_angle,gsl_rng *r){
	double coefficient;
    double k_spp_abs = 2.0 * M_PI /lambda_spp;
    double k_light_abs = 2.0 * M_PI /lambda_light;
    double phi = 1400e-10;
    coefficient = cexp(-1/4*phi*phi*(k_spp_abs*k_spp_abs + k_light_abs*k_light_abs -2*k_light_abs*k_spp_abs
    		*cos(radiation_angle))*(1-cos(radiation_angle)));
	return coefficient;
}

__inline double angle_between_scattering(scatterer_t scatt_prev, scatterer_t scatt_curr, scatterer_t scatt_next){
	double angle = 0.0;
	double sin_angle = 0.0;
	double vector_a_x = scatt_prev.x - scatt_curr.x;
	double vector_a_y = scatt_prev.y - scatt_curr.y;
	double vector_a_z = scatt_prev.z - scatt_curr.z;
    double vector_b_x = scatt_next.x - scatt_curr.x;
    double vector_b_y = scatt_next.y - scatt_curr.y;
    double vector_b_z = scatt_next.z - scatt_curr.z;
    sin_angle = (vector_a_x*vector_b_x + vector_a_y*vector_b_y +vector_a_z*vector_b_z)
    		/(sqrt(vector_a_x*vector_a_x + vector_a_y*vector_a_y +vector_a_z*vector_a_z)*sqrt(vector_b_x*vector_b_x + vector_b_y*vector_b_y +vector_b_z*vector_b_z));
    angle = asin(sin_angle);
    return angle;
}

__inline int next_scatter_index(double current_scatt_index, int scatterers_number, gsl_rng *r){

	int next_scatter_index = random_int(r,scatterers_number);
	while(next_scatter_index == scatterers_number){
		next_scatter_index = random_int(r,scatterers_number);
	}
	return next_scatter_index;
}

__inline complex double single_field_spp(int first_scatt_index, scatterer_t *scatts, gsl_rng *r){
	/*Constances definition
	 *
	 * Here we define the spp happens on the gold film surface.
	 * */
	double lambda_light = 632.8e-9; /* wavelength of light */
	double lambda_spp = 610.0e-9;/*wavelength of spp*/
	double mean_free_path = 18.0e-6;/*mean free path of spp*/
	int next_scatterer_index = next_scatter_index(first_scatt_index,sizeof(scatts),r);
	complex double field;
	scatterer_t scatt_cur = scatts[first_scatt_index];
	scatterer_t scatt_prev = scatts[first_scatt_index];
	scatterer_t scatt_next = scatts[next_scatterer_index];
	double pathlength = 0.0;
	double radi_angle = 0.0;
	double coefficient = spp_radiation_co(lambda_spp, lambda_light, radi_angle,r);
	while(is_radiation_out(pathlength,mean_free_path, r))
	{
		/*
		 * Radiates to the next scatterer.
		 */
		pathlength += distance_two_scatters(scatt_cur, scatt_next);
		scatt_prev = scatt_cur;
		scatt_cur = scatt_next;
		next_scatterer_index = next_scatter_index(first_scatt_index,sizeof(scatts),r);
		scatt_next = scatts[next_scatter_index];
		radi_angle = angle_between_scattering(scatt_prev, scatt_cur, scatt_next);
		coefficient = spp_radiation_co(lambda_spp, lambda_light, radi_angle,r);
		// TODO: the coefficient should be the only first scatter here.
		// it is not correct here
	}
	field = coefficient*cexp(2.0i*M_PI/lambda_spp*pathlength);

	return field;
}

int main(int argc, char **argv) {

	/* program variables */
	  unsigned int i; /* variables to iterate over */
	  double x,y,z;
	  double scanxy = 15e-6; /* the boundary of the region */
	  int NSCAT = 500; /*the scatterer number*/

	  /* seed the scatterers */
	  scatterer_t *scatt = malloc(NSCAT*sizeof(scatterer_t));
	  assert(scatt!=NULL);

	  /* To read the scatterers positions from the file into the scatterer array*/
	  FILE *fp = fopen("/home/jiapeng/Documents/Master thesis with Dr.Frank Vollmer/codes/random_scatterers_creater/Scatterers.txt","r");

	  if (fp == 0) {
		fprintf(stderr, "failed to open the file");
		exit(1);
	  }

	  for (i = 0; i < NSCAT; ++i) {
		fscanf(fp,"%lf %lf %lf",&scatt[i].x,&scatt[i].y,&scatt[i].z);
		fscanf(fp,"\n");
		printf("%12.12f %12.12f %12.12f\n",scatt[i].x,scatt[i].y,scatt[i].z);
	  }
	  fclose(fp);


	  return 0;
}


