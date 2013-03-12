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

/*a field with field strength and wave vector k*/
typedef struct {
	k_t k;
	double field_abs;
} field_t;

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

__inline k_t random_spp_wavevector(double k_spp_abs, gsl_rng *r){
	k_t k;
	double theta = random_double_range(r, 0.0, 2.0*M_PI);
	k.kz = 0.0;
	k.kx = k_spp_abs * sin(theta);
	k.ky = k_spp_abs * cos(theta);
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

__inline double spp_radiation_co(double k_spp_abs, double radiation_angle){
	double coefficient;
    double phi = 1400.0e-10;
    double delta_k_abs = 2.0*k_spp_abs*sin(radiation_angle/2.0);
    double exp_co = -1.0/4.0*phi*phi*(delta_k_abs)*(delta_k_abs);
    coefficient = cexp(exp_co)*(1-cos(radiation_angle))*(1-cos(radiation_angle));
	return coefficient;
}

/**
 *
 * This is method to calculate the angle between two scattering process, this method maybe useful after
 * introduce the scattering coefficiency between two scattering process.
 *
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
*/

/**
 * This is method to calculate the angle between two k vector.
 */
__inline double angle_between_vector(k_t k_1, k_t k_2){
	double angle = 0.0;
	double sin_angle = 0.0;
   	sin_angle = (k_1.kx*k_2.kx + k_1.ky*k_2.ky  +k_1.kz*k_2.kz)
   				/(sqrt(k_1.kx*k_1.kx + k_1.ky*k_1.ky  +k_1.kz*k_1.kz )*sqrt(k_2.kx*k_2.kx + k_2.ky*k_2.ky  + k_2.kz*k_2.kz));
   	angle = asin(sin_angle);
   	return angle;
}

/*
 * this is the method to create the index of the next scatter */
__inline int next_scatter_index(double current_scatt_index, int scatterers_number, gsl_rng *r){

	int next_scatter_index = random_int(r,scatterers_number);
	while(next_scatter_index == scatterers_number){
		next_scatter_index = random_int(r,scatterers_number);
	}
	return next_scatter_index;
}

__inline field_t single_field_spp(int first_scatt_index, scatterer_t *scatts, gsl_rng *r){
	/*Constances definition
	 *
	 * Here we define the spp happens on the gold film surface.
	 * */
	double lambda_light = 632.8e-9; /* wavelength of light */
	double k_spp_abs = 2.00329*2.0*M_PI/lambda_light*sin(32.0*M_PI/180.0);/*spp wave vector on LAH79 with angle of theta_sp
	  = 32.0 degree*/
	double mean_free_path = 18.0e-6;/*mean free path of spp*/
	int next_scatterer_index = next_scatter_index(first_scatt_index,sizeof(scatts),r);
	field_t field;
	scatterer_t scatt_cur = scatts[first_scatt_index];
	scatterer_t scatt_prev = scatts[first_scatt_index];
	scatterer_t scatt_next = scatts[next_scatterer_index];
	double pathlength = 0.0;
	double radi_angle = 0.0;
	k_t k_in;
	k_in.kx = 1.0;
	k_in.ky = 0.0;
	k_in.kz = 0.0;
	k_t k_out = random_spp_wavevector(k_spp_abs, r);
	double coefficient = 1.0;
	while(is_radiation_out(pathlength,mean_free_path, r))
	{
		/*
		 * Radiates to the next scatterer.
		 */
		pathlength += distance_two_scatters(scatt_cur, scatt_next);
		scatt_prev = scatt_cur;
		scatt_cur = scatt_next;
		next_scatterer_index = next_scatter_index(first_scatt_index,sizeof(scatts),r);
		//scatt_next = scatts[next_scatter_index];
		//radi_angle = angle_between_scattering(scatt_prev, scatt_cur, scatt_next);
		//coefficient = spp_radiation_co(lambda_spp, lambda_light, radi_angle,r);
		// TODO: the coefficient should be the only first scatter here.
		// it is not correct here
	}
	radi_angle = angle_between_vector(k_in, k_out);
	coefficient = spp_radiation_co(k_spp_abs,radi_angle);
	field.field_abs = sqrt(coefficient)*cexp(1.0i*k_spp_abs*pathlength);
	field.k = k_out;
	return field;
}

/**
 * This is the code for the FFT method for one set of all the scatterers.
 */
__inline int fft_transfer(field_t *field, )



int main(int argc, char **argv) {

	/* program variables */
	  unsigned int i,j; /* variables to iterate over */
	  double x,y,z;
	  double scanxy = 15e-6; /* the boundary of the region */
	  double lambda_light = 632.8e-9; /* wavelength of light */
	  int NSCAT = 500; /*the scatterer number*/
	  int NSA = 360; /* the simulation iteration for the angle for each scatter*/
	  field_t(*field)[NSCAT] = malloc((sizeof *field)*NSA);
	  assert(field!=NULL);

	/* seed the scatterers */
	  scatterer_t *scatts = malloc(NSCAT*sizeof(scatterer_t));
	  assert(scatts!=NULL);

	/* To read the scatterers positions from the file into the scatterer array*/
	  FILE *fp = fopen("/home/jiapeng/Documents/Master thesis with Dr.Frank Vollmer/codes/random_scatterers_creater/Scatterers.txt","r");

	  if (fp == 0) {
		fprintf(stderr, "failed to open the file");
		exit(1);
	  }

	  for (i = 0; i < NSCAT; ++i) {
		fscanf(fp,"%lf %lf %lf",&scatts[i].x,&scatts[i].y,&scatts[i].z);
		fscanf(fp,"\n");
		//printf("%12.12f %12.12f %12.12f\n",scatts[i].x,scatts[i].y,scatts[i].z);
	  }
	  fclose(fp);

	/* initialize the random number generator */
	  int rseed = (int)time(NULL);
	  const gsl_rng_type *T;
	  gsl_rng *r;
	  gsl_rng_env_setup();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc(T);
	  gsl_rng_set(r,rseed);

	/*
	 *
	 */
	  for(i=0;i<NSA;++i){
		  for (j = 0; j < NSCAT; ++j) {
			field[i][j] = single_field_spp(j, scatts, r);
		}
	  }

	  /* For the perpendicular component of the wave_vector*/
	  double k_per_abs = 2.00329*2.0*M_PI/lambda_light*cos(32.0*M_PI/180.0);/*wave vector on LAH79 with angle of theta_sp
	  = 32.0 degree*/
	  k_t k_per;
	  k_per.kx = 0.0;
	  k_per.ky = 0.0;
	  k_per.kz = k_per_abs;

	  /*
	   * Add the perpendicular wave-vector on the field.
	   */
	  for(i=0;i<NSA;++i){
		  for (j = 0; j < NSCAT; ++j) {
			field[i][j].k.kz += k_per.kz;
		}
	  }	  printf("%s\n","the program is finished");
	  return 0;
}


