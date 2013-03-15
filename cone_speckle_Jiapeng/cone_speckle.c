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
	complex double field_abs;
} field_t;

/* a struct stands for the index of lower and upper bounds*/
typedef struct{
	int index;
	double lower_bound;
	double upper_bound;
} index_bounds;


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
__inline double distance_two_scatters(scatterer_t scatt_a,scatterer_t scatt_b){
	double distance = 0.0;
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
    coefficient =(k_spp_abs*k_spp_abs*k_spp_abs*k_spp_abs)* cexp(exp_co);//*(1-cos(radiation_angle))*(1-cos(radiation_angle));
	return coefficient;
}

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
__inline int next_scatter_index(int current_index,int NSCAT, gsl_rng *r, index_bounds *matrix){
	int next_scatter_index,i;
	double random_probability = 0.0;
	random_probability = random_double_range(r, 0.0, 1.0);
	for (i = 0; i < NSCAT-1; ++i) {
		if (random_probability<matrix[i].upper_bound && random_probability > matrix[i].lower_bound) {
			next_scatter_index = matrix[i].index;
		}
	}
	return next_scatter_index;
}

/**
 * This is a method to create a conductance matrix
 *
 */
__inline int conduc_matrix(int NSCAT,scatterer_t *scatts,int scatt_index, index_bounds *matrix){
	double coeff[NSCAT-1];
	int i,j;
	double pathlength, coeff_sum;
	for (i = 0; i < NSCAT; ++i) {
		if (i!=scatt_index) {
			pathlength = distance_two_scatters(scatts[i],scatts[scatt_index]);
			pathlength *= 1.0e+6;
			coeff[i] = 1/(pathlength*pathlength*pathlength);
			coeff_sum += coeff[i];
		}
	}
	for (i = 0; i < NSCAT; ++i) {
		if (i!=scatt_index) {
			coeff[i] = coeff[i]/coeff_sum ;
		}
	}

	for (i = 0; i < NSCAT; ++i) {
		if (i!=scatt_index) {
			//printf("%f   %d\n", coeff[i],i) ;
		}
	}

	coeff_sum = 0.0;
	for (i = 0,j=0; i < NSCAT; ++i,++j) {
		if (i!=scatt_index) {

			matrix[j].index = i;
			matrix[j].lower_bound = coeff_sum;
			coeff_sum += coeff[i];
			matrix[j].upper_bound = coeff_sum;
		} else {
			--j;
		}
	}

	return 0;
}



__inline field_t single_field_spp(int first_scatt_index, scatterer_t *scatts,int NSCAT, gsl_rng *r, double *free_path, int free_path_index,index_bounds *matrix){
	/*Constances definition
	 *
	 * Here we define the spp happens on the gold film surface.
	 * */
	double lambda_light = 632.8e-9; /* wavelength of light */
	double k_spp_abs = 2.00329*2.0*M_PI/lambda_light*sin(32.0*M_PI/180.0);/*spp wave vector on LAH79 with angle of theta_sp
	  = 32.0 degree*/
	double mean_free_path = 18.0e-6;/*mean free path of spp*/
	int next_scatterer_index = next_scatter_index(first_scatt_index,NSCAT,r,matrix);
	int curr_scatterer_index = first_scatt_index;
	field_t field;
	scatterer_t scatt_cur = scatts[first_scatt_index];
	scatterer_t scatt_next = scatts[next_scatterer_index];
	double pathlength = 0.0;
	double radi_angle = 0.0;
	k_t k_in;
	k_in.kx = 1.0;
	k_in.ky = 0.0;
	k_in.kz = 0.0;
	k_t k_out = random_spp_wavevector(k_spp_abs, r);
	double coefficient = 1.0;
	while(!is_radiation_out(pathlength,mean_free_path, r))
	{
		/*
		 * Radiates to the next scatterer.
		 */
		pathlength += distance_two_scatters(scatt_cur, scatt_next);
		scatt_cur = scatt_next;
		curr_scatterer_index = next_scatterer_index;
		next_scatterer_index = next_scatter_index(curr_scatterer_index,NSCAT,r,matrix);
		scatt_next = scatts[next_scatterer_index];
		//radi_angle = angle_between_scattering(scatt_prev, scatt_cur, scatt_next);
		//coefficient = spp_radiation_co(lambda_spp, lambda_light, radi_angle,r);
		// TODO: the coefficient should be the only first scatter here.
		// it is not correct here
	}
	free_path[free_path_index] = pathlength;
	radi_angle = angle_between_vector(k_in, k_out);
	coefficient = spp_radiation_co(k_spp_abs,radi_angle);
	field.field_abs = sqrt(coefficient)*cexp(1.0i*k_spp_abs*pathlength);
	field.k = k_out;
	return field;
}

/**
 * This is the method for the free space transform  for one set of all the scatterers.
 */
__inline scatterer_t *fst_transfer(field_t *field, double z0, scatterer_t *scatts_orig,int NSCAT, scatterer_t *locations){
	int i;
	for (i = 0; i < NSCAT; ++i) {
		locations[i].x = scatts_orig[i].x + field[i].k.kx/field[i].k.kz*z0;
		locations[i].y = scatts_orig[i].y + field[i].k.ky/field[i].k.kz*z0;
		locations[i].z = scatts_orig[i].z + z0;
	}
	return locations;
}

/**
 * This is the method to create the field distribution on the ccd chip from a certain distance for the orignal
 * scattering plane. Here we put the ccd chip parallize the x-y plane and has a distance from the orignal point
 * of z0.
 */
__inline complex double* field_on_ccd(double distance, scatterer_t *locations, field_t *field, camera_t cam,int NSCAT, complex double *ccd){
	  double dx = cam.cam_lx/(cam.cam_sx-1);
	  double dy = cam.cam_ly/(cam.cam_sy-1);
	  int cam_sx = cam.cam_sx;
	  int cam_sy = cam.cam_sy;
	  assert(ccd!=NULL);
	  int i,pixel_index;
	  int size = NSCAT;
	  for (i = 0; i < size; ++i) {
		  // TODO: double check
		  pixel_index = ((int)((locations[i].x+cam.cam_lx/2)/dx)*cam_sy) + (int)((locations[i].y+cam.cam_ly/2)/dy);
		  ccd[pixel_index] += field[i].field_abs;
	}
	  return ccd;
}


int main(int argc, char **argv) {

	/* program variables */
	  unsigned int i,j,k; /* variables to iterate over */
	  double x,y,z;
	  double scanxy = 15e-6; /* the boundary of the region */
	  double lambda_light = 632.8e-9; /* wavelength of light */
	  int NSCAT = 100; /*the scatterer number*/
	  int NSA = 36; /* the simulation iteration for the angle for each scatter*/
	  camera_t cam = {0.20,0.20,512,512}; /*initialize the cam struct*/
	  double z0 = 0.13; /* the camera distance from the orginal plane*/
	  field_t(*field)[NSCAT] = malloc((sizeof *field)*NSA);
	  assert(field!=NULL);
	/* seed the scatterers */
	  scatterer_t *scatts = malloc(NSCAT*sizeof(scatterer_t));
	  assert(scatts!=NULL);
	/* the locations */
	  scatterer_t *locations = malloc(NSCAT*sizeof(scatterer_t));
	  assert(locations!=NULL);
	/*free path array*/
	  double *free_path = malloc(NSA*NSCAT*sizeof(double));

	  complex double *ccd_all = malloc(cam.cam_sx*cam.cam_sy*sizeof(complex double));
	  bzero(ccd_all,cam.cam_sx*cam.cam_sy*sizeof(complex double));
	  complex double *ccd_one = malloc(cam.cam_sx*cam.cam_sy*sizeof(complex double));
	  bzero(ccd_one,cam.cam_sx*cam.cam_sy*sizeof(complex double));


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

	  index_bounds (*Matrix)[NSCAT-1] = malloc((sizeof *Matrix)*NSCAT);
	  index_bounds *matrix;
	  for (i = 0; i < NSCAT; ++i) {
		  matrix = Matrix[i];
		  conduc_matrix(NSCAT, scatts, i, matrix);

		for (k = 0; k < NSCAT-1; ++k) {
			//printf("%12.12f\n",Matrix[i][k].lower_bound);
		}
		printf("%d\n",i);
	  }

	  printf("%s","conductive matrix finished");
	  //scatts[500].x = 5.0e-6;
	  //scatts[500].y = 4.0e-6;
	  //scatts[500].z = 0.2e-6;

	/* initialize the random number generator */
	  int rseed = (int)time(NULL);
	  const gsl_rng_type *T;
	  gsl_rng *r;
	  gsl_rng_env_setup();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc(T);
	  gsl_rng_set(r,rseed);

	/*
	 *Create the path length array
	 */
	  for(i=0;i<NSA;++i){
		  for (j = 0; j < NSCAT; ++j) {
			matrix = Matrix[j];
			field[i][j] = single_field_spp(j, scatts, NSCAT, r, free_path,(i*NSCAT + j),matrix);
		}

	  }

	  /* For the perpendicular component of the wave_vector*/
	  double k_per_abs = 2.00329*2.0*M_PI/lambda_light*cos(32.0*M_PI/180.0);/*wave vector on LAH79 with angle of theta_sp
	  = 32.0 degree*/
	  field_t *field_one = malloc(NSCAT*sizeof(field_t));
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
			field_one[j] = field[i][j];
		}

		  locations = fst_transfer(field_one,z0,scatts,NSCAT, locations);
		  ccd_one = field_on_ccd(z0,locations,field_one,cam,NSCAT,ccd_one);

		  for (k = 0; k < cam.cam_sx*cam.cam_sy; ++k) {
			  ccd_all[k] += ccd_one[k];
		}
	  }

	  printf("%s\n","the program is finished");

	  /* output file */
	    hid_t file,dataset,dataspace;
	    herr_t status;
	    hsize_t dims[2]; /* dimensionality of the set */
	    dims[0]=cam.cam_sx;
	    dims[1]=cam.cam_sy;
	    dataspace = H5Screate_simple(2,dims,NULL);
	    file = H5Fcreate("out_2.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
	    dataset = H5Dcreate1(file,"/e2",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);

	    double *tmp = malloc(cam.cam_sx*cam.cam_sy*sizeof(double));
	    assert(tmp!=NULL);

	  /* normalize tmp, first by finding max */
	    double max = 0;
	    for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
	      tmp[i]=cabs(ccd_all[i])*cabs(ccd_all[i]);
	      if(tmp[i]>max){ max=tmp[i]; }

	    }
	  //  printf("%g\n",max);
	  /* then dividing by max */
	    for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
	  	  tmp[i]/=max;
	    }

	    /**
	     * Write the tmp array as a way to do the cross-correlation function.
	     */
	    FILE *fp_tmp = fopen("tmp_500.txt","w");

	    if (fp_tmp == 0) {
	  	fprintf(stderr, "failed to open the file");
	  	exit(1);
	    }

	    int size_none_zero = 0;

	    for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
	    	if(tmp[i]!=0.0){
	    		++size_none_zero;
	    	}
	    }

	    double *tmp_N = malloc(size_none_zero*sizeof(double));
	    j = 0;
	    for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
	    	if(tmp[i]!=0.0){
	    		tmp_N[j] = tmp[i];
	    		j++;
	    	}
	    }
	    for (i = 0; i < size_none_zero; ++i) {
	  	fprintf(fp_tmp,"%-12.12f\n",tmp_N[i]);
	    }
	    fclose(fp_tmp);


	    /**
	     * Write the free path array as a way to do the cross-correlation function.
	     */
	    FILE *fp_path = fopen("free_path.txt","w");

	    if (fp_path == 0) {
	  	fprintf(stderr, "failed to open the file");
	  	exit(1);
	    }
	    for (i = 0; i < NSA*NSCAT; ++i) {
	  	fprintf(fp_path,"%-12.12f\n",free_path[i]);
	    }
	    fclose(fp_path);


	    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

	  /* write the real and imaginary parts as well */
	    dataset = H5Dcreate1(file,"/er",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);
	    for(i=0;i<cam.cam_sx*cam.cam_sy;++i){ tmp[i]=creal(ccd_all[i]); }
	    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

	    dataset = H5Dcreate1(file,"/ei",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);
	    for(i=0;i<cam.cam_sx*cam.cam_sy;++i){ tmp[i]=cimag(ccd_all[i]); }
	    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

	  /* clean up */
	    status = H5Dclose(dataset);
	    status = H5Sclose(dataspace);
	    status = H5Fclose(file);

	    free(tmp);
	    free(ccd_one);
	    free(ccd_all);
	    free(scatts);
	    free(locations);

	    printf("%s","finished");
	  return 0;
}


