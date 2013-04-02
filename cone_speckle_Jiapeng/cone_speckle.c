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
#include <hdf5.h>        /* hdf5 output */
#include <time.h>
#include <ctype.h>
#include <unistd.h>      /* command line arguments, get memory */
#include <getopt.h>      /* command line arguments */
#include <assert.h>
#include <complex.h>
#include <gsl/gsl_rng.h> /* random number generator using mt19937 */

#include <errno.h>
#include <error.h>
#include "zed-debug.h"   /* useful debug macros */

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
	int scatterer_index; // stands for which scatterer generate the field;
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
bool is_radiation_out(double pathlength, double mean_free_pathlength,gsl_rng *r){
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
 * This is the method to create a Gaussian distribution wave vector for the SPP.
 * In this case, the wave vector of the spp is depended on the postion of scatterer. Here, we define the (0,0) in x-,y- axis is the centeral of the Gaussian beam with the vaule of k_spp_abs.
 *
 */

k_t random_spp_wavevector(double k_spp_abs,scatterer_t scatt,double distance, gsl_rng *r){
	k_t k;
	double k_abs = 0.0;
	k_abs = k_spp_abs + (sqrt(scatt.x*scatt.x + scatt.y*scatt.y)/distance)*k_spp_abs/sin(32.0*M_PI/180.0);
	double theta = random_double_range(r, 0.0, 2.0*M_PI);
	k.kz = 0.0;
	k.kx = k_abs * sin(theta);
	k.ky = k_abs * cos(theta);
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

double spp_radiation_co(double k_spp_abs, double radiation_angle){
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
double angle_between_vector(k_t k_1, k_t k_2){
	double angle = 0.0;
	double sin_angle = 0.0;
   	sin_angle = (k_1.kx*k_2.kx + k_1.ky*k_2.ky  +k_1.kz*k_2.kz)
   				/(sqrt(k_1.kx*k_1.kx + k_1.ky*k_1.ky  +k_1.kz*k_1.kz )*sqrt(k_2.kx*k_2.kx + k_2.ky*k_2.ky  + k_2.kz*k_2.kz));
   	angle = asin(sin_angle);
   	return angle;
}

/*
 * this is the method to create the index of the next scatter */
int next_scatter_index(int current_index,int NSCAT, gsl_rng *r, index_bounds *matrix){
	int next_scatter_index,i;
	double random_probability = 0.0;
	random_probability = random_double_range(r, 0.0, 1.0);
	for (i = 0; i < NSCAT-1; ++i) {
		if (random_probability<matrix[i].upper_bound && random_probability > matrix[i].lower_bound) {
			next_scatter_index = matrix[i].index;
			break;
		}
	}
	return next_scatter_index;
}

/**
 * This is a method to create a conductance matrix
 *
 */
int conduc_matrix(int NSCAT,scatterer_t *scatts,int scatt_index, index_bounds *matrix){
	double coeff[NSCAT-1];
	int i,j;
	double pathlength, coeff_sum;
	for (i = 0; i < NSCAT; ++i) {
		if (i!=scatt_index) {
			pathlength = distance_two_scatters(scatts[i],scatts[scatt_index]);
			pathlength *= 1.0e+6;
			coeff[i] = 1/(pathlength*pathlength);
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



field_t single_field_spp(int first_scatt_index, scatterer_t *scatts,int NSCAT, gsl_rng *r, double *free_path, int free_path_index,index_bounds *matrix){
	/*Constances definition
	 *
	 * Here we define the spp happens on the gold film surface.
	 * */
	double lambda_light = 632.8e-9; /* wavelength of light */
	double k_spp_abs = 2.00329*2.0*M_PI/lambda_light*sin(32.0*M_PI/180.0);/*spp wave vector on LAH79 with angle of theta_sp
	  = 32.0 degree*/
	/* For the perpendicular component of the wave_vector*/
	/* wavevector on LAH79 with angle of theta_sp = 32.0 degree*/
	double k_per_abs = 2.00329*2.0*M_PI/lambda_light*cos(32.0*M_PI/180.0);
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
	// TODO: what the distance should be take here, test just 600e-6
	k_t k_out = random_spp_wavevector(k_spp_abs,scatt_cur,600e-6, r);
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
	field.k.kx = k_out.kx;
	field.k.ky = k_out.ky;
	field.k.kz = k_per_abs;
	field.scatterer_index = curr_scatterer_index;
	return field;
}
/**
 * This function create the phase shift between the scatterer to a position on the cone.
 * The equation is based on the equation from Bert's thesis page.57.
 */
complex double phase_scatt_to_cone(k_t k_out, double distance, scatterer_t scatt){
	complex double phase = 1.0i;
	double lambda_light = 632.8e-9; /* wavelength of light */
	double k_light_abs = 2.0*M_PI/lambda_light;
	double radi_angle = 0.0;
	double shift = 0.0;
	double tan_theta = tan(32.0*M_PI/180.0);
	k_t k_in;
	k_in.kx = 1.0;
	k_in.ky = 0.0;
	k_in.kz = 0.0;
	radi_angle = angle_between_vector(k_in, k_out);
	shift = sqrt((distance*tan_theta*cos(radi_angle) + scatt.x)*(distance*tan_theta*cos(radi_angle) + scatt.x)
			+ (distance*tan_theta*sin(radi_angle) + scatt.y)*(distance*tan_theta*sin(radi_angle) + scatt.y)
			+ (distance*distance));
	phase = cexp(1.0i*shift*k_light_abs);
	return phase;
}



/**
 * This is the method for the free space transform  for one set of all the scatterers.
 */
scatterer_t *fst_transfer(field_t *field, double z0, scatterer_t *scatts_orig,int NSCAT, scatterer_t *locations){
	int i;
	for (i = 0; i < NSCAT; ++i) {
		locations[i].x = scatts_orig[field[i].scatterer_index].x + field[i].k.kx/field[i].k.kz*z0;
		locations[i].y = scatts_orig[field[i].scatterer_index].y + field[i].k.ky/field[i].k.kz*z0;
		locations[i].z = scatts_orig[field[i].scatterer_index].z + z0;
		// Phase shift
		field[i].field_abs *= phase_scatt_to_cone(field[i].k,z0,scatts_orig[field[i].scatterer_index]);
	}
	return locations;
}

/**
 * This is the method to create the field distribution on the ccd chip from a certain distance for the orignal
 * scattering plane. Here we put the ccd chip parallize the x-y plane and has a distance from the orignal point
 * of z0.
 */
int field_on_ccd(double distance, scatterer_t *locations, field_t *field, camera_t cam,int NSCAT, complex double *ccd){
	  double dx = cam.cam_lx/(cam.cam_sx-1);
	  double dy = cam.cam_ly/(cam.cam_sy-1);
	  //int cam_sx = cam.cam_sx;  /* unused */
	  int cam_sy = cam.cam_sy;
		//check_mem(ccd);
	  int i,pixel_index;
	  for (i = 0; i < NSCAT; ++i) {
		  // TODO: double check
		  pixel_index =   (int)((locations[i].x+cam.cam_lx/2)/dx)*cam_sy 
				            + (int)((locations[i].y+cam.cam_ly/2)/dy);
		  ccd[pixel_index] += field[i].field_abs;
	}
	  return 0;
}


int main(int argc, char **argv) {

	/* program variables */
	  unsigned int i,j,m;                       /* variables to iterate over */
	  int NSCAT = 65;                         /* the scatterer number */
	  /*
	   * To avoid the memory limitation, here we use a LOW_NI and UP_NI to get enough iteration times.
	   * The total iteration times should be LOW_NI*UP_NI.
	   */
	  int LOW_NI = 5000;                       /* lower iteration times for each scatterer*/
	  int UP_NI = 5000;                        /*Upper iteration times for each scatterer*/

	  camera_t cam = {0.20,0.20,512,512};     /* initialize the cam struct*/
	  double z0 = 0.13;                       /* the camera distance from the orginal plane*/

		field_t **field = NULL;
	  field = malloc(LOW_NI*sizeof(field_t*));
		check_mem(field);
		for(i=0;i<LOW_NI;++i){
			field[i] = (field_t*)malloc(NSCAT*sizeof(field_t)); 
			check_mem(field[i]);
		}

	/* seed the scatterers */
	  scatterer_t *scatts = NULL;
	  scatts = malloc(NSCAT*sizeof(scatterer_t));
		check_mem(scatts);
	/* the locations */
	  scatterer_t *locations = NULL;
	  locations = malloc(NSCAT*sizeof(scatterer_t));
		check_mem(locations);
	/*free path array*/
	  double *free_path = NULL;
	  free_path = malloc(LOW_NI*NSCAT*sizeof(double));

	/* the radiation counting for NSCATS */
	  int *visiting_array = NULL;
	  visiting_array = malloc(NSCAT*sizeof(int));
	  bzero(visiting_array,NSCAT*sizeof(int));

	  complex double *ccd_all = NULL;
	  ccd_all = malloc(cam.cam_sx*cam.cam_sy*sizeof(complex double));
	  bzero(ccd_all,cam.cam_sx*cam.cam_sy*sizeof(complex double));
	/* To read the scatterers positions from the file into the scatterer array*/
		char filename[FILENAME_MAX];
		bzero(filename,FILENAME_MAX*sizeof(char));
		sprintf(filename,"%s","Scatterers.txt");
		log_info("Reading from %s.",filename);
	  FILE *fp = fopen(filename,"r");
		check(fp, "Failed to open %s for reading.", filename); 
	  for (i = 0; i < NSCAT; ++i) {
			fscanf(fp,"%lf %lf %lf",&scatts[i].x,&scatts[i].y,&scatts[i].z);
			fscanf(fp,"\n");
			//printf("%12.12f %12.12f %12.12f\n",scatts[i].x,scatts[i].y,scatts[i].z);
	  }
	  fclose(fp);

	  index_bounds **Matrix = NULL;
		Matrix = (index_bounds**)malloc(NSCAT*sizeof(index_bounds*));
		check_mem(Matrix);
		for(i=0;i<NSCAT;++i){ 
			Matrix[i] = (index_bounds*)malloc(NSCAT*sizeof(index_bounds)); 
			check_mem(Matrix[i]);
		  conduc_matrix(NSCAT, scatts, i, Matrix[i]);
		}
		log_info("Conductance matrix finished");

	  //index_bounds *matrix;
	  //for (i = 0; i < NSCAT; ++i) {
		 // matrix = Matrix[i];
		//for (k = 0; k < NSCAT-1; ++k) {
		//	printf("%12.12f   %d     ",Matrix[i][k].upper_bound,Matrix[i][k].index);
		//}
		//printf("%d\n",i);
	  //}
		//

	/* initialize the random number generator */
	  int rseed = (int)time(NULL);
	  const gsl_rng_type *T;
	  gsl_rng *r;
	  gsl_rng_env_setup();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc(T);
	  gsl_rng_set(r,rseed);

	  for (m = 0; m < UP_NI; ++m) {
		  for(i=0;i<LOW_NI;++i){
			  for (j=0;j<NSCAT;++j) {
				  field[i][j] = single_field_spp(j, scatts, NSCAT, r, free_path,(i*NSCAT + j),Matrix[j]);
				  ++visiting_array[field[i][j].scatterer_index];
			  }
		  }

		  for(i=0;i<LOW_NI;++i) {
			  locations = fst_transfer(field[i],z0,scatts,NSCAT, locations);
			  field_on_ccd(z0,locations,field[i],cam,NSCAT,ccd_all);
		  }
	  }
		/**
		 * Write the tmp array as a way to do the cross-correlation function.
		 */
		bzero(filename,FILENAME_MAX*sizeof(char));
		sprintf(filename,"%s","visiting_counter_65_2.txt");
		FILE *fp_counting = fopen(filename,"w");
		check(fp_counting, "Failed to open %s for writing.", filename);

		log_info("Writing to %s.",filename);
		for (i = 0; i < NSCAT; ++i) {
			fprintf(fp_counting,"%d\n",visiting_array[i]);
		}
		fclose(fp_counting);

		double *tmp = NULL;
		tmp = malloc(cam.cam_sx*cam.cam_sy*sizeof(double));
		check_mem(tmp);

	  /* normalize tmp, first by finding max */
		double max = 0;
		for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
			tmp[i]=cabs(ccd_all[i])*cabs(ccd_all[i]);
			if(tmp[i]>max){ max=tmp[i]; }

		}
		/* then dividing by max */
		for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
			tmp[i]/=max;
		}

		/**
		 * Write the tmp array as a way to do the cross-correlation function.
		 */
		bzero(filename,FILENAME_MAX*sizeof(char));
		sprintf(filename,"%s","tmp_500.txt");
		FILE *fp_tmp = fopen(filename,"w");
		check(fp_tmp, "Failed to open %s for writing.", filename); 

		int size_none_zero = 0;

		for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
			if(tmp[i]!=0.0){
				++size_none_zero;
			}
		}

	    double *tmp_N = NULL;
	    tmp_N = malloc(size_none_zero*sizeof(double));
			check_mem(tmp_N);
	    j = 0;
	    for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
	    	if(tmp[i]!=0.0){
	    		tmp_N[j] = tmp[i];
	    		j++;
	    	}
	    }
			log_info("Writing to %s.",filename);
	    for (i = 0; i < size_none_zero; ++i) {
				fprintf(fp_tmp,"%-12.12f\n",tmp_N[i]);
	    }
	    fclose(fp_tmp);


	    /**
	     * Write the free path array as a way to do the cross-correlation function.
	     */
		bzero(filename,FILENAME_MAX*sizeof(char));
		sprintf(filename,"%s","free_path.txt");
	    FILE *fp_path = fopen("free_path.txt","w");
		check(fp_path, "Failed to open %s for writing.", filename);

		log_info("Writing to %s.",filename);
	    for (i = 0; i < LOW_NI*NSCAT; ++i) {
				fprintf(fp_path,"%-12.12f\n",free_path[i]);
	    }
	    fclose(fp_path);


	  /* output file */
		hid_t file,dataset,dataspace;
		herr_t status = 0;
		hsize_t dims[2]; /* dimensionality of the set */
		dims[0]=cam.cam_sx;
		dims[1]=cam.cam_sy;
		dataspace = H5Screate_simple(2,dims,NULL);
		bzero(filename,FILENAME_MAX*sizeof(char));
		sprintf(filename,"%s","out2.h5");
		file = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
		dataset = H5Dcreate1(file,"/e2",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);

		log_info("Writing to HDF5 file %s.",filename);
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
		check(status>=0,"Can't write dataset /e2 to %s",filename);

	  /* write the real and imaginary parts as well */
		dataset = H5Dcreate1(file,"/er",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);
		for(i=0;i<cam.cam_sx*cam.cam_sy;++i){ tmp[i]=creal(ccd_all[i]); }
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
		check(status>=0,"Can't write dataset /er to %s",filename);

		dataset = H5Dcreate1(file,"/ei",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);
		for(i=0;i<cam.cam_sx*cam.cam_sy;++i){ tmp[i]=cimag(ccd_all[i]); }
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
		check(status>=0,"Can't write dataset /ei to %s",filename);

	  /* clean up */
		log_info("Cleaning up.");
		gsl_rng_free(r);
		status = H5Dclose(dataset);
		status = H5Sclose(dataspace);
		status = H5Fclose(file);

		free(tmp);
		free(scatts);
		free(locations);
		free(free_path);
		free(visiting_array);
		free(ccd_all);
		for(i=0;i<NSCAT;++i){ free(Matrix[i]); }
		free(Matrix);
		for(i=0;i<LOW_NI;++i){ free(field[i]); }
		free(field);

		log_info("Program finished.");

	  return 0;
		error:
			return 1;
}


