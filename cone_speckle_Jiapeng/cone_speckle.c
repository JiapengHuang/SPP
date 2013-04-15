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

int NSCAT = 300;                         /* the scatterer number */
double z0 = 0.13;                       /* the camera distance from the orginal plane*/
double waist = 10.0e-6;                   /*minimun 5.0e-6*/
double lambda_light = 632.8e-9; /* wavelength of light */
double scanx_edge = -(50.0e-6)/2;
double scany_edge = -(50.0e-6)/2;

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

/* Calculate the distance between two scatters*/
__inline double distance_two_scatters(scatterer_t scatt_a,scatterer_t scatt_b){
	double distance = 0.0;
	distance = sqrt((scatt_a.x - scatt_b.x)*(scatt_a.x - scatt_b.x) + ((scatt_a.y - scatt_b.y)*(scatt_a.y - scatt_b.y)) + ((scatt_a.z - scatt_b.z)*(scatt_a.z - scatt_b.z)));
	return distance;
}

/**
 * Gaussian Beam creator
 */
int gaussian_profile_creator(scatterer_t* scatts,double waist,unsigned long num_iter,int NSCAT,unsigned long *gaussian_matrix){
	double coeff[NSCAT];
	int i;
	double dis, coeff_sum;
	long error = 0;
	for (i = 0; i < NSCAT; ++i) {
			dis = sqrt((scatts[i].x*scatts[i].x)+(scatts[i].y*scatts[i].y));
			coeff[i] = exp(-2.0*dis*dis/(waist*waist));
			coeff_sum += coeff[i];
	}
	for (i = 0; i < NSCAT; ++i) {
			coeff[i] = coeff[i]/coeff_sum;
	}

	for (i = 0; i < NSCAT; ++i) {
		gaussian_matrix[i] = coeff[i]*num_iter;
		error += gaussian_matrix[i];
	}
	error -= num_iter;
	gaussian_matrix[0] -= error;

	error = 0;
	for (i = 0; i < NSCAT; ++i) {
		error += gaussian_matrix[i];
	}
	if (error == num_iter) {
		log_info("Gaussian profile created successfully.");
	} else {
		log_err("Gaussian profile created failure!");
	}
	return 0;
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

k_t random_spp_wavevector(double k_spp_abs, double k_light, double* gaussian_k, gsl_rng *r){
	k_t k;
	double k_abs = 0.0;
	int randn = random_int(r,9999);
	k_abs = gaussian_k[randn];
	double theta = random_double_range(r, 0.0, 2.0*M_PI);
	k.kz = sqrt(k_light*k_light - k_abs*k_abs);
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
	next_scatter_index = 0;
	double random_probability = 0.0;
	random_probability = random_double_range(r, 0.0, 1.0);
	for (i = 0; i < NSCAT-1; ++i) {
		if (random_probability < matrix[i].upper_bound && random_probability >= matrix[i].lower_bound) {
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

/**
 * This is a function to create the phase different introduce by the different positions of the scatterers.
 * This is very important to observe the Coherent Back Scattering.
 * The light is coming in in xoz plane, so the phase shift is accord to the scatterers' positions from the right
 * adge of the simulation area
 */
double phase_by_scatterer(scatterer_t scatt){
	double phase = 0.0;
	phase = 2.0*M_PI*scatt.x/lambda_light*sin(32.0*M_PI/180.0);
	return phase;
}

field_t single_field_spp(int first_scatt_index, scatterer_t *scatts,int NSCAT, gsl_rng *r,index_bounds *matrix, double*gaussian_k){
	/* Constances definition
	 *
	 * Here we define the spp happens on the gold film surface.
	 * */
	double k_spp_abs = 2.00329*2.0*M_PI/lambda_light*sin(32.0*M_PI/180.0);/*spp wave vector on LAH79 with angle of theta_sp
	  = 32.0 degree*/
	/* For wave_vector of light*/
	/* wavevector on LAH79 with angle of theta_sp = 32.0 degree*/
	double k_light = 2.00329*2.0*M_PI/lambda_light;
	double mean_free_path = 18.0e-6;/*mean free path of spp*/
	int next_scatterer_index = next_scatter_index(first_scatt_index,NSCAT,r,matrix);
	int curr_scatterer_index = first_scatt_index;
	field_t field;
	scatterer_t scatt_cur = scatts[first_scatt_index];
	scatterer_t scatt_next = scatts[next_scatterer_index];
	double pathlength = 0.0;
	// TODO: what the distance should be take here, test just 600e-6
	k_t k_out = random_spp_wavevector(k_spp_abs,k_light,gaussian_k, r);
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
	}
	// TODO: we first get the coefficient out
	// coefficient = spp_radiation_co(k_spp_abs,radi_angle);
	// field.field_abs = sqrt(coefficient)*cexp(1.0i*k_spp_abs*pathlength);
	field.field_abs = cexp(1.0i*k_spp_abs*pathlength)*cexp(1.0i*phase_by_scatterer(scatts[first_scatt_index]));
	field.k.kx = k_out.kx;
	field.k.ky = k_out.ky;
	field.k.kz = k_out.kz;
	field.scatterer_index = curr_scatterer_index;
	return field;
}

/**
 * This is the method for the free space transform  for one set of one scatterers.
 */
scatterer_t fst_transfer(field_t field, double z0, scatterer_t scatt_orig){
	scatterer_t location;
	// TODO: should I just get ride of the scatterers original position here.
	location.x = field.k.kx/field.k.kz*z0;
	location.y = field.k.ky/field.k.kz*z0;
	location.z = z0;
	return location;
}

/**
 * This is the method to create the field distribution on the ccd chip from a certain distance for the orignal
 * scattering plane. Here we put the ccd chip parallize the x-y plane and has a distance from the orignal point
 * of z0.
 */
int field_on_ccd(double distance, scatterer_t location,scatterer_t scatt, field_t field, camera_t cam, complex double *ccd){
	  double dx = cam.cam_lx/(cam.cam_sx-1);
	  double dy = cam.cam_ly/(cam.cam_sy-1);
	  double phase = 0.0;
	  //int cam_sx = cam.cam_sx;  /* unused */
	  int cam_sy = cam.cam_sy;
		//check_mem(ccd);
	  int pixel_index;
	  pixel_index =   (int)((location.x+cam.cam_lx/2)/dx)*cam_sy
				            + (int)((location.y+cam.cam_ly/2)/dy);

	  //phase caursed by the out light
	  phase = (field.k.kx*scatt.x + field.k.ky*scatt.y + field.k.kz*scatt.z)
			  /sqrt(field.k.kx*field.k.kx + field.k.ky*field.k.ky + field.k.kz*field.k.kz)*(2.0*M_PI/lambda_light);
	  field.field_abs *= cexp(1.0i*phase);
	  ccd[pixel_index] += field.field_abs;
	  return 0;
}


int main(int argc, char **argv) {

	/* program variables */
	  unsigned int i,j,m,n;                       /* variables to iterate over */
	  /*
	   * To avoid the memory limitation, here we use a LOW_NI and UP_NI to get enough iteration times.
	   * The total iteration times should be LOW_NI*UP_NI.
	   */
	  unsigned long LOW_NI = 30000;                       /* lower iteration times for each scatterer*/
	  unsigned long UP_NI = 50000;                        /*Upper iteration times for each scatterer*/ 
	  camera_t cam = {0.20,0.20,512,512};     /* initialize the cam struct*/

	  field_t field;

	/* seed the scatterers */
	  scatterer_t *scatts = NULL;
	  scatts = malloc(NSCAT*sizeof(scatterer_t));
	  check_mem(scatts);
	/* the locations */
	  scatterer_t location;
	/*free path array*/
	 // double *free_path = NULL;
	//  free_path = malloc(LOW_NI*UP_NI*sizeof(double));

	/* the radiation counting for NSCATS */
	  int *visiting_array = NULL;
	  visiting_array = malloc(NSCAT*sizeof(int));
	  bzero(visiting_array,NSCAT*sizeof(int));

	/*Gaussian profile in k space read from file */
	  double *gaussian_k = NULL;
	  gaussian_k = malloc(10000*sizeof(double));
	  bzero(gaussian_k,10000*sizeof(double));

	  complex double *ccd_all = NULL;
	  ccd_all = malloc(cam.cam_sx*cam.cam_sy*sizeof(complex double));
	  bzero(ccd_all,cam.cam_sx*cam.cam_sy*sizeof(complex double));

	  unsigned long *gaussian_beam = malloc(NSCAT*sizeof(unsigned long));
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
	  }
	  fclose(fp);
	/* To read the gaussian profile for k-vector*/
	  char filename2[FILENAME_MAX];
	  bzero(filename2,FILENAME_MAX*sizeof(char));
	  sprintf(filename2,"%s","Gaussian_k.txt");
	  log_info("Reading from %s.",filename2);
	  FILE *fp2 = fopen(filename2,"r");
	  check(fp2, "Failed to open %s for reading.", filename2);
	  for (i = 0; i < 10000; ++i) {
		  fscanf(fp2,"%lf",&gaussian_k[i]);
		  fscanf(fp2,"\n");
	  }
	  fclose(fp2);


	  index_bounds **Matrix = NULL;
	  Matrix = (index_bounds**)malloc(NSCAT*sizeof(index_bounds*));
	  check_mem(Matrix);
	  for(i=0;i<NSCAT;++i){
		  Matrix[i] = (index_bounds*)malloc(NSCAT*sizeof(index_bounds));
		  check_mem(Matrix[i]);
		  conduc_matrix(NSCAT, scatts, i, Matrix[i]);
		}

	  //for (i = 0; i < NSCAT; ++i) {
		 // for (m = 0; m < NSCAT-1; ++m) {
			//  printf("%12.12f %d ",Matrix[i][m].upper_bound,Matrix[i][m].index);
		//  }
		//  printf("%d\n",i);
	 // }

	  log_info("Conductance matrix finished");

	/* initialize the random number generator */
	  int rseed = (int)time(NULL);
	  const gsl_rng_type *T;
	  gsl_rng *r;
	  gsl_rng_env_setup();
	  T = gsl_rng_default;
	  r = gsl_rng_alloc(T);
	  gsl_rng_set(r,rseed);

	  int gaussian_sum = 0;
	  j = 0;
	  n = 1;
	  gaussian_profile_creator(scatts,waist, UP_NI*LOW_NI,NSCAT,gaussian_beam);
	  for (m = 0; m < UP_NI; ++m) {
		  for (i = 0; i < LOW_NI; ++i) {
			  if (gaussian_sum < gaussian_beam[j]) {
				++gaussian_sum;
				field = single_field_spp(j, scatts, NSCAT,r,Matrix[j],gaussian_k);
				++visiting_array[field.scatterer_index];
				location = fst_transfer(field,z0,scatts[field.scatterer_index]);
				field_on_ccd(z0,location,scatts[field.scatterer_index],field,cam,ccd_all);
			} else {
				++j;
				gaussian_sum = 0;
				if (j == (int)(0.1*n*NSCAT)) {
					log_info("%d0 percent has been done.",n);
					++n;
				}
			}
		  }
	  }
	  log_info("Field generated.");

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
	  //bzero(filename,FILENAME_MAX*sizeof(char));
	  //sprintf(filename,"%s","free_path.txt");
	  //FILE *fp_path = fopen("free_path.txt","w");
	  //check(fp_path, "Failed to open %s for writing.", filename);

	  //log_info("Writing to %s.",filename);
	  //for (i = 0; i < LOW_NI*NSCAT; ++i) {
		  //fprintf(fp_path,"%-12.12f\n",free_path[i]);
	  //}
	  //fclose(fp_path);

	  /* output file */
	  hid_t file,dataset,dataspace;
	  herr_t status = 0;
	  hsize_t dims[2]; /* dimensionality of the set */
	  dims[0]=cam.cam_sx;
	  dims[1]=cam.cam_sy;
	  dataspace = H5Screate_simple(2,dims,NULL);
	  bzero(filename,FILENAME_MAX*sizeof(char));

	  sprintf(filename,"%s","out_s300_w10_50000_30000.h5");
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
	  for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
		  tmp[i]=cimag(ccd_all[i]);
	  }
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
	  //free(free_path);
	  free(visiting_array);
	  free(ccd_all);
	  for(i=0;i<NSCAT;++i){
		  free(Matrix[i]);
	  }
	  free(Matrix);
	  free(gaussian_k);
	  log_info("Program finished.");

	  return 0;
	  error:
	  return 1;
}
