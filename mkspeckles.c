/* 
 * mkspeckles.c
 * This program makes specklegrams using single scattering statistics and
 * writes them to HDF5 files
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
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
__inline double random_double_range(gsl_rng *r,double min, double max) {
	 return (max-min)*(gsl_rng_uniform(r))+min;
}

/* return a random int */
__inline int random_int(gsl_rng *r,int max){
	return gsl_rng_uniform_int(r,max+1);
}

/* a scatterer */
typedef struct {
	double x,y,z; /* scatterer position */
} scatterer_t;


int main(int argc, char **argv) {
	/* initialize the random number generator */
	int rseed = (int)time(NULL);
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r,rseed);

	/* program variables */
	unsigned int i,j,k,l; /* variables to iterate over */
	double lambda = 632.8e-9; /* wavelength of light */ 
	double scanxy = 15e-6; /* the boundary of the region */
	unsigned int NSCAT = 50;
	unsigned int EVENTS = 10001;
	double ELLIPSE_A    = 5.0e-6; /* size of the illuminated area as an ellipse */
	double ELLIPSE_B    = 3.0e-6;
	double cam_lx = 0.20; /* size of camera ccd x */
	double cam_ly = 0.20; /* size of camera ccd y */
	unsigned int cam_sx = 512; /* camera pixels x */
	unsigned int cam_sy = 512; /* camera pixels y */

	/* seed the scatterers */
	scatterer_t *scatt = malloc(NSCAT*sizeof(scatterer_t));
	assert(scatt!=NULL);

	for(i=0;i<NSCAT;++i){
		scatt[i].x=random_double_range(r,-scanxy/2.0,scanxy/2.0);
		scatt[i].y=random_double_range(r,-scanxy/2.0,scanxy/2.0);
		scatt[i].z=random_double_range(r,0,1000e-9);
		printf("%g\t%g\n",scatt[i].x,scatt[i].y);
	}
	
	/* compute the path lengths */
	double **pathlengths = malloc(NSCAT*sizeof(double*));
	for(i=0;i<NSCAT;++i){ pathlengths[i]= malloc(NSCAT*sizeof(double)); }

	/* fill the pathlengths with the proper pathlengths */
	for(i=0;i<NSCAT;++i){
		for(j=0;j<NSCAT;++j){
			pathlengths[i][j]=sqrt( 
		 				(scatt[i].x-scatt[j].x)*(scatt[i].x-scatt[j].x)
					+ (scatt[i].y-scatt[j].y)*(scatt[i].y-scatt[j].y));
	}}


	/* this holds the amount of scattering canidates */
	int *iscanidates = calloc(NSCAT,sizeof(int)); 
	assert(iscanidates!=NULL);

	/* far field camera */
	complex double *ccd = malloc(cam_sx*cam_sy*sizeof(complex double));
	bzero(ccd,cam_sx*cam_sy*sizeof(complex double));
	assert(ccd!=NULL);

	double z = 0.02; /* distance to ccd */
	double x,y;
	double dx = cam_lx/(cam_sx-1);
	double dy = cam_ly/(cam_sy-1);
	for(l=0,x=-cam_lx/2,i=0;i<cam_sx;++i,x+=dx){
	for(y=-cam_lx/2,j=0;j<cam_sy;++j,y+=dy,++l){
		for(k=0;k<NSCAT;++k){

			ccd[l]+=cexp(2.0i*M_PI/lambda*
								sqrt( (x-scatt[k].x)*(x-scatt[k].x) 
											+ (y-scatt[k].y)*(y-scatt[k].y)
											+ (z-scatt[k].z)*(z-scatt[k].z)
											));
		}
	}}


	/* output file */
	hid_t file,dataset,dataspace;
	herr_t status;
	hsize_t dims[2]; /* dimensionality of the set */
	dims[0]=cam_sx;
	dims[1]=cam_sy;
	dataspace = H5Screate_simple(2,dims,NULL);
	file = H5Fcreate("out.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
	dataset = H5Dcreate1(file,"/e2",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);

	double *tmp = malloc(cam_sx*cam_sy*sizeof(double));
	assert(tmp!=NULL);

	/* normalize tmp, first by finding max */
	double max = 0;
	for(i=0;i<cam_sx*cam_sy;++i){ 
		tmp[i]=cabs(ccd[i])*cabs(ccd[i]); 
		if(tmp[i]>max){ max=tmp[i]; }
	
	}
	/* then dividing by max */
	for(i=0;i<cam_sx*cam_sy;++i){ 
		tmp[i]/=max;
	}


	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

	/* write the real and imaginary parts as well */
	dataset = H5Dcreate1(file,"/er",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);
	for(i=0;i<cam_sx*cam_sy;++i){ tmp[i]=creal(ccd[i]); }
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

	dataset = H5Dcreate1(file,"/ei",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);
	for(i=0;i<cam_sx*cam_sy;++i){ tmp[i]=cimag(ccd[i]); }
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

	/* clean up */
	status = H5Dclose(dataset);
	status = H5Sclose(dataspace);
	status = H5Fclose(file);

	free(tmp);
	free(ccd);
	free(scatt);
	gsl_rng_free(r);
	return 0;
}
