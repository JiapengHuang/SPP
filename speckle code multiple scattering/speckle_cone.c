/*
 * speckles_cone.c
 * This is program to run the simulation of the cone speckle structure from Bert's thesis.  * As the result of this program a HDF5 file will be write to. 
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
__inline double random_double_range(gsl_rng *r,double min, double max)
{
	return ((max-min)*gsl_rng_uniform(r)) + min;
}

/* return a random int */
__inline int random_int(gsl_rng *r,int max){
  return gsl_rng_uniform_int(r,max+1);
}

/* a scatterer */
typedef struct {
  double x,y,z; /* scatterer position */
} scatterer_t;

/* a camera */

typedef struct{
  double cam_lx; /* size of camera ccd x*/
  double cam_ly; /* size of camera ccd y*/
  unsigned int cam_sx; /*camera pixels x*/
  unsigned int cam_sy; /*camera pixels y*/
} camera_t;

/* calculate the scattering field for a single scatterer scattering*/ 
__inline complex double * field_single_scatterer(scatterer_t scatterer,camera_t cam,complex double *ccd,complex double phase_shift){
	assert(ccd!=NULL);

	double lambda_light = 632.8e-9;
	unsigned int i,j,k,l; /* variables to iterate over */	
	double z = 0.02; /* distance to ccd */
  	double x,y;
  	double dx =cam.cam_lx/(cam.cam_sx-1);
  	double dy = cam.cam_ly/(cam.cam_sy-1);
  	for(l=0,x=-cam.cam_lx/2,i=0;i<cam.cam_sx;++i,x+=dx){
  	for(y=-cam.cam_lx/2,j=0;j<cam.cam_sy;++j,y+=dy,++l){
      	ccd[l]+=cexp(2.0i*M_PI/lambda_light*
                sqrt( (x-scatterer.x)*(x-scatterer.x)
                      + (y-scatterer.y)*(y-scatterer.y)
                      + (z-scatterer.z)*(z-scatterer.z)
                      ))*phase_shift;
  }}
	return ccd;
}

/* find the next scatter*/
__inline scatterer_t next_scatter(scatterer_t current_scatter, int scatterers_number,scatterer_t *scatterers,gsl_rng *r){
    scatterer_t next_scatterer = current_scatter;
	int next_scatter_index = random_int(r,scatterers_number);
	while(scatterers[next_scatter_index].x == current_scatter.x && scatterers[next_scatter_index].y ==current_scatter.y)
	{
		next_scatter_index = random_int(r,scatterers_number);
	}
	//printf("%i\n",next_scatter_index);
	next_scatterer =scatterers[next_scatter_index];
	return next_scatterer;
}

/* Calculate the distance between two scatters*/
__inline double distance_two_scatters(scatterer_t current_scatterer,scatterer_t next_scatterer){
	double distance = 0.0;
	distance = sqrt((current_scatterer.x - next_scatterer.x)*(current_scatterer.x - next_scatterer.x) + ((current_scatterer.y - next_scatterer.y)*(current_scatterer.y - next_scatterer.y)));
	//printf("%g\n",distance);
	return distance;
}

/* Calculate the scattering field after multiple scattering process. Here we have the assumption that the scattering process stop when the path length is larger than the max path langth and when the scatterer comes back to the orignal one.*/

__inline complex double * field_multi_scattering(scatterer_t original_scatterer,int scatterers_number,scatterer_t *scatterers,camera_t cam,double maxPathLength,double lambda_spp,complex double * ccd_field,gsl_rng *r)
{
	double pathLength = 0.0;
	scatterer_t current_scatterer = original_scatterer;
	scatterer_t nextScatterer = next_scatter(current_scatterer,scatterers_number, scatterers,r);
	while(pathLength <= maxPathLength && (nextScatterer.x!= original_scatterer.x && nextScatterer.y!= original_scatterer.y))
	{
		pathLength = pathLength + distance_two_scatters(current_scatterer,nextScatterer);
		current_scatterer = nextScatterer;
		nextScatterer = next_scatter(nextScatterer,scatterers_number,scatterers,r);
	}
	//printf("%g\n",pathLength);
	//printf("%s\n","able to come to A");
	complex double phase_shift = cexp(2.0i*M_PI/lambda_spp*pathLength);
	field_single_scatterer(nextScatterer,cam,ccd_field,phase_shift);
	return ccd_field;
}

int main(int argc, char **argv){
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
  double lambda_spp = 610.0e-9;
  double scanxy = 15e-6; /* the boundary of the region */
  double MAX_PATH_LENGTH = 18.0e-6;
  unsigned int NSCAT = 50;
  unsigned int EVENTS = 10001;
  double ELLIPSE_A    = 5.0e-6; /* size of the illuminated area as an ellipse */
  double ELLIPSE_B    = 3.0e-6;
  camera_t cam = {0.20,0.20,512,512}; /*initialize the cam struct*/

  /* seed the scatterers */
  scatterer_t *scatt = malloc(NSCAT*sizeof(scatterer_t));
  assert(scatt!=NULL);

  for(i=0;i<NSCAT;++i){
    scatt[i].x=random_double_range(r,-scanxy/2.0,scanxy/2.0);
    scatt[i].y=random_double_range(r,-scanxy/2.0,scanxy/2.0);	   
    scatt[i].z=0.0;//random_double_range(r,0,1000e-9);
  }
    
  scatterer_t current_scatterer = scatt[0];
  complex double *ccd = malloc(cam.cam_sx*cam.cam_sy*sizeof(complex double));
  bzero(ccd,cam.cam_sx*cam.cam_sy*sizeof(complex double));
  assert(ccd!=NULL);	
  
  /*add all the multiple scattering from individual scatterers*/
  for(k=0;k<NSCAT;++k){
    field_multi_scattering(current_scatterer,NSCAT,scatt,cam,MAX_PATH_LENGTH,lambda_spp,ccd,r);
	current_scatterer = scatt[k+1];
  }

  /* output file */
  hid_t file,dataset,dataspace;
  herr_t status;
  hsize_t dims[2]; /* dimensionality of the set */
  dims[0]=cam.cam_sx;
  dims[1]=cam.cam_sy;
  dataspace = H5Screate_simple(2,dims,NULL);
  file = H5Fcreate("out.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  dataset = H5Dcreate1(file,"/e2",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);

  double *tmp = malloc(cam.cam_sx*cam.cam_sy*sizeof(double));
  assert(tmp!=NULL);

/* normalize tmp, first by finding max */
  double max = 0;
  for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
    tmp[i]=cabs(ccd[i])*cabs(ccd[i]);
    if(tmp[i]>max){ max=tmp[i]; }

  }
//  printf("%g\n",max);
/* then dividing by max */
  for(i=0;i<cam.cam_sx*cam.cam_sy;++i){
	  tmp[i]/=max;
  }

  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

/* write the real and imaginary parts as well */
  dataset = H5Dcreate1(file,"/er",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);
  for(i=0;i<cam.cam_sx*cam.cam_sy;++i){ tmp[i]=creal(ccd[i]); }
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

  dataset = H5Dcreate1(file,"/ei",H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT);
  for(i=0;i<cam.cam_sx*cam.cam_sy;++i){ tmp[i]=cimag(ccd[i]); }
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

/* clean up */
  status = H5Dclose(dataset);
  status = H5Sclose(dataspace);
  status = H5Fclose(file);

  free(tmp);
  free(ccd);
  free(scatt);
  return 0;
  
}

