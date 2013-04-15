/**
 * This a code to produce a random scatterers array in the space and store the 3-D data in to the file.
 *
 * @author - Jiapeng Huang - MPL - 2013.03.08
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

/* a scatterer */
typedef struct {
  double x,y,z; /* scatterer position */
} scatterer_t;

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
  unsigned int i; /* variables to iterate over */
  double scanxy = 50.0e-6; /* the boundary of the region */
  int NSCAT = 300; /*the scatterer number*/

  /* seed the scatterers */
  scatterer_t *scatt = malloc(NSCAT*sizeof(scatterer_t));
  assert(scatt!=NULL);

  for(i=0;i<NSCAT;++i){
    scatt[i].x=random_double_range(r,-scanxy/2.0,scanxy/2.0);
    scatt[i].y=random_double_range(r,-scanxy/2.0,scanxy/2.0);
    scatt[i].z=random_double_range(r,0,1000e-9);
  }

/*
 * Creat a txt file to store the scatterers position in a format looks like:
 *
 * x1 y2 z3
 * x2 y2 z3
 *    .
 *    .
 *    .
 * xn yn zn
 */

  FILE *fp = fopen("Scatterers.txt","w");

  if (fp == 0) {
	fprintf(stderr, "failed to open the file");
	exit(1);
  }

  for (i = 0; i < NSCAT; ++i) {
	fprintf(fp,"%-12.12f %-12.12f %-12.12f\n",scatt[i].x,scatt[i].y,scatt[i].z);
  }
  fclose(fp);

  FILE *fp_x = fopen("Scatterers_x.txt","w");

  if (fp_x == 0) {
	fprintf(stderr, "failed to open the file");
	exit(1);
  }

  for (i = 0; i < NSCAT; ++i) {
	fprintf(fp_x,"%-12.12f\n",scatt[i].x);
  }
  fclose(fp_x);

  FILE *fp_y = fopen("Scatterers_y.txt","w");

  if (fp_y == 0) {
	fprintf(stderr, "failed to open the file");
	exit(1);
  }

  for (i = 0; i < NSCAT; ++i) {
	fprintf(fp_y,"%-12.12f\n",scatt[i].y);
  }
  fclose(fp_y);

  return 0;
}
