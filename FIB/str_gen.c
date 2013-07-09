/** 

This a program to generate the stream file for the FIB litho
-Jiapeng Huang, 9.July.2013
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <stdbool.h>
#include <string.h>
#include <malloc.h>
#include <sys/types.h>
#include <float.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h> /* command line arguments, get memory */
#include <getopt.h> /* command line arguments */
#include <assert.h>
#include <complex.h>
#include <gsl/gsl_rng.h> /* random number generator using mt19937 */
#include <gsl/gsl_randist.h>
#include <errno.h>
#include <error.h>

int num_scatterers = 100;
double scat_radius = 0.25; /*in micron*/
double x_area = 100;
double y_area = 100;

/* the X has a range from 0-4095*/
int x_index_low = 500;
int x_index_up = 3500;
/* the Y has a range from 280-3816*/

int y_index_low = 500;
int y_index_up = 3500;

/*dwell time */
int dwell_t = 96; /* dewll_t * 0.1 mircsec
/*repeating time*/
int rep = 20;

/* return a random int */
int random_int(gsl_rng *r,int max){
  return gsl_rng_uniform_int(r,max+1);
}

/* a scatterer */
typedef struct {
  int x,y; /* scatterer index position */
} scatterer_t;

/*Bresenham Algorithm for a full circle*/
void rasterCircle(int x0, int y0, int radius)
{
  int f = 1 - radius;
  int ddF_x = 1;
  int ddF_y = -2 * radius;
  int x = 0;
  int y = radius;
 
  setPixel(x0, y0 + radius);
  setPixel(x0, y0 - radius);
  setPixel(x0 + radius, y0);
  setPixel(x0 - radius, y0);
 
  while(x < y)
  {
    // ddF_x == 2 * x + 1;
    // ddF_y == -2 * y;
    // f == x*x + y*y - radius*radius + 2*x - y + 1;
    if(f >= 0) 
    {
      y--;
      ddF_y += 2;
      f += ddF_y;
    }
    x++;
    ddF_x += 2;
    f += ddF_x;    
    setPixel(x0 + x, y0 + y);
    setPixel(x0 - x, y0 + y);
    setPixel(x0 + x, y0 - y);
    setPixel(x0 - x, y0 - y);
    setPixel(x0 + y, y0 + x);
    setPixel(x0 - y, y0 + x);
    setPixel(x0 + y, y0 - x);
    setPixel(x0 - y, y0 - x);
  }

}

void setPixel(int x,int y){
        FILE *fp_1 = fopen("stream.s","a");
        if (fp_1 == 0) {
                fprintf(stderr, "failed to open the file");
                exit(1);
        }
        fprintf(fp_1,"%d %d %d\n",dwell_t,x,y);
    	fclose(fp_1);
}

int main(int argc, char **argv) {
	
	/*index range for a circle radi*/
	int scat_radi_index;
	scat_radi_index = scat_radius/x_area*(x_index_up-x_index_low);
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
  	int scanx = x_index_up-x_index_low; /* the boundary of the region */
  	int scany = y_index_up-y_index_low; /* the boundary of the region */

	/* seed the scatterers */
  	scatterer_t *scatt = malloc(num_scatterers*sizeof(scatterer_t));
  	assert(scatt!=NULL);

  	for(i=0;i<num_scatterers;++i){
    		scatt[i].x=random_int(r,scanx);
    		scatt[i].y=random_int(r,scany);
  	}

	FILE *fp = fopen("stream.s","w");
       	if (fp == 0) {
        	fprintf(stderr, "failed to open the file");
        	exit(1);
  	}
	
	fprintf(fp,"s\n%d\n",rep);	
	fclose(fp);
	int radi_index ;
	for(i=0;i<num_scatterers;++i){
		radi_index = scat_radi_index;
		while(radi_index >= 0){
                rasterCircle(scatt[i].x,scatt[i].y,radi_index);
                radi_index = radi_index -1;
		}
	}
	return 0;
}


