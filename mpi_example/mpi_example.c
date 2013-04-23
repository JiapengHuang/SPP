/* mpi_example.c 
 * Example of using mpi gather to do monte carlo simulations in parallel */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h> /* distributed computing */ 
#include <malloc.h>
#include <string.h>
#include <unistd.h> /* parse command line arguments */
#include <getopt.h> /* command line arguments */
#include <sys/stat.h>
#include <ctype.h>
#include <time.h>
#include <float.h>
#include <gsl/gsl_rng.h> /* random number generator using mt19937 */
#include <errno.h>
#include <search.h>



int main( int argc, char *argv[]) {
	unsigned int x_i,y_i,o = 0;

	/* MPI initialization */
	/* The variables rank and nprocs are set by the MPI process.  rank is
	 * the number of the individual process, nprocs is the total number of
	 * processors */
	int rank, nprocs; 
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	/* each process should get a different seed */
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	int rseed = time(0);
	srand(rseed+rank);
	gsl_rng_set(r,rseed);

	unsigned int SX=10;
	unsigned int SY=10;
	double *buf = calloc(SX*SY,sizeof(double));			 /* a buf for each process */
	double *buf_root = calloc(SX*SY,sizeof(double)); /* what will be the sum */

	/* raster */
	fprintf(stdout,"Rank %i of %i started...\n",rank,nprocs);
	for(o=0,x_i=0;x_i<SX;++x_i){
		for(y_i=0;y_i<SY;++y_i,++o){
			buf[o]+=1;
		}
	}

	/* combine all the event runs and add them up */
	MPI_Reduce(buf,buf_root,SX*SY,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	/* only print out the rank 0 "buf_root" array */
	if(rank==0){
		for(o=0,x_i=0;x_i<SX;++x_i){
			for(y_i=0;y_i<SY;++y_i,++o){
				printf("%g\t",buf_root[o]);
			}
			printf("\n");
		}
	}



	free(buf); free(buf_root);
	gsl_rng_free(r);
	MPI_Finalize();

	return 0;
}
