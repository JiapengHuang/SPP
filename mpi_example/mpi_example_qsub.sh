#!/bin/bash
#PBS -N mpi_example
#PBS -o mpi_example.log
#PBS -l walltime=00:05:00
#PBS -l nodes=1:ppn=4

# go to proper location
cd $PBS_O_WORKDIR

# command to be executed
time mpirun mpi_example
