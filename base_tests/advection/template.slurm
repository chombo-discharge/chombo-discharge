#!/bin/bash
#SBATCH --account=nn9453k
#SBATCH --job-name=advection_gauss
#SBATCH --time=0-00:30:00
#SBATCH --qos=devel
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32

## Safety settings, I don't know why this is recommended in the documentation but fine
set -o errexit
set -o nounset

## Software modules
module restore system
module load intel/2017a
module load HDF5/1.10.1-intel-2017a

# Program to run and input script
executable=main2d.Linux.64.mpiicpc.ifort.OPT.MPI.ex
input_file=gauss_wave.inputs

# Create directory, copy executable and input script
workdir=$USERWORK/$SLURM_JOB_ACCOUNT/$SLURM_JOB_NAME
if [ ! -d $workdir ]; then mkdir -p $workdir; fi
cp -f $executable $workdir/$executable
cp -f inputs/$input_file $workdir/$input_file

# Navigate to directory and run
cd $workdir
mpirun ./$executable $input_file
