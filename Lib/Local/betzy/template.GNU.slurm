#!/bin/bash
#SBATCH --account=nn9887k
#SBATCH --job-name=job_template
#SBATCH --time=0-00:10:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --qos=devel

## Safety settings, I don't know why this is recommended in the documentation but fine
set -o errexit
set -o nounset

## Software modules
module restore system
module load foss/2023a
module load HDF5/1.14.0-foss-2023a

# Program to run and input script
executable=main3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex

# Input file
input_file=armadillo3d.inputs
mesh_file=armadillo.ply

# Create directory, copy executable and input script
workdir=$USERWORK/$SLURM_JOB_NAME
if [ ! -d $workdir ]; then mkdir -p $workdir; fi

# Only pout.0
export CH_OUTPUT_INTERVAL=9999999

# Lustre striping
lfs setstripe --stripe-count 4 --stripe-size 32M $workdir

# Copy executables to work directory
cp -f $executable $workdir/$executable
cp -f $input_file $workdir/$input_file
cp -f $mesh_file  $workdir/$mesh_file

# Navigate to directory and run
cd $workdir
mpirun ./$executable $input_file
