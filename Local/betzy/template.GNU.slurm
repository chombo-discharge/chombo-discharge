#!/bin/bash
#SBATCH --account=nn9887k
#SBATCH --job-name=job_template
#SBATCH --time=0-00:10:00
#SBATCH --nodes=4
#SBATCH --qos=devel
#SBATCH --ntasks-per-node=128

## Safety settings, I don't know why this is recommended in the documentation but fine
set -o errexit
set -o nounset

## Software modules
module restore system
module load foss/2020a
module load HDF5/1.10.6-gompi-2020a

# Program to run and input script
executable=main3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex

# Input file
mesh_file=armadillo.ply
input_file=armadillo3d.inputs

# Create directory, copy executable and input script
workdir=$USERWORK/$SLURM_JOB_NAME
if [ ! -d $workdir ]; then mkdir -p $workdir; fi

# Only pout.0
export CH_OUTPUT_INTERVAL=9999999

# Copy executables to work directory
cp -f $executable $workdir/$executable
cp -f $mesh_file $workdir/$mesh_file
cp -f $input_file $workdir/$input_file

# Lustre striping
lfs setstripe --stripe-count 8 --stripe-size 32M $workdir

# Navigate to directory and run
cd $workdir
mpirun ./$executable $input_file # AmrMesh.max_amr_depth=2 PhaseRealm.verbosity=10
