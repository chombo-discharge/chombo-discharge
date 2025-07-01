#!/bin/bash
#SBATCH --account=nn12041k
#SBATCH --job-name=cylinder
#SBATCH --time=0-00:10:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
##SBATCH --cpus-per-task=64
##SBATCH --qos=devel

## Safety settings, I don't know why this is recommended in the documentation but fine
set -o errexit
set -o nounset

## Software modules
module restore system
module load foss/2023a
module load HDF5/1.14.0-gompi-2023a

# Program to run and input script
executable=program3d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex

# Input file
input_file=example.inputs
chemistry_file=no_energy_sharing.json
transport_file=transport_data.dat

# Create directory, copy executable and input script
workdir=$USERWORK/$SLURM_JOB_ACCOUNT/DBD/$SLURM_JOB_NAME
if [ ! -d $workdir ]; then mkdir -p $workdir; fi

# Only pout.0
export CH_OUTPUT_INTERVAL=9999999

# Lustre striping
lfs setstripe --stripe-count 4 --stripe-size 32M $workdir

# Copy executables to work directory
cp -f $executable $workdir/$executable
cp -f $input_file $workdir/$input_file
cp -f $chemistry_file  $workdir/$chemistry_file
cp -f $transport_file  $workdir/$transport_file

# Navigate to directory and run
cd $workdir

mpirun ./$executable $input_file ItoKMCJSON.chemistry_file=$chemistry_file Driver.max_steps=100
