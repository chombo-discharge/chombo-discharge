#!/bin/bash -l
#SBATCH -q regular
#SBATCH -C knl,quad,cache
#SBATCH -N 1   # number of nodes
#SBATCH -S 4   # Cori KNL reservation of excess cores. we use 64, not 68
#SBATCH -L cscratch1
#SBATCH -J myJobName  # what you want your job called in the queue
#SBATCH -t 00:02:00  #time for reservation HH:MM:SS


export CH_TIMER=1  # environment variables used for Chombo
export CH_OUTPUT_INTERVAL=1
#export OMP_PROC_BIND=spread
export OMP_PROC_BIND=true
export OMP_PLACES=threads


ProgName=/global/cscratch1/sd/bvs/ChomboKNL/test1/test1.exe
input="run1.inputs"
rm -f time.table.* pout.*

# example of doing a hybrid tuning sweep
for nranks in 64 32 16 8 4 2; do
    for nthreads in 1 2 4 8 16 32; do
        ranks=$(( $nranks * SLURM_JOB_NUM_NODES))
        p=$(( $nranks * $nthreads ))
        if (( $p == 64 )); then
            ncpus=$(( 4 * $nthreads ))
            export OMP_NUM_THREADS=$nthreads
            echo srun -n $ranks -c $ncpus --cpu_bind=cores $ProgName $input >log
                 srun -n $ranks -c $ncpus --cpu_bind=cores $ProgName $input >log
            res_dir="_run_mpi${ranks}_omp${nthreads}"
            if [ ! -f $res_dir ]; then
                mkdir $res_dir
            fi
            mv time.table* pout* log "$res_dir"
        fi
    done
done
