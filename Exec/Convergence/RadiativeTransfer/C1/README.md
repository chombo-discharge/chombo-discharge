## Convergence/RadiativeTransfer/C1

This example computes temporal convergence rates for Physics/RadiativeTransfer using the Eddington SP1 approximation in a coaxial cable geometry. 
The solution errors are computed by coarsening a solution with finer (2x) grid resolution. 

# Compilation

To compile:

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 program```

# Running the example

To run with MPI:

```mpirun -np <num_proc> program2d.*ex convergence2d.inputs```

# Output

The output shows the Linf/L1/L2 errors. 
To get the output in a plot-friendly form, re-route the standard output to a file, e.g.

```mpirun -np <num_proc> program2d.*ex convergence2d.inputs >& convergence.dat```