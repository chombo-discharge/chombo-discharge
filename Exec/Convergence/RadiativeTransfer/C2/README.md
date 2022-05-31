## Convergence/RadiativeTransfer/C2.

This example computes temporal convergence rates for Physics/RadiativeTransfer using the Eddington SP1 approximation in a coaxial cable geometry. 
The solution errors are computed by subtracting a solution with a finer time step. 

# Compilation

To compile:

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 program```

# Running the example

To run with MPI:

```mpirun -np <num_proc> program2d.*ex convergence2d.inputs```

To use first/second order integration in time, set ```EddingtonSP1.use_tga``` to false/true. 

# Output

The output shows the Linf/L1/L2 errors. 
To get the output in a plot-friendly form, re-route the standard output to a file, e.g.

```mpirun -np <num_proc> program2d.*ex convergence2d.inputs >& convergence.dat```