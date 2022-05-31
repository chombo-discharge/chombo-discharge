## Convergence/Electrostatics/C2

This example computes convergence rates for the Electrostatics physics module, using a profiled surface geometry. 
The solution errors are computed by coarsening a numerical solution with a finer (2x) resolution. 

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
