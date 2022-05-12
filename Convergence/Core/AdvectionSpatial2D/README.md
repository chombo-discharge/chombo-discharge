## Convergence/Core/AdvectionSpatial2D

This example computes spatial convergence rates for CdrGodunov, a MUSCL-type advection solver.
It uses the initial data and velocity fields from Physics/AdvectionDiffusion.
If changes are made to the velocity/initial fields in that physics module, this needs needs to be modified. 

# Compilation

To compile:

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 program```

# Running the example

To run with MPI:

```mpirun -np <num_proc> program2d.*ex convergence2d.inputs```

# Output

The output shows the Linf/L1/L2 errors which are obtained with the help of an exact solution.
To get the output in a plot-friendly form, re-route the standard output to a file, e.g.

```mpirun -np <num_proc> program2d.*ex convergence2d.inputs >& convergence.dat```