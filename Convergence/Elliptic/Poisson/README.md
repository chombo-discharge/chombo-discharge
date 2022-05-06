## Convergence/Poisson

This example computes convergence rates for MFHelmholtzOp, a Helmholtz operator with support for discontinuous material coefficients.
The example uses a "coaxial cable" geometry with an embedded dielectric.

# Compilation

To compile:

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 program```

# Running the example

To run with MPI:

```mpirun -np <num_proc> program2d.*ex convergence2d.inputs```

# Output

The output shows the Linf/L1/L2 errors which are obtained with the help of an exact solution.
