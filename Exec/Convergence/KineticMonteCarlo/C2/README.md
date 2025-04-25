## Convergence/KineticMonteCarlo/C2

This example runs a Kinetic Monte Carlo solver for the Schlogl model.
The rate constants are the same as in [Cao. et. al in J. Chem. Phys. Vol. 124 044109 (2006)](https://doi.org/10.1063/1.2159468)
For each run of the model the population will converge to one of the bi-stable states. 

# Compilation

To compile:

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 program```

# Running the example

To run with MPI:

```mpirun -np <num_proc> program2d.*ex convergence.inputs```

The user can select between different algorithms and initial conditions in the input script. 

# Output

Output is given in the pout.* files.
The files contain data in the format

"Time" "KMC solution"
