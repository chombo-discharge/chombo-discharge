## Exec/Examples/KineticMonteCarlo/AirDischarge

This example runs a Kinetic Monte Carlo solver for a simple discharge model in air. 
The transport data is computed using BOLSIG+, and the user sets a uniform background field.

# Compilation

To compile:

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 program```

# Running the example

To run with MPI:

```mpirun -np <num_proc> program2d.*ex example.inputs```

The user can select between different algorithms and initial conditions in the input script. 

# Output

Output is given in what are normally the pout.* files.
To compare reactions algorithms, these are instead named by the reaction algorithm, e.g., instead of pout.0 one will see explicit_euler.0, ssa.0, etc.
The files contain data in the format

"Time" "Species population numbers (array)"

If running in parallel with MPI, each rank will use a different RNG seed when running the program and output a different state.
