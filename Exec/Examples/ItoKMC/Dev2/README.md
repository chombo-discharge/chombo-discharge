## Exec/Examples/ItoKMC/AirDBD

This example runs a model of 2D and 3D streamer discharges over a dielectric surface.
the input scripts and plasma chemistry are defined in

* Plasma chemistry: air_chemistry.json
* Input files: example.inputs

This model uses an Ito-KMC model.
For a further explanation to the parameters in the chemistry file, see the [chombo-discharge documentation](https://chombo-discharge.github.io/chombo-discharge/Applications/ItoKMC.html#json-0d-chemistry-interface).

# Compilation

To compile:

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 program```

or in 3D

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=3 program```

# Running the example

To run with MPI:

```mpirun -np <num_proc> program2d.*ex example.inputs

or in 3D

```mpirun -np <num_proc> program3d.*ex example.inputs

The user can select between different algorithms and initial conditions in the input script and chemistry file. 

# Output

Output is given to HDF5 files in the plt folder.
