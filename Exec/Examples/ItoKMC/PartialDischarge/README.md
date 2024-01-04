## Exec/Examples/ItoKMC/PartialDischarge

This example (which is a work-in-progress) runs a model of 2D and 3D streamer discharges in a void geometry.
The input scripts and plasma chemistry are defined in

* Plasma chemistry: simple_air_chemistry.json
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

# Caveats

This simulation is fully stochastic and starts from a single physical electron.
There is a substantial probability that the initial electron is attached before it ionizes, in which case the simulation will still proceed but discharges fail to initiate.
Another common event is that the initial electron undergoes ionization but the resulting ions do not generate secondary electron emission at the cathode side of the avoid.

# Output

Output is given to HDF5 files in the plt folder.
