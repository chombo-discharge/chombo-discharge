## Examples/StreamerInception/MechanicalShaft

This example solves a streamer inception problem for complex geometry consisting of an electrode and an insulator.
It was set up from $DISCHARGE_HOME/Physics/StreamerInception using

```./setup.py -base_dir=Exec/Examples/StreamerInception -app_name=MechanicalShaft -geometry=MechanicalShaft```

The program uses input data for atmospheric-pressure air, computed using BOLSIG+.

The example program is set up to run in stationary mode (no transient evolution).
To compile it in 2D

```make -s -j<num_proc> OPT=HIGH DIM=2 program```

To run it, type

```mpirun -np <num_proc> program2d.*.ex example.inputs```

Likewise in 3D, do

```make -s -j<num_proc> OPT=HIGH DIM=3 program```

To run it, type

```mpirun -np <num_proc> program3d.*.ex example.inputs```

Note that the 3D example requires lot of computational resources, and the example needs to be scaled down if running it on a laptop/workstation.