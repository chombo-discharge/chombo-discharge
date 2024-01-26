## Examples/StreamerInception/WireWire

This example solves a streamer inception problem for two contacting wires. 
It was set up from $DISCHARGE_HOME/Physics/StreamerInception using

```./setup.py -base_dir=Exec/Examples/StreamerInception -app_name=PartialDischarge -geometry=WireWire```

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

The 3D example might require substantial computational resources, and the example might need to be scaled down if running it on a laptop/workstation.