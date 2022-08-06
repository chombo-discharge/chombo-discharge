## Exec/Examples/StreamerInception/Cylinder

This example runs a streamer inception model (for air) for an cylinder electrode. 
It was set up from $DISCHARGE_HOME/Physics/StreamerInception using

> ./setup.py -base_dir=Exec/Examples/StreamerInception -app_name=Cylinder -geometry=Cylinder

To compile it in 2D/3D

> make -s -j<num_proc> OPT=HIGH DIM=3 program
> make -s -j<num_proc> OPT=HIGH DIM=3 program

To run it in 2D/3D, type

> mpirun -np <num_proc> program2d.*.ex example2d.inputs
> mpirun -np <num_proc> program3d.*.ex example3d.inputs
