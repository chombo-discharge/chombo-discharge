## Exec/Examples/AdvectionDiffusion/PipeFlow.

This example solve an advection-diffusion problem with a flow diagonal to the grid inside a pipe.
It was set up from $DISCHARGE_HOME/Physics/AdvectionDiffusion using

> ./setup.py -base_dir=Exec/Examples/AdvectionDiffusion -app_name=DiagonalFlowNoEB

To compile and run in 2D/3D:

```
make -s -j<num_proc> OPT=HIGH DIM=2 program
```

and to run:

```
mpirun -np <num_proc> program<bunch_of_options>.ex example.inputs
```