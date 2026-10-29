## Exec/Tests/ItoKMC/DomainEmission

Regression test for secondary emission on the faces of the computational domain in the Ito-KMC model.

The geometry (a rod in 2D, a needle in 3D) and the solver settings are those of `Exec/Tests/ItoKMC/JSON`. The
differences are that the electrode is at a negative potential, so the domain faces are anodes and the electrons
drift into them, and that the time step is not capped, so the length of the run is set by the step count in
`ItoKMC.ini`. The chemistry (`chemistry.json`) turns every electron, and half of the photons, that leave the domain
into a neutral tracer species `M`, released in the boundary cell they left through. Field emission on the low y face
keeps the domain emission pass running every step.

# Compilation

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=2 main```

or in 3D

```make -s -j<num_proc> OPT=HIGH DEBUG=FALSE DIM=3 main```

# Running the test

```mpirun -np <num_proc> main2d.*ex regression2d.inputs```

or in 3D

```mpirun -np <num_proc> main3d.*ex regression3d.inputs```
