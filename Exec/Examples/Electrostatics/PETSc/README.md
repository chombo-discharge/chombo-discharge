## Exec/Examples/Electrostatics/PETSc

Example of `FieldSolverAMG`, which assembles the electrostatic operator as a PETSc matrix and solves it with a PETSc
Krylov method preconditioned by algebraic multigrid.

The domain is the unit box without embedded boundaries. The x-low side is at 1 V and the x-high side is grounded; the
other sides have zero normal field. A Gaussian space charge blob sits in the middle of the box
(`FieldStepper.init_rho`, `rho_center`, `rho_radius`). `RefinementBoxTagger` refines two nested boxes around the blob,
giving three AMR levels.

The same problem can be solved with `FieldSolverGMG` by setting `PETScExample.field_solver = gmg`, so that the two
plot files can be compared.

PETSc must be installed and `PETSC_DIR` set. To compile, type

```make -s -j<num_proc> USE_PETSC=TRUE MPI=TRUE DIM=2```

```make -s -j<num_proc> USE_PETSC=TRUE MPI=TRUE DIM=3```

To run it, type

```mpirun -np <num_proc> main2d.*.ex example2d.inputs```

```mpirun -np <num_proc> main3d.*.ex example3d.inputs```

The solver settings are `FieldSolverAMG.*` in the inputs file. PETSc options can also be given on the command line with
the prefix `-fieldsolveramg_`, e.g. `-fieldsolveramg_ksp_monitor` to print the residual at every iteration.

# Output

Output is given to HDF5 files in the `plt` folder: the potential, the space charge, the electric field and the residual
on every level.
