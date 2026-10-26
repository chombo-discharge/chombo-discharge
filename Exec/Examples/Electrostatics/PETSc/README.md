## Exec/Examples/Electrostatics/PETSc

Example of `FieldSolverAMG`, which assembles the electrostatic operator as a PETSc matrix and solves it with a PETSc
Krylov method preconditioned by algebraic multigrid.

The problem is the same as in `Exec/Examples/Electrostatics/ProfiledSurface`: a live wheel electrode above a dielectric
plate with a profiled surface (`DiskProfiledPlane`). The mesh is refined around the electrode and the dielectric, and
the inputs files are those of `ProfiledSurface` with a `FieldSolverAMG` block added. The embedded-boundary stencils of
`FieldSolverAMG` (`eb_order = 1`, `eb_weight = 4`) are those of `FieldSolverGMG` with
`gmg_use_default_settings = true`.

The same problem can be solved with `FieldSolverGMG` by setting `PETScExample.field_solver = gmg`, so that the two
plot files can be compared.

PETSc must be installed and `PETSC_DIR` set. To compile, type

```make -s -j<num_proc> USE_PETSC=TRUE MPI=TRUE DIM=2```

```make -s -j<num_proc> USE_PETSC=TRUE MPI=TRUE DIM=3```

To run it, type

```mpirun -np <num_proc> main2d.*.ex example2d.inputs```

```mpirun -np <num_proc> main3d.*.ex example3d.inputs```

The `regression2d.inputs` and `regression3d.inputs` files are coarser versions of the same problem, for the test suite.

The solver settings are `FieldSolverAMG.*` in the inputs file. `FieldSolverAMG.petsc_verbosity` prints the residual
at every Krylov iteration (1), the true residual as well (2), and the solver and preconditioner (3). Other PETSc options can
be given with the prefix `-fieldsolveramg_` through the `PETSC_OPTIONS` environment variable, e.g.

```PETSC_OPTIONS="-fieldsolveramg_pc_gamg_threshold 0.02" mpirun -x PETSC_OPTIONS -np <num_proc> main2d.*.ex example2d.inputs```

Options without a value, e.g. `-fieldsolveramg_ksp_monitor`, can also be given on the command line. Options with a
value cannot, because `ParmParse` reads the command line too and aborts on them.

# Output

Output is given to HDF5 files in the `plt` folder: the potential, the space charge, the electric field and the surface
charge on every level.
