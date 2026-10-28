## Exec/Examples/Electrostatics/DielectricVoid

Electrostatics in the geometry of `Exec/Examples/ItoKMC/PartialDischarge`: a gas void (a sphere in 3D, a cylinder in
2D) inside a dielectric, in a uniform field set up by the potentials 0 and 1 on the low and high y sides of the domain.
The mesh and geometry inputs are those of the PartialDischarge example.

The void is small compared to the domain, so the exact solution is that of a void in an infinite dielectric: a uniform
field inside the void, and the uniform field plus a dipole outside it. After the solve, the example compares the
following with the exact solution, in the cells near the void:

* A: the computed potential at the cell centers;
* B: the electric field from the computed potential, at the cell centers;
* C: the electric field from the exact potential (the gradient operator alone), at the cell centers;
* D: B interpolated to the cell centroids, which is what the plot files show;
* E: the exact electric field interpolated from the cell centers to the centroids (the interpolation alone);
* F: C interpolated to the centroids.

The errors are written to `pout.0`, per phase, level and cell type (regular or cut), together with the cut cells with the
largest error D. Neglecting the images of the dipole in the domain sides gives the exact solution an error of the order
of (radius/domain size)^2 in 2D and (radius/domain size)^3 in 3D, which bounds how small the measured errors can become.

The field solver is selected with `DielectricVoid.field_solver`: `gmg` (`FieldSolverGMG`) or `amg` (`FieldSolverAMG`,
which needs `USE_PETSC=TRUE`). With `DielectricVoid.box_tagger = true` and `Driver.initial_regrids = 1`, the box of
`RefinementBoxTagger` is refined as well, e.g. to put the whole void on the finest level.

To compile, type

```make -s -j<num_proc> USE_PETSC=TRUE MPI=TRUE DIM=2```

```make -s -j<num_proc> USE_PETSC=TRUE MPI=TRUE DIM=3```

or leave out `USE_PETSC=TRUE` to build without `FieldSolverAMG`. To run it, type

```mpirun -np <num_proc> main2d.*.ex example2d.inputs```

```mpirun -np <num_proc> main3d.*.ex example3d.inputs```
