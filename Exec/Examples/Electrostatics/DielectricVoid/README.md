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

By default the interface condition is imposed at the centroids of the EB faces of the cut cells, which lie off the curved
void surface, and the electric field in the cut cells is first order. With `FieldSolverGMG.gmg_interface_points = true`
(or `FieldSolverAMG.interface_points = true`) it is imposed at the closest points on the void surface instead. With
second-order interface stencils, i.e.

```
FieldSolverGMG.gmg_use_default_settings = false
FieldSolverGMG.gmg_jump_order           = 2
FieldSolverGMG.gmg_jump_weight          = 1
FieldSolverGMG.gmg_bc_order             = 2
FieldSolverGMG.gmg_bc_weight            = 1
```

(or `FieldSolverAMG.eb_order = 2` and `FieldSolverAMG.eb_weight = 1`), the error D in the cut cells is then second order
in 2D. The refinement study below puts the void on level `N` of a larger domain, so that the images of the dipole matter
less:

```
mpirun -np <num_proc> main2d.*.ex example2d.inputs "AmrMesh.lo_corner=-8E-3 -8E-3" "AmrMesh.hi_corner=8E-3 8E-3" \
  AmrMesh.max_amr_depth=N "AmrMesh.ref_rat=2 2 2 2 2 2 2 2 2" Driver.refine_dielectrics=N \
  DielectricVoid.box_tagger=true Driver.initial_regrids=1 RefinementBoxTagger.ref_box1_lvl=N
```

To compile, type

```make -s -j<num_proc> USE_PETSC=TRUE MPI=TRUE DIM=2```

```make -s -j<num_proc> USE_PETSC=TRUE MPI=TRUE DIM=3```

or leave out `USE_PETSC=TRUE` to build without `FieldSolverAMG`. To run it, type

```mpirun -np <num_proc> main2d.*.ex example2d.inputs```

```mpirun -np <num_proc> main3d.*.ex example3d.inputs```
