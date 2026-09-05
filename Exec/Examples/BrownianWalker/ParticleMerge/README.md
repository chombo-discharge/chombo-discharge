## Exec/Examples/BrownianWalker/ParticleMerge

This example measures what `ItoSolver`'s super-particle merging does to a particle population that is
already correct, so that every deviation it reports is the merge's own doing.

The particles do not move: there is no drift, no diffusion, no chemistry and no geometry. Each round
seeds a fresh uniform population on top of whatever survived the previous merge and merges the union
back down to `ItoSolver.particles_per_cell`. Repeating that is what exposes any error that compounds,
which is the failure mode that matters -- a merge run once per step for thousands of steps.

The grid is two-level. `ParticleMergeTagger` refines a fixed slab, bounded in one direction and
spanning the domain in the others, so the coarse-fine boundary is a flat plane at a position the user
chose. Because the refined region does not depend on the solution it does not move, and the same
geometry is seen on every round. That matters because a merge's sub-cell artifacts are locked to the
*local* cell width: the two levels alias on lattices whose periods differ by the refinement ratio, and
they meet at that plane. The report is broken out per level so the two can be compared directly.

### What it reports, per level

- **Sub-cell position distribution.** Histograms of the position within the owning cell, plus the
  Fourier amplitudes `|Z_k| = |<exp(2*pi*i*k*u)>|` for `k = 1..8`. A correct merge leaves these at the
  shot floor (`~1/sqrt(N)`). A merge that places each group at its centroid drives the particles onto
  a sub-cell lattice, which shows up as a large `|Z_k|` at `k = ppc^(1/D)`.
- **Sub-cell spread.** `E[(u-0.5)^2]`, which is `1/12 = 0.0833333` for a uniform distribution.
  Centroid placement contracts it.
- **Per-cell density.** The per-cell weight before and after the merge. A merge that leaves the
  density alone reproduces it exactly; one that moves weight across cell faces does not, and the error
  compounds over rounds.
- **Block frequency.** Mean density against the cell's phase within its patch. A merge that treats
  patch-boundary cells differently from interior cells shows up here and nowhere else, which matters
  because that error moves when the decomposition does.
- **Weight distribution.** `N_eff = (sum w)^2 / (sum w^2)` per cell, the heaviest particle's share of
  its cell, and the weight range. `N_eff` equal to the target means the cell's mass is spread evenly
  over its particles.
- Total weight, population, and particles per cell.

Cells within `ParticleMerge.margin` of a domain face are excluded, so every analysed cell sits in a
fully populated neighbourhood and the domain edge cannot contaminate the result. Set `margin` to a
multiple of the block size to keep the block-frequency bins balanced.

### Running

```
make -s -j<num_proc> OPT=HIGH DIM=2
mpirun -np <num_proc> main2d*.ex example2d.inputs
```

The report is written to `pout.0`. Set `ParticleMerge.write_particles = 1` to also write the initial
and final populations to `parts_<tag>_initial.h5part` and `parts_<tag>_final.h5part`, which can be
opened directly in a particle viewer -- the sub-cell lattice is visible by eye.

Point `ItoSolver.merge_algorithm` at any of the merge algorithms to compare them. The most instructive
contrasts:

- `kd_cell` with `ItoSolver.kd_cell_placement` set to `centroid`, `sample` and `random` in turn -- the
  same partition, differing only in where each leaf's particle is put.
- `ItoSolver.kd_cell_partition` set to `weight` and `count` -- the same placement, differing only in
  how a node is divided. Watch `N_eff`.
- `ParticleMerge.stratified = 1` against `0` -- an exactly-`init_ppc`-per-cell fill makes a
  count-median tree bisect onto cell faces, which is a degenerate case worth knowing about but not
  representative of a running simulation.

### A note on the location

This example lives under `Exec/Examples/BrownianWalker` because merging is what a Brownian-walker
application spends its particle budget on, but it does not use `BrownianWalkerStepper` or any other
physics module -- it drives `AmrMesh` and `ItoSolver` directly, so the only dependency is
`discharge-lib`.
