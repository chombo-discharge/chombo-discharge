## Exec/Examples/BrownianWalker/ParticleMerge

This example measures what `ItoSolver`'s super-particle merging does to a particle population that is
already correct, so that every deviation it reports is the merge's own doing.

It is a `BrownianWalkerStepper` application, so it runs real time steps through the real machinery --
transport, remap, deposit, regrid, merge. `ParticleMergeStepper` seeds
`ParticleMergeStepper.new_particles_per_cell` fresh unit-weight particles per cell before every advance,
on top of whatever survived the previous merge, and the base stepper then transports the union and merges
it back down to `ItoSolver.particles_per_cell`. Repeating that is what exposes an error that compounds,
which is the failure mode that matters for an operator applied once per step for thousands of steps.

Run it with `ParticleMergeStepper.diffusion = false` so the only motion is the advective field. Under
`ParticleMergeStepper.velocity_field = constant` the population translates rigidly, which leaves it uniform and
keeps the merge the only thing acting on its shape. There is no geometry and no chemistry.

The grid is two-level. `ParticleMergeTagger` refines a fixed slab, bounded in one direction and
spanning the domain in the others, so the coarse-fine boundary is a flat plane at a position the user
chose. Because the refined region does not depend on the solution it does not move, and the same
geometry is seen on every round. That matters because a merge's sub-cell artifacts are locked to the
*local* cell width: the two levels alias on lattices whose periods differ by the refinement ratio, and
they meet at that plane. The report is broken out per level so the two can be compared directly.

Statistics are off by default. Set `ParticleMergeStepper.print_stats = true` to have the per-level report
written to `pout` on every step, and `ParticleMergeStepper.plot_particles = true` to write the particles
themselves to H5Part (be aware that at the shipped resolution each file is order a gigabyte).

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

Cells covered by a finer level are excluded from both the seeding and the statistics: they hold no
particles, because a particle lives on the finest level that covers it, so counting them would report
empty cells that are not empty. Note that the density ratio spans transport *and* merge -- advection
moves particles between cells too -- while the sub-cell, weight-spread and block-frequency numbers
describe the post-merge population directly.

### Running

```
make -s -j<num_proc> OPT=HIGH DIM=2
mpirun -np <num_proc> main2d*.ex example2d.inputs
```

The report is written to `pout.0`. Set `ParticleMergeStepper.plot_particles = true` to also write the
particles on every plot step to `particles/<solver name>/<solver name>.step*.<dim>d.h5part`, which can be
opened directly in a particle viewer -- the sub-cell lattice is visible by eye.

Point `ItoSolver.merge_method` at any of the merge algorithms to compare them. The most instructive
contrasts:

- `kd_cell` with `ItoSolver.kd_placement` set to `centroid`, `sample` and `random` in turn -- the
  same partition, differing only in where each leaf's particle is put.
- `ItoSolver.kd_partition` set to `weight` and `count` -- the same placement, differing only in
  how a node is divided. Watch `N_eff`.
- `kd_cell` against `kd_patch` and `kd_amr` -- the same two axes over a wider scope. Their leaves are
  not confined to one cell, so they can merge across a cell face; watch the per-cell density instead.
- `ParticleMergeTagger.ref_box1_*` -- move or resize the refined slab, or set `num_ref_boxes = 0` to run
  on a single level.

### A note on the option prefix

`ParticleMergeStepper` derives from `BrownianWalkerStepper`, and hands its own class name down to the
base class, which reads every one of its options under that prefix instead of the default
`BrownianWalker`. An input script therefore has one prefix per class, and it is the class the
application actually instantiates. `CD_ParticleMergeStepper.options` lists the full inherited set
alongside the two options this class adds.
