## Exec/Examples/ItoDiffusion/ParticleMerge

This example measures what `ItoSolver`'s super-particle merging does to a particle population that is
already correct, so that every deviation it reports is the merge's own doing.

It fills a single-level grid with a uniform particle population, merges it down to
`ItoSolver.particles_per_cell`, and reduces the result into a set of diagnostics. The fill can be
repeated: each round adds a fresh population on top of the survivors and merges the union back down,
which is what a running simulation does and what exposes any error that compounds.

There is no geometry, no transport and no chemistry -- only the merge.

### What it reports

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

Edge effects are excluded by analysing only cells at least `ParticleMerge.margin` cells from every
domain face, so every analysed cell sits in a fully populated neighbourhood. Set `margin` to a
multiple of the block size to keep the block-frequency bins balanced.

### Running

```
make -s -j<num_proc> OPT=HIGH DIM=2
mpirun -np <num_proc> main2d*.ex example2d.inputs
```

The report is written to `pout.0`. Set `ParticleMerge.write_particles = 1` to also write the initial
and final populations to `parts_<tag>_initial.h5part` and `parts_<tag>_final.h5part`, which can be
opened directly in a particle viewer -- the sub-cell lattice is visible by eye.

Point `ItoSolver.merge_algorithm` at any of the merge algorithms to compare them. Two useful
contrasts:

- `equal_weight_kd` against `equal_weight_kd_sampled` -- the same per-cell partition, differing only
  in whether each leaf is placed at its weighted centroid or at one of its own members.
- `ParticleMerge.stratified = 1` against `0` -- an exactly-`init_ppc`-per-cell fill makes a
  count-median tree bisect onto cell faces, which is a degenerate case worth knowing about but not
  representative of a running simulation.
