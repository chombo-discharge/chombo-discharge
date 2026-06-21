# Container architecture: SoA lives per-patch

Decision: the SoA storage granularity is **per-patch**. `ParticleSoA<P>` holds all
the particles of a single grid patch (box) — it is the direct replacement for the
`List<P>` that today lives inside one `BinFab<P>` / Chombo `ParticleData` box.

```
ParticleContainer<P>                         (AMR hierarchy, one per realm)
  └─ per level l:  LayoutData<ParticleSoA<P>>   (one ParticleSoA per box on this rank)
       └─ ParticleSoA<P>                        (one patch: tuple<vector<col>...>)
            └─ column<&P::field>()              (contiguous per-field array)
```

The parallel data holders that `ParticleContainer` already keeps — valid, buffer
(grown grids), mask (halo), cache (regrid) — each become a
`LayoutData<ParticleSoA<P>>` with the same per-patch granularity.

## Why per-patch (vs per-cell or per-level)

- **Per-cell storage** would shatter each column into thousands of tiny arrays, one
  per cell — destroying the contiguity that is the whole point of SoA, and adding
  huge allocation overhead. Cell structure is instead expressed as a lightweight
  *index* over the per-patch arrays (below).
- **Per-level storage** would fight the existing box/rank ownership model and the
  `#pragma omp parallel for` over boxes. Per-patch maps cleanly onto both.
- **Per-patch** keeps columns large and contiguous (vectorizable deposit/interp,
  cheap whole-column reductions) while still being the unit of MPI ownership and
  OpenMP parallelism.

## How the existing operations map

- **OpenMP**: the current `#pragma omp parallel for` over boxes parallelizes over
  per-patch `ParticleSoA`s; each thread owns its patch → no contention for
  deposit/interp/per-cell work.
- **Remap / transfer**: particles staying on a patch cost nothing; particles moving
  between patches on the same rank are appended to the destination SoA and
  swap-popped from the source. Across ranks: `linearizeParticle` the movers, scatter
  bytes, `delinearizeAndAppend` into the destination patch SoA. Whole-patch moves
  (regrid caching) are `std::move` of the column vectors.
- **MPI**: generic per-particle (or batched per-column) linearization already
  prototyped; column-major batching is a future speedup over today's per-list-item
  packing.

## LOCKED design decisions (converged)

- **Storage: arena-backed SoA per patch.** One aligned allocation holds all columns as
  offset slices (`CD_ParticleSoAArena.H`). Rationale: one-allocation build, *robust*
  no-intrinsics bulk pack/copy (a single `memcpy`, or zero-copy `MPI_Send(data())`),
  aligned SIMD. Requires `reserve`/container-reuse (growth without it is ~7x).
- **Vector fields are stored PER-COMPONENT as raw `ParticleReal` `x`/`y`/`z` scalar
  columns — NOT `RealVect` columns.** Each scalar member is one column (existing column
  machinery, no auto-split). Position is designated by `SpaceDim` member pointers; the
  `position(i)` accessor returns a *promoted* `double` `RealVect` by value for cell-lookup
  / Chombo interop. Rationale: per-component vectorizes the field-update kernels (Euler
  advance 4x, transform 5.6x; interleaved `RealVect` stayed scalar), is GPU-optimal, and
  raw scalars are what allow `ParticleReal != Real` (`RealVect` is hardwired `double`).
  Costs: dimension-independent struct declaration needs `D_DECL`/`#if CH_SPACEDIM`; more
  columns → slightly worse scatter/remap (~1.2-1.4x vs `vector<P>`, still ~2x over
  `List`). Deferred sugar: a `ParticleRealVect = std::array<ParticleReal,SpaceDim>` member
  with container auto-split. Supersedes the earlier per-field `vector<RealVect>` choice.
- **Particle precision: a compile-time `ParticleReal`, independent of the mesh `Real`.**
  Particle columns are `ParticleReal` (e.g. `float`) while the grid/`FArrayBox` stays
  `Real` (`double`). Benefits: ~half memory/bandwidth, 8-wide SIMD, smaller MPI buffers.
  Mixed-precision kernels: compute cell index + CIC weights in `double` (promote particle
  position), accumulate into the `double` grid, demote on interpolate. (If `float` global
  positions lose sub-cell resolution on large grids, switch to cell-relative coords.)
- **Cell binning is a property of the SoA's order + a CSR offsets array — NOT a separate
  container.** Drop today's second representation (`BinFab<List<P>>`) and the
  `organizeByCell/byPatch` copy-conversion entirely.

## Cell-sorting & KMC: one container, cell-sorted, with CSR offsets

Today `ParticleContainer` keeps TWO representations — patch-`List<P>` and a separate
cell-sorted `BinFab<List<P>>` — and physically copies particles between them. The
`BinFab<List>` is doubly allocation-heavy (a `List` per cell + a node per particle: the
~6x build cost measured in BENCHMARKS.md), plus linked-list cache misses per cell.

With the contiguous SoA this collapses to ONE container that is canonically cell-sorted:

1. compute a scratch cell key per particle (`floor((x-lo)/dx)` -> Fortran linear index),
2. **counting-sort** the SoA columns by that key (O(N), no allocation),
3. build per-cell **CSR offsets** `cellStart[c]` so cell `c` owns `[cellStart[c], cellStart[c+1])`.

KMC reactions, per-cell accounting, and merge/split then operate on **contiguous ranges**
of the single SoA (per-cell `Sum weight` is a contiguous reduction; merge works on a
cell's range). Deposition over the cell-sorted SoA also gets write-locality for free.
This removes the second container, the per-cell `List` allocations, the per-node `new`,
and the copy-conversion — replaced by a counting sort + an `int` offsets array.

Maintenance:
- **Motion:** CFL < 1 cell/step means after `advance` the array is *nearly* sorted, so an
  incremental re-bin of only the boundary-crossers keeps it ordered (cheaper than a full
  rebuild).
- **KMC create/destroy:** append new particles to an unsorted tail, swap-pop/tombstone
  removals during the per-cell pass, then one re-sort/compact at the end of the phase.
- Physically reorder the columns (not just an index permutation) so per-cell access is
  contiguous for the repeated KMC/merge sweeps (a one-time permutation amortizes).

## Open follow-on questions (for later)

1. **Per-component expression:** does the user particle struct declare `RealVect`
   members that the container auto-splits into `SpaceDim` columns (keeps the `RealVect`
   AoS view for `gather`/`scatter`), or declare raw `x`/`y`/`z` scalars? (Auto-split is
   nicer ergonomically, more traits machinery.)
2. **Particle metadata:** do particles still carry `particleID`/`rankID` (determinism,
   checkpoint-restart, debugging)? As columns?
3. **Checkpoint format + versioning:** column order/layout is the on-disk ABI, and
   per-component changes it — needs a checkpoint version key.
4. **Reserve/capacity lifecycle:** reserve-at-regrid heuristic (prev count x ~1.3),
   never shrink, pool/reuse containers across regrids so the arena self-sizes.
5. **AMR integration:** how `LayoutData<arena-SoA>` plugs into the box/`DataIndex`
   machinery, regrid caching, and the valid/buffer/mask/cache holders (same type or a
   lighter variant for holders that never bin or checkpoint).
6. **Remap protocol:** whole-patch transfer on ownership change (zero-copy `MPI_Send`)
   vs subset gather-pack for boundary-crossers — two code paths.
