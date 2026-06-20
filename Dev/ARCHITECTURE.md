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

## Cell-sorting & KMC on a per-patch SoA (the demanding case, §8 of USAGE_PATTERNS)

Because storage is per-patch, "sort by cell" becomes an **index over the patch
arrays**, not a separate per-cell container:

- **CSR-style binning** (recommended): permute the patch's columns into cell order
  once, and keep a `begin/end` offset per cell. KMC reactions, per-cell accounting,
  and merge/split then operate on **contiguous ranges** of the patch columns — cache
  friendly, and many containers (Ito species, CDR photo products, photons) can share
  the same cell-offset structure.
- Reactions and merge/split change per-cell counts: append new particles to the
  patch SoA, swap-pop removed ones, then rebuild the cell offsets (cheap relative to
  the chemistry). A per-split field-reconcile callback operates through the
  member-pointer column accessors.

## Open follow-on questions (for later)

1. Exact cell-index structure: CSR offsets vs per-cell index vectors vs a hybrid for
   cut cells. Interacts with how dynamic the per-cell counts are during KMC.
2. Whether `LayoutData<ParticleSoA<P>>` wraps Chombo's existing box/DataIndex
   machinery directly, or we introduce a thin owner that also tracks per-patch
   capacity for reuse across steps.
3. Buffer/mask/cache holders: same `ParticleSoA` type, or a lighter variant (some
   never need cell-sorting or HDF5).
