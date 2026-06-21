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
- **ALL mandatory fields are container-owned columns; the user struct is payload-only.**
  The container always allocates `position` (`SpaceDim` raw `ParticleReal` x/y/z columns),
  `weight` (`ParticleReal`), `particleID`, and `rankID`. The user declares ONLY the extra
  payload (e.g. `velocity`, `mobility`, `oldPosition`) via a payload struct + its
  `ParticleTraits`; the payload may be empty, so `ParticleSoA<>` is a ready-made
  point/tracer particle. Accessors `position(i)`/`weight(i)`/`particleID(i)`/`rankID(i)`.
  Rationale: the position/weight layout is *already* fixed by the locked decisions (raw
  per-component `ParticleReal`), so user control over them buys nothing but a chance to
  misdeclare — baking them in gives a uniform mandatory schema, deletes the
  `positionPtr`/`weightPtr` traits and the mandatory-field `static_assert`, and makes
  omission/reordering impossible. NO inheritance (re-adds BinItem-style coupling, breaks
  standard-layout). Trade-off: position/weight are reached via accessors rather than the
  payload struct, so `gather`/`scatter`/merge split into "mandatory accessors + payload
  columns" — at least as clear as today's whole-object access. (Supersedes the earlier
  "position/weight are user-declared physics fields enforced by traits" choice.)
- **No checkpoint versioning / capacity never-shrinks / full redesign (no facade).** H5 =
  the type's `h5Columns` (no cross-type restart compat expected); arenas grow
  geometrically, never shrink, and are reused across steps (self-sizing).

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

## The AMR / container layer (LOCKED)

Everything above describes the per-patch arena-SoA *leaf*. The `ParticleContainer`
equivalent owns particles across the AMR hierarchy and ties into Chombo's grid machinery:

1. **Per-level ownership:** `LayoutData<arena-SoA>` indexed by `DataIndex` over the
   level's `DisjointBoxLayout` — one arena-SoA per locally-owned box, iterated with
   `DataIterator` (the `#pragma omp for` over boxes). Mirrors today's `ParticleData<List>`.
   *Template constraints check (Chombo `LayoutData<T>`):* `LayoutData` stores `Vector<T*>`,
   `new T`-default-constructs each element, indexes by reference, and `delete`s them — it
   **never copies, moves, or assigns `T`** (its own copy ctor/assignment are private). So the
   only requirements on `T = ParticleSoA<P>` are **default-constructible + destructible**,
   both satisfied; **move-only is a non-issue**. (Cross-rank exchange will NOT go through
   `LevelData<T>::copyTo` — which would demand `T::define`/linearization — but through our own
   `linearizeParticle`-based MPI, so that heavier interface does not apply.)
2. **Halo/transfer holders use the SAME arena-SoA leaf type** — the buffer (grown grids),
   mask (halo), cache (regrid), and the remap pool are all `LayoutData<arena-SoA>` (the
   send side is a small set of per-destination-rank arena-SoAs). Same columns → zero
   conversion: `append`/`gather`/`linearizeParticle`/`delinearizeAndAppend` and the
   zero-copy `data()`/`byteSpan()` whole-pool send all work unchanged. The CSR cell-sort
   state simply stays empty in transient buffers — NOT worth a second "lighter" type.
3. **Regrid:** `preRegrid` caches particles off the old layout (into the cache holder);
   `regrid` rebuilds the `LayoutData` on the new `DisjointBoxLayout` and redistributes to
   the new owning boxes/ranks.
4. **Remap protocol (the pool model):** collect all movers into a pool, keep the ones that
   still belong to a box on this rank (local append into the destination arena-SoA), and
   scatter the rest by rank. The pool IS an arena-SoA (point 2). Two transfer shapes:
   **whole-patch ownership change** → zero-copy `MPI_Send(data())`; **individual
   boundary-crossers** → gather the leaving subset into a per-rank send arena-SoA.
5. **The leaf does NOT store its Box.** The `DisjointBoxLayout` (`dbl[dataIndex]`) is the
   single source of truth; geometry-dependent leaf methods (cell-sort key, deposition
   region) take `box`/`dx`/`probLo` as arguments — the container already holds the
   `DataIndex` when it iterates, so `leaf.binByCell(box, dx, probLo)` is trivial. Storing
   the box would duplicate state that desyncs on regrid and would make the leaf awkward as
   a halo/remap-pool buffer (where "the box" is absent). The arena stays columns-only
   regardless, so `data()`/`byteSpan()` are unaffected. (Fallback if too many call sites
   need it: store `m_box`/`m_dx` as plain members *outside* the arena — start pure.)

The mechanics (LayoutData/DataIterator/regrid hooks) mirror today's `ParticleContainer`.
