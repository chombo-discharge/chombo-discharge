# Particle library redesign — working notes (Dev/)

## Goal

Redesign the chombo-discharge particle library to get rid of Chombo's `List<P>` and
replace the per-bin/per-patch storage with a **Struct-of-Arrays (SoA)** layout.

This is a large job. The work in `Dev/` is a **staging area** — we do **not** modify
production code under `Source/`, `Physics/`, or `Geometries/` yet. The point of
`Dev/` is to:

1. Capture how particles are *actually* used across the codebase (see
   `USAGE_PATTERNS.md`).
2. Harden the new SoA particle core (`CD_ParticleSoA.H`, `CD_ParticleLoops.H`) and a
   particle type (`CD_DemoParticle.H`) by exercising them on real Chombo objects via
   a chombo-discharge test executable (`main.cpp`, built with the local `GNUmakefile`).
3. Decide the core methodology before touching production code.

## Hard constraints (discovered)

- **C++17** (decided — bumped up from the previous C++14 baseline). Enables the
  single-token member-pointer field selector (`template <auto>`), inline
  `static constexpr` data members, `if constexpr`, fold expressions.

  > **⚠️ REMINDER — project-wide C++17 migration still TODO.** `Make.defs.local` has
  > been bumped to C++17 locally, but the rest of the tree has **not** been updated
  > yet. Before/at integration we must:
  > - set `CXXSTD = 17` in **all** `Lib/Local/**/Make.defs.*` (and
  >   `Make.defs.local.template`), including the `GitHub/` CI defs;
  > - confirm the GCC and Intel oneAPI CI compilers build at C++17;
  > - update docs that state the standard: `CLAUDE.md` (§1 build), the Sphinx docs,
  >   and any README mentioning C++14.
  > Search hint: `grep -rn "CXXSTD" Lib/Local` and `grep -rni "c++14\|c++11" Docs CLAUDE.md`.
- Particle properties must remain **trivially copyable** (memcpy linearization for
  MPI and HDF5).
- Must keep working with the **AMR + embedded-boundary** machinery
  (`ParticleContainer`, remapping, regridding, deposition across refinement
  boundaries, cut cells).
- Must support **OpenMP + MPI** hybrid parallelism (current code declares OMP
  reductions over `List<P>` and scatters particles across ranks).

## Minimum requirements for the new particle core

- A particle is a user-defined **plain struct** (AoS view). The user names their own
  fields, e.g. `struct MyParticle { RealVect pos; Real weight; Real lightsaberDamage; };`.
- The container is templated on that struct and stores it column-wise (SoA).
- **Only two fields are mandatory**: a `RealVect` position and a `Real` weight. The
  SoA container must expose `position(i)` and `weight(i)`. Everything else is generic.
- Generic field access: get a whole column (contiguous, vectorizable) or a single
  particle's field.
- Obvious container ops: append/insert, remove (swap-and-pop), size, reserve, clear,
  gather a single particle back to the AoS view, scatter an AoS particle in.
- **Linearization** for MPI (and an HDF5-subset variant) — pack/unpack one particle
  to/from a byte buffer.

## Current status

- [x] Surveyed `Source/Particle` + all `ParticleContainer` users → `USAGE_PATTERNS.md`.
- [x] SoA container `CD_ParticleSoA.H` (traits-driven, member-pointer columns; derived
      indices; member-pointer field selector; MPI + HDF5-subset linearization;
      append/remove/gather/scatter/get/column; move-only).
- [x] `CD_ParticleLoops.H` — SIMD-decorated particle loop (BoxLoops analogue).
- [x] `CD_DemoParticle.H` — example particle type.
- [x] Exposed to real Chombo: types are `Real`/`RealVect`, namespace `ChomboDischarge`.
      `main.cpp` builds and runs as a chombo-discharge executable via the local
      `GNUmakefile`, exercising add / iterate+remove / SIMD kernel / center-of-mass
      merge. Next: deposition on FArrayBoxes, EB intersections, merging on real grids.

## Decisions (locked)

- **C++ standard: bump `CXXSTD` to 17** — enables the single-token member-pointer
  field selector.
- **Field selector: single-token member pointer, type deduced** —
  `soa.column<&P::velocity>()`, `soa.get<&P::weight>(i)`, `deposit<&P::weight>(...)`.
  Index derived from the member pointer (reorder-safe); return type deduced (no `Ret`
  arg like today).
- **SoA granularity: PER-COMPONENT** — vector fields are stored as separate `x`/`y`/`z`
  scalar columns, NOT `RealVect` columns. (Supersedes the earlier "per-field
  `vector<RealVect>`" choice; driven by the benchmarks — Euler advance 4x, transform
  5.6x, interleaved `RealVect` stayed scalar — and GPU-readiness. Cost: more columns →
  slightly worse scatter/remap, accepted.)
- **Arena backing: single allocation** — all columns are offset slices of one aligned
  buffer (`CD_ParticleSoAArena.H`). Gives one-allocation build, robust no-intrinsics
  bulk pack/copy (one `memcpy` / zero-copy send), aligned SIMD. Requires `reserve` /
  container-reuse (growth without it is ~7x). The per-column-`memcpy` SoA pack is
  bimodal (~40% of launches 4x slow at L3-scale); the arena removes that variance.
- **Particle precision: a compile-time `ParticleReal`, independent of `Real`** — particle
  columns use `ParticleReal` (e.g. `float`) which may differ from the mesh/`Real`
  precision (`double`). Synergistic with per-component (no `RealVect`, which is hardwired
  `double`). Benefits: ~half the memory, 8-wide (vs 4-wide) SIMD, faster bandwidth-bound
  ops + packing.
- **Per-component expression: raw `ParticleReal` x/y/z scalar columns** (NOT `RealVect`
  members, which would force `double`). Each scalar member is one column; the existing
  column machinery works unchanged (no auto-split). Position is designated by `SpaceDim`
  member pointers (`positionPtrs = std::make_tuple(D_DECL(&P::m_x, ...))`); `position(i)`
  returns a *promoted* `double` `RealVect` (by value) for cell-lookup / Chombo interop.
  Cost: dimension-independent struct declaration needs `D_DECL`/`#if CH_SPACEDIM` guards.
  Alternative deferred: a `ParticleRealVect = std::array<ParticleReal,SpaceDim>` member
  with container auto-split — pure ergonomic sugar, add later if raw decls get noisy.
- **Mixed-precision kernels:** compute cell index + CIC weights in `double` (promote the
  `float` particle position), accumulate into the `double` grid, demote on interpolate.
  (If float global positions lose sub-cell resolution on large grids, switch to
  cell-relative coordinates.)
- **Cell binning: order + CSR offsets, not a separate container** — drop the
  `BinFab<List<P>>` second representation; the per-patch SoA is canonically cell-sorted
  (counting sort + `cellStart[]`). See `ARCHITECTURE.md`.
- **All mandatory fields are container-owned; the user struct is payload-only** — the
  container always allocates `position` (SpaceDim raw `ParticleReal` x/y/z), `weight`
  (`ParticleReal`), `particleID`, and `rankID`. The user declares ONLY the extra payload
  via a payload struct + `ParticleTraits`; an empty payload is valid, so `ParticleSoA<>`
  is a ready-made point/tracer particle. Accessors `position(i)`/`weight(i)`/
  `particleID(i)`/`rankID(i)`; id/rank are in MPI linearization but NOT HDF5 (regenerated
  on restart). **No inheritance** (re-adds BinItem coupling, breaks standard-layout).
  Rationale: position/weight layout is already fixed by the locked decisions, so baking
  them in loses nothing, deletes the `positionPtr`/`weightPtr` traits + the mandatory-field
  `static_assert`, and makes misdeclaration impossible. Trade-off: position/weight reached
  via accessors, so gather/scatter/merge = "mandatory accessors + payload columns".
  (Supersedes the earlier "position/weight are user-declared, enforced by traits" choice.)
- **No checkpoint versioning** — no cross-particle-type restart compatibility is
  expected; H5 output = the type's `h5Columns`; write→read round-trips within a build.
  The only invariant preserved is restart on a *different MPI rank count* (same type),
  which is orthogonal to the type definition.
- **Migration: full redesign, no facade / no drop-in.**
- **Capacity: never-shrink + self-resizing** — geometric growth, `clear()` keeps
  capacity, containers reused across steps so each arena warms to its high-water mark
  once. No reserve prediction required.
- **Schema: explicit traits, NOT a macro** — `ParticleTraits<P>` with a `columns`
  member-pointer tuple + `positionPtrs` (SpaceDim) / `weightPtr` + optional `h5Columns`.
  Indices are all **derived** from member pointers (no literal indices → reordering is
  safe). A `CD_PARTICLE_LAYOUT` macro was evaluated and rejected; the file was removed.
- **Storage granularity: per-patch** — one container per patch; the AMR container holds
  `LayoutData<arena-SoA>` per level. See `ARCHITECTURE.md`.

## Benchmark findings (see `BENCHMARKS.md`)

`Dev/Benchmark/` compares `List<P>` vs `vector<P>` (AoS) vs `ParticleSoA<P>` on
deposition, interpolation, streaming transform, build, remap, and MPI packing.
Headline: **escaping the linked list is the big win, not the SoA column split** —
`vector<P>` ties or beats SoA everywhere except SIMD streaming transforms, where SoA
is 1.7x (per-component 4x) faster. So full SoA is justified mainly by hot SIMD
particle kernels and/or GPU offload; otherwise `vector<P>` is the sweet spot. This
reframes the "drop-in vs clean break" question below into "SoA vs vector<P> at all".

## Big open questions

**None — the design is fully locked.** The AMR/container layer is now resolved too (see
`ARCHITECTURE.md` "The AMR / container layer (LOCKED)"):
- per-level `LayoutData<arena-SoA>` over the `DisjointBoxLayout`, indexed by `DataIndex`;
- halo/buffer/mask/cache/remap-pool all use the **same arena-SoA leaf type** (CSR state
  just stays empty in transient buffers — no lighter variant);
- regrid: `preRegrid` cache + `regrid` redistribute onto the new `DisjointBoxLayout`;
- remap = **pool model**: collect movers, keep same-rank locally, scatter the rest;
  whole-patch transfer → zero-copy `MPI_Send(data())`, boundary-crossers → gather subset;
- the **leaf does NOT store its Box** — the `DisjointBoxLayout` is the single source of
  truth; geometry-dependent methods take `box`/`dx`/`probLo` as arguments.

Everything else decided — see `ARCHITECTURE.md` (LOCKED design decisions) and `BENCHMARKS.md`.
Next step is implementation, starting from the arena-SoA leaf (container-owned mandatory
columns + payload, `position(i)` promotion, cell-sort/CSR), then the AMR container layer.
