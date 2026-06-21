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
- **SoA granularity: per-field** — one `std::vector<RealVect>` per vector column
  (not split into per-component `x[]/y[]/z[]`).
- **Schema: explicit traits, NOT a macro** — `ParticleTraits<P>` with a `columns`
  member-pointer tuple + `positionPtr`/`weightPtr` + optional `h5Columns`. Indices
  for position/weight/HDF5 are all **derived** from member pointers (no literal
  indices anywhere → reordering `columns` is safe). A `CD_PARTICLE_LAYOUT` macro was
  evaluated and rejected (opaque to tooling); the file was removed.
- **Storage granularity: per-patch** — `ParticleSoA<P>` = one patch; the AMR
  container holds `LayoutData<ParticleSoA<P>>` per level. See `ARCHITECTURE.md`.

## Benchmark findings (see `BENCHMARKS.md`)

`Dev/Benchmark/` compares `List<P>` vs `vector<P>` (AoS) vs `ParticleSoA<P>` on
deposition, interpolation, streaming transform, build, remap, and MPI packing.
Headline: **escaping the linked list is the big win, not the SoA column split** —
`vector<P>` ties or beats SoA everywhere except SIMD streaming transforms, where SoA
is 1.7x (per-component 4x) faster. So full SoA is justified mainly by hot SIMD
particle kernels and/or GPU offload; otherwise `vector<P>` is the sweet spot. This
reframes the "drop-in vs clean break" question below into "SoA vs vector<P> at all".

## Big open questions (remaining)

1. **Drop-in vs. clean break.** Keep `GenericParticle<M,N>` / `PointParticle` /
   `Photon` / `ItoParticle` as a thin AoS facade over SoA, or rewrite call sites?
   (Plus: do we need the tree compilable at every commit, i.e. incremental?)
2. **Cell-index structure for KMC** (now that storage is per-patch): CSR offsets vs
   per-cell index vectors vs hybrid for cut cells (see `ARCHITECTURE.md` follow-ons
   and §8 of `USAGE_PATTERNS.md`).
3. **Container layering details**: how `LayoutData<ParticleSoA<P>>` wraps Chombo's
   box/`DataIndex` machinery; whether buffer/mask/cache holders use a lighter SoA.

See `USAGE_PATTERNS.md` (per-feature requirements) and `ARCHITECTURE.md` (per-patch
layering).
