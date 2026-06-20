# Particle library redesign — working notes (Dev/)

## Goal

Redesign the chombo-discharge particle library to get rid of Chombo's `List<P>` and
replace the per-bin/per-patch storage with a **Struct-of-Arrays (SoA)** layout.

This is a large job. The work in `Dev/` is **design and prototyping only** — we do
**not** modify any production code under `Source/`, `Physics/`, or `Geometries/` yet.
The point of `Dev/` is to:

1. Capture how particles are *actually* used across the codebase (see
   `USAGE_PATTERNS.md`).
2. Stand up a small, standalone, compilable **proxy** of the new SoA particle core
   (`CD_ParticleSoA.H` + `test_ParticleSoA.cpp`) that we can iterate on and test
   against the real usage patterns, with no Chombo build required.
3. Decide the core methodology before touching production code.

## Hard constraints (discovered)

- **C++14**. Every `Make.defs.*` sets `CXXSTD = 14`. No `if constexpr`, no C++17
  fold expressions, no inline `static constexpr` data members, no C++20 NTTP string
  literals. The prototype is written to compile under `-std=c++14 -Wall -Wextra`.
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
- [x] Prototype SoA container `CD_ParticleSoA.H` (traits-driven, member-pointer
      columns; derived indices; member-pointer field selector; MPI + HDF5-subset
      linearization; append/remove/gather/scatter/get/column).
- [x] `test_ParticleSoA.cpp` — basic ops, deposition loop, interpolation in-place
      update, MPI linearization, field selector, reorder safety, non-mandatory field
      update, HDF5 subset. `test_ParticleLayoutMacro.cpp` — macro vs explicit traits.
      Both build + pass under **C++17** (`make check`).

## Decisions (locked)

- **C++ standard: bump `CXXSTD` to 17** — enables the single-token member-pointer
  field selector.
- **Field selector: single-token member pointer, type deduced** —
  `soa.column<&P::velocity>()`, `soa.get<&P::weight>(i)`, `deposit<&P::weight>(...)`.
  Index derived from the member pointer (reorder-safe); return type deduced (no `Ret`
  arg like today).
- **SoA granularity: per-field** — one `std::vector<RealVect>` per vector column
  (not split into per-component `x[]/y[]/z[]`).
- **Schema: explicit traits, NOT the macro** — `ParticleTraits<P>` with a `columns`
  member-pointer tuple + `positionPtr`/`weightPtr` + optional `h5Columns`. Indices
  for position/weight/HDF5 are all **derived** from member pointers (no literal
  indices anywhere → reordering `columns` is safe). `CD_PARTICLE_LAYOUT` exists in
  `CD_ParticleLayoutMacro.H` for reference but is not the chosen path.
- **Storage granularity: per-patch** — `ParticleSoA<P>` = one patch; the AMR
  container holds `LayoutData<ParticleSoA<P>>` per level. See `ARCHITECTURE.md`.

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
