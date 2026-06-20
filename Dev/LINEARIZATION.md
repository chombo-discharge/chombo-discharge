# Linearization: can it be SIMPLE and GENERIC?

**Short answer: yes — one generic implementation covers every particle type, and
the prototype already does it. You only need type-specific code if (a) particles
grow variable-length/runtime fields again, or (b) you want a *subset* of columns
for HDF5 — and even (b) stays generic with column tags.**

## Why generic works here

The `ParticleTraits<P>` column descriptor (a tuple of pointer-to-member) hands the
compiler the *complete, ordered list of fields*. Every field is required to be
**trivially copyable** (enforced by `static_assert`). So packing one particle is
just "walk the columns, `memcpy` each", and the per-particle byte size is a
compile-time `sumSizes(sizeof(col)...)`. No per-type linearizer:

```cpp
// generic, in CD_ParticleSoA.H — works for PointParticle, Photon, ItoParticle, MyParticle...
soa.linearizeParticle(buffer, i);     // pack column[i] of every column
soa.delinearizeAndAppend(buffer);     // unpack one particle, push onto every column
static constexpr size_t s_bytesPerParticle = /* compile-time sum of sizeof(columns) */;
```

The prototype test ships an `ItoLikeParticle` (104 bytes) through a byte buffer and
round-trips it exactly, with zero type-specific code.

### Bonus: SoA makes *batch* linearization faster than today

Because each column is contiguous, you can `memcpy` a whole run of particles per
column in one shot (column-major buffer) instead of field-by-field per particle.
That is strictly more cache-friendly than the current per-`List`-element `linearOut`.
Worth doing for the MPI remap path where we move many particles at once.

## When you'd need specialization (and how to avoid it)

1. **Variable-length / runtime-sized fields.** The old `ItoParticle` once carried
   runtime scalar buffers; that needed a dynamic `size()`. The *current* particles
   are all fixed-size, so generic fixed-size packing is valid. If runtime fields
   come back, generic-by-memcpy breaks and that type needs a custom pack/size.
   → Recommendation: keep particles fixed-size; forbid runtime buffers in the SoA core.

2. **HDF5 = a subset of columns** (Reals only, no `particleID`/`rankID`, often no
   scratch/`tmp` fields). This is *column selection*, not really per-type code.
   Keep it generic by tagging columns in the traits, e.g.

   ```cpp
   // index_sequence of columns to checkpoint to HDF5
   using H5Columns = std::index_sequence<0, 1, 2>;  // pos, weight, energy ...
   ```
   Then a single generic `h5LinearizeParticle` walks `H5Columns` instead of all
   columns. Same machinery, different index list — still one implementation.

   **Implemented in the prototype.** `ParticleSoA` exposes
   `h5BytesPerParticle()`, `h5LinearizeParticle()`, `h5DelinearizeAndAppend()`,
   parameterized by `Traits::H5Columns`. If the traits don't declare `H5Columns`
   it **defaults to all columns** (detected via SFINAE). On restart, columns not in
   the subset are left default-constructed (they get recomputed), exactly like
   `GenericParticle::H5linearIn`. Tested in `test_ParticleSoA.cpp::testH5Subset`:
   `ItoLikeParticle` ships 104 B for MPI but only 40 B (pos+weight+energy) to HDF5,
   and `velocity`/`mobility` come back default after a checkpoint round-trip.

3. **Endianness / heterogeneous clusters.** `memcpy` assumes homogeneous byte
   order. This matches Chombo's current behavior (it also blits raw bytes), so it's
   not a regression. Flag only if we ever target mixed-endian runs.

## Recommendation

- One **generic** `linearize/delinearize` over all columns for MPI.
- One **generic** variant parameterized by a compile-time **column index list** for
  HDF5 (and for any "trim these fields before sending" use case).
- No per-type linearizers. The only per-type input is the `ParticleTraits<P>`
  column list (which the user writes anyway), plus an optional `H5Columns` subset.
