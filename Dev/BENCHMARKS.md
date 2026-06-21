# Particle storage benchmarks — results & decision record

Decision-support data comparing three per-patch particle storages for chombo-discharge:

- **`List<P>`** — Chombo's linked list of Array-of-Structs particles (today's layout,
  inside `BinFab`). Benchmarked both *compact* (freshly built, contiguous nodes — best
  case) and *fragmented* (heap churned between node allocations — the realistic case).
- **`vector<P>`** — a plain contiguous Array-of-Structs `std::vector<BenchParticle>`.
- **`ParticleSoA<P>`** — the Struct-of-Arrays container (`Dev/CD_ParticleSoA.H`), one
  contiguous column per field. Also a *per-component* variant (`x[]`,`y[]`,`z[]` as
  separate `Real` arrays instead of one interleaved `vector<RealVect>`).

## Re-validation on the merged ParticleSoA (2026-06-21)

The two prototypes (`ParticleSoA` + `ParticleSoAArena`) were merged into one arena-backed
`ParticleSoA` with **container-owned mandatory columns** (per-component `ParticleReal`
position, weight, `particleID`, `rankID`) + a user payload. The benchmark was ported to it
(SoA uses `ParticleSoA<>` / `<MoverPayload>` / `<MergePayload>`; the AoS baselines carry a
matching `particleID`/`rankID` so bytes/particle are comparable — 48 B AoS vs 44 B SoA
columns). **Every earlier decision reproduced** (64K particles, 16^3+2, GCC `-O3 -march=native`):

| Section | Result | vs `List` | vs `vector<P>` | Verdict |
|---|---|---|---|---|
| (1) deposition | SoA 11.4 ns/p | 1.06x | 1.05x | tie |
| (2) interpolation | SoA 7.9 ns/p | 1.03x | 0.89x | tie (vector slightly ahead) |
| (3) transform | **SoA 0.28 ns/p** | **4.26x** | **2.30x** | **SoA wins (per-component AVX)** |
| (4) build | SoA-reserve 1.46 ns/p | 2.17x | 0.92x | needs reserve to match vector |
| (5) remap | SoA-2pass 4.72 ns/p | 1.07x | 0.93x | needs two-pass to match vector |
| (6) MPI pack | arena 0.77 ns/p | 1.50x | ~vector | arena one-memcpy robust |
| (7) vec interp | SoA 16.5 ns/p | 1.08x | 1.05x | tie |
| (8) Euler advance | **SoA 0.45 ns/p** | **4.88x** | **4.83x** | **SoA wins (per-component v)** |
| (9) KD merge | SoA 177 ns/p | 1.65x | 1.03x | tie with vector |
| (10) container copy | arena 0.78 ns/p | 8.2x | 1.17x | arena one-memcpy wins |

Notes from the port:
- **Per-component position helps the streaming kernels and hurts nothing measurable**, but
  random-access whole-particle work (KD merge) must read position via the per-component
  columns (`positionColumn(d)[i]`), NOT `position(i)` — the latter builds a promoted
  `RealVect` by value and, called per comparison in `nth_element`, was ~2x slower. The
  production merge path operates on cell-sorted contiguous ranges, where this is a non-issue.
- The advance win grew (4.9x vs the earlier ~4x) because velocity is now ALSO per-component
  (`MoverPayload{vx,vy,vz}`), so position and velocity updates are all unit-stride.
- The arena removes the bimodal per-column-`memcpy` pack variance: `data()` is one contiguous
  span, so pack/copy are a single `memcpy` (and a real whole-patch MPI send skips even that).
- The separate "N-vector vs arena" comparison rows are gone: the merged container IS the
  arena, so there is one SoA path per section.

The full corrected table for the current run lives in `Dev/Benchmark/pout.0`. The detailed
per-op analysis below is from the original prototype runs and still applies.

## Setup

- Source: `Dev/Benchmark/main.cpp`, `Dev/Benchmark/GNUmakefile`.
- Build/run (3D, optimized, matches the prebuilt OPTHIGH 3d libraries):
  ```
  cd Dev/Benchmark
  make DIM=3 OPT=HIGH DEBUG=FALSE MPI=TRUE CXXSTD=17
  mpirun -np 1 ./main3d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex   # results in pout.0
  ```
  (Build serially — concurrent `make` invocations race while generating `.d`
  dependency files and corrupt them. `CXXSTD=17` is a one-off override; the tree is
  still C++14 — see `TODO.md`.)
- Grid: all-regular `FArrayBox`, 16 valid + 2 ghost cells/coord (20^3 allocated),
  single component. Particles: 65,536 (16^3 valid cells x 16 per cell).
- Every storage carries identical data (position + weight) and calls the **identical**
  per-particle kernel, so differences are storage layout only. Deposition verified
  mass-conserving and bitwise-identical across storages for the same particle order.
- Hardware/compiler: GCC, `-O3 -march=native` (AVX, 4-wide double). Numbers are
  averages over 200 reps (deposition/interp/transform/packing) or 50 reps
  (build/remap); single-threaded. Absolute ns/particle will vary by machine — read
  the **ratios**, not the absolute values.

## Results — baseline grid (65,536 particles; 16^3 valid + 2 ghost; ns/particle, lower is better)

| Operation | List (fragmented) | vector\<P\> AoS | SoA | best |
|---|---:|---:|---:|---|
| (1) Deposition (scatter -> grid)    | 11.93 | 11.98 | 10.89 | ~tie |
| (2) Interpolation (gather <- grid)  |  8.22 |  7.06 |  6.96 | contiguity (vec≈SoA, 1.18x vs List) |
| (3) Streaming transform `w=½|x|²`   |  1.24*| 0.78  | 0.45 / **0.31**† | **SoA / per-component** |
| (4) Build (allocation)              |  4.43 | **0.65** | 1.25 | **vector\<P\>** |
| (5) Remap (scatter -> 64 buckets)   |  6.36 | **3.34** | 5.65 | **vector\<P\>** |
| (6) MPI packing (serialize 32 B/p)  |  1.23 | 0.59  | 0.57 | contiguity (vec≈SoA, ~2.1x vs List) |
| (7) Vector-field interpolation (3-comp) | 17.88 | 17.94 | 17.07 | ~tie (all within 5%) |
| (8) Euler advance `x += v·dt`       |  2.23 | 1.98  | 2.20 / **0.46**‡ | **SoA per-component only** |

‡ Euler advance: SoA `vector<RealVect>` = 2.20 (no better than List!), SoA *per-component* = 0.46.

\* transform List figure is the *compact* list (best case). † SoA `vector<RealVect>`
= 0.45, SoA *per-component* = 0.31.

### Speedups vs today's `List<P>`

| Operation | vector\<P\> | SoA |
|---|---:|---:|
| Build       | **6.8x** | 3.5x |
| Remap       | **1.9x** | 1.1x |
| MPI packing | 2.1x | 2.1x |
| Interpolation | 1.2x | 1.2x |
| Deposition  | ~1.0x | ~1.1x |
| Transform   | 1.6x | 2.8x (per-component **4.0x**) |

## Results — large grid (1,048,576 particles; 32^3 valid + 4 ghost; 40^3 allocated)

Same benchmark, config set via the `constexpr`s at the top of `main()`
(`nCell=32, nGhost=4, ppc=32, fastReps=50, buildReps=20`). The working sets now exceed
cache (single-comp field 512 KB, 3-comp field 1.5 MB, columns ~24–56 MB), so several
results shift relative to the cache-resident baseline.

| Operation | List (fragmented) | vector\<P\> AoS | SoA | SoA vs List |
|---|---:|---:|---:|---:|
| (1) Deposition (sorted)            | 12.48 | 12.14 | 10.95 | 1.14x |
| (2) Interpolation (sorted)         |  8.91 |  7.00 |  6.94 | 1.28x (vs vector\<P\> 1.01x) |
| (3) Streaming transform `w=½|x|²`  |  2.88*| 1.46  | 0.88 / **0.51**† | 3.3x (per-comp **5.6x**) |
| (4) Build (allocation)             |  6.33 | **2.10** | 2.65 | vector\<P\> 3.0x; SoA 2.4x |
| (5) Remap (-> 64 buckets)          | 10.57 | **4.44** | 5.40 | vector\<P\> 2.4x; SoA 2.0x |
| (6) MPI packing                    |  3.22 | **1.47** | 5.55 | **SoA 0.58x — slower than List!** |
| (7) Vector-field interp (3-comp)   | 17.92 | 17.86 | 17.08 | 1.05x (vs vector\<P\> 1.00x) |
| (8) Euler advance `x += v·dt`      |  4.78 | 3.10  | 2.40 / **1.20**‡ | SoA(RealVect) 2.0x; per-comp 4.0x |

\* compact List. † SoA `vector<RealVect>`=0.88, per-component=0.51. ‡ SoA
`vector<RealVect>`=2.40, per-component=1.20.

### What changes at scale (vs the cache-resident baseline)

- **Memory-bound kernels pull further ahead of `List`.** Once data exceeds cache the
  cost is DRAM bandwidth, so contiguity matters more. The clearest case is the **Euler
  advance**: SoA `vector<RealVect>` was *no better than List* at 64K (1.0x, cache-bound)
  but is **2.0x at 1M** — contiguous columns beat pointer-chased nodes on bandwidth even
  though the kernel stays scalar. Per-component (SIMD) widens to 4.0x; the streaming
  transform widens to 5.6x (per-component).
- **Deposition and vector-field interpolation stay layout-insensitive** (1.05–1.14x) —
  still scatter- / grid-gather-bound, exactly as at 64K. The larger grid does not change
  this.
- **Build/remap**: `vector<P>` still wins and still beats SoA; the build ratio compressed
  (6.8x -> 3.0x) as everything becomes allocation/bandwidth-bound, while SoA's remap
  improved (1.1x -> 2.0x).
- **MPI packing reverses — and this is the surprise.** At 64K, SoA packing tied
  `vector<P>` (~2.1x over List). At 1M, **SoA packing is 0.58x — slower than `List`**,
  and 3.8x slower than `vector<P>`. Both do bulk `memcpy` of the same total bytes, so the
  most likely cause is the glibc `memcpy` **non-temporal-store threshold**: the single
  ~32 MB AoS copy (`vector<P>`) gets streaming NT stores (no read-for-ownership), while
  SoA's split per-column copies (~24 MB + ~8 MB) fall in a regime that uses ordinary
  stores and pays read-for-ownership traffic.

  **Confirmed (non-temporal stores).** Re-running the packing with explicit AVX
  non-temporal (streaming) stores for the SoA columns closes the gap completely:

  | pack variant (1M particles) | ns/p | vs List |
  |---|---:|---:|
  | List (per-particle)          | 3.32 | 1.0x |
  | vector\<P\> (1 bulk memcpy)  | 1.42 | 2.34x |
  | vector\<P\> (split 2 memcpy) | 1.66 | 2.0x |
  | SoA (per-column memcpy)      | 5.99 | **0.55x** |
  | **SoA (per-column, NT stores)** | **1.45** | **2.29x** |

  Forcing NT stores makes SoA packing **4.1x faster** and brings it to parity with
  `vector<P>`. The control — splitting the *contiguous AoS* copy into the same two sizes
  (24 MB + 8 MB) — stays fast (1.66), so the slowdown is **not** "two `memcpy` calls": it
  is that `glibc`'s `memcpy` did not take the non-temporal path for the SoA per-column
  copies, so they paid RFO write traffic. Explicit NT stores recover full bandwidth.

  **Reproducibility — characterized over many launches (decisive).** The per-column SoA
  pack is **bimodal**: across 8+ process launches it is ~1.4–1.5 ns/p (fast, ~2.3x) on
  ~60% of launches and ~5.5–6.0 ns/p (slow, ~0.55x) on ~40% — it really does hit the
  slow path, just not every time (earlier "always slow" / "always fast" impressions were
  small samples of this distribution). `vector<P>`, SoA-NT, and the arena pack are
  **stable** at ~1.4 ns/p every launch.

  **Assembly + mechanism.** `packSoA` disassembles to exactly two `call memcpy@plt` (the
  `R_X86_64_PLT32 memcpy` relocations confirm it) — our code just calls `glibc memcpy`;
  all the variance is inside libc. This machine's L3 is **32 MB**, and the copies are
  24 MB (pos) + 8 MB (weight) with a 32 MB destination — i.e. **right at glibc's
  non-temporal-store threshold** (a fraction of L3). At that size the path choice for the
  *separate-source* per-column copies flips with per-launch address/alignment (ASLR/
  allocator), hence bimodal. A single 32 MB contiguous copy (arena / `vector<P>`)
  reliably crosses the threshold → always NT → stable. So **the arena (or explicit NT)
  is a robustness fix, not just an optimization** — it removes a ~40%-of-launches 4x
  packing/copy slowdown. `vector<P>` is stable for free.
## Streaming transform, in detail (the one place SoA beats vector\<P\>)

| Layout | ns/p | vs List | deinterleave shuffles in the vectorized loop |
|---|---:|---:|---|
| List compact            | 1.24 | 1.0x | — (pointer-chased, not vectorizable) |
| vector\<P\> AoS         | 0.78 | 1.6x | **8 `vunpcklpd`** (every field strided in the 32 B struct) |
| SoA `vector<RealVect>`  | 0.45 | 2.8x | **2 `vunpcklpd`** (only position interleaved) |
| SoA per-component       | 0.31 | 4.0x | **0** — clean 4-wide AVX (`vmulpd`/`vfmadd231pd` on `ymm`) |

Shuffle count tracks performance exactly. SoA is **1.74x** faster than `vector<P>` here,
and per-component a further **1.43x** — this is the only operation where the column
split pays off beyond plain contiguity.

## Euler advance `x += v·dt`, in detail (the canonical particle update)

This reads *two* vector fields (x, v) and writes *one* (x) — the common case (velocity
advance, position advance). The result is the most important refinement of the whole
suite:

| Layout | ns/p | vs List | vectorized? (assembly) |
|---|---:|---:|---|
| List                   | 2.23 | 1.0x | no — scalar, pointer-chased |
| vector\<P\> AoS        | 1.98 | 1.1x | partially (shuffle-heavy) |
| **SoA `vector<RealVect>`** | 2.20 | **1.0x** | **NO — fully scalar** (`vmulsd`/`vfmadd*sd`, zero packed) |
| **SoA per-component**  | 0.46 | **4.8x** | yes — packed FMA |

**Interleaved `vector<RealVect>` SoA gives essentially zero benefit over `List` for the
canonical particle update.** `advanceSoA` compiles to entirely scalar FP: `pos[i] +=
vel[i]*dt` goes through `RealVect`'s user-defined `operator+=`/`operator*` on interleaved
storage, and the compiler does not vectorize it across particles. Only the
**per-component** layout (plain `Real` arrays) vectorizes — and then it is **4.8x**.

Contrast with the `½|x|²` transform (which writes a *scalar*): there the `RealVect`
column did vectorize (with shuffles). So the SoA-over-`vector<P>` advantage is **highly
workload-dependent**, and for genuine vector-in/vector-out updates it materializes
**only with per-component columns**, not with the per-field `vector<RealVect>` default.

## Vector-field interpolation, in detail (3-component FArrayBox -> particle RealVect)

A 3-component `FArrayBox` (component-major: each component is a contiguous plane, i.e.
SoA-for-the-grid) interpolated to a particle `RealVect` via CIC. Per particle this is a
2x2x2 stencil x 3 components = **24 grid reads** + 24 FMA + 8 CIC weights ("side A"),
plus one position read and one `RealVect` write ("side B").

| Layout | ns/p | vs List | vs vector\<P\> |
|---|---:|---:|---:|
| List      | 17.88 | 1.0x  | — |
| vector\<P\> | 17.94 | 1.0x | 1.0x |
| SoA       | 17.07 | 1.05x | 1.00x |

**All three layouts tie (within 5%).** This confirms the static prediction: side A (the
gather) operates on the `FArrayBox`, which is passed **`const` and is identical for every
particle storage**, and it dominates. The 3-component field is `20^3 x 3 x 8 B ≈ 192 KB`
→ L2-resident, so the gather is fast but still the bulk of the work, and **independent of
the particle container**. Side B (the particle read/write) is a small fraction, so the
storage layout barely matters.

Note the absolute cost (~18 ns/p) is ~2.5x the scalar interpolation (~7 ns/p) — exactly
the 3x grid-gather growth (8 -> 24 reads) dominating, while the SoA *relative* advantage
**shrinks** (1.05x here vs 1.19x for scalar interpolation). Interpolation is a
data-dependent gather and does not vectorize across particles, so per-component storage
would not help either. **SoA gives essentially nothing for vector-field interpolation —
the benefit is hidden behind the const, cache-resident grid gather.**

## No-intrinsics SoA variants for build / remap / copy / pack (1M particles)

The per-column SoA loses to `vector<P>` on bulk operations because it does N separate
operations (N allocations, N memcpies) where `vector<P>` does one. Three no-intrinsics
mitigations were tested at 1M:

| Operation | List | vector\<P\> | SoA (naive) | SoA fix | fix vs vector\<P\> |
|---|---:|---:|---:|---:|---:|
| Remap (-> 64 buckets) | 10.79 | 4.40 | 5.45 | **5.33** (two-pass reserve) | 0.82x |
| MPI packing           |  3.29 | 1.35 | 5.97 | **1.35** (arena, 1 memcpy) | 1.00x |
| Container copy        | 13.11 | 2.27 | 7.65 | **2.66** (arena, 1 memcpy) | 0.86x |

(ns/particle. "arena" = all columns in one contiguous allocation, so the bulk op is a
single `std::memcpy` that `glibc` auto-streams; no intrinsics.)

**Findings:**
- **Packing and copy are fully fixed by an *arena* layout** (single backing allocation):
  one `memcpy` → `glibc` auto-NT → `vector<P>` parity (1.00x / 0.86x), **no intrinsics**.
  This is the answer to "fast SoA pack/copy without intrinsics": store the columns in
  one buffer so the bulk op is one big copy. (Equivalent to the explicit-NT result, but
  with plain `memcpy`.) For a *full-container* send the arena even enables zero-copy
  (the arena *is* the wire layout).
- **Remap is *not* fixed by two-pass `reserve`** (5.45 -> 5.33, ~2%). The cost is not
  reallocation churn — it is the inherent **per-particle scatter into N columns** (N
  stores vs one struct store). So SoA scatter-remap stays ~1.2x slower than `vector<P>`.
  **But it is still 2.0x FASTER than today's `List`** — so adopting SoA is a remap
  *improvement* over the status quo, not a regression; it is only slower than the
  hypothetical `vector<P>` alternative.

## Arena-backed SoA (single allocation) — end-to-end vs N-vector SoA (1M particles)

`CD_ParticleSoAArena.H` stores all columns in ONE 64-byte-aligned allocation (column
base pointers cached). Measured against the N-vector `ParticleSoA` and `vector<P>`:

| Operation | vector\<P\> | SoA (N-vector) | **SoA arena** | note |
|---|---:|---:|---:|---|
| Build (reserve)        | 2.11 | 2.60 | **2.29** | arena ≈ one allocation, beats N-vector |
| Build (NO reserve)     |  —   |  —   | **16.2** | **7x worse** — growth reallocs whole buffer |
| Deposition (sorted)    | 12.2 | 10.8 | 12.2 | scatter-bound; offset access is free (≈ noise) |
| Transform `½|x|²`      | 1.18 | 0.93 | **0.72** | arena 1.29x over N-vec SoA — 64B-aligned SIMD |
| Remap (two-pass)       | 4.37 | 5.23 | 6.15 | arena slightly *worse*; doesn't fix scatter |
| MPI pack               | 1.35 | 1.5–6.0† | **1.36** | arena reliably = vector\<P\>, no intrinsics |

(ns/particle.) † per-column SoA pack measured ~1.5 (2.1x) in this session but ~6.0
(0.55x) in an earlier session — see below.

**Findings:**
- **Build:** arena with `reserve` (2.29) is one allocation and beats N-vector SoA
  (2.60), near `vector<P>` (2.11). **Without `reserve` it is catastrophic (16.2, 7x
  worse)** — every doubling reallocates and moves the *whole* buffer. So **arena makes
  `reserve` mandatory**; with it, build is competitive.
- **Element access is free:** arena deposition/transform read through cached column
  pointers and match (transform: beat) the N-vector version — the offset indirection
  costs nothing measurable.
- **Transform is *faster* with arena (1.29x):** 64-byte-aligned columns let the
  compiler use aligned SIMD (no shuffle/peel that the 16-byte `std::vector` data forces).
  A free bonus of the arena layout (still below per-component's 5.6x).
- **Pack/copy: arena is reliably at `vector<P>` parity, no intrinsics** — one `memcpy`
  of the contiguous arena (or zero-copy `MPI_Send` of `data()` for a whole-container
  transfer). Importantly, the per-column SoA pack is **not reliable**: it was fast
  (2.1x) in this session but slow (0.55x, the "reversal") in an earlier one, because
  glibc's NT decision for the split copy is borderline at this size. **The arena removes
  that variance.**
- **Remap is not helped** (slightly worse, 6.15 vs 5.23): scatter-into-buckets is
  per-particle and the arena's `append` indirection adds a touch of overhead. Remap
  stays the SoA soft spot (~1.2–1.4x vs `vector<P>`, still ~2x faster than `List`).

**Net:** an arena-backed SoA gives one-allocation builds, *faster* aligned-SIMD
transforms, and *robust* no-intrinsics pack/copy — at the cost of a hard `reserve`
requirement (growth is 7x) and no remap improvement. It directly addresses the
copy/pack-without-intrinsics goal; it does not change the remap trade.

## KD-tree merge, in detail (whole-particle access — the AoS-favouring case?)

Superparticle-style merge: 131,072 **7-field** particles (3 RealVect + 4 Real, ≈
`ItoParticle`; 32/cell on the 16^3 grid). A KD-tree recursively median-splits the
particle **indices by position** (this reads positions only), then merges each
nearest-pair leaf into one particle conserving total weight + center of mass (this
**gathers whole particles**). Verified weight- and COM-conserving.

| Storage | ns/p | vs List | vs vector\<P\> |
|---|---:|---:|---:|
| List (copy to vector + KD) | 275 | 1.0x | — |
| vector\<P\> AoS            | 184 | 1.49x | — |
| SoA                        | 179 | 1.54x | **1.03x** |

**SoA ties `vector<P>` (1.03x) — it does *not* lose, contrary to the "whole-particle
access favours AoS" expectation.** Reason: the cost is dominated by the O(N log N)
median-partitioning (`nth_element`), which reads **position only** → SoA streams just
the 3 MB position column, while `vector<P>` drags the whole 104-byte particle through
cache to use 24 bytes of it. That position-only phase *favours* SoA and outweighs the
leaf-merge phase (the O(N) whole-particle `gather`, which does favour AoS). With fine
pairwise merging (leaf=2) the partition dominates ~17:1, so SoA stays competitive.

Caveat: this is specific to a **sort/partition-dominated** merge. A merge that is pure
whole-particle random access (no position-only sort — e.g. fixed-bucket gather/merge)
would shift toward AoS, since SoA `gather` touches all 7 columns (7 cache lines) vs one
struct. But for the canonical KD merge, **AoS has no advantage**.

## Why deposition/interpolation don't benefit (assembly)

`depositSoA`/`depositList` compile to **entirely scalar** FP (`vmulsd`, `vfmadd*sd`,
`vroundsd` for `floor`, `vcvttsd2si`); **no** `vmulpd`/`vgatherqpd`/`vscatter`. The 2^D
corner updates write `rho(iv,0) += …` to **data-dependent, potentially-aliasing
addresses** (two corners or two particles can hit the same cell → a real WAW
dependency). The compiler cannot vectorize that scatter, and `CD_PRAGMA_SIMD`/`ivdep`
cannot override a *true* dependency. So SIMD is genuinely inapplicable, and the kernel
is bound by the per-particle compute + grid scatter — which the particle storage layout
does not change. (Fragmenting the list cost only ~2%: the CPU prefetches `next` and
hides the latency behind the per-particle work.)

Deposition could only be sped up by **cell-sorting** (write locality; ~5–10%, observed)
or an **algorithmic change** (per-cell binning/accumulation → turn the scatter into a
local reduction) — not by layout or SIMD.

## Key conclusions

1. **The big, cheap win is escaping the linked list (`List` → contiguous), not the SoA
   column split.** `vector<P>` matches or *beats* `ParticleSoA` on the non-SIMD
   operations. The only place `vector<P>` *structurally* beats SoA is the
   **allocation-churn** ops (build 6.8x and beats SoA; remap 1.9x and beats SoA — SoA
   pays for N column buffers instead of one).
1b. **The whole-particle-access concern did not materialize for KD merge.** SoA *ties*
   `vector<P>` (1.03x) on the KD-tree merge, because the dominant cost is the
   position-only median-partition (which favours SoA's single column) — not the
   whole-particle leaf gather. AoS only wins whole-particle access if the algorithm is
   gather-dominated *without* a position-only sort.
2. **`ParticleSoA` uniquely wins only on SIMD per-particle field kernels — and only
   with per-component columns.** For the canonical Euler advance `x += v·dt`, the
   per-field `vector<RealVect>` SoA is **no faster than `List`** (it stays scalar);
   only the **per-component** layout vectorizes (4.8x). The `½|x|²` transform (scalar
   output) is the friendlier case where `vector<RealVect>` already gets ~2.8x. So the
   SIMD win is real but **conditional on per-component storage** for vector-field
   updates, which is a heavier design (3 columns per vector field, RealVect view
   rebuilt on gather, more complex linearization).
3. **Deposition, interpolation**: grid-bound and layout-insensitive (`vector<P>` ≈ SoA,
   both ~1.05–1.3x over `List`), at both 64K and 1M.
4. **MPI packing / container copy**: `vector<P>` is fast and **stable** for free. SoA's
   per-column `memcpy` is **bimodal** — confirmed over many launches it hits a 4x slow
   path ~40% of the time (glibc's NT threshold is borderline at these ~L3-sized copies;
   the assembly is just two `memcpy` calls, so the variance is libc-internal). **Two
   fixes, both stable at `vector<P>` parity:** explicit NT stores (intrinsics), OR — no
   intrinsics — an **arena layout** (one allocation → one big `memcpy` that glibc
   reliably streams). So for SoA, an arena (or NT pack) is a *robustness* requirement,
   not just an optimization.
4b. **Remap is the one residual SoA cost.** Two-pass `reserve` did **not** help (~2%):
   the cost is the inherent per-particle scatter into N columns, not reallocation. SoA
   scatter-remap stays ~1.2x slower than `vector<P>` — **but 2.0x FASTER than today's
   `List`.** So SoA is a remap *improvement* over the status quo; it is only slower than
   the `vector<P>` alternative.
5. **Scale amplifies the contiguity win.** When data exceeds cache (1M particles), the
   memory-bound kernels (transform, Euler advance, scalar interpolation) pull *further*
   ahead of `List`, and SoA `vector<RealVect>` starts beating `List` on the Euler advance
   (1.0x → 2.0x) purely via bandwidth. Real per-patch particle counts in dense
   avalanches are large, so the scale-amplified numbers are the more representative ones.

## Recommendation

- **Highest ROI / lowest risk: replace `List<P>` with `std::vector<P>` (AoS) per
  patch.** It captures the allocation (6.8x), remap (1.9x), MPI-packing (2x) and
  interpolation (1.2x) wins, keeps `p.field()` ergonomics, and is far simpler than SoA.
- **Adopt full `ParticleSoA` only if** (a) per-particle field-update / SIMD streaming
  kernels are hot in the real pipeline AND you go **per-component** for the vector
  fields (per-field `vector<RealVect>` gives ~nothing for `x += v·dt` — see §Euler
  advance), or (b) **GPU offload** is on the roadmap (SoA, naturally per-component, is
  the standard device layout — then it is close to mandatory regardless of these CPU
  numbers). Note per-component is the heavier variant (more columns, RealVect view
  reconstructed on gather, more involved linearization).
- **If you choose SoA, use a single-allocation "arena" backing** (`CD_ParticleSoAArena.H`,
  columns as offset slices of one buffer). **Measured** (1M, end-to-end): build = one
  allocation (beats N-vector SoA); pack/copy reliably at `vector<P>` parity with **no
  intrinsics** (one `memcpy`, or zero-copy `MPI_Send` of `data()` for whole-container
  sends); and transforms are **1.29x faster** than N-vector SoA via 64-byte-aligned SIMD.
  Costs: growth without `reserve` is **7x catastrophic** (whole-buffer realloc) → reserve
  to capacity is mandatory; and remap is *not* improved (still ~1.2–1.4x vs `vector<P>`,
  ~2x faster than `List`).
- **Load imbalance** is orthogonal: a uniform per-rank speedup still helps the
  bottleneck rank; the lever for imbalance is load balancing, which `vector<P>`/SoA both
  make cheaper (vector more so via faster remap).

## Open follow-ups (not yet measured)

- Allocation under **OpenMP** — `List`'s per-node `new` contends across threads, so the
  build/remap gaps for `vector<P>`/SoA likely widen in the threaded regions (KMC,
  merge/split, regrid).
- **Profile a real ItoKMC run** for the particle-time fraction and per-op breakdown —
  the missing input that turns "which layout" into "how much wall-clock."
- ~~Confirm the MPI-packing reversal~~ **DONE**: it is a non-temporal-store effect;
  explicit NT stores make SoA packing 4.1x faster (parity with `vector<P>`). See the MPI
  packing detail section.
