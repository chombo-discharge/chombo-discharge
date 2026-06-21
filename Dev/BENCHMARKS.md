# Particle storage benchmarks — results & decision record

Decision-support data comparing three per-patch particle storages for chombo-discharge:

- **`List<P>`** — Chombo's linked list of Array-of-Structs particles (today's layout,
  inside `BinFab`). Benchmarked both *compact* (freshly built, contiguous nodes — best
  case) and *fragmented* (heap churned between node allocations — the realistic case).
- **`vector<P>`** — a plain contiguous Array-of-Structs `std::vector<BenchParticle>`.
- **`ParticleSoA<P>`** — the Struct-of-Arrays container (`Dev/CD_ParticleSoA.H`), one
  contiguous column per field. Also a *per-component* variant (`x[]`,`y[]`,`z[]` as
  separate `Real` arrays instead of one interleaved `vector<RealVect>`).

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
  Conclusion: **SoA packing is not inherently slow — it just needs an NT-aware pack**
  (or a tuned `memcpy`); with that, SoA packs at `vector<P>` speed. (`vector<P>`'s plain
  bulk copy gets there for free, which is still a point in its favour.)
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
   column split.** `vector<P>` matches or *beats* `ParticleSoA` on 5 of 6 operations,
   including the allocation-heavy ones (build 6.8x and beats SoA; remap 1.9x and beats
   SoA — SoA loses there to multi-column reallocation churn).
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
4. **MPI packing**: `vector<P>` (single bulk `memcpy`) is the robust winner for free.
   SoA's per-column `memcpy` *tied* `vector<P>` at 64K but **regressed to slower-than
   `List` (0.55x) at 1M** because `glibc` did not take the non-temporal store path —
   **confirmed**: forcing NT stores makes SoA packing 4.1x faster, back to `vector<P>`
   parity. So SoA can pack fast, but only with an NT-aware pack; `vector<P>` gets it for
   free.
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
