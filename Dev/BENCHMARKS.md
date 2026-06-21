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

## Results (ns/particle; lower is better)

| Operation | List (fragmented) | vector\<P\> AoS | SoA | best |
|---|---:|---:|---:|---|
| (1) Deposition (scatter -> grid)    | 11.93 | 11.98 | 10.89 | ~tie |
| (2) Interpolation (gather <- grid)  |  8.22 |  7.06 |  6.96 | contiguity (vec≈SoA, 1.18x vs List) |
| (3) Streaming transform `w=½|x|²`   |  1.24*| 0.78  | 0.45 / **0.31**† | **SoA / per-component** |
| (4) Build (allocation)              |  4.43 | **0.65** | 1.25 | **vector\<P\>** |
| (5) Remap (scatter -> 64 buckets)   |  6.36 | **3.34** | 5.65 | **vector\<P\>** |
| (6) MPI packing (serialize 32 B/p)  |  1.23 | 0.59  | 0.57 | contiguity (vec≈SoA, ~2.1x vs List) |

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
2. **`ParticleSoA` uniquely wins only the SIMD streaming transform** (1.7x over
   `vector<P>`; 4x over List with per-component). That is the kernel class
   `ParticleLoops` targets.
3. **Deposition, interpolation, MPI packing**: contiguity is what matters; `vector<P>`
   and SoA tie, both ~1.2–2.1x over `List`.

## Recommendation

- **Highest ROI / lowest risk: replace `List<P>` with `std::vector<P>` (AoS) per
  patch.** It captures the allocation (6.8x), remap (1.9x), MPI-packing (2x) and
  interpolation (1.2x) wins, keeps `p.field()` ergonomics, and is far simpler than SoA.
- **Adopt full `ParticleSoA` only if** (a) per-particle field-update / SIMD streaming
  kernels are hot in the real pipeline (then SoA is 1.7–4x there, per-component best),
  or (b) **GPU offload** is on the roadmap (SoA is the standard device layout — then it
  is close to mandatory regardless of these CPU numbers).
- **Load imbalance** is orthogonal: a uniform per-rank speedup still helps the
  bottleneck rank; the lever for imbalance is load balancing, which `vector<P>`/SoA both
  make cheaper (vector more so via faster remap).

## Open follow-ups (not yet measured)

- Allocation under **OpenMP** — `List`'s per-node `new` contends across threads, so the
  build/remap gaps for `vector<P>`/SoA likely widen in the threaded regions (KMC,
  merge/split, regrid).
- **Profile a real ItoKMC run** for the particle-time fraction and per-op breakdown —
  the missing input that turns "which layout" into "how much wall-clock."
