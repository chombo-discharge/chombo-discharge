/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   main.cpp
  @brief  Benchmarks comparing Chombo List<P> with ParticleSoA<P>.
  @author Robert Marskar
  @details
  Micro-benchmarks on a single grid patch, all using the identical particle data and
  identical math, so the measured difference is the storage layout (linked-list
  Array-of-Structs vs contiguous Struct-of-Arrays). Three storages are compared:
  List<P>, vector<P> (AoS), and ParticleSoA<P> (with a per-component variant):

    1. CIC deposition   (particle weight -> FArrayBox; a scatter)
    2. CIC interpolation(FArrayBox -> particle scalar; a gather)
    3. streaming transform (weight = 0.5|x|^2; SoA path uses ParticleLoops SIMD)
    4. build            (construct the container; allocation cost)
    5. remap            (scatter particles into K destination containers)
    6. MPI packing      (serialize position + weight into a byte buffer)
    7. vector-field interpolation (3-component FArrayBox -> particle RealVect)
    8. Euler advance    (x += v*dt; the canonical particle update)
    9. KD-tree merge    (median-split by position, merge nearest pairs; whole-particle access)

  Grid: all-regular FArrayBox, 16 valid + 2 ghost cells per coordinate (20^3
  allocated), single component. Particles: 16 per valid cell (32/cell for the merge).
  All builds are 3D (make DIM=3 ...; the output prints SpaceDim=3).
*/

// Std includes
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <immintrin.h>
#include <numeric>
#include <random>
#include <vector>

// Chombo includes
#include <FArrayBox.H>
#include <List.H>
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_ParticleSoA.H>
#include <CD_ParticleLoops.H>

namespace ChomboDischarge {

  /** @brief Benchmark particle: position + weight. */
  struct BenchParticle
  {
    RealVect m_position = RealVect::Zero;
    Real     m_weight   = 0.0;

    inline const RealVect&
    position() const noexcept
    {
      return m_position;
    }
    inline Real
    weight() const noexcept
    {
      return m_weight;
    }
  };

  /** @brief SoA layout for BenchParticle. */
  template <>
  struct ParticleTraits<BenchParticle>
  {
    static constexpr auto columns     = std::make_tuple(&BenchParticle::m_position, &BenchParticle::m_weight);
    static constexpr auto positionPtr = &BenchParticle::m_position;
    static constexpr auto weightPtr   = &BenchParticle::m_weight;
  };

  /** @brief Particle carrying a velocity, for the Euler advance x += v*dt. */
  struct MoverParticle
  {
    RealVect m_position = RealVect::Zero;
    RealVect m_velocity = RealVect::Zero;
    Real     m_weight   = 0.0; // required by ParticleSoA (mandatory weight column); unused by the advance
  };

  /** @brief SoA layout for MoverParticle. */
  template <>
  struct ParticleTraits<MoverParticle>
  {
    static constexpr auto columns     = std::make_tuple(&MoverParticle::m_position,
                                                    &MoverParticle::m_velocity,
                                                    &MoverParticle::m_weight);
    static constexpr auto positionPtr = &MoverParticle::m_position;
    static constexpr auto weightPtr   = &MoverParticle::m_weight;
  };

  /** @brief A richer particle (~ItoParticle: 3 RealVect + 4 Real) for the merge benchmark. */
  struct MergeParticle
  {
    RealVect m_position  = RealVect::Zero;
    RealVect m_velocity  = RealVect::Zero;
    RealVect m_oldPos    = RealVect::Zero;
    Real     m_weight    = 0.0;
    Real     m_mobility  = 0.0;
    Real     m_diffusion = 0.0;
    Real     m_energy    = 0.0;
  };

  /** @brief SoA layout for MergeParticle (7 columns). */
  template <>
  struct ParticleTraits<MergeParticle>
  {
    static constexpr auto columns     = std::make_tuple(&MergeParticle::m_position,
                                                    &MergeParticle::m_velocity,
                                                    &MergeParticle::m_oldPos,
                                                    &MergeParticle::m_weight,
                                                    &MergeParticle::m_mobility,
                                                    &MergeParticle::m_diffusion,
                                                    &MergeParticle::m_energy);
    static constexpr auto positionPtr = &MergeParticle::m_position;
    static constexpr auto weightPtr   = &MergeParticle::m_weight;
  };

} // namespace ChomboDischarge

using namespace ChomboDischarge;

namespace {

  using SoA = ParticleSoA<BenchParticle>;

  volatile std::uint64_t g_sink = 0; // defeats dead-code elimination of build/remap

  // ---------------------------------------------------------------------------
  // Shared per-particle kernels.
  // ---------------------------------------------------------------------------

  /** @brief CIC deposition of one particle weight onto a single-component FArrayBox. */
  inline void
  depositOneCIC(FArrayBox&      a_rho,
                const RealVect& a_x,
                const Real      a_w,
                const RealVect& a_probLo,
                const RealVect& a_invDx,
                const Real      a_invVol) noexcept
  {
    IntVect  i0;
    RealVect f;
    for (int d = 0; d < SpaceDim; d++) {
      const Real g = (a_x[d] - a_probLo[d]) * a_invDx[d] - 0.5;
      const int  i = static_cast<int>(std::floor(g));
      i0[d]        = i;
      f[d]         = g - static_cast<Real>(i);
    }

    constexpr int nCorners = 1 << SpaceDim;
    for (int c = 0; c < nCorners; c++) {
      Real    wcic = 1.0;
      IntVect iv   = i0;
      for (int d = 0; d < SpaceDim; d++) {
        if (c & (1 << d)) {
          wcic *= f[d];
          iv[d] += 1;
        }
        else {
          wcic *= (1.0 - f[d]);
        }
      }
      a_rho(iv, 0) += a_w * a_invVol * wcic;
    }
  }

  /** @brief CIC interpolation of a single-component FArrayBox to a scalar at the particle. */
  inline Real
  interpolateOneCIC(const FArrayBox& a_rho,
                    const RealVect&  a_x,
                    const RealVect&  a_probLo,
                    const RealVect&  a_invDx) noexcept
  {
    IntVect  i0;
    RealVect f;
    for (int d = 0; d < SpaceDim; d++) {
      const Real g = (a_x[d] - a_probLo[d]) * a_invDx[d] - 0.5;
      const int  i = static_cast<int>(std::floor(g));
      i0[d]        = i;
      f[d]         = g - static_cast<Real>(i);
    }

    constexpr int nCorners = 1 << SpaceDim;
    Real          val      = 0.0;
    for (int c = 0; c < nCorners; c++) {
      Real    wcic = 1.0;
      IntVect iv   = i0;
      for (int d = 0; d < SpaceDim; d++) {
        if (c & (1 << d)) {
          wcic *= f[d];
          iv[d] += 1;
        }
        else {
          wcic *= (1.0 - f[d]);
        }
      }
      val += wcic * a_rho(iv, 0);
    }
    return val;
  }

  /** @brief Streaming per-particle update: weight = 0.5 * |position|^2 (reads vector, writes scalar). */
  inline Real
  transformOne(const RealVect& a_x) noexcept
  {
    return 0.5 * (D_TERM(a_x[0] * a_x[0], +a_x[1] * a_x[1], +a_x[2] * a_x[2]));
  }

  // ---------------------------------------------------------------------------
  // Per-storage loops.
  // ---------------------------------------------------------------------------

  void
  depositList(FArrayBox&                 a_rho,
              const List<BenchParticle>& a_p,
              const RealVect&            a_lo,
              const RealVect&            a_iDx,
              const Real                 a_iVol)
  {
    for (ListIterator<BenchParticle> lit(a_p); lit.ok(); ++lit) {
      depositOneCIC(a_rho, lit().position(), lit().weight(), a_lo, a_iDx, a_iVol);
    }
  }

  void
  depositSoA(FArrayBox& a_rho, const SoA& a_soa, const RealVect& a_lo, const RealVect& a_iDx, const Real a_iVol)
  {
    const std::vector<RealVect>& pos = a_soa.positions();
    const std::vector<Real>&     w   = a_soa.weights();
    const std::size_t            n   = a_soa.size();
    for (std::size_t i = 0; i < n; i++) {
      depositOneCIC(a_rho, pos[i], w[i], a_lo, a_iDx, a_iVol);
    }
  }

  void
  depositVector(FArrayBox&                        a_rho,
                const std::vector<BenchParticle>& a_v,
                const RealVect&                   a_lo,
                const RealVect&                   a_iDx,
                const Real                        a_iVol)
  {
    const std::size_t n = a_v.size();
    for (std::size_t i = 0; i < n; i++) {
      depositOneCIC(a_rho, a_v[i].m_position, a_v[i].m_weight, a_lo, a_iDx, a_iVol);
    }
  }

  void
  interpolateList(const FArrayBox& a_rho, List<BenchParticle>& a_p, const RealVect& a_lo, const RealVect& a_iDx)
  {
    for (ListIterator<BenchParticle> lit(a_p); lit.ok(); ++lit) {
      lit().m_weight = interpolateOneCIC(a_rho, lit().position(), a_lo, a_iDx);
    }
  }

  void
  interpolateVector(const FArrayBox&            a_rho,
                    std::vector<BenchParticle>& a_v,
                    const RealVect&             a_lo,
                    const RealVect&             a_iDx)
  {
    const std::size_t n = a_v.size();
    for (std::size_t i = 0; i < n; i++) {
      a_v[i].m_weight = interpolateOneCIC(a_rho, a_v[i].m_position, a_lo, a_iDx);
    }
  }

  void
  interpolateSoA(const FArrayBox& a_rho, SoA& a_soa, const RealVect& a_lo, const RealVect& a_iDx)
  {
    const std::vector<RealVect>& pos = a_soa.positions();
    std::vector<Real>&           w   = a_soa.weights();
    const std::size_t            n   = a_soa.size();
    for (std::size_t i = 0; i < n; i++) {
      w[i] = interpolateOneCIC(a_rho, pos[i], a_lo, a_iDx);
    }
  }

  /**
    @brief CIC interpolation of a SpaceDim-component (vector) FArrayBox to a RealVect at the particle.
    @details The FArrayBox is component-major (component varies slowest), so the SpaceDim values of a
    given cell are numPts() apart. This gather is identical for every particle storage (the FArrayBox
    is shared and const), so it is the layout-independent part of the cost.
  */
  inline RealVect
  interpolateVecFieldCIC(const FArrayBox& a_rho,
                         const RealVect&  a_x,
                         const RealVect&  a_probLo,
                         const RealVect&  a_invDx) noexcept
  {
    IntVect  i0;
    RealVect f;
    for (int d = 0; d < SpaceDim; d++) {
      const Real g = (a_x[d] - a_probLo[d]) * a_invDx[d] - 0.5;
      const int  i = static_cast<int>(std::floor(g));
      i0[d]        = i;
      f[d]         = g - static_cast<Real>(i);
    }

    constexpr int nCorners = 1 << SpaceDim;
    RealVect      val      = RealVect::Zero;
    for (int c = 0; c < nCorners; c++) {
      Real    wcic = 1.0;
      IntVect iv   = i0;
      for (int d = 0; d < SpaceDim; d++) {
        if (c & (1 << d)) {
          wcic *= f[d];
          iv[d] += 1;
        }
        else {
          wcic *= (1.0 - f[d]);
        }
      }
      for (int comp = 0; comp < SpaceDim; comp++) {
        val[comp] += wcic * a_rho(iv, comp);
      }
    }
    return val;
  }

  void
  interpVecList(const FArrayBox& a_rho, List<MoverParticle>& a_p, const RealVect& a_lo, const RealVect& a_iDx)
  {
    for (ListIterator<MoverParticle> lit(a_p); lit.ok(); ++lit) {
      lit().m_velocity = interpolateVecFieldCIC(a_rho, lit().m_position, a_lo, a_iDx);
    }
  }

  void
  interpVecVector(const FArrayBox& a_rho, std::vector<MoverParticle>& a_v, const RealVect& a_lo, const RealVect& a_iDx)
  {
    const std::size_t n = a_v.size();
    for (std::size_t i = 0; i < n; i++) {
      a_v[i].m_velocity = interpolateVecFieldCIC(a_rho, a_v[i].m_position, a_lo, a_iDx);
    }
  }

  void
  interpVecSoA(const FArrayBox& a_rho, ParticleSoA<MoverParticle>& a_soa, const RealVect& a_lo, const RealVect& a_iDx)
  {
    const std::vector<RealVect>& pos = a_soa.column<&MoverParticle::m_position>();
    std::vector<RealVect>&       vel = a_soa.column<&MoverParticle::m_velocity>();
    const std::size_t            n   = a_soa.size();
    for (std::size_t i = 0; i < n; i++) {
      vel[i] = interpolateVecFieldCIC(a_rho, pos[i], a_lo, a_iDx);
    }
  }

  void
  transformList(List<BenchParticle>& a_p)
  {
    for (ListIterator<BenchParticle> lit(a_p); lit.ok(); ++lit) {
      lit().m_weight = transformOne(lit().position());
    }
  }

  // noinline so the loop gets its own symbol for assembly inspection; the single
  // call per 65k-iteration loop is negligible for timing.
  __attribute__((noinline)) void
  transformSoA(SoA& a_soa)
  {
    const std::vector<RealVect>& pos = a_soa.column<&BenchParticle::m_position>();
    std::vector<Real>&           w   = a_soa.column<&BenchParticle::m_weight>();
    ParticleLoops::loop(a_soa, [&](std::size_t i) {
      w[i] = transformOne(pos[i]);
    });
  }

  // Per-component (x[], y[], z[]) layout: the same transform over three SEPARATE
  // contiguous Real arrays. Unit-stride, no interleaving -> should hit wide AVX.
  __attribute__((noinline)) void
  transformPerComponent(const std::vector<Real>& a_x,
                        const std::vector<Real>& a_y,
                        const std::vector<Real>& a_z,
                        std::vector<Real>&       a_w)
  {
    const Real* x = a_x.data();
    const Real* y = a_y.data();
    const Real* z = a_z.data();
    Real*       w = a_w.data();
    ParticleLoops::loop(a_w.size(), [&](std::size_t i) {
      w[i] = 0.5 * (x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
    });
  }

  // Contiguous AoS: the same transform over a vector<BenchParticle>. Memory is
  // contiguous (unlike List) but each field is strided by sizeof(BenchParticle),
  // so SIMD must gather/scatter within the struct.
  __attribute__((noinline)) void
  transformVector(std::vector<BenchParticle>& a_v)
  {
    BenchParticle* p = a_v.data();
    ParticleLoops::loop(a_v.size(), [&](std::size_t i) {
      p[i].m_weight = transformOne(p[i].m_position);
    });
  }

  // ---------------------------------------------------------------------------
  // Euler advance: x += v*dt (the canonical particle update; reads two vector
  // fields, writes one). Same math in all four layouts.
  // ---------------------------------------------------------------------------
  using MoverSoA = ParticleSoA<MoverParticle>;

  __attribute__((noinline)) void
  advanceList(List<MoverParticle>& a_p, const Real a_dt)
  {
    for (ListIterator<MoverParticle> lit(a_p); lit.ok(); ++lit) {
      lit().m_position += lit().m_velocity * a_dt;
    }
  }

  __attribute__((noinline)) void
  advanceVector(std::vector<MoverParticle>& a_v, const Real a_dt)
  {
    MoverParticle* p = a_v.data();
    ParticleLoops::loop(a_v.size(), [&](std::size_t i) {
      p[i].m_position += p[i].m_velocity * a_dt;
    });
  }

  __attribute__((noinline)) void
  advanceSoA(MoverSoA& a_soa, const Real a_dt)
  {
    std::vector<RealVect>&       pos = a_soa.column<&MoverParticle::m_position>();
    const std::vector<RealVect>& vel = a_soa.column<&MoverParticle::m_velocity>();
    ParticleLoops::loop(a_soa, [&](std::size_t i) {
      pos[i] += vel[i] * a_dt;
    });
  }

  __attribute__((noinline)) void
  advancePerComponent(Real*             a_px,
                      Real*             a_py,
                      Real*             a_pz,
                      const Real*       a_vx,
                      const Real*       a_vy,
                      const Real*       a_vz,
                      const std::size_t a_n,
                      const Real        a_dt)
  {
    ParticleLoops::loop(a_n, [&](std::size_t i) {
      a_px[i] += a_vx[i] * a_dt;
      a_py[i] += a_vy[i] * a_dt;
      a_pz[i] += a_vz[i] * a_dt;
    });
  }

  // ---------------------------------------------------------------------------
  // KD-tree merge: recursively median-split particle INDICES by position (this
  // reads positions only), then merge each leaf's particles into one conserving
  // total weight + center of mass (this gathers WHOLE particles). The split phase
  // favours column storage (position-only); the leaf merge favours AoS (whole
  // particle). a_get(k) returns the full MergeParticle for index entry k.
  // ---------------------------------------------------------------------------
  template <typename Get>
  MergeParticle
  mergeRange(const int a_lo, const int a_hi, Get&& a_get)
  {
    Real     W    = 0.0;
    RealVect pos  = RealVect::Zero;
    RealVect vel  = RealVect::Zero;
    RealVect old  = RealVect::Zero;
    Real     mob  = 0.0;
    Real     dif  = 0.0;
    Real     ener = 0.0;
    for (int k = a_lo; k < a_hi; k++) {
      const MergeParticle p = a_get(k);
      const Real          w = p.m_weight;
      W += w;
      pos += w * p.m_position;
      vel += w * p.m_velocity;
      old += w * p.m_oldPos;
      mob += w * p.m_mobility;
      dif += w * p.m_diffusion;
      ener += w * p.m_energy;
    }
    const Real    inv = 1.0 / W;
    MergeParticle m;
    m.m_weight    = W;
    m.m_position  = pos * inv;
    m.m_velocity  = vel * inv;
    m.m_oldPos    = old * inv;
    m.m_mobility  = mob * inv;
    m.m_diffusion = dif * inv;
    m.m_energy    = ener * inv;
    return m;
  }

  template <typename PosDim, typename EmitLeaf>
  void
  kdMerge(std::vector<int>& a_idx,
          const int         a_lo,
          const int         a_hi,
          const int         a_depth,
          const int         a_leafSize,
          PosDim&&          a_posDim,
          EmitLeaf&&        a_emitLeaf)
  {
    if (a_hi - a_lo <= a_leafSize) {
      a_emitLeaf(a_lo, a_hi);
      return;
    }
    const int dim = a_depth % SpaceDim;
    const int mid = (a_lo + a_hi) / 2;
    std::nth_element(a_idx.begin() + a_lo, a_idx.begin() + mid, a_idx.begin() + a_hi, [&](int a, int b) {
      return a_posDim(a, dim) < a_posDim(b, dim);
    });
    kdMerge(a_idx, a_lo, mid, a_depth + 1, a_leafSize, a_posDim, a_emitLeaf);
    kdMerge(a_idx, mid, a_hi, a_depth + 1, a_leafSize, a_posDim, a_emitLeaf);
  }

  std::vector<MergeParticle>
  mergeVectorKD(const std::vector<MergeParticle>& a_v, const int a_leafSize)
  {
    const int        n = static_cast<int>(a_v.size());
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<MergeParticle> out;
    out.reserve(n / a_leafSize + 1);
    kdMerge(
      idx,
      0,
      n,
      0,
      a_leafSize,
      [&](int a, int d) {
        return a_v[a].m_position[d];
      },
      [&](int lo, int hi) {
        out.push_back(mergeRange(lo, hi, [&](int k) {
          return a_v[idx[k]];
        }));
      });
    return out;
  }

  ParticleSoA<MergeParticle>
  mergeSoAKD(ParticleSoA<MergeParticle>& a_soa, const int a_leafSize)
  {
    const int        n = static_cast<int>(a_soa.size());
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    ParticleSoA<MergeParticle> out;
    out.reserve(n / a_leafSize + 1);
    kdMerge(
      idx,
      0,
      n,
      0,
      a_leafSize,
      [&](int a, int d) {
        return a_soa.position(a)[d];
      }, // position column only
      [&](int lo, int hi) {
        out.append(mergeRange(lo, hi, [&](int k) {
          return a_soa.gather(idx[k]);
        })); // gathers all 7 columns
      });
    return out;
  }

  std::vector<MergeParticle>
  mergeListKD(const List<MergeParticle>& a_p, const int a_leafSize)
  {
    std::vector<MergeParticle> tmp; // realistic: List is copied into a contiguous buffer for KD work
    for (ListIterator<MergeParticle> lit(a_p); lit.ok(); ++lit) {
      tmp.push_back(lit());
    }
    return mergeVectorKD(tmp, a_leafSize);
  }

  // ---------------------------------------------------------------------------
  // Helpers.
  // ---------------------------------------------------------------------------

  Real
  totalMass(const FArrayBox& a_rho)
  {
    const Real*       p = a_rho.dataPtr(0);
    const std::size_t n = a_rho.box().numPts();
    Real              s = 0.0;
    for (std::size_t i = 0; i < n; i++) {
      s += p[i];
    }
    return s;
  }

  Real
  maxDiff(const FArrayBox& a_a, const FArrayBox& a_b)
  {
    const Real*       pa = a_a.dataPtr(0);
    const Real*       pb = a_b.dataPtr(0);
    const std::size_t n  = a_a.box().numPts();
    Real              m  = 0.0;
    for (std::size_t i = 0; i < n; i++) {
      m = std::max(m, std::abs(pa[i] - pb[i]));
    }
    return m;
  }

  // ---------------------------------------------------------------------------
  // MPI buffer packing: serialize (position + weight) into a byte send buffer.
  // ---------------------------------------------------------------------------
  static_assert(sizeof(BenchParticle) == sizeof(RealVect) + sizeof(Real),
                "BenchParticle must be padding-free for the AoS bulk memcpy");
  constexpr std::size_t s_bytesPerParticle = sizeof(RealVect) + sizeof(Real);

  /** @brief Pack a List per particle (pointer-chased; the realistic Chombo path). */
  void
  packList(unsigned char* a_buf, const List<BenchParticle>& a_p)
  {
    unsigned char* q = a_buf;
    for (ListIterator<BenchParticle> lit(a_p); lit.ok(); ++lit) {
      std::memcpy(q, &lit().m_position, sizeof(RealVect));
      q += sizeof(RealVect);
      std::memcpy(q, &lit().m_weight, sizeof(Real));
      q += sizeof(Real);
    }
  }

  /** @brief Pack a contiguous AoS vector in a single bulk memcpy (buffer is AoS-formatted). */
  void
  packVector(unsigned char* a_buf, const std::vector<BenchParticle>& a_v)
  {
    std::memcpy(a_buf, a_v.data(), a_v.size() * sizeof(BenchParticle));
  }

  /** @brief Pack SoA as per-column bulk memcpies (buffer is column-major). */
  void
  packSoA(unsigned char* a_buf, const SoA& a_soa)
  {
    const std::vector<RealVect>& pos = a_soa.positions();
    const std::vector<Real>&     w   = a_soa.weights();
    const std::size_t            n   = a_soa.size();
    std::memcpy(a_buf, pos.data(), n * sizeof(RealVect));
    std::memcpy(a_buf + n * sizeof(RealVect), w.data(), n * sizeof(Real));
  }

  /** @brief Copy a_nd doubles src->dst with non-temporal (streaming) stores. dst must be 32-byte aligned. */
  inline void
  ntCopyDoubles(double* a_dst, const double* a_src, const std::size_t a_nd) noexcept
  {
    std::size_t i = 0;
    for (; i + 4 <= a_nd; i += 4) {
      _mm256_stream_pd(a_dst + i, _mm256_loadu_pd(a_src + i)); // a_dst+i is 32-aligned (i multiple of 4)
    }
    for (; i < a_nd; i++) {
      a_dst[i] = a_src[i];
    }
    _mm_sfence();
  }

  /** @brief Pack SoA columns with NON-TEMPORAL stores (tests the memcpy NT-threshold hypothesis). */
  void
  packSoANT(unsigned char* a_buf, const SoA& a_soa)
  {
    const std::vector<RealVect>& pos = a_soa.positions();
    const std::vector<Real>&     w   = a_soa.weights();
    const std::size_t            n   = a_soa.size();
    double* const                d   = reinterpret_cast<double*>(a_buf);
    ntCopyDoubles(d, reinterpret_cast<const double*>(pos.data()), n * SpaceDim);
    ntCopyDoubles(d + n * SpaceDim, reinterpret_cast<const double*>(w.data()), n);
  }

  /** @brief Control: pack the contiguous AoS in TWO memcpy calls (same sizes as the SoA columns). */
  void
  packVectorSplit(unsigned char* a_buf, const std::vector<BenchParticle>& a_v)
  {
    const std::size_t    n   = a_v.size();
    const unsigned char* src = reinterpret_cast<const unsigned char*>(a_v.data());
    std::memcpy(a_buf, src, n * sizeof(RealVect));
    std::memcpy(a_buf + n * sizeof(RealVect), src + n * sizeof(RealVect), n * sizeof(Real));
  }

  /** @brief Time the average wall-time (ns) of a_op over a_reps repetitions (one warm-up). */
  template <typename Op>
  double
  timeOp(const int a_reps, Op&& a_op)
  {
    using clock = std::chrono::steady_clock;
    a_op();
    double ns = 0.0;
    for (int r = 0; r < a_reps; r++) {
      const auto t0 = clock::now();
      a_op();
      const auto t1 = clock::now();
      ns += std::chrono::duration<double, std::nano>(t1 - t0).count();
    }
    return ns / static_cast<double>(a_reps);
  }

} // namespace

int
main(int argc, char* argv[])
{
  initialize(argc, argv);

  {
    // ----- Configuration --------------------------------------------------------
    constexpr int nCell     = 16;
    constexpr int nGhost    = 2;
    constexpr int ppc       = 16;
    constexpr int fastReps  = 200; // deposition / interpolation / transform
    constexpr int buildReps = 50;  // build / remap (each rep allocates)
    constexpr int nBuckets  = 64;  // remap destination patches

    const Box         valid(IntVect::Zero, (nCell - 1) * IntVect::Unit);
    const Box         grown  = grow(valid, nGhost);
    const RealVect    probLo = RealVect::Zero;
    const RealVect    invDx  = RealVect::Unit; // dx = 1
    const Real        invVol = 1.0;
    const std::size_t nPart  = static_cast<std::size_t>(nCell) * nCell * nCell * ppc;

    // ----- Particles in FArrayBox (Fortran) order, plus a randomized copy --------
    std::mt19937                         rng(20260620u);
    std::uniform_real_distribution<Real> u(0.0, 1.0);

    std::vector<BenchParticle> sortedParticles;
    sortedParticles.reserve(nPart);
    for (int k = 0; k < nCell; k++) {
      for (int j = 0; j < nCell; j++) {
        for (int i = 0; i < nCell; i++) {
          for (int p = 0; p < ppc; p++) {
            const Real    ux = u(rng);
            const Real    uy = u(rng);
            const Real    uz = u(rng);
            BenchParticle bp;
            bp.m_position = probLo + RealVect(D_DECL(i + ux, j + uy, k + uz));
            bp.m_weight   = 1.0;
            sortedParticles.push_back(bp);
          }
        }
      }
    }
    std::vector<BenchParticle> randomParticles = sortedParticles;
    std::shuffle(randomParticles.begin(), randomParticles.end(), rng);

    // ----- Container builders ---------------------------------------------------
    auto buildCompactList = [](const std::vector<BenchParticle>& a_src) {
      List<BenchParticle> lst;
      for (const BenchParticle& p : a_src) {
        lst.append(p);
      }
      return lst;
    };
    std::vector<std::vector<char>> junk; // keeps the heap fragmented for the next build
    junk.reserve(nPart);
    auto buildFragmentedList = [&](const std::vector<BenchParticle>& a_src) {
      List<BenchParticle> lst;
      for (const BenchParticle& p : a_src) {
        lst.append(p);
        junk.emplace_back(64 + (rng() % 256));
      }
      return lst;
    };
    auto buildSoA = [](const std::vector<BenchParticle>& a_src) {
      SoA soa;
      soa.reserve(a_src.size());
      for (const BenchParticle& p : a_src) {
        soa.append(p);
      }
      return soa;
    };

    List<BenchParticle> listSorted     = buildCompactList(sortedParticles);
    List<BenchParticle> listRandom     = buildCompactList(randomParticles);
    List<BenchParticle> listFragSorted = buildFragmentedList(sortedParticles);
    List<BenchParticle> listFragRandom = buildFragmentedList(randomParticles);
    SoA                 soaSorted      = buildSoA(sortedParticles);
    SoA                 soaRandom      = buildSoA(randomParticles);

    FArrayBox rho(grown, 1);

    // =====================================================================
    // (1) DEPOSITION
    // =====================================================================
    rho.setVal(0.0);
    depositList(rho, listSorted, probLo, invDx, invVol);
    FArrayBox refList(grown, 1);
    refList.copy(rho);
    const Real massList = totalMass(rho);

    rho.setVal(0.0);
    depositSoA(rho, soaSorted, probLo, invDx, invVol);
    const Real massSoA    = totalMass(rho);
    const Real listVsSoA  = maxDiff(rho, refList);
    const Real expectMass = static_cast<Real>(nPart) * invVol;
    if (std::abs(massList - expectMass) > 1.E-6 * expectMass || std::abs(massSoA - expectMass) > 1.E-6 * expectMass ||
        listVsSoA > 1.E-12 * expectMass) {
      MayDay::Abort("benchmark: deposition correctness check failed");
    }

    auto timeDeposit = [&](auto&& a_deposit) {
      return timeOp(fastReps, [&]() {
        rho.setVal(0.0);
        a_deposit();
      });
    };
    const double dListS = timeDeposit([&]() {
      depositList(rho, listSorted, probLo, invDx, invVol);
    });
    const double dListR = timeDeposit([&]() {
      depositList(rho, listRandom, probLo, invDx, invVol);
    });
    const double dFragS = timeDeposit([&]() {
      depositList(rho, listFragSorted, probLo, invDx, invVol);
    });
    const double dFragR = timeDeposit([&]() {
      depositList(rho, listFragRandom, probLo, invDx, invVol);
    });
    const double dSoAS  = timeDeposit([&]() {
      depositSoA(rho, soaSorted, probLo, invDx, invVol);
    });
    const double dSoAR  = timeDeposit([&]() {
      depositSoA(rho, soaRandom, probLo, invDx, invVol);
    });
    const double dVecS  = timeDeposit([&]() {
      depositVector(rho, sortedParticles, probLo, invDx, invVol);
    });
    const double dVecR  = timeDeposit([&]() {
      depositVector(rho, randomParticles, probLo, invDx, invVol);
    });

    // =====================================================================
    // (2) INTERPOLATION  (fill rho once, then gather to particle weights)
    // =====================================================================
    rho.setVal(0.0);
    depositSoA(rho, soaSorted, probLo, invDx, invVol); // some non-trivial field to gather

    interpolateList(rho, listSorted, probLo, invDx);
    interpolateSoA(rho, soaSorted, probLo, invDx);

    const double iListS = timeOp(fastReps, [&]() {
      interpolateList(rho, listSorted, probLo, invDx);
    });
    const double iListR = timeOp(fastReps, [&]() {
      interpolateList(rho, listRandom, probLo, invDx);
    });
    const double iFragS = timeOp(fastReps, [&]() {
      interpolateList(rho, listFragSorted, probLo, invDx);
    });
    const double iFragR = timeOp(fastReps, [&]() {
      interpolateList(rho, listFragRandom, probLo, invDx);
    });
    const double iSoAS  = timeOp(fastReps, [&]() {
      interpolateSoA(rho, soaSorted, probLo, invDx);
    });
    const double iSoAR  = timeOp(fastReps, [&]() {
      interpolateSoA(rho, soaRandom, probLo, invDx);
    });
    const double iVecS  = timeOp(fastReps, [&]() {
      interpolateVector(rho, sortedParticles, probLo, invDx);
    });
    const double iVecR  = timeOp(fastReps, [&]() {
      interpolateVector(rho, randomParticles, probLo, invDx);
    });

    // =====================================================================
    // (3) STREAMING TRANSFORM  (SoA path via ParticleLoops SIMD)
    // =====================================================================
    // Per-component (x[], y[], z[]) arrays for the same particles, to test whether a
    // non-interleaved position layout vectorizes wider than vector<RealVect>.
    std::vector<Real> xs(nPart), ys(nPart), zs(nPart), ws(nPart);
    for (std::size_t i = 0; i < nPart; i++) {
      const RealVect& x = sortedParticles[i].m_position;
      xs[i]             = x[0];
      ys[i]             = x[1];
      zs[i]             = x[2];
    }

    const double tListS = timeOp(fastReps, [&]() {
      transformList(listSorted);
    });
    const double tFragS = timeOp(fastReps, [&]() {
      transformList(listFragSorted);
    });
    const double tVecS  = timeOp(fastReps, [&]() {
      transformVector(sortedParticles);
    });
    const double tSoAS  = timeOp(fastReps, [&]() {
      transformSoA(soaSorted);
    });
    const double tPCS   = timeOp(fastReps, [&]() {
      transformPerComponent(xs, ys, zs, ws);
    });

    // =====================================================================
    // (4) BUILD  (construct container from the source array; allocation cost)
    // =====================================================================
    const double bList = timeOp(buildReps, [&]() {
      List<BenchParticle> l;
      for (const BenchParticle& p : sortedParticles) {
        l.append(p);
      }
      g_sink += l.isEmpty() ? 0u : 1u;
    });
    const double bSoA  = timeOp(buildReps, [&]() {
      SoA s;
      s.reserve(nPart);
      for (const BenchParticle& p : sortedParticles) {
        s.append(p);
      }
      g_sink += s.size();
    });
    const double bVec  = timeOp(buildReps, [&]() {
      std::vector<BenchParticle> v;
      v.reserve(nPart);
      for (const BenchParticle& p : sortedParticles) {
        v.push_back(p);
      }
      g_sink += v.size();
    });

    // =====================================================================
    // (5) REMAP  (scatter particles into nBuckets destination containers)
    // =====================================================================
    std::vector<int> bucketId(nPart);
    for (std::size_t i = 0; i < nPart; i++) {
      const RealVect& x = sortedParticles[i].m_position;
      const int       b = (static_cast<int>(x[0]) + nCell * static_cast<int>(x[1])) % nBuckets;
      bucketId[i]       = (b + nBuckets) % nBuckets;
    }
    const double rList = timeOp(buildReps, [&]() {
      std::vector<List<BenchParticle>> dst(nBuckets);
      for (std::size_t i = 0; i < nPart; i++) {
        dst[bucketId[i]].append(sortedParticles[i]);
      }
      for (const List<BenchParticle>& d : dst) {
        g_sink += d.isEmpty() ? 0u : 1u;
      }
    });
    const double rSoA  = timeOp(buildReps, [&]() {
      std::vector<SoA> dst(nBuckets);
      for (std::size_t i = 0; i < nPart; i++) {
        dst[bucketId[i]].append(sortedParticles[i]);
      }
      for (const SoA& d : dst) {
        g_sink += d.size();
      }
    });
    const double rVec  = timeOp(buildReps, [&]() {
      std::vector<std::vector<BenchParticle>> dst(nBuckets);
      for (std::size_t i = 0; i < nPart; i++) {
        dst[bucketId[i]].push_back(sortedParticles[i]);
      }
      for (const std::vector<BenchParticle>& d : dst) {
        g_sink += d.size();
      }
    });

    // =====================================================================
    // (6) MPI PACKING  (serialize position + weight into a byte send buffer)
    // =====================================================================
    const std::size_t bufBytes = nPart * s_bytesPerParticle;
    const std::size_t bufAlloc = ((bufBytes + 63) / 64) * 64; // round up for aligned_alloc
    unsigned char*    packBuf  = static_cast<unsigned char*>(std::aligned_alloc(64, bufAlloc));

    // Correctness: the split AoS pack and the NT SoA pack must produce identical bytes
    // to their plain counterparts.
    {
      std::vector<unsigned char> refV(bufBytes), refS(bufBytes);
      packVector(refV.data(), sortedParticles);
      packVectorSplit(packBuf, sortedParticles);
      if (std::memcmp(refV.data(), packBuf, bufBytes) != 0) {
        MayDay::Abort("benchmark: packVectorSplit produced different bytes");
      }
      packSoA(refS.data(), soaSorted);
      packSoANT(packBuf, soaSorted);
      if (std::memcmp(refS.data(), packBuf, bufBytes) != 0) {
        MayDay::Abort("benchmark: packSoANT produced different bytes");
      }
    }

    const double pList     = timeOp(fastReps, [&]() {
      packList(packBuf, listSorted);
      g_sink += packBuf[0];
    });
    const double pVec      = timeOp(fastReps, [&]() {
      packVector(packBuf, sortedParticles);
      g_sink += packBuf[0];
    });
    const double pVecSplit = timeOp(fastReps, [&]() {
      packVectorSplit(packBuf, sortedParticles);
      g_sink += packBuf[0];
    });
    const double pSoA      = timeOp(fastReps, [&]() {
      packSoA(packBuf, soaSorted);
      g_sink += packBuf[0];
    });
    const double pSoANT    = timeOp(fastReps, [&]() {
      packSoANT(packBuf, soaSorted);
      g_sink += packBuf[0];
    });
    std::free(packBuf);

    // =====================================================================
    // Mover containers (position + velocity), shared by (7) and (8).
    // =====================================================================
    const Real dt = 0.1;

    std::vector<MoverParticle> moverSrc(nPart);
    for (std::size_t i = 0; i < nPart; i++) {
      moverSrc[i].m_position = sortedParticles[i].m_position;
      moverSrc[i].m_velocity = RealVect(D_DECL(1.0, 0.5, 0.25)); // arbitrary drift
    }

    List<MoverParticle> moverList;
    for (const MoverParticle& p : moverSrc) {
      moverList.append(p);
    }
    std::vector<MoverParticle> moverVec = moverSrc;
    MoverSoA                   moverSoA;
    moverSoA.reserve(nPart);
    for (const MoverParticle& p : moverSrc) {
      moverSoA.append(p);
    }
    std::vector<Real> px(nPart), py(nPart), pz(nPart), vx(nPart), vy(nPart), vz(nPart);
    for (std::size_t i = 0; i < nPart; i++) {
      px[i] = moverSrc[i].m_position[0];
      py[i] = moverSrc[i].m_position[1];
      pz[i] = moverSrc[i].m_position[2];
      vx[i] = moverSrc[i].m_velocity[0];
      vy[i] = moverSrc[i].m_velocity[1];
      vz[i] = moverSrc[i].m_velocity[2];
    }

    // =====================================================================
    // (7) VECTOR-FIELD INTERPOLATION  (3-component FArrayBox -> particle RealVect)
    // Runs BEFORE the advance, while positions are still inside the box.
    // =====================================================================
    FArrayBox vfield(grown, SpaceDim);
    for (int comp = 0; comp < SpaceDim; comp++) {
      vfield.setVal(static_cast<Real>(comp + 1), comp); // constant per component: CIC must return (1,2,..)
    }
    interpVecSoA(vfield, moverSoA, probLo, invDx);
    {
      const MoverParticle q = moverSoA.gather(123);
      for (int comp = 0; comp < SpaceDim; comp++) {
        if (std::abs(q.m_velocity[comp] - static_cast<Real>(comp + 1)) > 1.E-9) {
          MayDay::Abort("benchmark: vector interpolation of a constant field is wrong");
        }
      }
    }

    const double viList = timeOp(fastReps, [&]() {
      interpVecList(vfield, moverList, probLo, invDx);
    });
    const double viVec  = timeOp(fastReps, [&]() {
      interpVecVector(vfield, moverVec, probLo, invDx);
    });
    const double viSoA  = timeOp(fastReps, [&]() {
      interpVecSoA(vfield, moverSoA, probLo, invDx);
    });

    // =====================================================================
    // (8) EULER ADVANCE  x += v*dt  (canonical particle update)
    // =====================================================================
    const double aList = timeOp(fastReps, [&]() {
      advanceList(moverList, dt);
    });
    const double aVec  = timeOp(fastReps, [&]() {
      advanceVector(moverVec, dt);
    });
    const double aSoA  = timeOp(fastReps, [&]() {
      advanceSoA(moverSoA, dt);
    });
    const double aPC   = timeOp(fastReps, [&]() {
      advancePerComponent(px.data(), py.data(), pz.data(), vx.data(), vy.data(), vz.data(), nPart, dt);
    });

    // =====================================================================
    // (9) KD-TREE MERGE  (32 particles/cell on the 16^3 grid; whole-particle access)
    // =====================================================================
    constexpr int     mergePpc  = 32;
    constexpr int     leafSize  = 2; // merge nearest pairs
    constexpr int     mergeReps = 30;
    const std::size_t nMerge    = static_cast<std::size_t>(nCell) * nCell * nCell * mergePpc;

    std::vector<MergeParticle> mergeSrc;
    mergeSrc.reserve(nMerge);
    for (int k = 0; k < nCell; k++) {
      for (int j = 0; j < nCell; j++) {
        for (int i = 0; i < nCell; i++) {
          for (int p = 0; p < mergePpc; p++) {
            const Real    ux = u(rng);
            const Real    uy = u(rng);
            const Real    uz = u(rng);
            MergeParticle mp;
            mp.m_position = probLo + RealVect(D_DECL(i + ux, j + uy, k + uz));
            mp.m_weight   = 1.0 + u(rng); // non-uniform weights
            mergeSrc.push_back(mp);
          }
        }
      }
    }

    List<MergeParticle> mergeList;
    for (const MergeParticle& p : mergeSrc) {
      mergeList.append(p);
    }
    ParticleSoA<MergeParticle> mergeSoa;
    mergeSoa.reserve(nMerge);
    for (const MergeParticle& p : mergeSrc) {
      mergeSoa.append(p);
    }

    // Correctness: merging conserves total weight and center of mass.
    auto weightAndCom = [](const std::vector<MergeParticle>& a_v) {
      Real     W   = 0.0;
      RealVect com = RealVect::Zero;
      for (const MergeParticle& p : a_v) {
        W += p.m_weight;
        com += p.m_weight * p.m_position;
      }
      com *= (1.0 / W);
      return std::make_pair(W, com);
    };
    {
      const auto                       before = weightAndCom(mergeSrc);
      const std::vector<MergeParticle> merged = mergeVectorKD(mergeSrc, leafSize);
      const auto                       after  = weightAndCom(merged);
      if (std::abs(after.first - before.first) > 1.E-6 * before.first) {
        MayDay::Abort("benchmark: KD merge did not conserve weight");
      }
      for (int d = 0; d < SpaceDim; d++) {
        if (std::abs(after.second[d] - before.second[d]) > 1.E-9 * (1.0 + std::abs(before.second[d]))) {
          MayDay::Abort("benchmark: KD merge did not conserve center of mass");
        }
      }
    }

    const double mList = timeOp(mergeReps, [&]() {
      const std::vector<MergeParticle> out = mergeListKD(mergeList, leafSize);
      g_sink += out.size();
    });
    const double mVec  = timeOp(mergeReps, [&]() {
      const std::vector<MergeParticle> out = mergeVectorKD(mergeSrc, leafSize);
      g_sink += out.size();
    });
    const double mSoA  = timeOp(mergeReps, [&]() {
      const ParticleSoA<MergeParticle> out = mergeSoAKD(mergeSoa, leafSize);
      g_sink += out.size();
    });

    // =====================================================================
    // Report
    // =====================================================================
    auto row = [&](const char* a_name, const double a_ns) {
      pout() << "    " << a_name << ": " << a_ns / 1.0e6 << " ms, " << a_ns / static_cast<double>(nPart)
             << " ns/particle" << endl;
    };

    pout() << "Particle benchmarks (SpaceDim=" << SpaceDim << ")" << endl;
    pout() << "  grid       : " << nCell << "^3 valid + " << nGhost << " ghost (" << grown.size(0) << "^3)" << endl;
    pout() << "  particles  : " << nPart << " (" << ppc << "/cell), fastReps=" << fastReps
           << ", buildReps=" << buildReps << endl;
    pout() << "  mass check : expected " << expectMass << ", List " << massList << ", SoA " << massSoA
           << " (List-vs-SoA " << listVsSoA << ")" << endl;

    pout() << "  (1) DEPOSITION (scatter to grid)" << endl;
    row("List fragmented sorted", dFragS);
    row("List fragmented random", dFragR);
    row("vector<P> AoS   sorted", dVecS);
    row("vector<P> AoS   random", dVecR);
    row("SoA             sorted", dSoAS);
    row("SoA             random", dSoAR);
    pout() << "    SoA vs fragmented List / vs vector<P> (sorted): " << dFragS / dSoAS << "x / " << dVecS / dSoAS << "x"
           << endl;

    pout() << "  (2) INTERPOLATION (gather from grid)" << endl;
    row("List fragmented sorted", iFragS);
    row("List fragmented random", iFragR);
    row("vector<P> AoS   sorted", iVecS);
    row("vector<P> AoS   random", iVecR);
    row("SoA             sorted", iSoAS);
    row("SoA             random", iSoAR);
    pout() << "    SoA vs fragmented List / vs vector<P> (sorted): " << iFragS / iSoAS << "x / " << iVecS / iSoAS << "x"
           << endl;

    pout() << "  (3) STREAMING TRANSFORM (weight = 0.5|x|^2)" << endl;
    row("List compact      sorted", tListS);
    row("List fragmented   sorted", tFragS);
    row("vector<P> AoS           ", tVecS);
    row("SoA vector<RealVect>    ", tSoAS);
    row("SoA per-component x/y/z ", tPCS);
    pout() << "    vs compact List: vector<P> " << tListS / tVecS << "x, SoA(RealVect) " << tListS / tSoAS
           << "x, per-component " << tListS / tPCS << "x" << endl;
    pout() << "    SoA(RealVect) vs vector<P>: " << tVecS / tSoAS
           << "x; per-component vs SoA(RealVect): " << tSoAS / tPCS << "x" << endl;

    pout() << "  (4) BUILD (allocation)" << endl;
    row("List     ", bList);
    row("vector<P>", bVec);
    row("SoA      ", bSoA);
    pout() << "    vs List: vector<P> " << bList / bVec << "x, SoA " << bList / bSoA << "x" << endl;

    pout() << "  (5) REMAP (scatter into " << nBuckets << " buckets)" << endl;
    row("List     ", rList);
    row("vector<P>", rVec);
    row("SoA      ", rSoA);
    pout() << "    vs List: vector<P> " << rList / rVec << "x, SoA " << rList / rSoA << "x" << endl;

    pout() << "  (6) MPI PACKING (serialize " << s_bytesPerParticle << " B/particle to send buffer)" << endl;
    row("List          (per-particle)   ", pList);
    row("vector<P>     (1 bulk memcpy)  ", pVec);
    row("vector<P>     (split 2 memcpy) ", pVecSplit);
    row("SoA           (per-column memcpy)", pSoA);
    row("SoA           (per-column NT)  ", pSoANT);
    pout() << "    vs List: vector<P> " << pList / pVec << "x, SoA " << pList / pSoA << "x, SoA-NT " << pList / pSoANT
           << "x" << endl;
    pout() << "    split-vs-single vector<P>: " << pVec / pVecSplit << "x; NT-vs-plain SoA: " << pSoA / pSoANT << "x"
           << endl;

    pout() << "  (7) VECTOR-FIELD INTERPOLATION (3-comp FArrayBox -> particle RealVect)" << endl;
    row("List     ", viList);
    row("vector<P>", viVec);
    row("SoA      ", viSoA);
    pout() << "    vs List: vector<P> " << viList / viVec << "x, SoA " << viList / viSoA << "x" << endl;

    pout() << "  (8) EULER ADVANCE (x += v*dt)" << endl;
    row("List                   ", aList);
    row("vector<P> AoS          ", aVec);
    row("SoA vector<RealVect>   ", aSoA);
    row("SoA per-component      ", aPC);
    pout() << "    vs List: vector<P> " << aList / aVec << "x, SoA(RealVect) " << aList / aSoA << "x, per-component "
           << aList / aPC << "x" << endl;
    pout() << "    SoA(RealVect) vs vector<P>: " << aVec / aSoA << "x; per-component vs SoA(RealVect): " << aSoA / aPC
           << "x" << endl;

    pout() << "  (9) KD-TREE MERGE (" << nMerge << " 7-field particles, " << mergePpc << "/cell, leaf=" << leafSize
           << ", reps=" << mergeReps << ")" << endl;
    auto mrow = [&](const char* a_name, const double a_ns) {
      pout() << "    " << a_name << ": " << a_ns / 1.0e6 << " ms, " << a_ns / static_cast<double>(nMerge)
             << " ns/particle" << endl;
    };
    mrow("List (copy+KD) ", mList);
    mrow("vector<P> AoS  ", mVec);
    mrow("SoA            ", mSoA);
    pout() << "    vs List: vector<P> " << mList / mVec << "x, SoA " << mList / mSoA
           << "x; SoA vs vector<P>: " << mVec / mSoA << "x" << endl;

    pout() << "  (sink=" << g_sink << ")" << endl;
  }

  return finalize();
}
