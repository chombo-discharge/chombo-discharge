/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   main.cpp
  @brief  Benchmarks comparing Chombo List<P> and vector<P> (AoS) with ParticleSoA.
  @author Robert Marskar
  @details
  Micro-benchmarks on a single grid patch, all using the identical particle data and math,
  so the measured difference is the storage layout. Three storages are compared: List<P>
  (linked-list AoS), vector<P> (contiguous AoS), and ParticleSoA (arena-backed SoA, the
  merged production container).

  ParticleSoA OWNS the position (per-component ParticleReal columns), weight, particleID and
  rankID; the user supplies only the extra payload. So for a byte-fair comparison the AoS
  baseline structs carry a matching particleID + rankID, and the SoA stores position+weight
  via accessors + per-component columns, with payloads:
    - BenchParticle (position+weight)            -> ParticleSoA<>            (empty payload)
    - MoverParticle (+velocity)                  -> ParticleSoA<MoverPayload> (per-component v)
    - MergeParticle (~ItoParticle, 7 fields)     -> ParticleSoA<MergePayload>

  Sections:
    1. CIC deposition          (particle weight -> FArrayBox; a scatter)
    2. CIC interpolation       (FArrayBox -> particle scalar; a gather)
    3. streaming transform     (weight = 0.5|x|^2; SoA path is per-component via ParticleLoops)
    4. build                   (construct the container; allocation cost)
    5. remap                   (scatter particles into K destination containers)
    6. MPI packing             (serialize a container to a byte buffer)
    7. vector-field interp     (3-component FArrayBox -> particle velocity)
    8. Euler advance           (x += v*dt; SoA path is per-component)
    9. KD-tree merge           (median-split by position, merge nearest; whole-particle access)
   10. container copy          (copy whole container -> fresh container)

  Because ParticleSoA is arena-backed, the SoA bulk paths (pack/copy) are a single memcpy of
  data() (glibc auto-streams; no intrinsics), and growth needs reserve(). All builds are 3D.
*/

// Std includes
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <random>
#include <utility>
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

  /** @brief Benchmark particle (AoS baseline): position + weight + container-owned metadata. */
  struct BenchParticle
  {
    RealVect   m_position = RealVect::Zero;
    Real       m_weight   = 0.0;
    ParticleID m_id       = -1;
    RankID     m_rank     = -1;

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

  /** @brief AoS baseline with an (interleaved) velocity, for the Euler advance x += v*dt. */
  struct MoverParticle
  {
    RealVect   m_position = RealVect::Zero;
    RealVect   m_velocity = RealVect::Zero;
    Real       m_weight   = 0.0;
    ParticleID m_id       = -1;
    RankID     m_rank     = -1;
  };

  /** @brief SoA payload for a mover: PER-COMPONENT velocity (unit-stride -> wide SIMD). */
  struct MoverPayload
  {
    ParticleReal vx = 0.0;
    ParticleReal vy = 0.0;
    ParticleReal vz = 0.0;
  };
  template <>
  struct ParticleTraits<MoverPayload>
  {
    static constexpr auto columns = std::make_tuple(&MoverPayload::vx, &MoverPayload::vy, &MoverPayload::vz);
  };

  /** @brief A richer AoS baseline (~ItoParticle: 3 RealVect + 4 Real + metadata) for the merge. */
  struct MergeParticle
  {
    RealVect   m_position  = RealVect::Zero;
    RealVect   m_velocity  = RealVect::Zero;
    RealVect   m_oldPos    = RealVect::Zero;
    Real       m_weight    = 0.0;
    Real       m_mobility  = 0.0;
    Real       m_diffusion = 0.0;
    Real       m_energy    = 0.0;
    ParticleID m_id        = -1;
    RankID     m_rank      = -1;
  };

  /** @brief SoA payload for the merge particle (position + weight are container-owned). */
  struct MergePayload
  {
    RealVect velocity  = RealVect::Zero;
    RealVect oldPos    = RealVect::Zero;
    Real     mobility  = 0.0;
    Real     diffusion = 0.0;
    Real     energy    = 0.0;
  };
  template <>
  struct ParticleTraits<MergePayload>
  {
    static constexpr auto columns = std::make_tuple(&MergePayload::velocity,
                                                    &MergePayload::oldPos,
                                                    &MergePayload::mobility,
                                                    &MergePayload::diffusion,
                                                    &MergePayload::energy);
  };

} // namespace ChomboDischarge

using namespace ChomboDischarge;

namespace {

  using BenchSoA = ParticleSoA<>;
  using MoverSoA = ParticleSoA<MoverPayload>;
  using MergeSoA = ParticleSoA<MergePayload>;

  volatile std::uint64_t g_sink = 0; // defeats dead-code elimination of build/remap

  // ---------------------------------------------------------------------------
  // Shared per-particle kernels (unchanged math across all storages).
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

  /**
    @brief CIC interpolation of a SpaceDim-component (vector) FArrayBox to a RealVect.
    @details The FArrayBox is component-major (component varies slowest); this gather is
    identical for every storage (the box is shared and const).
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

  /** @brief Streaming per-particle update: 0.5 * |position|^2. */
  inline Real
  transformOne(const RealVect& a_x) noexcept
  {
    return 0.5 * (D_TERM(a_x[0] * a_x[0], +a_x[1] * a_x[1], +a_x[2] * a_x[2]));
  }

  // ---------------------------------------------------------------------------
  // (1) DEPOSITION
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
  depositSoA(FArrayBox& a_rho, const BenchSoA& a_soa, const RealVect& a_lo, const RealVect& a_iDx, const Real a_iVol)
  {
    const ParticleReal* px = a_soa.positionColumn(0);
    const ParticleReal* py = a_soa.positionColumn(1);
    const ParticleReal* pz = a_soa.positionColumn(2);
    const ParticleReal* w  = a_soa.weightColumn();
    const std::size_t   n  = a_soa.size();
    for (std::size_t i = 0; i < n; i++) {
      depositOneCIC(a_rho, RealVect(D_DECL(px[i], py[i], pz[i])), w[i], a_lo, a_iDx, a_iVol);
    }
  }

  // ---------------------------------------------------------------------------
  // (2) INTERPOLATION
  // ---------------------------------------------------------------------------
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
  interpolateSoA(const FArrayBox& a_rho, BenchSoA& a_soa, const RealVect& a_lo, const RealVect& a_iDx)
  {
    const ParticleReal* px = a_soa.positionColumn(0);
    const ParticleReal* py = a_soa.positionColumn(1);
    const ParticleReal* pz = a_soa.positionColumn(2);
    ParticleReal*       w  = a_soa.weightColumn();
    const std::size_t   n  = a_soa.size();
    for (std::size_t i = 0; i < n; i++) {
      w[i] = interpolateOneCIC(a_rho, RealVect(D_DECL(px[i], py[i], pz[i])), a_lo, a_iDx);
    }
  }

  // ---------------------------------------------------------------------------
  // (3) STREAMING TRANSFORM (SoA path is per-component -> wide SIMD)
  // ---------------------------------------------------------------------------
  void
  transformList(List<BenchParticle>& a_p)
  {
    for (ListIterator<BenchParticle> lit(a_p); lit.ok(); ++lit) {
      lit().m_weight = transformOne(lit().position());
    }
  }

  __attribute__((noinline)) void
  transformVector(std::vector<BenchParticle>& a_v)
  {
    BenchParticle* p = a_v.data();
    ParticleLoops::loop(a_v.size(), [&](std::size_t i) {
      p[i].m_weight = transformOne(p[i].m_position);
    });
  }

  __attribute__((noinline)) void
  transformSoA(BenchSoA& a_soa)
  {
    const ParticleReal* x = a_soa.positionColumn(0);
    const ParticleReal* y = a_soa.positionColumn(1);
    const ParticleReal* z = a_soa.positionColumn(2);
    ParticleReal*       w = a_soa.weightColumn();
    ParticleLoops::loop(a_soa, [&](std::size_t i) {
      w[i] = 0.5 * (x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
    });
  }

  // ---------------------------------------------------------------------------
  // (7) VECTOR-FIELD INTERPOLATION
  // ---------------------------------------------------------------------------
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
  interpVecSoA(const FArrayBox& a_rho, MoverSoA& a_soa, const RealVect& a_lo, const RealVect& a_iDx)
  {
    const ParticleReal* px = a_soa.positionColumn(0);
    const ParticleReal* py = a_soa.positionColumn(1);
    const ParticleReal* pz = a_soa.positionColumn(2);
    ParticleReal*       vx = a_soa.column<&MoverPayload::vx>();
    ParticleReal*       vy = a_soa.column<&MoverPayload::vy>();
    ParticleReal*       vz = a_soa.column<&MoverPayload::vz>();
    const std::size_t   n  = a_soa.size();
    for (std::size_t i = 0; i < n; i++) {
      const RealVect v = interpolateVecFieldCIC(a_rho, RealVect(D_DECL(px[i], py[i], pz[i])), a_lo, a_iDx);
      D_TERM(vx[i] = v[0];, vy[i] = v[1];, vz[i] = v[2];);
    }
  }

  // ---------------------------------------------------------------------------
  // (8) EULER ADVANCE  x += v*dt
  // ---------------------------------------------------------------------------
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
    ParticleReal*       px = a_soa.positionColumn(0);
    ParticleReal*       py = a_soa.positionColumn(1);
    ParticleReal*       pz = a_soa.positionColumn(2);
    const ParticleReal* vx = a_soa.column<&MoverPayload::vx>();
    const ParticleReal* vy = a_soa.column<&MoverPayload::vy>();
    const ParticleReal* vz = a_soa.column<&MoverPayload::vz>();
    ParticleLoops::loop(a_soa, [&](std::size_t i) {
      px[i] += vx[i] * a_dt;
      py[i] += vy[i] * a_dt;
      pz[i] += vz[i] * a_dt;
    });
  }

  // ---------------------------------------------------------------------------
  // (9) KD-tree merge (split by position; leaf merge gathers whole particles).
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

  /** @brief Reconstruct a whole MergeParticle from the SoA (accessors + payload gather). */
  inline MergeParticle
  soaToMerge(const MergeSoA& a_soa, const int a_idx)
  {
    const MergePayload pl = a_soa.gather(a_idx);
    MergeParticle      p;
    p.m_position  = a_soa.position(a_idx);
    p.m_weight    = a_soa.weight(a_idx);
    p.m_velocity  = pl.velocity;
    p.m_oldPos    = pl.oldPos;
    p.m_mobility  = pl.mobility;
    p.m_diffusion = pl.diffusion;
    p.m_energy    = pl.energy;
    return p;
  }

  MergeSoA
  mergeSoAKD(MergeSoA& a_soa, const int a_leafSize)
  {
    const int        n = static_cast<int>(a_soa.size());
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    MergeSoA out;
    out.reserve(n / a_leafSize + 1);
    kdMerge(
      idx,
      0,
      n,
      0,
      a_leafSize,
      [&](int a, int d) {
        return a_soa.positionColumn(d)[a];
      }, // per-component position only (no RealVect construction)
      [&](int lo, int hi) {
        const MergeParticle m = mergeRange(lo, hi, [&](int k) {
          return soaToMerge(a_soa, idx[k]);
        });
        const MergePayload  pl{m.m_velocity, m.m_oldPos, m.m_mobility, m.m_diffusion, m.m_energy};
        out.append(m.m_position, m.m_weight, pl);
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
  // MPI packing helpers.
  // ---------------------------------------------------------------------------
  /** @brief Pack a List per particle (pointer-chased; the realistic Chombo path). */
  void
  packList(unsigned char* a_buf, const List<BenchParticle>& a_p)
  {
    unsigned char* q = a_buf;
    for (ListIterator<BenchParticle> lit(a_p); lit.ok(); ++lit) {
      std::memcpy(q, &lit(), sizeof(BenchParticle));
      q += sizeof(BenchParticle);
    }
  }

  /** @brief Pack a contiguous AoS vector in a single bulk memcpy. */
  void
  packVector(unsigned char* a_buf, const std::vector<BenchParticle>& a_v)
  {
    std::memcpy(a_buf, a_v.data(), a_v.size() * sizeof(BenchParticle));
  }

  /** @brief Pack one SoA column [colPtr, colPtr + size*sizeof) and advance q. */
  template <typename SoAT, std::size_t K>
  inline void
  packOneColumn(unsigned char*& a_q, const SoAT& a_soa)
  {
    using T                 = typename SoAT::template ColumnType<K>;
    const std::size_t bytes = a_soa.size() * sizeof(T);
    std::memcpy(a_q, a_soa.template columnByIndex<K>(), bytes);
    a_q += bytes;
  }
  template <typename SoAT, std::size_t... K>
  inline void
  packPerColumnImpl(unsigned char* a_buf, const SoAT& a_soa, std::index_sequence<K...>)
  {
    unsigned char* q = a_buf;
    using expander   = int[];
    (void)expander{0, ((void)packOneColumn<SoAT, K>(q, a_soa), 0)...};
  }
  /** @brief Pack the SoA as N per-column bulk memcpies (the historically bimodal-slow path). */
  void
  packSoAPerColumn(unsigned char* a_buf, const BenchSoA& a_soa)
  {
    packPerColumnImpl(a_buf, a_soa, std::make_index_sequence<BenchSoA::s_numColumns>{});
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
    // ----- Configuration ------------------------------------------------------
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

    // ----- Particles in FArrayBox (Fortran) order, plus a randomized copy ------
    std::mt19937                         rng(20260620u);
    std::uniform_real_distribution<Real> u(0.0, 1.0);

    std::vector<BenchParticle> sortedParticles;
    sortedParticles.reserve(nPart);
    for (int k = 0; k < nCell; k++) {
      for (int j = 0; j < nCell; j++) {
        for (int i = 0; i < nCell; i++) {
          for (int p = 0; p < ppc; p++) {
            BenchParticle bp;
            bp.m_position = probLo + RealVect(D_DECL(i + u(rng), j + u(rng), k + u(rng)));
            bp.m_weight   = 1.0;
            sortedParticles.push_back(bp);
          }
        }
      }
    }
    std::vector<BenchParticle> randomParticles = sortedParticles;
    std::shuffle(randomParticles.begin(), randomParticles.end(), rng);

    // ----- Container builders --------------------------------------------------
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
      BenchSoA soa;
      soa.reserve(a_src.size());
      for (const BenchParticle& p : a_src) {
        soa.append(p.m_position, p.m_weight);
      }
      return soa;
    };

    List<BenchParticle> listSorted     = buildCompactList(sortedParticles);
    List<BenchParticle> listRandom     = buildCompactList(randomParticles);
    List<BenchParticle> listFragSorted = buildFragmentedList(sortedParticles);
    List<BenchParticle> listFragRandom = buildFragmentedList(randomParticles);
    BenchSoA            soaSorted      = buildSoA(sortedParticles);
    BenchSoA            soaRandom      = buildSoA(randomParticles);

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
    const double dFragS = timeDeposit([&]() {
      depositList(rho, listFragSorted, probLo, invDx, invVol);
    });
    const double dFragR = timeDeposit([&]() {
      depositList(rho, listFragRandom, probLo, invDx, invVol);
    });
    const double dVecS  = timeDeposit([&]() {
      depositVector(rho, sortedParticles, probLo, invDx, invVol);
    });
    const double dVecR  = timeDeposit([&]() {
      depositVector(rho, randomParticles, probLo, invDx, invVol);
    });
    const double dSoAS  = timeDeposit([&]() {
      depositSoA(rho, soaSorted, probLo, invDx, invVol);
    });
    const double dSoAR  = timeDeposit([&]() {
      depositSoA(rho, soaRandom, probLo, invDx, invVol);
    });

    // =====================================================================
    // (2) INTERPOLATION
    // =====================================================================
    rho.setVal(0.0);
    depositSoA(rho, soaSorted, probLo, invDx, invVol); // some non-trivial field to gather

    interpolateList(rho, listSorted, probLo, invDx);
    interpolateSoA(rho, soaSorted, probLo, invDx);

    const double iFragS = timeOp(fastReps, [&]() {
      interpolateList(rho, listFragSorted, probLo, invDx);
    });
    const double iFragR = timeOp(fastReps, [&]() {
      interpolateList(rho, listFragRandom, probLo, invDx);
    });
    const double iVecS  = timeOp(fastReps, [&]() {
      interpolateVector(rho, sortedParticles, probLo, invDx);
    });
    const double iVecR  = timeOp(fastReps, [&]() {
      interpolateVector(rho, randomParticles, probLo, invDx);
    });
    const double iSoAS  = timeOp(fastReps, [&]() {
      interpolateSoA(rho, soaSorted, probLo, invDx);
    });
    const double iSoAR  = timeOp(fastReps, [&]() {
      interpolateSoA(rho, soaRandom, probLo, invDx);
    });

    // =====================================================================
    // (3) STREAMING TRANSFORM
    // =====================================================================
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

    // =====================================================================
    // (4) BUILD
    // =====================================================================
    const double bList    = timeOp(buildReps, [&]() {
      List<BenchParticle> l;
      for (const BenchParticle& p : sortedParticles) {
        l.append(p);
      }
      g_sink += l.isEmpty() ? 0u : 1u;
    });
    const double bVec     = timeOp(buildReps, [&]() {
      std::vector<BenchParticle> v;
      v.reserve(nPart);
      for (const BenchParticle& p : sortedParticles) {
        v.push_back(p);
      }
      g_sink += v.size();
    });
    const double bSoA     = timeOp(buildReps, [&]() {
      BenchSoA s;
      s.reserve(nPart);
      for (const BenchParticle& p : sortedParticles) {
        s.append(p.m_position, p.m_weight);
      }
      g_sink += s.size();
    });
    const double bSoAGrow = timeOp(buildReps, [&]() {
      BenchSoA s; // no reserve -> exercises arena growth (whole-buffer realloc per doubling)
      for (const BenchParticle& p : sortedParticles) {
        s.append(p.m_position, p.m_weight);
      }
      g_sink += s.size();
    });

    // =====================================================================
    // (5) REMAP (scatter particles into nBuckets destination containers)
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
    const double rVec  = timeOp(buildReps, [&]() {
      std::vector<std::vector<BenchParticle>> dst(nBuckets);
      for (std::size_t i = 0; i < nPart; i++) {
        dst[bucketId[i]].push_back(sortedParticles[i]);
      }
      for (const std::vector<BenchParticle>& d : dst) {
        g_sink += d.size();
      }
    });
    const double rSoA  = timeOp(buildReps, [&]() {
      std::vector<BenchSoA> dst(nBuckets);
      for (std::size_t i = 0; i < nPart; i++) {
        dst[bucketId[i]].append(sortedParticles[i].m_position, sortedParticles[i].m_weight);
      }
      for (const BenchSoA& d : dst) {
        g_sink += d.size();
      }
    });
    // Two-pass remap: count per bucket, reserve, then fill (no reallocation churn).
    const double rSoA2 = timeOp(buildReps, [&]() {
      std::vector<std::size_t> counts(nBuckets, 0);
      for (std::size_t i = 0; i < nPart; i++) {
        counts[bucketId[i]]++;
      }
      std::vector<BenchSoA> dst(nBuckets);
      for (int b = 0; b < nBuckets; b++) {
        dst[b].reserve(counts[b]);
      }
      for (std::size_t i = 0; i < nPart; i++) {
        dst[bucketId[i]].append(sortedParticles[i].m_position, sortedParticles[i].m_weight);
      }
      for (const BenchSoA& d : dst) {
        g_sink += d.size();
      }
    });

    // =====================================================================
    // (6) MPI PACKING (serialize a whole container to a byte send buffer)
    // =====================================================================
    const std::size_t bufBytes = std::max(nPart * sizeof(BenchParticle), soaSorted.byteSpan());
    const std::size_t bufAlloc = ((bufBytes + 63) / 64) * 64;
    unsigned char*    packBuf  = static_cast<unsigned char*>(std::aligned_alloc(64, bufAlloc));

    const double pList   = timeOp(fastReps, [&]() {
      packList(packBuf, listSorted);
      g_sink += packBuf[0];
    });
    const double pVec    = timeOp(fastReps, [&]() {
      packVector(packBuf, sortedParticles);
      g_sink += packBuf[0];
    });
    const double pSoACol = timeOp(fastReps, [&]() {
      packSoAPerColumn(packBuf, soaSorted);
      g_sink += packBuf[0];
    });
    // Arena: the whole compact container is ONE contiguous span -> a single memcpy of data()
    // (glibc auto-streams; no intrinsics). A real whole-container MPI send skips even this copy.
    const std::size_t arenaSpan = soaSorted.byteSpan();
    const double      pArena    = timeOp(fastReps, [&]() {
      std::memcpy(packBuf, soaSorted.data(), arenaSpan);
      g_sink += packBuf[0];
    });

    // =====================================================================
    // (10) CONTAINER COPY (copy whole container -> fresh container)
    // =====================================================================
    const double cList = timeOp(buildReps, [&]() {
      List<BenchParticle> d;
      for (ListIterator<BenchParticle> lit(listSorted); lit.ok(); ++lit) {
        d.append(lit());
      }
      g_sink += d.isEmpty() ? 0u : 1u;
    });
    const double cVec  = timeOp(buildReps, [&]() {
      const std::vector<BenchParticle> d(sortedParticles); // bulk copy ctor
      g_sink += d.size();
    });
    const double cSoA  = timeOp(buildReps, [&]() {
      BenchSoA d;
      d.resize(nPart);                                               // allocate the destination arena
      std::memcpy(d.data(), soaSorted.data(), soaSorted.byteSpan()); // one bulk copy (same offsets)
      g_sink += static_cast<const unsigned char*>(d.data())[0];      // force the copy (defeat DCE)
    });

    std::free(packBuf);

    // =====================================================================
    // Mover containers (position + velocity), shared by (7) and (8).
    // =====================================================================
    const Real dt = 0.1;

    std::vector<MoverParticle> moverSrc(nPart);
    for (std::size_t i = 0; i < nPart; i++) {
      moverSrc[i].m_position = sortedParticles[i].m_position;
      moverSrc[i].m_velocity = RealVect(D_DECL(1.0, 0.5, 0.25));
    }

    List<MoverParticle> moverList;
    for (const MoverParticle& p : moverSrc) {
      moverList.append(p);
    }
    std::vector<MoverParticle> moverVec = moverSrc;
    MoverSoA                   moverSoA;
    moverSoA.reserve(nPart);
    for (const MoverParticle& p : moverSrc) {
      MoverPayload pl;
      D_TERM(pl.vx = p.m_velocity[0];, pl.vy = p.m_velocity[1];, pl.vz = p.m_velocity[2];);
      moverSoA.append(p.m_position, p.m_weight, pl);
    }

    // =====================================================================
    // (7) VECTOR-FIELD INTERPOLATION (runs while positions are inside the box)
    // =====================================================================
    FArrayBox vfield(grown, SpaceDim);
    for (int comp = 0; comp < SpaceDim; comp++) {
      vfield.setVal(static_cast<Real>(comp + 1), comp); // constant per component: CIC returns (1,2,..)
    }
    interpVecSoA(vfield, moverSoA, probLo, invDx);
    {
      const ParticleReal* vx = moverSoA.column<&MoverPayload::vx>();
      if (std::abs(vx[123] - 1.0) > 1.E-6) {
        MayDay::Abort("benchmark: vector interpolation of a constant field is wrong");
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
    // (8) EULER ADVANCE
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

    // =====================================================================
    // (9) KD-TREE MERGE
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
            MergeParticle mp;
            mp.m_position = probLo + RealVect(D_DECL(i + u(rng), j + u(rng), k + u(rng)));
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
    MergeSoA mergeSoa;
    mergeSoa.reserve(nMerge);
    for (const MergeParticle& p : mergeSrc) {
      const MergePayload pl{p.m_velocity, p.m_oldPos, p.m_mobility, p.m_diffusion, p.m_energy};
      mergeSoa.append(p.m_position, p.m_weight, pl);
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
      const MergeSoA out = mergeSoAKD(mergeSoa, leafSize);
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
    pout() << "  bytes/part : AoS struct " << sizeof(BenchParticle) << " B, SoA cols " << BenchSoA::bytesPerParticle()
           << " B" << endl;
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
    row("SoA per-component       ", tSoAS);
    pout() << "    vs compact List: vector<P> " << tListS / tVecS << "x, SoA " << tListS / tSoAS
           << "x; SoA vs vector<P>: " << tVecS / tSoAS << "x" << endl;

    pout() << "  (4) BUILD (allocation)" << endl;
    row("List              ", bList);
    row("vector<P>         ", bVec);
    row("SoA arena (reserve)", bSoA);
    row("SoA arena (grow)  ", bSoAGrow);
    pout() << "    vs List: vector<P> " << bList / bVec << "x, SoA " << bList / bSoA << "x, SoA-grow "
           << bList / bSoAGrow << "x" << endl;

    pout() << "  (5) REMAP (scatter into " << nBuckets << " buckets)" << endl;
    row("List              ", rList);
    row("vector<P>         ", rVec);
    row("SoA (naive append)", rSoA);
    row("SoA (two-pass)    ", rSoA2);
    pout() << "    vs List: vector<P> " << rList / rVec << "x, SoA-2pass " << rList / rSoA2
           << "x; SoA-2pass vs vector<P>: " << rVec / rSoA2 << "x" << endl;

    pout() << "  (6) MPI PACKING (serialize whole container to send buffer)" << endl;
    row("List      (per-particle)     ", pList);
    row("vector<P> (1 bulk memcpy)    ", pVec);
    row("SoA       (per-column memcpy)", pSoACol);
    row("SoA arena (1 memcpy data())  ", pArena);
    pout() << "    vs List: vector<P> " << pList / pVec << "x, SoA-percol " << pList / pSoACol << "x, arena "
           << pList / pArena << "x; arena vs per-column SoA: " << pSoACol / pArena << "x" << endl;

    pout() << "  (7) VECTOR-FIELD INTERPOLATION (3-comp FArrayBox -> particle velocity)" << endl;
    row("List     ", viList);
    row("vector<P>", viVec);
    row("SoA      ", viSoA);
    pout() << "    vs List: vector<P> " << viList / viVec << "x, SoA " << viList / viSoA << "x" << endl;

    pout() << "  (8) EULER ADVANCE (x += v*dt)" << endl;
    row("List              ", aList);
    row("vector<P> AoS     ", aVec);
    row("SoA per-component ", aSoA);
    pout() << "    vs List: vector<P> " << aList / aVec << "x, SoA " << aList / aSoA
           << "x; SoA vs vector<P>: " << aVec / aSoA << "x" << endl;

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

    pout() << "  (10) CONTAINER COPY (copy whole container -> fresh container)" << endl;
    row("List               ", cList);
    row("vector<P>          ", cVec);
    row("SoA arena (1 memcpy)", cSoA);
    pout() << "    vs List: vector<P> " << cList / cVec << "x, SoA-arena " << cList / cSoA
           << "x; SoA-arena vs vector<P>: " << cVec / cSoA << "x" << endl;

    pout() << "  (sink=" << g_sink << ")" << endl;
  }

  return finalize();
}
