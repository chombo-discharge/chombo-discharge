/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_BucketTreeCarveParticleMergeImplem.H
 * @brief  Implementation of CD_BucketTreeCarveParticleMerge.H
 * @author Robert Marskar
 */

#ifndef CD_BUCKETTREECARVEPARTICLEMERGEIMPLEM_H
#define CD_BUCKETTREECARVEPARTICLEMERGEIMPLEM_H

// Std includes
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <utility>
#include <vector>

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>
#include <SPMD.H>
#include <DataIterator.H>

// Our includes
#include <CD_AmrMesh.H>
#include <CD_ParticleGhostMask.H>
#include <CD_BucketTreeCarveParticleMerge.H>
#include <CD_NamespaceHeader.H>

namespace ParticleManagement {

namespace detail {

inline bool
operator<(const BucketBoxKey& a_lhs, const BucketBoxKey& a_rhs) noexcept
{
  if (a_lhs.volume != a_rhs.volume) {
    return a_lhs.volume < a_rhs.volume;
  }

  return a_lhs.anchor < a_rhs.anchor;
}

template <typename T>
inline std::vector<T>
bucketTreeExchangeByRank(const std::vector<std::vector<T>>& a_sendByRank)
{
  const int numRanks = static_cast<int>(a_sendByRank.size());

#ifdef CH_MPI
  if (numRanks <= 1) {
    return a_sendByRank.empty() ? std::vector<T>() : a_sendByRank[0];
  }

  std::vector<int> sendCounts(numRanks, 0);

  for (int r = 0; r < numRanks; r++) {
    sendCounts[r] = static_cast<int>(a_sendByRank[r].size());
  }

  std::vector<int> recvCounts(numRanks, 0);
  MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, Chombo_MPI::comm);

  std::vector<int> sdispl(numRanks, 0);
  std::vector<int> rdispl(numRanks, 0);

  long stot = 0;
  long rtot = 0;

  for (int r = 0; r < numRanks; r++) {
    sdispl[r] = static_cast<int>(stot);
    stot += sendCounts[r];

    rdispl[r] = static_cast<int>(rtot);
    rtot += recvCounts[r];
  }

  std::vector<T> sflat(stot);

  for (int r = 0; r < numRanks; r++) {
    if (sendCounts[r] > 0) {
      std::copy(a_sendByRank[r].begin(), a_sendByRank[r].end(), sflat.begin() + sdispl[r]);
    }
  }

  std::vector<T> rflat(rtot);

  // Scale record counts/displacements up to bytes -- MPI only knows bytes.
  std::vector<int> sendBytes(numRanks);
  std::vector<int> recvBytes(numRanks);
  std::vector<int> sdisplBytes(numRanks);
  std::vector<int> rdisplBytes(numRanks);

  for (int r = 0; r < numRanks; r++) {
    sendBytes[r]   = sendCounts[r] * static_cast<int>(sizeof(T));
    recvBytes[r]   = recvCounts[r] * static_cast<int>(sizeof(T));
    sdisplBytes[r] = sdispl[r] * static_cast<int>(sizeof(T));
    rdisplBytes[r] = rdispl[r] * static_cast<int>(sizeof(T));
  }

  MPI_Alltoallv(sflat.data(),
                sendBytes.data(),
                sdisplBytes.data(),
                MPI_BYTE,
                rflat.data(),
                recvBytes.data(),
                rdisplBytes.data(),
                MPI_BYTE,
                Chombo_MPI::comm);

  return rflat;
#else
  return a_sendByRank.empty() ? std::vector<T>() : a_sendByRank[0];
#endif
}

/**
 * @brief Axis-aligned bounding box of particles[a_lo, a_hi).
 * @param[in] a_particles Particle buffer.
 * @param[in] a_lo Start of the range (inclusive).
 * @param[in] a_hi End of the range (exclusive).
 * @param[out] a_boxLo Low corner.
 * @param[out] a_boxHi High corner.
 */
template <typename Packed>
inline void
bucketTreeBBox(const std::vector<BucketParticle<Packed>>& a_particles,
               const std::size_t                          a_lo,
               const std::size_t                          a_hi,
               RealVect&                                  a_boxLo,
               RealVect&                                  a_boxHi) noexcept
{
  a_boxLo = a_particles[a_lo].position;
  a_boxHi = a_particles[a_lo].position;

  for (std::size_t i = a_lo + 1; i < a_hi; i++) {
    const RealVect& pos = a_particles[i].position;

    for (int dir = 0; dir < SpaceDim; dir++) {
      a_boxLo[dir] = std::min(a_boxLo[dir], pos[dir]);
      a_boxHi[dir] = std::max(a_boxHi[dir], pos[dir]);
    }
  }
}

/**
 * @brief Largest per-axis leaf extent, in cell widths, that may still be merged.
 * @details buildBucketTreeLeaves() and mergeBucketTreeCarve() must use the same value: the former
 * sets each cell's leaf budget assuming a leaf merges to one particle, the latter vetoes leaves
 * exceeding this bound. A disagreement spends budget on leaves that are then rejected.
 */
static constexpr Real s_bucketTreeMaxLeafExtent = 1.0;

/**
 * @brief Largest per-axis extent of [a_boxLo,a_boxHi], expressed as a fraction of that axis's own
 * cell width.
 * @details The correctness backstop buildBucketTreeLeaves()/mergeBucketTreeCarve() rely on to
 * decide whether a leaf is compact enough to merge: unlike an aggregate VOLUME ratio (a product
 * across axes), this cannot be fooled by an anisotropic (long-and-thin) bounding box -- e.g. a box
 * 5 cells long in one axis and 0.2 cells wide in another has volume ratio 1.0 (deceptively "one
 * cell's worth") but a maxAxisSpan of 5.0 (correctly flagged as still spanning several cells).
 * Particles hugging a curved EB surface, or funneled along a narrow L-shaped coarse-fine strip,
 * produce exactly this kind of elongated cluster.
 * @param[in] a_boxLo Low corner.
 * @param[in] a_boxHi High corner.
 * @param[in] a_dx    Cell spacing.
 * @return max_dir (a_boxHi[dir] - a_boxLo[dir]) / a_dx[dir].
 */
inline Real
bucketTreeMaxAxisSpan(const RealVect& a_boxLo, const RealVect& a_boxHi, const RealVect& a_dx) noexcept
{
  Real maxAxisSpan = 0.0;

  for (int dir = 0; dir < SpaceDim; dir++) {
    maxAxisSpan = std::max(maxAxisSpan, (a_boxHi[dir] - a_boxLo[dir]) / a_dx[dir]);
  }

  return maxAxisSpan;
}

/**
 * @brief The unclamped, position-derived cell index a physical position falls in.
 * @details Same convention as the nearest-neighbor merge's cellKeyOf() -- not clamped to any
 * patch's box, so positions in different patches' local frames stay directly comparable.
 * @param[in] a_position Physical position.
 * @param[in] a_probLo Physical position of the domain's low corner.
 * @param[in] a_dx Cell spacing.
 * @return floor((a_position - a_probLo) / a_dx), component-wise.
 */
inline IntVect
bucketTreeCellKeyOf(const RealVect& a_position, const RealVect& a_probLo, const RealVect& a_dx) noexcept
{
  IntVect iv;

  for (int dir = 0; dir < SpaceDim; dir++) {
    iv[dir] = static_cast<int>(std::floor((a_position[dir] - a_probLo[dir]) / a_dx[dir]));
  }

  return iv;
}

/**
 * @brief Strict weak ordering over IntVect, lexicographic component-wise.
 * @details Used to sort/binary-search the per-patch, per-cell TRUE particle-count histogram (see
 * buildBucketTreeLeaves()'s own docs) -- IntVect has no built-in total order.
 * @param[in] a_lhs Left-hand cell index.
 * @param[in] a_rhs Right-hand cell index.
 * @return True iff a_lhs sorts strictly before a_rhs.
 */
inline bool
bucketTreeIntVectLess(const IntVect& a_lhs, const IntVect& a_rhs) noexcept
{
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (a_lhs[dir] != a_rhs[dir]) {
      return a_lhs[dir] < a_rhs[dir];
    }
  }

  return false;
}

/**
 * @brief Split particles[a_lo,a_hi) in place by the longest axis at the count-median -- the sole
 * split rule used by buildBucketTreeLeaves().
 * @details The split point is the position-sorted median particle (via nth_element), balancing
 * COUNT as evenly as possible -- mirrors equal_weight_kd's own single-cell logic, but by count,
 * never by weight, and with NO synthetic particle splitting (a particle is never divided across
 * the two children). nth_element divides the range in half by count however the positions are
 * distributed, so a tightly clustered node still splits evenly.
 * @param[in,out] a_particles Particle buffer; reordered in [a_lo,a_hi).
 * @param[in] a_lo Start of the range (inclusive).
 * @param[in] a_hi End of the range (exclusive); a_hi - a_lo must be >= 2.
 * @return The split index: [a_lo,mid) and [mid,a_hi) are the two children, each non-empty.
 */
template <typename Packed>
inline std::size_t
bucketTreeSplitCountMedian(std::vector<BucketParticle<Packed>>& a_particles,
                           const std::size_t                    a_lo,
                           const std::size_t                    a_hi) noexcept
{
  RealVect boxLo, boxHi;
  bucketTreeBBox(a_particles, a_lo, a_hi, boxLo, boxHi);

  const int axis = (boxHi - boxLo).maxDir(true);

  const std::size_t mid = a_lo + (a_hi - a_lo) / 2;

  std::nth_element(a_particles.begin() + a_lo,
                   a_particles.begin() + mid,
                   a_particles.begin() + a_hi,
                   [axis](const BucketParticle<Packed>& a_p1, const BucketParticle<Packed>& a_p2) noexcept -> bool {
                     return a_p1.position[axis] < a_p2.position[axis];
                   });

  return mid;
}

/**
 * @brief Split particles[a_lo,a_hi) in place by the longest axis at the WEIGHT median -- the plane
 * with half this node's total particle weight on either side.
 * @details Used instead of the count median so that the resulting leaves carry roughly equal total
 * weight, and therefore so do the super-particles they collapse into. Uniform super-particle weight
 * is what keeps statistical noise down, and is the defining property of the equal_weight_kd merge.
 *
 * Unlike equal_weight_kd this NEVER divides a particle synthetically, so exact equality is not
 * attainable: once a single particle's weight exceeds a leaf's target weight, that particle alone
 * sets a floor no partition can get under. The bound equal_weight_kd guarantees (max/min <= 2) is
 * bought precisely by that synthetic division and is therefore out of reach here by construction.
 *
 * A weight median can land very unevenly BY COUNT -- one dominant particle carrying half the node's
 * weight forces a 1-vs-(n-1) split. That is harmless here only because buildBucketTreeLeaves()
 * enforces a hard per-cell leaf quota independently of how the split plane behaves: an uneven split
 * can cost merge quality but cannot leave a cell above target. Under any of the earlier,
 * quota-less formulations the same degeneracy was fatal.
 * @param[in,out] a_particles Particle buffer; reordered in [a_lo,a_hi).
 * @param[in] a_lo Start of the range (inclusive).
 * @param[in] a_hi End of the range (exclusive); a_hi - a_lo must be >= 2.
 * @return The split index: [a_lo,mid) and [mid,a_hi) are the two children, each non-empty.
 */
template <typename Packed>
inline std::size_t
bucketTreeSplitWeightMedian(std::vector<BucketParticle<Packed>>& a_particles,
                            const std::size_t                    a_lo,
                            const std::size_t                    a_hi) noexcept
{
  RealVect boxLo, boxHi;
  bucketTreeBBox(a_particles, a_lo, a_hi, boxLo, boxHi);

  const int axis = (boxHi - boxLo).maxDir(true);

  std::sort(a_particles.begin() + a_lo,
            a_particles.begin() + a_hi,
            [axis](const BucketParticle<Packed>& a_p1, const BucketParticle<Packed>& a_p2) noexcept -> bool {
              return a_p1.position[axis] < a_p2.position[axis];
            });

  Real totalWeight = 0.0;

  for (std::size_t idx = a_lo; idx < a_hi; idx++) {
    totalWeight += a_particles[idx].weight;
  }

  const Real half = 0.5 * totalWeight;

  Real        acc = 0.0;
  std::size_t mid = a_lo + 1;

  for (std::size_t idx = a_lo; idx < a_hi; idx++) {
    acc += a_particles[idx].weight;

    if (acc >= half) {
      mid = idx + 1;

      break;
    }
  }

  return std::max(a_lo + 1, std::min(a_hi - 1, mid));
}

template <typename Packed, typename CellTrueCount>
inline void
buildBucketTreeLeaves(std::vector<BucketParticle<Packed>>& a_particles,
                      const int                            a_ppc,
                      const Real                           a_weightSplitScale,
                      const RealVect&                      a_dx,
                      const RealVect&                      a_probLo,
                      const CellTrueCount&                 a_cellTrueCount,
                      std::vector<BucketLeaf>&             a_leaves) noexcept
{
  CH_TIME("ParticleManagement::buildBucketTreeLeaves");

  a_leaves.clear();

  if (a_particles.empty() || a_ppc <= 0) {
    return;
  }

  // Partition into the particles that may be merged (those in TRUE-crowded cells) and those that
  // may not, moving the mergeable ones to the front so every leaf below is still a contiguous range.
  // Particles in uncrowded cells are emitted as singleton leaves: the merge threshold is 2 members,
  // so a singleton is never merged, and an uncrowded cell can therefore never be drained below
  // target by a neighbour's merge.
  const auto firstUnmergeable = std::stable_partition(
    a_particles.begin(),
    a_particles.end(),
    [&](const BucketParticle<Packed>& a_p) {
      return a_cellTrueCount(bucketTreeCellKeyOf(a_p.position, a_probLo, a_dx)) > a_ppc;
    });

  const std::size_t numMergeable = static_cast<std::size_t>(firstUnmergeable - a_particles.begin());

  for (std::size_t idx = numMergeable; idx < a_particles.size(); idx++) {
    RealVect boxLo, boxHi;
    bucketTreeBBox(a_particles, idx, idx + 1, boxLo, boxHi);

    a_leaves.push_back(BucketLeaf{idx, idx + 1, boxLo, boxHi});
  }

  if (numMergeable == 0) {
    return;
  }

  // LIVE PER-CELL LEAF BUDGET, enforced as a hard constraint at split time rather than computed
  // per node and hoped for afterwards.
  //
  // Every leaf becomes exactly ONE super-particle at its own weighted centroid, so the number of
  // leaves whose centroid falls in a cell is that cell's post-merge population. Track that count
  // live and refuse any split that would push a cell past a_ppc. The per-cell ceiling is a shared
  // resource and cannot be expressed by any test a node applies to itself alone.
  //
  // Splits are handed out heaviest-first: the quota is claimed first-come-first-served, so a
  // lighter ordering lets light leaves take a cell's slots before a heavy leaf is considered, and
  // the heavy one is then refused and retired whole as a single very heavy super-particle.
  auto centroidCellOf = [&](const std::size_t a_lo, const std::size_t a_hi) -> IntVect {
    RealVect centroid = RealVect::Zero;
    Real     totalW   = 0.0;

    for (std::size_t idx = a_lo; idx < a_hi; idx++) {
      centroid += a_particles[idx].weight * a_particles[idx].position;
      totalW += a_particles[idx].weight;
    }

    if (totalW > 0.0) {
      centroid /= totalW;
    }
    else {
      centroid = a_particles[a_lo].position;
    }

    return bucketTreeCellKeyOf(centroid, a_probLo, a_dx);
  };

  struct BucketCand
  {
    Real        weight;
    std::size_t count;
    std::size_t lo;
    std::size_t hi;
    IntVect     cell;
  };

  // Hand out each cell's finite quota to its HEAVIEST candidate leaf first, not its most populous.
  // The quota is a scarce resource claimed first-come-first-served, so ordering by count lets light
  // leaves consume a cell's slots before a heavy leaf is ever considered -- the heavy one is then
  // refused and retired whole, becoming a single very heavy super-particle sitting next to the light
  // ones that took its slots. Ordering by weight makes the leaves that most need subdividing claim
  // slots first.
  const auto byWeight = [](const BucketCand& a_lhs, const BucketCand& a_rhs) noexcept -> bool {
    return a_lhs.weight < a_rhs.weight;
  };

  auto rangeWeight = [&](const std::size_t a_lo, const std::size_t a_hi) -> Real {
    Real w = 0.0;

    for (std::size_t idx = a_lo; idx < a_hi; idx++) {
      w += a_particles[idx].weight;
    }

    return w;
  };

  // Whether [a_lo,a_hi) is still too wide to be merged at all. mergeBucketTreeCarve() vetoes any
  // final leaf whose largest per-axis extent exceeds one cell width and leaves every one of its
  // members behind as an individual particle, so such a node does NOT cost its cell one particle --
  // it costs its full member count.
  auto maxAxisSpanOf = [&](const std::size_t a_lo, const std::size_t a_hi) -> Real {
    RealVect boxLo, boxHi;
    bucketTreeBBox(a_particles, a_lo, a_hi, boxLo, boxHi);

    return bucketTreeMaxAxisSpan(boxLo, boxHi, a_dx);
  };

  // Whether [a_lo,a_hi) has shrunk to the scale at which the weight median takes over from the
  // count median (see the split site below for why the crossover exists and where it sits).
  const auto cellLess = [](const IntVect& a_lhs, const IntVect& a_rhs) noexcept -> bool {
    return bucketTreeIntVectLess(a_lhs, a_rhs);
  };

  std::map<IntVect, int, decltype(cellLess)> used(cellLess);

  std::vector<BucketCand> heap;
  std::vector<BucketCand> finalLeaves;

  {
    const IntVect rootCell = centroidCellOf(0, numMergeable);

    used[rootCell]++;
    heap.push_back(BucketCand{rangeWeight(0, numMergeable), numMergeable, 0, numMergeable, rootCell});
  }

  while (!heap.empty()) {
    std::pop_heap(heap.begin(), heap.end(), byWeight);

    const BucketCand node = heap.back();
    heap.pop_back();

    if (node.count < 2) {
      finalLeaves.push_back(node);

      continue;
    }

    // The split itself only reorders within [lo,hi), which is harmless even if the split is then
    // refused -- the range stays one contiguous leaf either way.
    // While a node spans several cells a split apportions leaves BETWEEN cells, which is a question
    // of counts; once it is narrow the weight median is what drives the super-particles toward equal
    // weight. This is an extent test, not a "both corners in one cell" test, so a narrow node
    // straddling a face still qualifies.
    const Real nodeSpan = maxAxisSpanOf(node.lo, node.hi);

    const std::size_t mid = (nodeSpan <= a_weightSplitScale)
                              ? bucketTreeSplitWeightMedian(a_particles, node.lo, node.hi)
                              : bucketTreeSplitCountMedian(a_particles, node.lo, node.hi);

    const IntVect cellLeft  = centroidCellOf(node.lo, mid);
    const IntVect cellRight = centroidCellOf(mid, node.hi);

    used[node.cell]--;
    used[cellLeft]++;
    used[cellRight]++;

    // The quota constrains only leaves that can actually merge. A node wider than
    // s_bucketTreeMaxLeafExtent is vetoed by the caller and leaves every member behind
    // individually, so it costs its cell its full member count rather than one particle; splitting
    // it further lowers that cell's final population even though it raises the leaf count.
    if (nodeSpan <= s_bucketTreeMaxLeafExtent && (used[cellLeft] > a_ppc || used[cellRight] > a_ppc)) {
      used[cellLeft]--;
      used[cellRight]--;
      used[node.cell]++;

      finalLeaves.push_back(node);

      continue;
    }

    heap.push_back(BucketCand{rangeWeight(node.lo, mid), mid - node.lo, node.lo, mid, cellLeft});
    std::push_heap(heap.begin(), heap.end(), byWeight);

    heap.push_back(BucketCand{rangeWeight(mid, node.hi), node.hi - mid, mid, node.hi, cellRight});
    std::push_heap(heap.begin(), heap.end(), byWeight);
  }

  for (const BucketCand& node : finalLeaves) {
    RealVect boxLo, boxHi;
    bucketTreeBBox(a_particles, node.lo, node.hi, boxLo, boxHi);

    a_leaves.push_back(BucketLeaf{node.lo, node.hi, boxLo, boxHi});
  }
}

/**
 * @brief Whether a cell has any outgoing ghost-shipment target -- a duplicate of the nearest-
 * neighbor merge's isBoundaryExposed() (a one-line function, not worth an EBGeometry-pulling
 * dependency to reuse).
 * @param[in] a_cellIndex Cell index in the mask's own frame.
 * @param[in] a_ghostMask Outgoing particle ghost mask.
 * @return True iff the cell has at least one outgoing ghost-shipment target.
 */
inline bool
bucketTreeIsBoundaryExposed(const IntVect& a_cellIndex, const ParticleGhostMask& a_ghostMask) noexcept
{
  return a_ghostMask.numTargets(a_cellIndex) > 0;
}

/**
 * @brief Debug-only sanity check on a merged particle's position: a weighted centroid of positions
 * that all lie within [a_boxLo,a_boxHi], using strictly positive weights, is mathematically
 * guaranteed to itself lie within that same box.
 * @details CH_assert this explicitly (rather than trust the arithmetic silently) so that a stray
 * sign error, a mis-set/negative weight, or a wrong id feeding the wrong position into the
 * accumulation surfaces immediately as an assertion in a DEBUG build, instead of silently producing
 * a merged superparticle placed somewhere non-physical. A relative tolerance accounts for floating-
 * point round-off in the summation, not for any real violation.
 * @param[in] a_centroid Computed weighted centroid.
 * @param[in] a_totalWeight Sum of weights used to compute a_centroid.
 * @param[in] a_boxLo Low corner of the bounding box every contributing position lies within.
 * @param[in] a_boxHi High corner of the bounding box every contributing position lies within.
 */
inline void
bucketTreeCheckCentroid(const RealVect& a_centroid,
                        const Real      a_totalWeight,
                        const RealVect& a_boxLo,
                        const RealVect& a_boxHi) noexcept
{
#ifndef NDEBUG
  CH_assert(a_totalWeight > 0.0);

  for (int dir = 0; dir < SpaceDim; dir++) {
    const Real tol = 1.0e-9 * (std::abs(a_boxLo[dir]) + std::abs(a_boxHi[dir]) + 1.0);

    CH_assert(a_centroid[dir] >= a_boxLo[dir] - tol);
    CH_assert(a_centroid[dir] <= a_boxHi[dir] + tol);
  }
#else
  (void)a_centroid;
  (void)a_totalWeight;
  (void)a_boxLo;
  (void)a_boxHi;
#endif
}

} // namespace detail

template <typename P,
          typename Packed,
          typename Traits,
          typename Gather,
          typename Combine,
          typename Scatter,
          typename Allocator,
          typename PosValid>
inline void
mergeBucketTreeCarve(AmrMesh&                      a_amr,
                     ParticleContainer<P, Traits>& a_particles,
                     const int                     a_ppc,
                     const Real                    a_weightSplitScale,
                     const Gather&                 a_gather,
                     const Combine&                a_combine,
                     const Scatter&                a_scatter,
                     const Allocator&              a_allocateID,
                     const PosValid&               a_isPositionValid,
                     const bool                    a_enableBoundaryCarve)
{
  using namespace detail;

  CH_TIME("ParticleManagement::mergeBucketTreeCarve");

  const std::string   realm       = a_particles.getRealm();
  const int           finestLevel = a_amr.getFinestLevel();
  const RealVect      probLo      = a_amr.getProbLo();
  const Vector<Real>& dxScalar    = a_amr.getDx();

  std::vector<RealVect> dxByLevel(finestLevel + 1);

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    dxByLevel[lvl] = dxScalar[lvl] * RealVect::Unit;
  }

  const int myRank   = procID();
  const int numRanks = numProc();

  // ---- Per-patch bookkeeping shared across the whole rank ----
  struct PatchWork
  {
    int       level;
    DataIndex din;
    RealVect  dx;
    Box       box;
  };

  struct RuntimeBox
  {
    BucketBoxKey              key;
    int                       patchIdx;
    std::vector<BucketMember> members;
  };

  struct MergedResult
  {
    BucketParticle<Packed> particle;
    int                    patchIdx; // which patch's leaf produced this -- the cheap placement guess
  };

  // One entry per (locally-owned, boundary-exposed particle) x (one of MY OWN boxes that lists it),
  // appended during STEP 1 and sorted by id once STEP 1 finishes -- see selfClaims' declaration
  // below for why an id can have more than one entry, and boxIdx < 0 for what a sentinel entry means.
  struct SelfClaimEntry
  {
    ParticleID   id;
    BucketBoxKey key;    // meaningless when boxIdx < 0.
    int          boxIdx; // -1 => listen-only sentinel: no actual claim from any of my own boxes, just
                         // marks that this rank must still process this id (an unmergeable leaf's
                         // exposed member, or a boundary leaf with < 2 members) when combining with
                         // foreign claims in STEP 3.
  };

  std::vector<PatchWork> patchWork;

  // Pooled: locals + ghosts, rank-wide. One entry per (id, patch-occurrence) is appended during
  // STEP 1 -- deduplicated (sorted + collapsed, preferring the non-ghost copy) once STEP 1 finishes,
  // rather than via a per-particle hash lookup while gathering. Point lookups after that go through
  // findParticle() (binary search) below.
  std::vector<std::pair<ParticleID, BucketParticle<Packed>>> particlesByID;
  std::vector<RuntimeBox>                                    myBoxes;
  std::vector<ParticleID>                                    consumedIDs; // sorted once, before STEP 7
  std::vector<MergedResult>                                  mergedResults;

  // A particle can be listed by more than one of MY OWN boxes at once -- it is gathered once into
  // its own home patch's build (as a local member there) but can ALSO be gathered as a ghost into a
  // NEIGHBORING patch's build on this SAME rank, and that patch's build is entirely independent -- it
  // may propose its own, different box that also lists this particle. Treating "owner == myRank" as
  // "no contest possible" (collapsing to a single entry per id) is exactly the bug this vector-of-
  // entries shape exists to prevent: two of my own boxes can genuinely compete for the same particle
  // with no other rank involved at all, and collapsing the earlier entry double-counts the particle's
  // weight if both boxes go on to commit. Sorted by id once STEP 1 finishes (see below); STEP 3
  // processes it as grouped runs of equal id rather than hashing.
  std::vector<SelfClaimEntry> selfClaims;

  // Cheap pre-pass (patch sizes only, no per-particle work) so the containers below can be sized
  // once up front instead of growing repeatedly as STEP 1 fills them. combinedCountUpperBound
  // includes ghosts, so it's also a safe (if slightly generous) upper bound for the local-only ones.
  std::size_t combinedCountUpperBound = 0;

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    const DisjointBoxLayout& dbl  = a_amr.getGrids(realm)[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

    for (int mybox = 0; mybox < nbox; mybox++) {
      combinedCountUpperBound += a_particles[lvl][dit[mybox]].size();
    }
  }

  particlesByID.reserve(combinedCountUpperBound);
  selfClaims.reserve(combinedCountUpperBound);
  consumedIDs.reserve(combinedCountUpperBound);

  // ==== STEP 1: gather + build + classify + commit interior, per patch ====
  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = a_amr.getGrids(realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const AMRParticleGhostMask& maskSame = a_amr.getParticleGhostMask(realm, 1);
    const AMRParticleGhostMask& maskC2F  = a_amr.getParticleGhostMaskCoarToFine(realm, 1);
    const AMRParticleGhostMask& maskF2C  = a_amr.getParticleGhostMaskFineToCoar(realm, 1);

    const int nbox = dit.size();

    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      const ParticleGhostMask* maskSameDin = maskSame[lvl].isNull() ? nullptr : &(*maskSame[lvl])[din];
      const ParticleGhostMask* maskC2FDin  = (lvl < finestLevel && !maskC2F[lvl].isNull()) ? &(*maskC2F[lvl])[din]
                                                                                           : nullptr;
      const ParticleGhostMask* maskF2CDin  = (lvl > 0 && !maskF2C[lvl].isNull()) ? &(*maskF2C[lvl])[din] : nullptr;

      std::vector<BucketParticle<Packed>> combined;
      combined.reserve(leaf.size());

      for (std::size_t i = 0; i < leaf.size(); i++) {
        BucketParticle<Packed> p;

        p.position  = leaf.position(i);
        p.weight    = leaf.weight(i);
        p.globalID  = leaf.particleID(i);
        p.ownerRank = leaf.rankID(i);
        p.isGhost   = leaf.isGhost(i);
        p.payload   = a_gather(leaf, i);

        combined.push_back(p);

        // Append unconditionally -- deduplicated by id (preferring the non-ghost copy, per
        // BucketParticle::isGhost's own docs) in the one pass below once STEP 1 finishes gathering,
        // rather than via a per-particle hash lookup here.
        particlesByID.emplace_back(p.globalID, p);
      }

      patchWork.push_back(PatchWork{lvl, din, dxByLevel[lvl], dbl[din]});
      const int patchIdx = static_cast<int>(patchWork.size()) - 1;

      if (combined.empty()) {
        continue;
      }

      // Per-patch, per-cell TRUE particle-count histogram -- computed once, directly from every
      // particle's own position, so it is unaffected by how buildBucketTreeLeaves() later happens to
      // carve up the space. This is what makes cellTrueCount (below) a reliable ground truth rather
      // than a node-local approximation -- see buildBucketTreeLeaves()'s own docs for why the
      // difference matters. Built the same sort-and-collapse-runs way as particlesByID above.
      std::vector<std::pair<IntVect, int>> cellCounts;
      cellCounts.reserve(combined.size());

      for (const BucketParticle<Packed>& p : combined) {
        cellCounts.emplace_back(bucketTreeCellKeyOf(p.position, probLo, dxByLevel[lvl]), 1);
      }

      std::sort(cellCounts.begin(), cellCounts.end(), [](const auto& a_lhs, const auto& a_rhs) {
        return bucketTreeIntVectLess(a_lhs.first, a_rhs.first);
      });

      {
        std::size_t writeIdx = 0;

        for (std::size_t readIdx = 0; readIdx < cellCounts.size();) {
          std::size_t runEnd = readIdx + 1;
          int         total  = cellCounts[readIdx].second;

          while (runEnd < cellCounts.size() && cellCounts[runEnd].first == cellCounts[readIdx].first) {
            total += cellCounts[runEnd].second;
            runEnd++;
          }

          cellCounts[writeIdx] = std::make_pair(cellCounts[readIdx].first, total);
          writeIdx++;
          readIdx = runEnd;
        }

        cellCounts.resize(writeIdx);
      }

      auto cellTrueCount = [&cellCounts](const IntVect& a_cell) -> int {
        const auto it = std::lower_bound(cellCounts.begin(),
                                         cellCounts.end(),
                                         a_cell,
                                         [](const std::pair<IntVect, int>& a_entry, const IntVect& a_key) {
                                           return bucketTreeIntVectLess(a_entry.first, a_key);
                                         });

        return (it != cellCounts.end() && it->first == a_cell) ? it->second : 0;
      };

      std::vector<BucketLeaf> leaves;
      buildBucketTreeLeaves(combined, a_ppc, a_weightSplitScale, dxByLevel[lvl], probLo, cellTrueCount, leaves);

      // Per-patch exposure cache: leaves are about a cell in extent and not grid-aligned, so many leaves
      // (and many members within one leaf) end up querying the SAME actual cell's exposure. A
      // resident (non-ghost) member's cell always lies within this patch's OWN box (that's what
      // "resident" means), so a dense array indexed by cell offset within the box is safe -- compute
      // each distinct cell's exposure (up to 3 mask lookups) once, then reuse the cached bit for
      // every subsequent member/leaf that touches the same cell.
      const IntVect patchLo   = dbl[din].smallEnd();
      const IntVect patchSize = dbl[din].size();
      std::size_t   cacheSize = 1;
      for (int dir = 0; dir < SpaceDim; dir++) {
        cacheSize *= static_cast<std::size_t>(patchSize[dir]);
      }
      std::vector<signed char> exposureCache(cacheSize, -1); // -1 = uncached, 0 = not exposed, 1 = exposed

      auto isExposedCached = [&](const IntVect& a_iv) -> bool {
        std::size_t cellIdx = 0;
        std::size_t stride  = 1;

        for (int dir = 0; dir < SpaceDim; dir++) {
          cellIdx += static_cast<std::size_t>(a_iv[dir] - patchLo[dir]) * stride;
          stride *= static_cast<std::size_t>(patchSize[dir]);
        }

        if (exposureCache[cellIdx] < 0) {
          bool isExposed = (maskSameDin != nullptr) && bucketTreeIsBoundaryExposed(a_iv, *maskSameDin);

          if (!isExposed && maskC2FDin != nullptr) {
            isExposed = bucketTreeIsBoundaryExposed(a_iv, *maskC2FDin);
          }
          if (!isExposed && maskF2CDin != nullptr) {
            isExposed = bucketTreeIsBoundaryExposed(a_iv, *maskF2CDin);
          }

          exposureCache[cellIdx] = isExposed ? 1 : 0;
        }

        return exposureCache[cellIdx] != 0;
      };

      for (const BucketLeaf& bl : leaves) {
        std::vector<BucketMember> members;
        members.reserve(bl.hi - bl.lo);

        bool       anyLocal = false;
        ParticleID anchor   = combined[bl.lo].globalID; // running min over ALL members (bl.hi > bl.lo
                                                        // always holds for a leaf); folded into this
                                                        // loop so the boundary branch below doesn't
                                                        // need its own separate pass just to compute it.

        for (std::size_t idx = bl.lo; idx < bl.hi; idx++) {
          members.push_back(BucketMember{combined[idx].globalID, combined[idx].ownerRank});
          anchor = std::min(anchor, combined[idx].globalID);

          // "Local" here means physically resident in THIS patch, not merely owned by this rank --
          // a particle owned by this rank via a DIFFERENT one of its own patches, appearing here
          // only as a ghost, gives this leaf no standing (see BucketParticle::isGhost's own docs).
          if (!combined[idx].isGhost) {
            anyLocal = true;
          }
        }

        if (!anyLocal) {
          continue; // Nothing physically resident here -- no standing to act on it.
        }

        Real leafVolume = 1.0;

        for (int dir = 0; dir < SpaceDim; dir++) {
          leafVolume *= (bl.boxHi[dir] - bl.boxLo[dir]);
        }

        // Per-axis extent, not aggregate volume -- see bucketTreeMaxAxisSpan()'s own docs for why a
        // volume ratio alone cannot be trusted here (a long, thin leaf can have volume <=
        // one cell while still spanning several cells along one axis). leafVolume itself
        // is still needed below, unrelated to this check -- it's the box key's tie-break priority.
        const bool unmergeable = bucketTreeMaxAxisSpan(bl.boxLo, bl.boxHi, dxByLevel[lvl]) > s_bucketTreeMaxLeafExtent;

        // Per-resident-member exposure -- needed either way (leaf-wide OR for the interior/
        // boundary decision below; per-particle for the unmergeable-leaf listen-only case). Ghosts
        // are skipped entirely, not just because they have no standing (above) but because their
        // position lies outside this patch's own box -- testing a ghost's cell against THIS
        // patch's own ghost mask would violate ParticleGhostMask::numTargets()'s own precondition
        // that the cell lie within the mask's box(). Filtering on isGhost (not ownerRank) is what
        // makes this safe: a particle owned by this rank via a different patch is exactly as much
        // a ghost here, position-wise, as one owned by a different rank entirely.
        std::vector<bool> exposed(bl.hi - bl.lo, false);
        bool              anyExposed = false;
        bool              anyGhost   = false;

        for (std::size_t idx = bl.lo; idx < bl.hi; idx++) {
          if (combined[idx].isGhost) {
            anyGhost = true;

            continue;
          }

          const IntVect cell      = bucketTreeCellKeyOf(combined[idx].position, probLo, dxByLevel[lvl]);
          const bool    isExposed = isExposedCached(cell);

          exposed[idx - bl.lo] = isExposed;
          anyExposed           = anyExposed || isExposed;
        }

        if (unmergeable) {
          for (std::size_t idx = bl.lo; idx < bl.hi; idx++) {
            // exposed[] is already false for every ghost by construction (the exposure loop above
            // skips them), so this is equivalent to also checking !isGhost -- kept explicit anyway
            // as the actual invariant being relied on, not an implicit one.
            if (!combined[idx].isGhost && exposed[idx - bl.lo]) {
              selfClaims.push_back(SelfClaimEntry{combined[idx].globalID, BucketBoxKey{}, -1}); // listen only.
            }
          }

          continue;
        }

        if (!anyExposed && !anyGhost) {
          // Interior: commit immediately, zero communication, direct insert -- safe ONLY because
          // anyGhost is false, i.e. this leaf genuinely contains no foreign data at all. Explicitly
          // requiring !anyGhost (not just !anyExposed) matters whenever a leaf spans cells: the old
          // argument ("a leaf whose members are all non-exposed can never contain a ghost, since
          // ghost-fill is symmetric at width 1") only holds when a leaf's own extent is <= 1 cell --
          // every member is then automatically within ghost-reach of every other member, so a ghost's
          // presence would force at least one resident's own cell to be exposed too. Once a leaf can
          // span SEVERAL cells, a ghost near one edge and a resident several
          // cells away at the other edge can coexist in the SAME leaf without that resident itself
          // ever being within ghost-width-1 of anything -- anyExposed alone would then be false even
          // though the leaf holds a foreign particle, and committing it here would merge that
          // particle's weight locally while its true home patch independently wins the SAME
          // particle's argmin via its own boundary box, double-counting it. Falling through to the
          // boundary-candidate branch below instead (rather than committing) is what lets a real
          // self-claim be raised for that foreign member, so the argmin actually arbitrates it.
          if (bl.hi - bl.lo >= 2) {
            std::vector<Packed> payloads;
            std::vector<Real>   weights;
            RealVect            centroid = RealVect::Zero;
            Real                totalW   = 0.0;

            payloads.reserve(bl.hi - bl.lo);
            weights.reserve(bl.hi - bl.lo);

            for (std::size_t idx = bl.lo; idx < bl.hi; idx++) {
              const BucketParticle<Packed>& p = combined[idx];

              CH_assert(!p.isGhost);

              payloads.push_back(p.payload);
              weights.push_back(p.weight);
              centroid += p.weight * p.position;
              totalW += p.weight;
            }

            centroid /= totalW;

            bucketTreeCheckCentroid(centroid, totalW, bl.boxLo, bl.boxHi);

            if (a_isPositionValid(centroid)) {
              BucketParticle<Packed> merged;

              merged.position  = centroid;
              merged.weight    = totalW;
              merged.globalID  = a_allocateID();
              merged.ownerRank = myRank;
              merged.payload   = a_combine(payloads.data(), weights.data(), payloads.size());

              mergedResults.push_back(MergedResult{merged, patchIdx});

              for (std::size_t idx = bl.lo; idx < bl.hi; idx++) {
                consumedIDs.push_back(combined[idx].globalID);
              }
            }
          }

          continue;
        }

        // Boundary candidate. Debug knob: when a_enableBoundaryCarve is false, no box is ever
        // proposed here -- these particles are simply left untouched this pass (see
        // mergeBucketTreeCarve()'s own @param docs), so a suspected bug can be isolated to the carve
        // protocol below versus the interior-tier commit path above.
        if (a_enableBoundaryCarve && members.size() >= 2) {
          const BucketBoxKey key{leafVolume, anchor};

          const int boxIdx = static_cast<int>(myBoxes.size());
          myBoxes.push_back(RuntimeBox{key, patchIdx, members});

          for (const BucketMember& m : members) {
            if (m.owner == myRank) {
              // Append, never collapse -- see selfClaims' own docs: this same particle can already
              // carry an entry from a DIFFERENT one of my own boxes (a different patch's independent
              // build), and both must survive to be argmin'd together in step 3.
              selfClaims.push_back(SelfClaimEntry{m.id, key, boxIdx});
            }
          }
        }
        else if (a_enableBoundaryCarve) {
          // Exactly one member, and it must be the resident, exposed one (members.size()>=1 and
          // anyLocal both hold) -- can never reach the merge threshold alone, listen only.
          selfClaims.push_back(SelfClaimEntry{members[0].id, BucketBoxKey{}, -1});
        }
      }
    }
  }

  // particlesByID was appended to once per (id, patch-occurrence) above -- an id can appear more
  // than once (its home patch, plus any neighboring patch that gathered it as a ghost). Sort once
  // and collapse each run to a single entry, preferring the non-ghost copy when the run has one
  // (ghost-vs-ghost or non-ghost-vs-non-ghost duplicates carry identical data, so it doesn't matter
  // which of those survives) -- this is the dedup rule the old per-particle hash lookup used to
  // enforce inline.
  std::sort(particlesByID.begin(), particlesByID.end(), [](const auto& a_lhs, const auto& a_rhs) {
    return a_lhs.first < a_rhs.first;
  });
  {
    std::size_t writeIdx = 0;

    for (std::size_t readIdx = 0; readIdx < particlesByID.size();) {
      std::size_t runEnd = readIdx + 1;

      while (runEnd < particlesByID.size() && particlesByID[runEnd].first == particlesByID[readIdx].first) {
        runEnd++;
      }

      std::size_t chosen = readIdx;
      for (std::size_t k = readIdx; k < runEnd; k++) {
        if (!particlesByID[k].second.isGhost) {
          chosen = k;

          break;
        }
      }

      if (writeIdx != chosen) {
        particlesByID[writeIdx] = particlesByID[chosen];
      }
      writeIdx++;
      readIdx = runEnd;
    }

    particlesByID.resize(writeIdx);
  }

  // Binary-search lookup into the now-sorted, deduplicated particlesByID -- replaces the old
  // unordered_map's .at().
  auto findParticle = [&particlesByID](const ParticleID a_id) -> const BucketParticle<Packed>& {
    const auto it = std::lower_bound(
      particlesByID.begin(),
      particlesByID.end(),
      a_id,
      [](const std::pair<ParticleID, BucketParticle<Packed>>& a_entry, const ParticleID a_key) {
        return a_entry.first < a_key;
      });

    CH_assert(it != particlesByID.end() && it->first == a_id);

    return it->second;
  };

  // selfClaims was appended to during STEP 1 -- sort by id once so STEP 3 can process it as grouped
  // runs of equal id (a particle's possibly-several self-claim entries, see SelfClaimEntry's own
  // docs) instead of hashing.
  std::sort(selfClaims.begin(), selfClaims.end(), [](const SelfClaimEntry& a_lhs, const SelfClaimEntry& a_rhs) {
    return a_lhs.id < a_rhs.id;
  });

  // ==== STEP 2: Phase 1 -- claims to owners ====
  std::vector<std::vector<BucketClaim>> claimSendByRank(numRanks);

  for (std::size_t boxIdx = 0; boxIdx < myBoxes.size(); boxIdx++) {
    const RuntimeBox& box = myBoxes[boxIdx];

    for (const BucketMember& m : box.members) {
      if (m.owner != myRank) {
        claimSendByRank[m.owner].push_back(BucketClaim{m.id, box.key, myRank, static_cast<int>(boxIdx)});
      }
    }
  }

  // Sorted by memberID once received, so STEP 3's point lookups (claimsRange, below) and STEP 4's
  // full pass (grouped runs) don't need a hash map either.
  std::vector<BucketClaim> incomingClaims = bucketTreeExchangeByRank(claimSendByRank);
  std::sort(incomingClaims.begin(), incomingClaims.end(), [](const BucketClaim& a_lhs, const BucketClaim& a_rhs) {
    return a_lhs.memberID < a_rhs.memberID;
  });

  auto claimsRange = [&incomingClaims](const ParticleID a_id) {
    const auto lo = std::lower_bound(incomingClaims.begin(),
                                     incomingClaims.end(),
                                     a_id,
                                     [](const BucketClaim& a_c, const ParticleID a_key) {
                                       return a_c.memberID < a_key;
                                     });
    const auto hi = std::upper_bound(incomingClaims.begin(),
                                     incomingClaims.end(),
                                     a_id,
                                     [](const ParticleID a_key, const BucketClaim& a_c) {
                                       return a_key < a_c.memberID;
                                     });

    return std::make_pair(lo, hi);
  };

  // ==== STEP 3: local argmin -- nominal winner only, no deletion yet ====
  // The winner of a particle is either one of MY OWN boxes (rank == myRank, boxIdx identifies
  // WHICH one -- see selfClaims' own docs for why this can't just be "myRank", plural competing
  // boxes on this same rank are exactly the bug that distinction exists to prevent) or a foreign
  // rank (boxIdx meaningless).
  struct Winner
  {
    RankID       rank;
    int          boxIdx;
    BucketBoxKey key; // the winning key itself -- needed in step 4 to answer EACH incoming claim
                      // individually (won iff that claim's own key is not strictly worse than this).
  };

  struct WinnerEntry
  {
    ParticleID id;
    Winner     winner;
  };

  // Built in ascending-id order for free -- selfClaims is sorted by id, and its grouped runs are
  // visited in that same order below -- so no separate sort is needed before findWinner's binary
  // search.
  std::vector<WinnerEntry> nominalWinner;
  nominalWinner.reserve(selfClaims.size());

  for (std::size_t idx = 0; idx < selfClaims.size();) {
    const ParticleID id = selfClaims[idx].id;

    std::size_t idxEnd = idx + 1;

    while (idxEnd < selfClaims.size() && selfClaims[idxEnd].id == id) {
      idxEnd++;
    }

    bool         haveWinner = false;
    BucketBoxKey bestKey{};
    RankID       bestRank   = myRank;
    int          bestBoxIdx = -1;

    for (std::size_t k = idx; k < idxEnd; k++) {
      if (selfClaims[k].boxIdx < 0) {
        continue; // Listen-only sentinel -- not an actual self-claim.
      }
      if (!haveWinner || selfClaims[k].key < bestKey) {
        haveWinner = true;
        bestKey    = selfClaims[k].key;
        bestRank   = myRank;
        bestBoxIdx = selfClaims[k].boxIdx;
      }
    }

    const auto range = claimsRange(id);

    for (auto cit = range.first; cit != range.second; ++cit) {
      const BucketClaim& c = *cit;

      if (!haveWinner || c.key < bestKey) {
        haveWinner = true;
        bestKey    = c.key;
        bestRank   = c.proposerRank;
        // The winning claim's OWN box index (on ITS rank), not a sentinel -- step 4 compares
        // (rank, boxIdx) as the winning claim's IDENTITY, not just its key value, precisely so
        // that an exact key tie between two INDEPENDENT claims doesn't tell both of them "you
        // won" (see step 4's own comment). This must be the real index even for a foreign
        // winner, or that identity comparison could never match a foreign claim at all.
        bestBoxIdx = c.proposerBoxIdx;
      }
    }

    if (haveWinner) {
      nominalWinner.push_back(WinnerEntry{id, Winner{bestRank, bestBoxIdx, bestKey}});
    }

    idx = idxEnd;
  }

  // Binary-search lookup into the now-sorted nominalWinner -- replaces the old unordered_map's
  // .find()/.at().
  auto findWinner = [&nominalWinner](const ParticleID a_id) -> const Winner* {
    const auto it = std::lower_bound(nominalWinner.begin(),
                                     nominalWinner.end(),
                                     a_id,
                                     [](const WinnerEntry& a_entry, const ParticleID a_key) {
                                       return a_entry.id < a_key;
                                     });

    return (it != nominalWinner.end() && it->id == a_id) ? &it->winner : nullptr;
  };

  // ==== STEP 4: Phase 2 -- verdicts to proposers, foreign members only ====
  // One verdict PER INCOMING CLAIM, not per rank: two of the SAME rank's own boxes can both claim
  // this particle (from two different patches that rank owns), and each needs its own, correct
  // answer -- a rank-level "your rank won" cannot tell them apart. See BucketClaim::proposerBoxIdx.
  std::vector<std::vector<BucketVerdict>> verdictSendByRank(numRanks);

  for (std::size_t idx = 0; idx < incomingClaims.size();) {
    const ParticleID id = incomingClaims[idx].memberID;

    std::size_t idxEnd = idx + 1;

    while (idxEnd < incomingClaims.size() && incomingClaims[idxEnd].memberID == id) {
      idxEnd++;
    }

    const Winner* winner = findWinner(id);

    // Every id that shows up in incomingClaims was claimed on one of MY OWN, locally-owned,
    // boundary-exposed particles -- STEP 1's classification is exhaustive over exposed members, so
    // it always left a (possibly listen-only) selfClaims entry for it, which always yields a
    // nominalWinner entry (haveWinner requires only ONE candidate, and the incoming claim itself
    // is one). This is a correctness invariant, not a defensive fallback.
    CH_assert(winner != nullptr);

    for (std::size_t k = idx; k < idxEnd; k++) {
      const BucketClaim& c = incomingClaims[k];

      // Compare IDENTITY (which specific claim was chosen), not key VALUE: an exact key tie
      // between two INDEPENDENT claims is real and not vanishingly rare here (BrownianWalker's
      // regular initial distribution can make two independently-built boxes end up with the exact
      // same volume), and comparing key value alone would tell BOTH tied claimants "you won" --
      // exactly the kind of double-count this whole protocol exists to prevent. Step 3 already
      // deterministically picked exactly one winner among any tie; only that specific claim,
      // identified by (rank, box index), is told so.
      const bool won = (c.proposerRank == winner->rank) && (c.proposerBoxIdx == winner->boxIdx);

      verdictSendByRank[c.proposerRank].push_back(BucketVerdict{id, c.proposerBoxIdx, won});
    }

    idx = idxEnd;
  }

  const std::vector<BucketVerdict> incomingVerdicts = bucketTreeExchangeByRank(verdictSendByRank);

  // wonForeignByBox[boxIdx] = every foreign member THAT SPECIFIC one of my own boxes was told it won.
  // Indexed directly by boxIdx (a dense range over myBoxes) rather than hashed -- each box's own
  // list is bounded by its own membership (at most ppc), so sorting it once and using binary_search
  // at the read site below avoids both a per-box hash-table allocation (unordered_set would add one
  // per box, of which there can be thousands) and, for large ppc, the O(box size^2) cost a linear
  // scan would have across all of a box's foreign-member lookups.
  std::vector<std::vector<ParticleID>> wonForeignByBox(myBoxes.size());

  for (const BucketVerdict& v : incomingVerdicts) {
    if (v.won) {
      wonForeignByBox[v.claimantBoxIdx].push_back(v.memberID);
    }
  }

  for (std::vector<ParticleID>& won : wonForeignByBox) {
    std::sort(won.begin(), won.end());
  }

  // ==== STEP 5: assemble each of my boxes' provisional membership, commit or not ====
  std::vector<std::vector<BucketCommit>> commitSendByRank(numRanks);

  for (std::size_t boxIdx = 0; boxIdx < myBoxes.size(); boxIdx++) {
    const RuntimeBox& box = myBoxes[boxIdx];

    std::vector<ParticleID> survivingLocal;
    std::vector<ParticleID> survivingForeign;

    for (const BucketMember& m : box.members) {
      if (m.owner == myRank) {
        const Winner* winner = findWinner(m.id);

        if (winner != nullptr && winner->rank == myRank && winner->boxIdx == static_cast<int>(boxIdx)) {
          survivingLocal.push_back(m.id);
        }
      }
      else {
        const std::vector<ParticleID>& won = wonForeignByBox[boxIdx];

        if (std::binary_search(won.begin(), won.end(), m.id)) {
          survivingForeign.push_back(m.id);
        }
      }
    }

    const std::size_t totalSurvivors = survivingLocal.size() + survivingForeign.size();

    bool committed = false;

    if (totalSurvivors >= 2) {
      std::vector<Packed> payloads;
      std::vector<Real>   weights;
      RealVect            centroid = RealVect::Zero;
      Real                totalW   = 0.0;
      RealVect            surviveBoxLo, surviveBoxHi;
      bool                haveSurviveBox = false;

      payloads.reserve(totalSurvivors);
      weights.reserve(totalSurvivors);

      auto accumulate = [&](const ParticleID a_id) {
        const BucketParticle<Packed>& p = findParticle(a_id);

        payloads.push_back(p.payload);
        weights.push_back(p.weight);
        centroid += p.weight * p.position;
        totalW += p.weight;

        if (!haveSurviveBox) {
          surviveBoxLo   = p.position;
          surviveBoxHi   = p.position;
          haveSurviveBox = true;
        }
        else {
          for (int dir = 0; dir < SpaceDim; dir++) {
            surviveBoxLo[dir] = std::min(surviveBoxLo[dir], p.position[dir]);
            surviveBoxHi[dir] = std::max(surviveBoxHi[dir], p.position[dir]);
          }
        }
      };

      for (const ParticleID id : survivingLocal) {
        accumulate(id);
      }

      for (const ParticleID id : survivingForeign) {
        accumulate(id);
      }

      centroid /= totalW;

      bucketTreeCheckCentroid(centroid, totalW, surviveBoxLo, surviveBoxHi);

      if (a_isPositionValid(centroid)) {
        committed = true;

        BucketParticle<Packed> merged;

        merged.position  = centroid;
        merged.weight    = totalW;
        merged.globalID  = a_allocateID();
        merged.ownerRank = myRank;
        merged.payload   = a_combine(payloads.data(), weights.data(), payloads.size());

        mergedResults.push_back(MergedResult{merged, box.patchIdx});

        for (const ParticleID id : survivingLocal) {
          consumedIDs.push_back(id); // I own these -- confirmed now, no need to wait for anything.
        }
      }
    }

    // Phase 3 outgoing: every foreign member this box nominally won needs to hear the outcome,
    // whether this box committed or not -- its owner is otherwise stuck not knowing whether to
    // delete it (see mergeBucketTreeCarve()'s own docs for why this can't be skipped).
    for (const ParticleID id : survivingForeign) {
      commitSendByRank[findParticle(id).ownerRank].push_back(BucketCommit{id, committed});
    }
  }

  // ==== STEP 6: Phase 3 -- commit/release, foreign members only ====
  const std::vector<BucketCommit> incomingCommits = bucketTreeExchangeByRank(commitSendByRank);

  for (const BucketCommit& c : incomingCommits) {
    if (c.committed) {
      consumedIDs.push_back(c.memberID);
    }
  }

  std::sort(consumedIDs.begin(), consumedIDs.end());

  // ==== STEP 7: remove consumed particles ====
  for (const PatchWork& pw : patchWork) {
    ParticleSoA<P, Traits>& leaf = a_particles[pw.level][pw.din];

    std::size_t i = 0;

    while (i < leaf.size()) {
      if (!leaf.isGhost(i) && std::binary_search(consumedIDs.begin(), consumedIDs.end(), leaf.particleID(i))) {
        leaf.remove(i);
      }
      else {
        i++;
      }
    }
  }

  a_particles.clearGhostParticles();

  // ==== STEP 8: place merged results ====
  std::vector<std::vector<BucketParticle<Packed>>> scatterByDestRank(numRanks);

  auto insertHere = [&](const BucketParticle<Packed>& a_p, const int a_level, const DataIndex& a_din) {
    ParticleSoA<P, Traits>& leaf = a_particles[a_level][a_din];

    a_scatter(leaf, a_p);
  };

  for (const MergedResult& mr : mergedResults) {
    const PatchWork& pw = patchWork[mr.patchIdx];

    RealVect boxRealLo, boxRealHi;

    for (int dir = 0; dir < SpaceDim; dir++) {
      boxRealLo[dir] = probLo[dir] + pw.box.smallEnd(dir) * pw.dx[dir];
      boxRealHi[dir] = probLo[dir] + (pw.box.bigEnd(dir) + 1) * pw.dx[dir];
    }

    bool inOwnBox = true;

    for (int dir = 0; dir < SpaceDim && inOwnBox; dir++) {
      if (mr.particle.position[dir] < boxRealLo[dir] || mr.particle.position[dir] >= boxRealHi[dir]) {
        inOwnBox = false;
      }
    }

    if (inOwnBox) {
      insertHere(mr.particle, pw.level, pw.din);

      continue;
    }

    const auto dst = a_particles.findDestination(mr.particle.position);

    if (!dst.valid) {
      MayDay::Error("ParticleManagement::mergeBucketTreeCarve -- merged particle not found in any patch");
    }

    if (dst.rank == myRank) {
      const DataIndex din = a_amr.getLevelTiles(realm)[dst.level]->getMyGrids().at(dst.gridIndex);

      insertHere(mr.particle, dst.level, din);
    }
    else {
      BucketParticle<Packed> corrected = mr.particle;
      corrected.ownerRank              = dst.rank;

      scatterByDestRank[dst.rank].push_back(corrected);
    }
  }

  const std::vector<BucketParticle<Packed>> incomingScattered = bucketTreeExchangeByRank(scatterByDestRank);

  for (const BucketParticle<Packed>& p : incomingScattered) {
    const auto dst = a_particles.findDestination(p.position);

    if (!dst.valid || dst.rank != myRank) {
      MayDay::Error("ParticleManagement::mergeBucketTreeCarve -- incoming scattered particle not "
                    "found in any of this rank's own patches");
    }

    const DataIndex din = a_amr.getLevelTiles(realm)[dst.level]->getMyGrids().at(dst.gridIndex);

    BucketParticle<Packed> corrected = p;
    corrected.ownerRank              = myRank;

    insertHere(corrected, dst.level, din);
  }
}

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#endif
