/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_KDParticleMergeImplem.H
 * @brief  Implementation of CD_KDParticleMerge.H
 * @author Robert Marskar
 */

#ifndef CD_KDPARTICLEMERGEIMPLEM_H
#define CD_KDPARTICLEMERGEIMPLEM_H

// Std includes
#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>
#include <vector>

// Chombo includes
#include <BaseFab.H>
#include <Box.H>
#include <CH_Timer.H>
#include <MayDay.H>
#include <SPMD.H>
#include <DataIterator.H>

// Our includes
#include <CD_LevelTiles.H>
#include <CD_KDParticleMerge.H>
#include <CD_NamespaceHeader.H>

namespace ParticleManagement {

namespace detail {

template <typename T>
inline std::vector<T>
kdExchangeByRank(const std::vector<std::vector<T>>& a_sendByRank)
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
 * @param[out] a_boxLo Low corner.
 * @param[out] a_boxHi High corner.
 * @param[in] a_particles Particle buffer.
 * @param[in] a_lo Start of the range (inclusive).
 * @param[in] a_hi End of the range (exclusive).
 */
template <typename Packed>
inline void
kdBBox(RealVect&                                 a_boxLo,
       RealVect&                                 a_boxHi,
       const std::vector<MergeParticle<Packed>>& a_particles,
       const std::size_t                         a_lo,
       const std::size_t                         a_hi) noexcept
{
  CH_assert(a_hi > a_lo);
  CH_assert(a_hi <= a_particles.size());

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
 * @details buildKDQuotaLeaves() and both merge entry points must use the same value: the former
 * sets each cell's leaf budget assuming a leaf merges to one particle, the latter vetoes leaves
 * exceeding this bound. A disagreement spends budget on leaves that are then rejected.
 */
static constexpr Real s_kdMaxLeafExtent = 1.0;

/**
 * @brief Largest per-axis extent of [a_boxLo,a_boxHi], expressed as a fraction of that axis's own
 * cell width.
 * @details The correctness backstop buildKDQuotaLeaves(), mergeKDCarve() and mergeKDPatch() rely on to
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
kdMaxAxisSpan(const RealVect& a_boxLo, const RealVect& a_boxHi, const RealVect& a_dx) noexcept
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
kdCellKeyOf(const RealVect& a_position, const RealVect& a_probLo, const RealVect& a_dx) noexcept
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
 * buildKDQuotaLeaves()'s own docs) -- IntVect has no built-in total order.
 * @param[in] a_lhs Left-hand cell index.
 * @param[in] a_rhs Right-hand cell index.
 * @return True iff a_lhs sorts strictly before a_rhs.
 */
inline bool
kdIntVectLess(const IntVect& a_lhs, const IntVect& a_rhs) noexcept
{
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (a_lhs[dir] != a_rhs[dir]) {
      return a_lhs[dir] < a_rhs[dir];
    }
  }

  return false;
}

/**
 * @brief Tally a per-cell particle-count histogram over one patch's gathered particles.
 * @details Built from every particle's own position, so it is unaffected by how buildKDQuotaLeaves()
 * later partitions the space -- a node-local count is NOT a safe substitute, see that function's own
 * docs for why. Shared by mergeKDCarve()/mergeKDPatch()/mergeKDInterior(): all three need the same
 * ground truth, and the invariant that it comes from true positions rather than from the partition is
 * easy to break if each keeps its own copy.
 *
 * @a a_counts must include at least one ghost cell, since a gathered ghost sits up to one cell
 * outside the patch's own box.
 *
 * The count is stored in @a a_counts' @c Real elements but is always a non-negative integer. A double
 * represents every integer exactly up to 2^53, so consumers -- e.g. the crowding test in
 * buildKDQuotaLeaves() -- may @c static_cast<int> it losslessly for any realistic per-cell count.
 * @tparam Packed See MergeParticle.
 * @param[out] a_counts Per-cell counts, zeroed here before tallying.
 * @param[in] a_particles Gathered particles for this patch.
 * @param[in] a_probLo Physical position of the domain's low corner.
 * @param[in] a_dx This patch's own cell spacing.
 */
template <typename Packed>
inline void
kdFillCellHistogram(FArrayBox&                                a_counts,
                    const std::vector<MergeParticle<Packed>>& a_particles,
                    const RealVect&                           a_probLo,
                    const RealVect&                           a_dx) noexcept
{
  a_counts.setVal(0);

  for (const MergeParticle<Packed>& p : a_particles) {
    const IntVect cell = kdCellKeyOf(p.position, a_probLo, a_dx);

    // Every gathered particle is resident in this patch or a ghost within one cell of it, so its
    // cell lies inside a_counts' box. The previous map tolerated any key silently; a cell-indexed
    // holder does not, which is the point -- the invariant is now checked instead of assumed.
    CH_assert(a_counts.box().contains(cell));

    a_counts(cell, 0)++;
  }
}

/**
 * @brief Split particles[a_lo,a_hi) in place by the longest axis at the count-median -- the sole
 * split rule used by buildKDQuotaLeaves().
 * @details The split point is the position-sorted median particle (via nth_element), balancing
 * COUNT as evenly as possible -- mirrors equal_weight_kd's own single-cell logic, but by count,
 * never by weight, and with NO synthetic particle splitting (a particle is never divided across
 * the two children). nth_element divides the range in half by count however the positions are
 * distributed, so a tightly clustered node still splits evenly.
 * @param[in,out] a_particles Particle buffer; reordered in [a_lo,a_hi).
 * @param[in] a_lo Start of the range (inclusive).
 * @param[in] a_hi End of the range (exclusive); a_hi - a_lo must be >= 2.
 * @param[in] a_axis Split axis -- the node's longest box axis, precomputed once by the caller (which
 * already needs the box for the span test) rather than re-derived here from a second bounding box.
 * @return The split index: [a_lo,mid) and [mid,a_hi) are the two children, each non-empty.
 */
template <typename Packed>
inline std::size_t
kdSplitCountMedian(std::vector<MergeParticle<Packed>>& a_particles,
                   const std::size_t                   a_lo,
                   const std::size_t                   a_hi,
                   const int                           a_axis) noexcept
{
  CH_assert(a_hi - a_lo >= 2);

  const std::size_t mid = a_lo + (a_hi - a_lo) / 2;

  std::nth_element(a_particles.begin() + a_lo,
                   a_particles.begin() + mid,
                   a_particles.begin() + a_hi,
                   [a_axis](const MergeParticle<Packed>& a_p1, const MergeParticle<Packed>& a_p2) noexcept -> bool {
                     return a_p1.position[a_axis] < a_p2.position[a_axis];
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
 * weight forces a 1-vs-(n-1) split. That is harmless here only because buildKDQuotaLeaves()
 * enforces a hard per-cell leaf quota independently of how the split plane behaves: an uneven split
 * can cost merge quality but cannot leave a cell above target. Under any of the earlier,
 * quota-less formulations the same degeneracy was fatal.
 * @param[in,out] a_particles Particle buffer; reordered in [a_lo,a_hi).
 * @param[in] a_lo Start of the range (inclusive).
 * @param[in] a_hi End of the range (exclusive); a_hi - a_lo must be >= 2.
 * @param[in] a_axis Split axis -- the node's longest box axis, precomputed once by the caller (which
 * already needs the box for the span test) rather than re-derived here from a second bounding box.
 * @return The split index: [a_lo,mid) and [mid,a_hi) are the two children, each non-empty.
 */
template <typename Packed>
inline std::size_t
kdSplitWeightMedian(std::vector<MergeParticle<Packed>>& a_particles,
                    const std::size_t                   a_lo,
                    const std::size_t                   a_hi,
                    const int                           a_axis) noexcept
{
  CH_assert(a_hi - a_lo >= 2);

  std::sort(a_particles.begin() + a_lo,
            a_particles.begin() + a_hi,
            [a_axis](const MergeParticle<Packed>& a_p1, const MergeParticle<Packed>& a_p2) noexcept -> bool {
              return a_p1.position[a_axis] < a_p2.position[a_axis];
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

  const std::size_t clamped = std::max(a_lo + 1, std::min(a_hi - 1, mid));

  // Both children must be non-empty or the caller's recursion cannot terminate.
  CH_assert(clamped > a_lo && clamped < a_hi);

  return clamped;
}

template <typename Packed>
inline void
buildKDQuotaLeaves(std::vector<MergeParticle<Packed>>& a_particles,
                   FArrayBox&                          a_used,
                   std::vector<KDLeaf>&                a_leaves,
                   const int                           a_ppc,
                   const Real                          a_splitWeightLeafDx,
                   const RealVect&                     a_dx,
                   const RealVect&                     a_probLo,
                   const FArrayBox&                    a_cellCounts) noexcept
{
  CH_TIME("ParticleManagement::buildKDQuotaLeaves");

  a_leaves.clear();

  if (a_particles.empty() || a_ppc <= 0) {
    return;
  }

  // Partition into the particles that may be merged (those in TRUE-crowded cells) and those that
  // may not, moving the mergeable ones to the front so every leaf below is still a contiguous range.
  // Particles in uncrowded cells are emitted as singleton leaves: the merge threshold is 2 members,
  // so a singleton is never merged, and an uncrowded cell can therefore never be drained below
  // target by a neighbour's merge.
  const auto firstUnmergeable = std::stable_partition(a_particles.begin(),
                                                      a_particles.end(),
                                                      [&](const MergeParticle<Packed>& a_p) {
                                                        const IntVect cell = kdCellKeyOf(a_p.position, a_probLo, a_dx);

                                                        // A cell outside the histogram box reads as 0
                                                        // (empty) and is therefore never crowded.
                                                        return a_cellCounts.box().contains(cell)
                                                                 ? static_cast<int>(a_cellCounts(cell, 0)) > a_ppc
                                                                 : false;
                                                      });

  const std::size_t numMergeable = static_cast<std::size_t>(firstUnmergeable - a_particles.begin());

  for (std::size_t idx = numMergeable; idx < a_particles.size(); idx++) {
    RealVect boxLo, boxHi;
    kdBBox(boxLo, boxHi, a_particles, idx, idx + 1);

    a_leaves.push_back(KDLeaf{idx, idx + 1, boxLo, boxHi});
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
  // Returns both the centroid's cell AND the range's total weight. The weight sum is a byproduct of
  // the centroid accumulation, so returning it here removes the separate per-child rangeWeight() pass
  // the heap loop used to make (same summation order, so the Real is bit-identical).
  auto centroidCellOf = [&](const std::size_t a_lo, const std::size_t a_hi) -> std::pair<IntVect, Real> {
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

    return {kdCellKeyOf(centroid, a_probLo, a_dx), totalW};
  };

  struct KDCand
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
  const auto byWeight = [](const KDCand& a_lhs, const KDCand& a_rhs) noexcept -> bool {
    return a_lhs.weight < a_rhs.weight;
  };

  // Live per-cell leaf count: how many leaves are currently centred in each cell. Caller-owned mesh
  // data, zeroed here; its ghost cell must cover every cell a leaf centroid can fall in, since
  // gathered ghosts sit up to one cell outside this patch.
  FArrayBox& used = a_used;
  used.setVal(0.0);

  // Reserved to numMergeable, a safe upper bound for both: leaves partition [0,numMergeable) with at
  // least one member each, so #leaves <= numMergeable, and the split frontier can never exceed the
  // leaf count. One allocation apiece instead of the push_back reallocation ladder, every patch.
  std::vector<KDCand> heap;
  std::vector<KDLeaf> finalLeaves;

  heap.reserve(numMergeable);
  finalLeaves.reserve(numMergeable);

  {
    const auto [rootCell, rootWeight] = centroidCellOf(0, numMergeable);

    // A centroid of positions gathered for this patch cannot fall outside the patch grown by the
    // ghost width. Asserted rather than assumed: the map this replaced would have accepted any cell
    // silently, which is exactly the class of bug a bounded holder is meant to surface.
    CH_assert(used.box().contains(rootCell));

    used(rootCell, 0)++;
    heap.push_back(KDCand{rootWeight, numMergeable, 0, numMergeable, rootCell});
  }

  while (!heap.empty()) {
    std::pop_heap(heap.begin(), heap.end(), byWeight);

    const KDCand node = heap.back();
    heap.pop_back();

    if (node.count < 2) {
      RealVect boxLo, boxHi;
      kdBBox(boxLo, boxHi, a_particles, node.lo, node.hi);

      finalLeaves.push_back(KDLeaf{node.lo, node.hi, boxLo, boxHi});

      continue;
    }

    // One bounding box per node, computed once here and reused three ways: the span test, the split
    // axis (formerly re-derived from a second box inside each split helper), and -- if the split is
    // refused below -- the emitted leaf (formerly re-derived a third time in the emission pass). A
    // split only reorders within [lo,hi), so the box is unchanged by it either way.
    RealVect nodeBoxLo, nodeBoxHi;
    kdBBox(nodeBoxLo, nodeBoxHi, a_particles, node.lo, node.hi);

    const Real nodeSpan  = kdMaxAxisSpan(nodeBoxLo, nodeBoxHi, a_dx);
    const int  splitAxis = (nodeBoxHi - nodeBoxLo).maxDir(true);

    // While a node spans several cells a split apportions leaves BETWEEN cells, which is a question
    // of counts; once it is narrow the weight median is what drives the super-particles toward equal
    // weight. This is an extent test, not a "both corners in one cell" test, so a narrow node
    // straddling a face still qualifies.
    const std::size_t mid = (nodeSpan <= a_splitWeightLeafDx)
                              ? kdSplitWeightMedian(a_particles, node.lo, node.hi, splitAxis)
                              : kdSplitCountMedian(a_particles, node.lo, node.hi, splitAxis);

    const auto [cellLeft, weightLeft]   = centroidCellOf(node.lo, mid);
    const auto [cellRight, weightRight] = centroidCellOf(mid, node.hi);

    CH_assert(used.box().contains(cellLeft));
    CH_assert(used.box().contains(cellRight));

    used(node.cell, 0)--;
    used(cellLeft, 0)++;
    used(cellRight, 0)++;

    // The quota constrains only leaves that can actually merge. A node wider than
    // s_kdMaxLeafExtent is vetoed by the caller and leaves every member behind
    // individually, so it costs its cell its full member count rather than one particle; splitting
    // it further lowers that cell's final population even though it raises the leaf count.
    if (nodeSpan <= s_kdMaxLeafExtent && (used(cellLeft, 0) > a_ppc || used(cellRight, 0) > a_ppc)) {
      used(cellLeft, 0)--;
      used(cellRight, 0)--;
      used(node.cell, 0)++;

      finalLeaves.push_back(KDLeaf{node.lo, node.hi, nodeBoxLo, nodeBoxHi});

      continue;
    }

    heap.push_back(KDCand{weightLeft, mid - node.lo, node.lo, mid, cellLeft});
    std::push_heap(heap.begin(), heap.end(), byWeight);

    heap.push_back(KDCand{weightRight, node.hi - mid, mid, node.hi, cellRight});
    std::push_heap(heap.begin(), heap.end(), byWeight);
  }

  for (const KDLeaf& leaf : finalLeaves) {
    a_leaves.push_back(leaf);
  }

  // The leaves must partition the input exactly: every particle in exactly one leaf, no overlap and
  // nothing dropped. Callers rely on this both for the merge itself and for the invariant that
  // total weight is conserved.
#ifndef NDEBUG
  std::size_t coveredMembers = 0;

  for (const KDLeaf& bl : a_leaves) {
    CH_assert(bl.hi > bl.lo);
    CH_assert(bl.hi <= a_particles.size());

    coveredMembers += bl.hi - bl.lo;
  }

  CH_assert(coveredMembers == a_particles.size());
#endif
}

/**
 * @brief Real-space bounds of a box: the half-open region [lo, hi) that its cells cover.
 * @details Half-open on purpose -- the high corner is the OUTER face of the last cell, so testing
 * pos >= lo && pos < hi per direction partitions space between abutting boxes with no gap and no
 * double-cover, which is exactly the "is this particle still in its own patch" test both callers make.
 * @param[out] a_boxLo  Low corner (inclusive).
 * @param[out] a_boxHi  High corner (exclusive).
 * @param[in]  a_box    Box whose bounds are wanted.
 * @param[in]  a_dx     Grid spacing at the box's level.
 * @param[in]  a_probLo Lower-left corner of the physical domain.
 */
inline void
kdBoxRealBounds(RealVect&       a_boxLo,
                RealVect&       a_boxHi,
                const Box&      a_box,
                const RealVect& a_dx,
                const RealVect& a_probLo) noexcept
{
  for (int dir = 0; dir < SpaceDim; dir++) {
    a_boxLo[dir] = a_probLo[dir] + a_box.smallEnd(dir) * a_dx[dir];
    a_boxHi[dir] = a_probLo[dir] + (a_box.bigEnd(dir) + 1) * a_dx[dir];
  }
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
kdCheckCentroid(const RealVect& a_centroid,
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
mergeKDCarve(ParticleContainer<P, Traits>& a_particles,
             EBAMRFAB&                     a_cellHistogram,
             EBAMRFAB&                     a_leafQuota,
             const AmrMesh&                a_amr,
             const int                     a_ppc,
             const Real                    a_splitWeightLeafDx,
             const Gather&                 a_gather,
             const Combine&                a_combine,
             const Scatter&                a_scatter,
             const Allocator&              a_allocateID,
             const PosValid&               a_isPositionValid)
{
  using namespace detail;

  // Carve-protocol types. Local to this function because nothing else constructs or consumes them:
  // they describe one pass of the claim/verdict/commit exchange and have no meaning outside it.

  // Box ranking key. Tightest box wins a contested particle; the anchor id is the deterministic
  // tiebreak, since AABB volume alone is not a strict total order.
  struct KDBoxKey
  {
    Real       volume;
    ParticleID anchor;

    bool
    operator<(const KDBoxKey& a_rhs) const noexcept
    {
      return (volume != a_rhs.volume) ? (volume < a_rhs.volume) : (anchor < a_rhs.anchor);
    }
  };

  // One member of a boundary box: the particle and the rank that owns it.
  struct KDMember
  {
    ParticleID id;
    RankID     owner;
  };

  // Phase 1 wire type: a box's claim on one member, routed to that member's owner.
  struct KDClaim
  {
    ParticleID memberID;
    KDBoxKey   key;
    RankID     proposerRank;
    int        proposerBoxIdx; // which of the proposer's own boxes claimed it
  };

  // Phase 2 wire type: the owner's nominal argmin winner, answered per claiming box.
  struct KDVerdict
  {
    ParticleID memberID;
    int        claimantBoxIdx;
    bool       won;
  };

  // Phase 3 wire type: the proposer's actual outcome. Only this authorizes a deletion.
  struct KDCommit
  {
    ParticleID memberID;
    bool       committed;
  };

  CH_TIMERS("ParticleManagement::mergeKDCarve");
  CH_TIMER("ParticleManagement::mergeKDCarve::build_classify", t_build);
  CH_TIMER("ParticleManagement::mergeKDCarve::carve_exchange", t_carve);
  CH_TIMER("ParticleManagement::mergeKDCarve::remove_consumed", t_remove);
  CH_TIMER("ParticleManagement::mergeKDCarve::place_results", t_place);

  const std::string realm       = a_particles.getRealm();
  const int         finestLevel = a_amr.getFinestLevel();
  const RealVect    probLo      = a_amr.getProbLo();

  // Per-cell scratch, as mesh data on this realm rather than per-patch buffers: both quantities are
  // cell data, and this is what the rest of the code uses for cell data. The caller owns the holders
  // so they are allocated once per regrid rather than once per call -- AmrMesh::allocate() builds a
  // Copier per level, which is cheap in isolation but adds up at scale. Their one ghost cell is what
  // the gathered ghosts and any just-outside leaf centroid need.
  EBAMRFAB& histogram = a_cellHistogram;
  EBAMRFAB& leafQuota = a_leafQuota;

  CH_assert(histogram[0]->nComp() == 1);
  CH_assert(leafQuota[0]->nComp() == 1);
  CH_assert(histogram[0]->ghostVect() >= IntVect::Unit);
  CH_assert(leafQuota[0]->ghostVect() >= IntVect::Unit);

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
    KDBoxKey              key;
    int                   patchIdx;
    std::vector<KDMember> members;
  };

  // particle: the merged super-particle itself.
  // patchIdx: which patch's leaf produced it, used as the cheap first guess when placing it.
  struct MergedResult
  {
    MergeParticle<Packed> particle;
    int                   patchIdx;
  };

  // One entry per (locally-owned, boundary-exposed particle) x (one of MY OWN boxes that lists it),
  // appended during STEP 1 and sorted by id once STEP 1 finishes -- see selfClaims' declaration
  // below for why an id can have more than one entry, and boxIdx < 0 for what a sentinel entry means.
  // id:     the locally-owned, boundary-exposed particle this entry concerns.
  // key:    the claiming box's key. Meaningless when boxIdx is negative.
  // boxIdx: index into myBoxes, or -1 for a listen-only sentinel. A sentinel carries no claim from
  //         any of this rank's own boxes; it records that this id must still be processed against
  //         foreign claims in STEP 3, which is the case for an unmergeable leaf's exposed member or
  //         a boundary leaf holding fewer than two members.
  struct SelfClaimEntry
  {
    ParticleID id;
    KDBoxKey   key;
    int        boxIdx;
  };

  std::vector<PatchWork> patchWork;

  // gatheredParticles: every particle this rank gathered, locals and ghosts alike, in gather order.
  //                 Append-only and never reordered, so an index into it is a stable handle. Holds the
  //                 full payloads particlesByID used to carry inline.
  // particlesByID:  one (id, slot-into-gatheredParticles) entry per (id, patch-occurrence), appended
  //                 during STEP 1, then deduplicated once STEP 1 finishes by sorting and collapsing
  //                 runs, preferring the non-ghost copy. Sorting these light (id, slot) keys instead of
  //                 whole particles is what keeps the dedup off the fat-struct memory path. Point
  //                 lookups after that go through findParticle() below, which binary-searches it and
  //                 dereferences the surviving slot.
  // myBoxes:        this rank's own boundary boxes, indexed by the boxIdx carried in claims.
  // consumedIDs:    ids whose particle has been merged away and must be deleted. Sorted once,
  //                 before STEP 7 consumes it by binary search.
  // mergedResults:  super-particles produced this pass, placed in STEP 8.
  std::vector<MergeParticle<Packed>>              gatheredParticles;
  std::vector<std::pair<ParticleID, std::size_t>> particlesByID;
  std::vector<RuntimeBox>                         myBoxes;
  std::vector<ParticleID>                         consumedIDs;
  std::vector<MergedResult>                       mergedResults;

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

  gatheredParticles.reserve(combinedCountUpperBound);
  particlesByID.reserve(combinedCountUpperBound);
  selfClaims.reserve(combinedCountUpperBound);
  consumedIDs.reserve(combinedCountUpperBound);

  // Per-patch leaf/member scratch, declared once and reused across every patch below rather than
  // reallocated inside the loop -- the same reuse pattern as the histogram/leafQuota holders above.
  // buildKDQuotaLeaves() clear()s leaves on entry and members is cleared per leaf, so each reuse
  // retains the heap capacity grown by earlier patches instead of starting from nothing.
  std::vector<KDLeaf>   leaves;
  std::vector<KDMember> members;

  // Which cells hold a particle that some other box can also see. Read from the realm, which derives it
  // from the three ghost masks when they are built; see Realm::m_particleGhostExposure.
  const AMRMask& exposure = a_amr.getParticleGhostExposure(realm, 1);

  // ==== STEP 1: gather + build + classify + commit interior, per patch ====
  CH_START(t_build);

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = a_amr.getGrids(realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const RealVect dx = a_amr.getDx()[lvl] * RealVect::Unit;

    const int nbox = dit.size();

    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      const BaseFab<bool>& exposureDin = (*exposure[lvl])[din];

      std::vector<MergeParticle<Packed>> combined;
      combined.reserve(leaf.size());

      for (std::size_t i = 0; i < leaf.size(); i++) {
        MergeParticle<Packed> p;

        p.position  = leaf.position(i);
        p.weight    = leaf.weight(i);
        p.globalID  = leaf.particleID(i);
        p.ownerRank = leaf.rankID(i);
        p.isGhost   = leaf.isGhost(i);
        p.payload   = a_gather(leaf, i);

        combined.push_back(p);

        // Append unconditionally -- deduplicated by id (preferring the non-ghost copy, per
        // MergeParticle::isGhost's own docs) in the one pass below once STEP 1 finishes gathering,
        // rather than via a per-particle hash lookup here. The full payload lives in gatheredParticles;
        // particlesByID carries only (id, slot) so the dedup sort stays off the fat-struct path.
        gatheredParticles.push_back(p);
        particlesByID.emplace_back(p.globalID, gatheredParticles.size() - 1);
      }

      patchWork.push_back(PatchWork{lvl, din, dx, dbl[din]});
      const int patchIdx = static_cast<int>(patchWork.size()) - 1;

      if (combined.empty()) {
        continue;
      }

      // Ground truth for how crowded each cell really is -- see kdFillCellHistogram(). Both this and
      // the live quota are cell data, so they live in the mesh holders allocated once above rather
      // than in per-patch scratch. Their ghost cell covers the gathered ghosts, which sit up to one
      // cell outside this patch, and any leaf centroid that lands just outside it.
      FArrayBox& cellCounts = (*histogram[lvl])[din];
      FArrayBox& quota      = (*leafQuota[lvl])[din];

      kdFillCellHistogram(cellCounts, combined, probLo, dx);

      buildKDQuotaLeaves(combined, quota, leaves, a_ppc, a_splitWeightLeafDx, dx, probLo, cellCounts);

      for (const KDLeaf& bl : leaves) {
        members.clear();
        members.reserve(bl.hi - bl.lo);

        bool anyLocal = false;

        // Running minimum over ALL members; bl.hi > bl.lo always holds for a leaf. Folded into the
        // loop below so the boundary branch does not need a separate pass just to compute it.
        ParticleID anchor = combined[bl.lo].globalID;

        for (std::size_t idx = bl.lo; idx < bl.hi; idx++) {
          members.push_back(KDMember{combined[idx].globalID, combined[idx].ownerRank});
          anchor = std::min(anchor, combined[idx].globalID);

          // "Local" here means physically resident in THIS patch, not merely owned by this rank --
          // a particle owned by this rank via a DIFFERENT one of its own patches, appearing here
          // only as a ghost, gives this leaf no standing (see MergeParticle::isGhost's own docs).
          if (!combined[idx].isGhost) {
            anyLocal = true;
          }
        }

        if (!anyLocal) {
          // Nothing physically resident here, so this rank has no standing to act on it.
          continue;
        }

        Real leafVolume = 1.0;

        for (int dir = 0; dir < SpaceDim; dir++) {
          leafVolume *= (bl.boxHi[dir] - bl.boxLo[dir]);
        }

        // Per-axis extent, not aggregate volume -- see kdMaxAxisSpan()'s own docs for why a
        // volume ratio alone cannot be trusted here (a long, thin leaf can have volume <=
        // one cell while still spanning several cells along one axis). leafVolume itself
        // is still needed below, unrelated to this check -- it's the box key's tie-break priority.
        const bool unmergeable = kdMaxAxisSpan(bl.boxLo, bl.boxHi, dx) > s_kdMaxLeafExtent;

        // Per-resident-member exposure -- needed either way (leaf-wide OR for the interior/
        // boundary decision below; per-particle for the unmergeable-leaf listen-only case). Ghosts
        // are skipped entirely, not just because they have no standing (above) but because their
        // position lies outside this patch's own box -- the exposure mask carries no ghost cells, so
        // a ghost's cell is not addressable in it at all. Filtering on isGhost (not ownerRank) is what
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

          const IntVect cell = kdCellKeyOf(combined[idx].position, probLo, dx);

          CH_assert(exposureDin.box().contains(cell));

          const bool isExposed = exposureDin(cell, 0);

          exposed[idx - bl.lo] = isExposed;
          anyExposed           = anyExposed || isExposed;
        }

        if (unmergeable) {
          for (std::size_t idx = bl.lo; idx < bl.hi; idx++) {
            // exposed[] is already false for every ghost by construction (the exposure loop above
            // skips them), so this is equivalent to also checking !isGhost -- kept explicit anyway
            // as the actual invariant being relied on, not an implicit one.
            if (!combined[idx].isGhost && exposed[idx - bl.lo]) {
              // Listen-only sentinel; see SelfClaimEntry.
              selfClaims.push_back(SelfClaimEntry{combined[idx].globalID, KDBoxKey{}, -1});
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
              const MergeParticle<Packed>& p = combined[idx];

              CH_assert(!p.isGhost);

              payloads.push_back(p.payload);
              weights.push_back(p.weight);
              centroid += p.weight * p.position;
              totalW += p.weight;
            }

            centroid /= totalW;

            kdCheckCentroid(centroid, totalW, bl.boxLo, bl.boxHi);

            if (a_isPositionValid(centroid)) {
              MergeParticle<Packed> merged;

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

        // Boundary candidate.
        if (members.size() >= 2) {
          const KDBoxKey key{leafVolume, anchor};

          const int boxIdx = static_cast<int>(myBoxes.size());
          myBoxes.push_back(RuntimeBox{key, patchIdx, members});

          for (const KDMember& m : members) {
            if (m.owner == myRank) {
              // Append, never collapse -- see selfClaims' own docs: this same particle can already
              // carry an entry from a DIFFERENT one of my own boxes (a different patch's independent
              // build), and both must survive to be argmin'd together in step 3.
              selfClaims.push_back(SelfClaimEntry{m.id, key, boxIdx});
            }
          }
        }
        else {
          // Exactly one member, and it must be the resident, exposed one (members.size()>=1 and
          // anyLocal both hold) -- can never reach the merge threshold alone, listen only.
          selfClaims.push_back(SelfClaimEntry{members[0].id, KDBoxKey{}, -1});
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
        if (!gatheredParticles[particlesByID[k].second].isGhost) {
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
  auto findParticle = [&particlesByID, &gatheredParticles](const ParticleID a_id) -> const MergeParticle<Packed>& {
    const auto it = std::lower_bound(particlesByID.begin(),
                                     particlesByID.end(),
                                     a_id,
                                     [](const std::pair<ParticleID, std::size_t>& a_entry, const ParticleID a_key) {
                                       return a_entry.first < a_key;
                                     });

    CH_assert(it != particlesByID.end() && it->first == a_id);

    return gatheredParticles[it->second];
  };

  // selfClaims was appended to during STEP 1 -- sort by id once so STEP 3 can process it as grouped
  // runs of equal id (a particle's possibly-several self-claim entries, see SelfClaimEntry's own
  // docs) instead of hashing.
  std::sort(selfClaims.begin(), selfClaims.end(), [](const SelfClaimEntry& a_lhs, const SelfClaimEntry& a_rhs) {
    return a_lhs.id < a_rhs.id;
  });

  CH_STOP(t_build);

  // ==== STEP 2: Phase 1 -- claims to owners ====
  CH_START(t_carve);
  std::vector<std::vector<KDClaim>> claimSendByRank(numRanks);

  for (std::size_t boxIdx = 0; boxIdx < myBoxes.size(); boxIdx++) {
    const RuntimeBox& box = myBoxes[boxIdx];

    for (const KDMember& m : box.members) {
      if (m.owner != myRank) {
        claimSendByRank[m.owner].push_back(KDClaim{m.id, box.key, myRank, static_cast<int>(boxIdx)});
      }
    }
  }

  // Sorted by memberID once received, so STEP 3's point lookups (claimsRange, below) and STEP 4's
  // full pass (grouped runs) don't need a hash map either.
  std::vector<KDClaim> incomingClaims = kdExchangeByRank(claimSendByRank);
  std::sort(incomingClaims.begin(), incomingClaims.end(), [](const KDClaim& a_lhs, const KDClaim& a_rhs) {
    return a_lhs.memberID < a_rhs.memberID;
  });

  auto claimsRange = [&incomingClaims](const ParticleID a_id) {
    const auto lo = std::lower_bound(incomingClaims.begin(),
                                     incomingClaims.end(),
                                     a_id,
                                     [](const KDClaim& a_c, const ParticleID a_key) {
                                       return a_c.memberID < a_key;
                                     });
    const auto hi = std::upper_bound(incomingClaims.begin(),
                                     incomingClaims.end(),
                                     a_id,
                                     [](const ParticleID a_key, const KDClaim& a_c) {
                                       return a_key < a_c.memberID;
                                     });

    return std::make_pair(lo, hi);
  };

  // ==== STEP 3: local argmin -- nominal winner only, no deletion yet ====
  // The winner of a particle is either one of MY OWN boxes (rank == myRank, boxIdx identifies
  // WHICH one -- see selfClaims' own docs for why this can't just be "myRank", plural competing
  // boxes on this same rank are exactly the bug that distinction exists to prevent) or a foreign
  // rank (boxIdx meaningless).
  // rank:   the rank owning the winning box.
  // boxIdx: that box's index within the winning rank's own box list.
  // key:    the winning key itself, needed in STEP 4 to answer each incoming claim individually --
  //         a claim won iff its own key is not strictly worse than this one.
  struct Winner
  {
    RankID   rank;
    int      boxIdx;
    KDBoxKey key;
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

    bool     haveWinner = false;
    KDBoxKey bestKey{};
    RankID   bestRank   = myRank;
    int      bestBoxIdx = -1;

    for (std::size_t k = idx; k < idxEnd; k++) {
      if (selfClaims[k].boxIdx < 0) {
        // Listen-only sentinel, not an actual self-claim.
        continue;
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
      const KDClaim& c = *cit;

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
  // answer -- a rank-level "your rank won" cannot tell them apart. See KDClaim::proposerBoxIdx.

  std::vector<std::vector<KDVerdict>> verdictSendByRank(numRanks);

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
      const KDClaim& c = incomingClaims[k];

      // Compare IDENTITY (which specific claim was chosen), not key VALUE: an exact key tie
      // between two INDEPENDENT claims is real and not vanishingly rare here (BrownianWalker's
      // regular initial distribution can make two independently-built boxes end up with the exact
      // same volume), and comparing key value alone would tell BOTH tied claimants "you won" --
      // exactly the kind of double-count this whole protocol exists to prevent. Step 3 already
      // deterministically picked exactly one winner among any tie; only that specific claim,
      // identified by (rank, box index), is told so.
      const bool won = (c.proposerRank == winner->rank) && (c.proposerBoxIdx == winner->boxIdx);

      verdictSendByRank[c.proposerRank].push_back(KDVerdict{id, c.proposerBoxIdx, won});
    }

    idx = idxEnd;
  }

  const std::vector<KDVerdict> incomingVerdicts = kdExchangeByRank(verdictSendByRank);

  // wonForeignByBox[boxIdx] = every foreign member THAT SPECIFIC one of my own boxes was told it won.
  // Indexed directly by boxIdx (a dense range over myBoxes) rather than hashed -- each box's own
  // list is bounded by its own membership (at most ppc), so sorting it once and using binary_search
  // at the read site below avoids both a per-box hash-table allocation (unordered_set would add one
  // per box, of which there can be thousands) and, for large ppc, the O(box size^2) cost a linear
  // scan would have across all of a box's foreign-member lookups.

  std::vector<std::vector<ParticleID>> wonForeignByBox(myBoxes.size());

  for (const KDVerdict& v : incomingVerdicts) {
    if (v.won) {
      wonForeignByBox[v.claimantBoxIdx].push_back(v.memberID);
    }
  }

  for (std::vector<ParticleID>& won : wonForeignByBox) {
    std::sort(won.begin(), won.end());
  }

  // ==== STEP 5: assemble each of my boxes' provisional membership, commit or not ====
  std::vector<std::vector<KDCommit>> commitSendByRank(numRanks);

  for (std::size_t boxIdx = 0; boxIdx < myBoxes.size(); boxIdx++) {
    const RuntimeBox& box = myBoxes[boxIdx];

    std::vector<ParticleID> survivingLocal;
    std::vector<ParticleID> survivingForeign;

    for (const KDMember& m : box.members) {
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
        const MergeParticle<Packed>& p = findParticle(a_id);

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

      kdCheckCentroid(centroid, totalW, surviveBoxLo, surviveBoxHi);

      if (a_isPositionValid(centroid)) {
        committed = true;

        MergeParticle<Packed> merged;

        merged.position  = centroid;
        merged.weight    = totalW;
        merged.globalID  = a_allocateID();
        merged.ownerRank = myRank;
        merged.payload   = a_combine(payloads.data(), weights.data(), payloads.size());

        mergedResults.push_back(MergedResult{merged, box.patchIdx});

        for (const ParticleID id : survivingLocal) {
          // This rank owns these, and the outcome is already known, so nothing need be awaited.
          consumedIDs.push_back(id);
        }
      }
    }

    // Phase 3 outgoing: every foreign member this box nominally won needs to hear the outcome,
    // whether this box committed or not, since its owner is otherwise left not knowing whether to
    // delete it.
    //
    // This round cannot be collapsed into Phase 2. Winning a member's nominal argmin does not mean
    // the winning box will actually commit: it can independently lose OTHER members to OTHER,
    // unrelated lower-keyed boxes and end up under the merge threshold itself. Deleting on the
    // Phase 2 verdict alone then destroys that particle's weight, because no merged particle
    // anywhere ends up holding it. Concretely, with three ranks and box keys Kx < Ky < Kz:
    //
    //   Bx{p1(X), p2(Y)},  By{p2(Y), p3(Z)},  Bz{p3(Z), p4(X)}
    //
    // Argmin gives p1,p2 -> Bx (Bx beats By on p2) and p3 -> By (By beats Bz). Assembling actual
    // memberships: Bx = {p1,p2} commits; By has lost p2 to Bx and is left with {p3} alone, so it
    // does NOT commit; Bz has lost p3 to By and is left with {p4} alone, so it does not either.
    // Z learns from Phase 2 only that By nominally won p3 -- deleting p3 there would lose it
    // outright. Phase 3 is what tells Z to release it instead. The same holds for chains of any
    // length.
    for (const ParticleID id : survivingForeign) {
      commitSendByRank[findParticle(id).ownerRank].push_back(KDCommit{id, committed});
    }
  }

  // ==== STEP 6: Phase 3 -- commit/release, foreign members only ====
  const std::vector<KDCommit> incomingCommits = kdExchangeByRank(commitSendByRank);

  for (const KDCommit& c : incomingCommits) {
    if (c.committed) {
      consumedIDs.push_back(c.memberID);
    }
  }

  std::sort(consumedIDs.begin(), consumedIDs.end());

  // STEP 7 below finds entries by binary search, which silently misses on an unsorted range.
  CH_assert(std::is_sorted(consumedIDs.begin(), consumedIDs.end()));

  CH_STOP(t_carve);

  // ==== STEP 7: remove consumed particles ====
  CH_START(t_remove);
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

  CH_STOP(t_remove);

  // ==== STEP 8: place merged results ====
  CH_START(t_place);
  std::vector<std::vector<MergeParticle<Packed>>> scatterByDestRank(numRanks);

  auto insertHere = [&](const MergeParticle<Packed>& a_p, const int a_level, const DataIndex& a_din) {
    ParticleSoA<P, Traits>& leaf = a_particles[a_level][a_din];

    a_scatter(leaf, a_p);
  };

  for (const MergedResult& mr : mergedResults) {
    const PatchWork& pw = patchWork[mr.patchIdx];

    RealVect boxRealLo, boxRealHi;

    kdBoxRealBounds(boxRealLo, boxRealHi, pw.box, pw.dx, probLo);

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
      MayDay::Error("ParticleManagement::mergeKDCarve -- merged particle not found in any patch");
    }

    if (dst.rank == myRank) {
      const DataIndex din = a_amr.getLevelTiles(realm)[dst.level]->getMyGrids().at(dst.gridIndex);

      insertHere(mr.particle, dst.level, din);
    }
    else {
      MergeParticle<Packed> corrected = mr.particle;
      corrected.ownerRank             = dst.rank;

      scatterByDestRank[dst.rank].push_back(corrected);
    }
  }

  const std::vector<MergeParticle<Packed>> incomingScattered = kdExchangeByRank(scatterByDestRank);

  for (const MergeParticle<Packed>& p : incomingScattered) {
    const auto dst = a_particles.findDestination(p.position);

    if (!dst.valid || dst.rank != myRank) {
      MayDay::Error("ParticleManagement::mergeKDCarve -- incoming scattered particle not "
                    "found in any of this rank's own patches");
    }

    const DataIndex din = a_amr.getLevelTiles(realm)[dst.level]->getMyGrids().at(dst.gridIndex);

    MergeParticle<Packed> corrected = p;
    corrected.ownerRank             = myRank;

    insertHere(corrected, dst.level, din);
  }

  CH_STOP(t_place);
}

template <typename P,
          typename Packed,
          typename Traits,
          typename Gather,
          typename Combine,
          typename Scatter,
          typename Allocator,
          typename PosValid>
inline void
mergeKDPatch(ParticleContainer<P, Traits>& a_particles,
             EBAMRFAB&                     a_cellHistogram,
             EBAMRFAB&                     a_leafQuota,
             const AmrMesh&                a_amr,
             const int                     a_ppc,
             const Real                    a_splitWeightLeafDx,
             const Gather&                 a_gather,
             const Combine&                a_combine,
             const Scatter&                a_scatter,
             const Allocator&              a_allocateID,
             const PosValid&               a_isPositionValid)
{
  // The patch-local merge is the ghost-free case of mergeKDInterior() writing back into its own
  // container, so it delegates rather than repeating the body.
  //
  // mergeKDInterior() commits a leaf only when the leaf holds no ghost. The caller here must not
  // have filled a ghost halo -- a ghost would be merged locally while its true owner merges it too --
  // so every gathered particle is resident, every leaf is ghost-free, and that test admits them all.
  // Passing a_particles as its own destination puts the super-particles back where the patch merge
  // has always put them.
  mergeKDInterior<P, Packed, Traits>(a_particles,
                                     a_particles,
                                     a_cellHistogram,
                                     a_leafQuota,
                                     a_amr,
                                     a_ppc,
                                     a_splitWeightLeafDx,
                                     a_gather,
                                     a_combine,
                                     a_scatter,
                                     a_allocateID,
                                     a_isPositionValid);
}

template <typename P,
          typename Packed,
          typename Traits,
          typename Gather,
          typename Combine,
          typename Scatter,
          typename Allocator,
          typename PosValid>
inline void
mergeKDInterior(ParticleContainer<P, Traits>& a_particles,
                ParticleContainer<P, Traits>& a_interior,
                EBAMRFAB&                     a_cellHistogram,
                EBAMRFAB&                     a_leafQuota,
                const AmrMesh&                a_amr,
                const int                     a_ppc,
                const Real                    a_splitWeightLeafDx,
                const Gather&                 a_gather,
                const Combine&                a_combine,
                const Scatter&                a_scatter,
                const Allocator&              a_allocateID,
                const PosValid&               a_isPositionValid)
{
  using namespace detail;

  CH_TIMERS("ParticleManagement::mergeKDInterior");
  CH_TIMER("ParticleManagement::mergeKDInterior::build", t_build);
  CH_TIMER("ParticleManagement::mergeKDInterior::commit", t_commit);
  CH_TIMER("ParticleManagement::mergeKDInterior::remove_place", t_place);

  const std::string realm       = a_particles.getRealm();
  const int         finestLevel = a_amr.getFinestLevel();
  const RealVect    probLo      = a_amr.getProbLo();

  // A leaf committed from patch (lvl,din) is written to the SAME patch of a_interior, which is only
  // meaningful when the two containers share a layout.
  CH_assert(a_interior.getRealm() == realm);

  // Caller-owned per-cell scratch; see mergeKDCarve() for why the holders are not allocated here.
  EBAMRFAB& histogram = a_cellHistogram;
  EBAMRFAB& leafQuota = a_leafQuota;

  CH_assert(histogram[0]->nComp() == 1);
  CH_assert(leafQuota[0]->nComp() == 1);
  CH_assert(histogram[0]->ghostVect() >= IntVect::Unit);
  CH_assert(leafQuota[0]->ghostVect() >= IntVect::Unit);

  const int myRank = procID();

  // Per-patch leaf scratch, declared once and reused across every patch below rather than reallocated
  // inside the loop -- the same reuse pattern as the histogram/leafQuota holders above.
  // buildKDQuotaLeaves() clear()s it on entry, so each reuse retains the heap capacity grown by
  // earlier patches instead of starting from nothing.
  std::vector<KDLeaf> leaves;

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = a_amr.getGrids(realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const RealVect dx = a_amr.getDx()[lvl] * RealVect::Unit;

    const int nbox = dit.size();

    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      // Ghosts ARE gathered, unlike mergeKDPatch(). They are never merged here, but they must take
      // part in the partition: a ghost is what marks the leaf it lands in as contested, and the
      // per-cell quota has to see the cell's true occupancy to hand out the right number of slots.
      std::vector<MergeParticle<Packed>> combined;
      combined.reserve(leaf.size());

      for (std::size_t i = 0; i < leaf.size(); i++) {
        MergeParticle<Packed> p;

        p.position  = leaf.position(i);
        p.weight    = leaf.weight(i);
        p.globalID  = leaf.particleID(i);
        p.ownerRank = leaf.rankID(i);
        p.isGhost   = leaf.isGhost(i);
        p.payload   = a_gather(leaf, i);

        combined.push_back(p);
      }

      if (combined.empty()) {
        continue;
      }

      CH_START(t_build);

      // Ground truth for how crowded each cell really is -- see kdFillCellHistogram(). Both this and
      // the live quota are cell data, so they live in the mesh holders allocated once above rather
      // than in per-patch scratch. Their ghost cell covers the gathered ghosts, which sit up to one
      // cell outside this patch, and any leaf centroid that lands just outside it.
      FArrayBox& cellCounts = (*histogram[lvl])[din];
      FArrayBox& quota      = (*leafQuota[lvl])[din];

      kdFillCellHistogram(cellCounts, combined, probLo, dx);

      buildKDQuotaLeaves(combined, quota, leaves, a_ppc, a_splitWeightLeafDx, dx, probLo, cellCounts);

      CH_STOP(t_build);
      CH_START(t_commit);

      std::vector<ParticleID>            consumedIDs;
      std::vector<MergeParticle<Packed>> mergedResults;

      for (const KDLeaf& bl : leaves) {
        if (bl.hi - bl.lo < 2) {
          continue;
        }

        // Cheap geometric test before the O(members) ghost scan below.
        if (kdMaxAxisSpan(bl.boxLo, bl.boxHi, dx) > s_kdMaxLeafExtent) {
          continue;
        }

        // THE rule for this tier: commit iff the leaf holds no ghost. Every member is then a
        // particle physically resident in this patch, so no other patch can commit any of them --
        // a patch that draws one of them into a leaf of its own necessarily sees it as a ghost, and
        // that leaf is disqualified here by this very test. Committed leaves therefore never share
        // a member across patches, and no weight is counted twice. Boundary exposure is deliberately
        // NOT consulted: an exposed leaf that happens to hold no ghost is uncontested and merges
        // here, which is exactly what shrinks the skin relative to the carve.
        bool anyGhost = false;

        for (std::size_t idx = bl.lo; idx < bl.hi && !anyGhost; idx++) {
          anyGhost = combined[idx].isGhost;
        }

        if (anyGhost) {
          continue;
        }

        std::vector<Packed> payloads;
        std::vector<Real>   weights;
        RealVect            centroid = RealVect::Zero;
        Real                totalW   = 0.0;

        payloads.reserve(bl.hi - bl.lo);
        weights.reserve(bl.hi - bl.lo);

        for (std::size_t idx = bl.lo; idx < bl.hi; idx++) {
          const MergeParticle<Packed>& p = combined[idx];

          // Guaranteed by the anyGhost test above; merging a ghost here would double-count its
          // weight against its true owner.
          CH_assert(!p.isGhost);

          payloads.push_back(p.payload);
          weights.push_back(p.weight);
          centroid += p.weight * p.position;
          totalW += p.weight;
        }

        centroid /= totalW;

        kdCheckCentroid(centroid, totalW, bl.boxLo, bl.boxHi);

        if (!a_isPositionValid(centroid)) {
          continue;
        }

        MergeParticle<Packed> merged;

        merged.position  = centroid;
        merged.weight    = totalW;
        merged.globalID  = a_allocateID();
        merged.ownerRank = myRank;
        merged.isGhost   = false;
        merged.payload   = a_combine(payloads.data(), weights.data(), payloads.size());

        mergedResults.push_back(merged);

        for (std::size_t idx = bl.lo; idx < bl.hi; idx++) {
          consumedIDs.push_back(combined[idx].globalID);
        }
      }

      CH_STOP(t_commit);

      if (mergedResults.empty()) {
        continue;
      }

      CH_START(t_place);

      // Consumed originals leave a_particles; the super-particles they became are written to
      // a_interior instead of back here. What is left behind in a_particles is precisely the skin.
      std::sort(consumedIDs.begin(), consumedIDs.end());

      // Removal below finds entries by binary search, which silently misses on an unsorted range.
      CH_assert(std::is_sorted(consumedIDs.begin(), consumedIDs.end()));

      {
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

      RealVect boxRealLo, boxRealHi;

      kdBoxRealBounds(boxRealLo, boxRealHi, dbl[din], dx, probLo);

      for (const MergeParticle<Packed>& mr : mergedResults) {
        bool inOwnBox = true;

        for (int dir = 0; dir < SpaceDim && inOwnBox; dir++) {
          if (mr.position[dir] < boxRealLo[dir] || mr.position[dir] >= boxRealHi[dir]) {
            inOwnBox = false;
          }
        }

        if (inOwnBox) {
          a_scatter(a_interior[lvl][din], mr);

          continue;
        }

        // Every member was resident in this patch, so their weighted centroid is too -- this branch
        // is unreachable in exact arithmetic and exists only for round-off at a patch face. It must
        // still land on this rank: there is no exchange here to carry a particle anywhere else.
        const auto dst = a_particles.findDestination(mr.position);

        if (!dst.valid || dst.rank != myRank) {
          MayDay::Error("ParticleManagement::mergeKDInterior -- merged particle left this rank's own patches");
        }

        const DataIndex dstDin = a_amr.getLevelTiles(realm)[dst.level]->getMyGrids().at(dst.gridIndex);

        a_scatter(a_interior[dst.level][dstDin], mr);
      }

      CH_STOP(t_place);
    }
  }
}

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#endif
