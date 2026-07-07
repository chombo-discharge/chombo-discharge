/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_NearestNeighborParticleMergeImplem.H
  @brief  Implementation of CD_NearestNeighborParticleMerge.H
  @author Robert Marskar
*/

#ifndef CD_NEARESTNEIGHBORPARTICLEMERGEIMPLEM_H
#define CD_NEARESTNEIGHBORPARTICLEMERGEIMPLEM_H

// Std includes
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>
#include <SPMD.H>
#include <DataIterator.H>
#include <BoxIterator.H>

// Our includes
#include <CD_NearestNeighborParticleMerge.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace ParticleManagement {

  // TEMPORARY DIAGNOSTIC -- do not commit.
  static unsigned long long g_nnDiagQueries        = 0;
  static unsigned long long g_nnDiagNodeVisits     = 0;
  static unsigned long long g_nnDiagLeafCandidates = 0;
  static unsigned long long g_nnDiagTreeSizeSum    = 0;
  static unsigned long long g_nnDiagTreeCount      = 0;
  static unsigned long long g_nnDiagLeavesVisited  = 0;
  static unsigned long long g_nnDiagLeafSizeSum    = 0;
  static unsigned long long g_nnDiagLeafSizeCount  = 0;
  static unsigned long long g_nnDiagMaxLeafSize    = 0;

  // ---------------------------------------------------------------------------------------------
  // Private helpers (not declared in CD_NearestNeighborParticleMerge.H -- implementation details
  // only, consumed exclusively within this file).
  // ---------------------------------------------------------------------------------------------

  /**
    @brief Squared Euclidean distance between two physical positions.
    @param[in] a_a First position.
    @param[in] a_b Second position.
    @return Squared distance between a_a and a_b.
  */
  inline Real
  nnMergeDistSquared(const RealVect& a_a, const RealVect& a_b) noexcept
  {
    Real d2 = 0.0;
    for (int dir = 0; dir < SpaceDim; dir++) {
      const Real delta = a_a[dir] - a_b[dir];
      d2 += delta * delta;
    }
    return d2;
  }

  /**
    @brief Chebyshev (max-per-axis) distance between two cell-index vectors -- see
    CD_NearestNeighborParticleMerge.H's a_maxCellDistance docs.
    @param[in] a_a First cell index.
    @param[in] a_b Second cell index.
    @return Max-per-axis distance between a_a and a_b.
  */
  inline int
  nnMergeChebyshevCellDistance(const IntVect& a_a, const IntVect& a_b) noexcept
  {
    int dist = 0;
    for (int dir = 0; dir < SpaceDim; dir++) {
      dist = std::max(dist, std::abs(a_a[dir] - a_b[dir]));
    }
    return dist;
  }

  /**
    @brief a_maxCellDistance check for two positions that may be on DIFFERENT AMR levels.
    @details When both participants share a level, this is IDENTICAL to the original,
    already-validated definition: Chebyshev distance between each participant's own cell index,
    both computed with that level's dx. Comparing raw cell indices computed with two DIFFERENT dx
    values has no consistent meaning on its own, so for a genuinely cross-level pair, both
    positions are re-expressed in the FINER of the two levels' own cell-size units first (higher
    AMR level number = finer, by this codebase's convention) -- the more conservative/precise
    common reference -- before comparing. Reduces exactly to the same-level case when
    a_levelA == a_levelB.
    @param[in] a_posA         First participant's position.
    @param[in] a_levelA       First participant's AMR level.
    @param[in] a_posB         Second participant's position.
    @param[in] a_levelB       Second participant's AMR level.
    @param[in] a_probLo       Physical position of the domain's low corner.
    @param[in] a_dxByLevel    Cell spacing per AMR level.
    @param[in] a_maxCellDistance The cap (see CD_NearestNeighborParticleMerge.H's a_maxCellDistance
    docs).
    @return True iff the pair exceeds a_maxCellDistance.
  */
  inline bool
  nnMergeCrossLevelTooFar(const RealVect&              a_posA,
                          const int                    a_levelA,
                          const RealVect&              a_posB,
                          const int                    a_levelB,
                          const RealVect&              a_probLo,
                          const std::vector<RealVect>& a_dxByLevel,
                          const int                    a_maxCellDistance) noexcept
  {
    const int       refLevel = std::max(a_levelA, a_levelB);
    const RealVect& refDx    = a_dxByLevel[refLevel];
    const IntVect   keyA     = cellKeyOf(a_posA, a_probLo, refDx);
    const IntVect   keyB     = cellKeyOf(a_posB, a_probLo, refDx);
    return nnMergeChebyshevCellDistance(keyA, keyB) > a_maxCellDistance;
  }

  /**
    @brief Weighted-centroid merge of two positions/weights -- shared by resolveTrivialTier() and
    judgeProposals() so the formula is written exactly once.
    @param[in]  a_posA        First particle's position.
    @param[in]  a_weightA     First particle's weight.
    @param[in]  a_posB        Second particle's position.
    @param[in]  a_weightB     Second particle's weight.
    @param[out] a_mergedPos   Weighted-centroid position.
    @param[out] a_mergedWeight Summed weight.
  */
  inline void
  nnMergeWeightedCentroid(const RealVect& a_posA,
                          const Real      a_weightA,
                          const RealVect& a_posB,
                          const Real      a_weightB,
                          RealVect&       a_mergedPos,
                          Real&           a_mergedWeight) noexcept
  {
    a_mergedWeight = a_weightA + a_weightB;
    for (int dir = 0; dir < SpaceDim; dir++) {
      a_mergedPos[dir] = (a_posA[dir] * a_weightA + a_posB[dir] * a_weightB) / a_mergedWeight;
    }
  }

  /**
    @brief Generic alltoallv-style exchange: send a per-destination-rank bucket of trivially-
    copyable records, return the flat concatenation of everything received.
    @details This is the "new general-purpose alltoallv-style utility" the header's own docs
    flagged as an open implementation question -- a minimal, self-contained one, matching the
    packing/MPI_Alltoall(byte counts)/MPI_Alltoallv(payload) style already used by
    ParticleContainer's own cross-rank particle exchange (CD_ParticleContainerImplem.H). In a
    serial (non-MPI) build, or when there is only one rank, this degenerates to "my own bucket,
    unchanged".
    @tparam T Trivially-copyable record type (NNMergeProposal<Packed>, NNMergeVerdict, or
    NNMergeParticle<Packed> -- all satisfy this given Packed itself is required to).
    @param[in] a_sendByRank One bucket per destination rank (size == number of ranks).
    @return Flat vector of every record received from every rank (including this rank's own
    bucket addressed to itself).
  */
  template <typename T>
  inline std::vector<T>
  nnMergeExchangeByRank(const std::vector<std::vector<T>>& a_sendByRank)
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

  // ---------------------------------------------------------------------------------------------
  // Public API (declared in CD_NearestNeighborParticleMerge.H).
  // ---------------------------------------------------------------------------------------------

  inline IntVect
  cellKeyOf(const RealVect& a_position, const RealVect& a_probLo, const RealVect& a_dx) noexcept
  {
    IntVect iv;
    for (int dir = 0; dir < SpaceDim; dir++) {
      iv[dir] = static_cast<int>(std::floor((a_position[dir] - a_probLo[dir]) / a_dx[dir]));
    }
    return iv;
  }

  inline bool
  isBoundaryExposed(const IntVect& a_cellIndex, const ParticleGhostMask& a_ghostMask) noexcept
  {
    return a_ghostMask.numTargets(a_cellIndex) > 0;
  }

  /**
    @brief Component-wise floor-division block index for a cell -- pure arithmetic, no lookup.
    @details Deliberately duplicated here rather than depending on MergeBlockGrid (which lives in
    Source/AmrMesh): this header stays AmrMesh-agnostic by design (see its own docs), so only plain
    primitives (blockSize, a Box) cross that boundary -- never the MergeBlockGrid type itself. This
    one-line computation is translation-invariant and stateless, so duplicating it here rather than
    introducing a dependency is the cheaper and more honest option.
    @param[in] a_cell      Cell index.
    @param[in] a_blockSize Block size in cells.
    @return floor(a_cell / a_blockSize), component-wise.
  */
  inline IntVect
  nnMergeBlockIndexOf(const IntVect& a_cell, const int a_blockSize) noexcept
  {
    IntVect block;
    for (int dir = 0; dir < SpaceDim; dir++) {
      block[dir] = (a_cell[dir] >= 0) ? (a_cell[dir] / a_blockSize) : -(((-a_cell[dir] - 1) / a_blockSize) + 1);
    }
    return block;
  }

  template <typename Packed>
  inline void
  buildNNCellBuckets(const std::vector<NNMergeParticle<Packed>>& a_localValid,
                     const std::vector<NNMergeParticle<Packed>>& a_ghosts,
                     const RealVect&                             a_probLo,
                     const RealVect&                             a_dx,
                     const int                                   a_blockSize,
                     const Box&                                  a_haloBlockBox,
                     NNCellBuckets<Packed>&                      a_buckets) noexcept
  {
    CH_TIME("ParticleManagement::buildNNCellBuckets");

    a_buckets.sorted.clear();
    a_buckets.cellOffsets.clear();
    a_buckets.blockCounts.assign(a_haloBlockBox.numPts(), 0);
    a_buckets.haloBlockBox    = a_haloBlockBox;
    a_buckets.blockSize       = a_blockSize;
    a_buckets.observedCellBox = Box();

    const std::size_t nLocal = a_localValid.size();
    const std::size_t nGhost = a_ghosts.size();
    const std::size_t total  = nLocal + nGhost;

    if (total == 0) {
      return;
    }

    // Pass 1: tally per-cell counts (as a temporary map -- number of DISTINCT cells is typically
    // far smaller than a dense range would require, hence a hash map rather than a dense array)
    // and roll counts up into blockCounts at O(1) extra cost per particle. Also tracks the
    // observed cell range (min/max per direction), for observedCellBox -- see its own docs.
    std::unordered_map<NNWalkCell, std::size_t, LevelTiles::TileHasher> tally;
    tally.reserve(total);

    IntVect obsLo   = IntVect::Zero;
    IntVect obsHi   = IntVect::Zero;
    bool    haveObs = false;

    auto tallyOne = [&](const RealVect& a_pos) {
      const IntVect cell  = cellKeyOf(a_pos, a_probLo, a_dx);
      const IntVect block = nnMergeBlockIndexOf(cell, a_blockSize);

      tally[cell]++;

      if (!haveObs) {
        obsLo   = cell;
        obsHi   = cell;
        haveObs = true;
      }
      else {
        for (int dir = 0; dir < SpaceDim; dir++) {
          obsLo[dir] = std::min(obsLo[dir], cell[dir]);
          obsHi[dir] = std::max(obsHi[dir], cell[dir]);
        }
      }

      if (a_haloBlockBox.contains(block)) {
        a_buckets.blockCounts[static_cast<std::size_t>(a_haloBlockBox.index(block))]++;
      }
      // else: this cell's block lies outside the expected 1-block halo -- see MergeBlockGrid's own
      // docs for the geometric assumption this relies on. Not fatal here (the cell is still
      // correctly tallied/bucketed below, and observedCellBox still correctly bounds it, so ring
      // expansion in findNearestNeighborCandidates() remains correct even if this assumption is
      // ever violated); the caller is responsible for asserting this never happens in debug builds.
    };

    for (const auto& p : a_localValid) {
      tallyOne(p.position);
    }
    for (const auto& p : a_ghosts) {
      tallyOne(p.position);
    }

    if (haveObs) {
      a_buckets.observedCellBox = Box(obsLo, obsHi);
    }

    // Convert tallies into contiguous [begin, end) ranges.
    a_buckets.sorted.resize(total);
    std::unordered_map<NNWalkCell, std::size_t, LevelTiles::TileHasher> cursor;
    cursor.reserve(tally.size());

    std::size_t running = 0;
    for (const auto& kv : tally) {
      a_buckets.cellOffsets[kv.first] = std::make_pair(running, running + kv.second);
      cursor[kv.first]                = running;
      running += kv.second;
    }

    // Pass 2: scatter every particle's indexed entry into its cell's next open slot.
    for (std::size_t i = 0; i < nLocal; i++) {
      const IntVect cell               = cellKeyOf(a_localValid[i].position, a_probLo, a_dx);
      a_buckets.sorted[cursor[cell]++] = NNIndexedParticle{i, true};
    }
    for (std::size_t i = 0; i < nGhost; i++) {
      const IntVect cell               = cellKeyOf(a_ghosts[i].position, a_probLo, a_dx);
      a_buckets.sorted[cursor[cell]++] = NNIndexedParticle{i, false};
    }
  }

  /**
    @brief Minimum possible squared distance from a point to any point inside an axis-aligned box.
    @param[in] a_point  Query point.
    @param[in] a_aabbLo Box's low corner.
    @param[in] a_aabbHi Box's high corner.
    @return 0 if a_point is inside the box; otherwise the squared distance to the nearest face.
  */
  inline Real
  nnMergePointToAABBDistSquared(const RealVect& a_point, const RealVect& a_aabbLo, const RealVect& a_aabbHi) noexcept
  {
    Real d2 = 0.0;
    for (int dir = 0; dir < SpaceDim; dir++) {
      Real delta = 0.0;
      if (a_point[dir] < a_aabbLo[dir]) {
        delta = a_aabbLo[dir] - a_point[dir];
      }
      else if (a_point[dir] > a_aabbHi[dir]) {
        delta = a_point[dir] - a_aabbHi[dir];
      }
      d2 += delta * delta;
    }
    return d2;
  }

  /**
    @brief Recursively build a node's count-median-split subtree in place, over the range given.
    @details Only NNIndexedParticle entries in a_sorted[a_lo, a_hi) are ever reordered (via
    std::nth_element against the longest-axis COUNT median of the current node's own particles) --
    never NNMergeParticle data. See NNSubNode/buildNNSpatialIndex()'s own docs for why a count
    median, not a pure geometric midpoint (this function's own original design): a midpoint split
    chooses its boundary from the node's bounding box alone, blind to where the particles actually
    are, and for data clustered unevenly across a wide bounding box (crowded cells beside sparse
    ones -- confirmed the norm, not an edge case, for a whole-patch tree) it can carve off only a
    couple of outliers per level while leaving a dense cluster nearly intact, producing wildly
    unbalanced leaves (empirically: leaf sizes from 1 up to 3211 particles in one test, against a
    target of 32). std::nth_element is O(n) per level (same order as std::partition, just a larger
    constant) and guarantees an exact (within 1) 50/50 split by COUNT regardless of how the
    underlying positions are distributed, so it cannot degenerate this way. Recursion stops at
    a_subPartitionThresh or a_maxDepth.
    @tparam Packed See NNMergeParticle.
    @param[in,out] a_sorted             The patch's full index array (see NNCellBuckets::sorted);
                                         only the sub-range [a_lo, a_hi) is ever reordered.
    @param[in]     a_localValid         This patch's own full valid-particle set.
    @param[in]     a_ghosts             Every ghost shipped to this patch this round.
    @param[in]     a_lo                 Start of this node's half-open range into a_sorted.
    @param[in]     a_hi                 End of this node's half-open range into a_sorted.
    @param[in]     a_aabbLo             This node's own bounding box, low corner.
    @param[in]     a_aabbHi             This node's own bounding box, high corner.
    @param[in]     a_subPartitionThresh Occupancy at/below which a node becomes a leaf.
    @param[in]     a_depth              This node's depth (0 for the cell's own root node).
    @param[in]     a_maxDepth           Hard recursion depth cap (defensive; see
                                         buildNNSpatialIndex()'s docs).
    @param[in,out] a_nodes              The flat node vector this tree is being built into; node 0
                                         is the tree's root.
    @return Index of the node just created (in a_nodes) for this range.
  */
  template <typename Packed>
  inline int
  buildNNSubNodeRecursive(std::vector<NNIndexedParticle>&             a_sorted,
                          const std::vector<NNMergeParticle<Packed>>& a_localValid,
                          const std::vector<NNMergeParticle<Packed>>& a_ghosts,
                          const std::size_t                           a_lo,
                          const std::size_t                           a_hi,
                          const RealVect&                             a_aabbLo,
                          const RealVect&                             a_aabbHi,
                          const int                                   a_subPartitionThresh,
                          const int                                   a_depth,
                          const int                                   a_maxDepth,
                          std::vector<NNSubNode>&                     a_nodes) noexcept
  {
    const std::size_t nodeIdx = a_nodes.size();
    a_nodes.push_back(NNSubNode{a_lo, a_hi, a_aabbLo, a_aabbHi, -1, -1});

    const std::size_t count = a_hi - a_lo;
    if (static_cast<int>(count) <= a_subPartitionThresh || a_depth >= a_maxDepth) {
      g_nnDiagLeafSizeSum += count;
      g_nnDiagLeafSizeCount++;
      g_nnDiagMaxLeafSize = std::max(g_nnDiagMaxLeafSize, static_cast<unsigned long long>(count));
      return static_cast<int>(nodeIdx); // leaf
    }

    int  axis = 0;
    Real span = a_aabbHi[0] - a_aabbLo[0];
    for (int dir = 1; dir < SpaceDim; dir++) {
      const Real s = a_aabbHi[dir] - a_aabbLo[dir];
      if (s > span) {
        span = s;
        axis = dir;
      }
    }
    auto getPos = [&](const NNIndexedParticle& a_entry) -> const RealVect& {
      return a_entry.isLocal ? a_localValid[a_entry.idx].position : a_ghosts[a_entry.idx].position;
    };

    // Count-median split via nth_element: after this call, every entry in [a_lo, splitPos) has an
    // axis-value <= the entry now AT splitPos, and every entry in [splitPos, a_hi) has one >= it --
    // an exact (within 1) 50/50 split by COUNT, regardless of how the underlying positions are
    // distributed (see this function's own docs for why this replaced a pure geometric midpoint).
    const std::size_t splitPos = a_lo + count / 2;

    std::nth_element(a_sorted.begin() + static_cast<std::ptrdiff_t>(a_lo),
                     a_sorted.begin() + static_cast<std::ptrdiff_t>(splitPos),
                     a_sorted.begin() + static_cast<std::ptrdiff_t>(a_hi),
                     [&](const NNIndexedParticle& a_a, const NNIndexedParticle& a_b) {
                       return getPos(a_a)[axis] < getPos(a_b)[axis];
                     });

    // Defensive only -- nth_element's split is exact by construction (within 1), so this can only
    // trigger for a degenerate a_subPartitionThresh of 0 with count == 1, not a realistic setting.
    if (splitPos <= a_lo || splitPos >= a_hi) {
      g_nnDiagLeafSizeSum += count;
      g_nnDiagLeafSizeCount++;
      g_nnDiagMaxLeafSize = std::max(g_nnDiagMaxLeafSize, static_cast<unsigned long long>(count));
      return static_cast<int>(nodeIdx); // leaf
    }

    // Safe (never too loose in the wrong direction) child AABBs: nth_element's own guarantee (see
    // above) means using the split entry's own axis-value as the shared boundary correctly bounds
    // BOTH children, even though it is not necessarily either child's tightest possible extent
    // (ties at the split value can land on either side).
    const Real splitValue = getPos(a_sorted[splitPos])[axis];

    RealVect leftHi  = a_aabbHi;
    leftHi[axis]     = splitValue;
    RealVect rightLo = a_aabbLo;
    rightLo[axis]    = splitValue;

    const int left  = buildNNSubNodeRecursive(a_sorted,
                                             a_localValid,
                                             a_ghosts,
                                             a_lo,
                                             splitPos,
                                             a_aabbLo,
                                             leftHi,
                                             a_subPartitionThresh,
                                             a_depth + 1,
                                             a_maxDepth,
                                             a_nodes);
    const int right = buildNNSubNodeRecursive(a_sorted,
                                              a_localValid,
                                              a_ghosts,
                                              splitPos,
                                              a_hi,
                                              rightLo,
                                              a_aabbHi,
                                              a_subPartitionThresh,
                                              a_depth + 1,
                                              a_maxDepth,
                                              a_nodes);

    a_nodes[nodeIdx].leftChild  = left;
    a_nodes[nodeIdx].rightChild = right;

    return static_cast<int>(nodeIdx);
  }

  template <typename Packed>
  inline void
  buildNNSpatialIndex(const std::vector<NNMergeParticle<Packed>>& a_localValid,
                      const std::vector<NNMergeParticle<Packed>>& a_ghosts,
                      const RealVect&                             a_probLo,
                      const RealVect&                             a_dx,
                      const int                                   a_blockSize,
                      const Box&                                  a_haloBlockBox,
                      const int                                   a_subPartitionThresh,
                      const int                                   a_numParticlesPerCellThresh,
                      const long                                  a_numCellsInPatch,
                      NNSpatialIndex<Packed>&                     a_index) noexcept
  {
    CH_TIME("ParticleManagement::buildNNSpatialIndex");

    buildNNCellBuckets(a_localValid, a_ghosts, a_probLo, a_dx, a_blockSize, a_haloBlockBox, a_index.buckets);

    a_index.tree.clear();

    const std::size_t total = a_index.buckets.sorted.size();
    if (total == 0) {
      return;
    }

    // Average-density pre-filter -- see this function's own docs. A cheap, O(1) "is this even
    // worth it" check mirroring what equal_weight_kd effectively gets for free via its per-cell
    // dispatch: if this patch's own particles, spread evenly over its own cells, wouldn't even
    // reach the crowding threshold on average, building a tree that (overwhelmingly likely) never
    // gets walked by a single query would be pure waste.
    if (static_cast<long>(total) <= a_numCellsInPatch * static_cast<long>(a_numParticlesPerCellThresh)) {
      return;
    }

    // Root AABB: the physical extent of every particle actually present (via observedCellBox,
    // which already tracks this in cell-index space -- see its own docs), NOT just this patch's
    // own box, since ghosts can extend beyond it.
    const IntVect& obsLo = a_index.buckets.observedCellBox.smallEnd();
    const IntVect& obsHi = a_index.buckets.observedCellBox.bigEnd();

    RealVect aabbLo, aabbHi;
    for (int dir = 0; dir < SpaceDim; dir++) {
      aabbLo[dir] = a_probLo[dir] + obsLo[dir] * a_dx[dir];
      aabbHi[dir] = a_probLo[dir] + (obsHi[dir] + 1) * a_dx[dir];
    }

    constexpr int maxDepth = 24; // generous -- reaching this depth is not a realistic occupancy

    buildNNSubNodeRecursive(a_index.buckets.sorted,
                            a_localValid,
                            a_ghosts,
                            0,
                            total,
                            aabbLo,
                            aabbHi,
                            a_subPartitionThresh,
                            0,
                            maxDepth,
                            a_index.tree);
    g_nnDiagTreeSizeSum += a_index.tree.size();
    g_nnDiagTreeCount++;
  }

  /**
    @brief Recursively walk the whole-patch tree, pruning subtrees by point-to-AABB distance and
    brute-forcing leaves.
    @details a_shouldPrune is re-evaluated fresh at every node (it closes over the caller's own,
    continuously-improving candidate buffer), so pruning gets more aggressive as better candidates
    are found deeper in the tree -- never a stale snapshot from before the walk started.
    @tparam Packed See NNMergeParticle.
    @tparam ConsiderCallback Deliberately a template parameter (see PayloadCombine's own docs for
    the inlining rationale) -- called once per non-consumed, non-self candidate encountered.
    @tparam PruneCheck Deliberately a template parameter, same reason -- given a node's AABB,
    returns true iff the whole subtree can be skipped.
    @param[in] a_nodes             The patch's whole-patch tree (see NNSpatialIndex::tree).
    @param[in] a_nodeIdx           Index of the node to visit (0 for the tree's root).
    @param[in] a_queryPos          The query particle's own position.
    @param[in] a_qi                The query particle's own index into a_localValid (for
                                    self-exclusion).
    @param[in] a_localValid        This patch's own full valid-particle set.
    @param[in] a_ghosts            Every ghost shipped to this patch this round.
    @param[in] a_sorted            The patch's full index array (see NNCellBuckets::sorted).
    @param[in] a_consumedIDs       Lazy-deletion skip set (see findNearestNeighborCandidates()'s
                                    own docs).
    @param[in] a_considerCandidate Called once per non-consumed, non-self candidate found in a leaf.
    @param[in] a_shouldPrune       Given a node's AABB, returns true iff its whole subtree can be
                                    skipped.
  */
  template <typename Packed, typename ConsiderCallback, typename PruneCheck>
  inline void
  nnMergeWalkSubNodes(const std::vector<NNSubNode>&               a_nodes,
                      const int                                   a_nodeIdx,
                      const RealVect&                             a_queryPos,
                      const std::size_t                           a_qi,
                      const std::vector<NNMergeParticle<Packed>>& a_localValid,
                      const std::vector<NNMergeParticle<Packed>>& a_ghosts,
                      const std::vector<NNIndexedParticle>&       a_sorted,
                      const std::set<ParticleID>&                 a_consumedIDs,
                      const ConsiderCallback&                     a_considerCandidate,
                      const PruneCheck&                           a_shouldPrune) noexcept
  {
    g_nnDiagNodeVisits++;
    const NNSubNode& node = a_nodes[static_cast<std::size_t>(a_nodeIdx)];

    if (a_shouldPrune(node.aabbLo, node.aabbHi)) {
      return;
    }

    if (node.leftChild < 0) {
      g_nnDiagLeavesVisited++;
      for (std::size_t s = node.lo; s < node.hi; s++) {
        g_nnDiagLeafCandidates++;
        const NNIndexedParticle& entry = a_sorted[s];

        if (entry.isLocal && entry.idx == a_qi) {
          continue; // self
        }

        const NNMergeParticle<Packed>& cand = entry.isLocal ? a_localValid[entry.idx] : a_ghosts[entry.idx];

        if (a_consumedIDs.count(cand.globalID) > 0) {
          continue;
        }

        a_considerCandidate(nnMergeDistSquared(a_queryPos, cand.position), entry.idx, entry.isLocal);
      }
      return;
    }

    // Visit whichever child is nearer to the query FIRST -- standard nearest-neighbor KD-tree
    // practice, and NECESSARY here (not just a minor optimization): visiting an arbitrary,
    // unconditional left-then-right order meant a query could descend deep into a subtree nowhere
    // near it before finding any candidate at all, leaving a_considerCandidate's own "worst-of-K"
    // baseline too weak to prune anything for a long time afterward. This mattered little for the
    // small, single-cell trees this walk originally served (a handful of leaves at most), but
    // becomes the dominant cost once the same walk serves a whole-patch tree spanning thousands of
    // particles (confirmed empirically: average node visits per query matched the ENTIRE tree size
    // before this fix). Visiting the near child first means its own results are already available
    // -- via a_considerCandidate's closure over the caller's "best" buffer -- by the time the far
    // child's own a_shouldPrune check runs, at the START of its own recursive call, letting it
    // correctly skip the far side whenever the near side already found something good enough.
    const NNSubNode& left  = a_nodes[static_cast<std::size_t>(node.leftChild)];
    const NNSubNode& right = a_nodes[static_cast<std::size_t>(node.rightChild)];

    const Real leftD2  = nnMergePointToAABBDistSquared(a_queryPos, left.aabbLo, left.aabbHi);
    const Real rightD2 = nnMergePointToAABBDistSquared(a_queryPos, right.aabbLo, right.aabbHi);

    const int nearChild = (leftD2 <= rightD2) ? node.leftChild : node.rightChild;
    const int farChild  = (leftD2 <= rightD2) ? node.rightChild : node.leftChild;

    nnMergeWalkSubNodes(a_nodes,
                        nearChild,
                        a_queryPos,
                        a_qi,
                        a_localValid,
                        a_ghosts,
                        a_sorted,
                        a_consumedIDs,
                        a_considerCandidate,
                        a_shouldPrune);
    nnMergeWalkSubNodes(a_nodes,
                        farChild,
                        a_queryPos,
                        a_qi,
                        a_localValid,
                        a_ghosts,
                        a_sorted,
                        a_consumedIDs,
                        a_considerCandidate,
                        a_shouldPrune);
  }

  template <typename Packed>
  inline void
  findNearestNeighborCandidates(const std::vector<NNMergeParticle<Packed>>&                a_localValid,
                                const std::vector<NNMergeParticle<Packed>>&                a_ghosts,
                                const std::set<ParticleID>&                                a_consumedIDs,
                                const NNSpatialIndex<Packed>&                              a_spatialIndex,
                                const std::unordered_map<NNCellKey, int, NNCellKeyHasher>& a_liveCellCount,
                                const RealVect&                                            a_probLo,
                                const RealVect&                                            a_dx,
                                const int                                                  a_numParticlesPerCellThresh,
                                std::vector<NNMergeEdge<Packed>>&                          a_edges,
                                const int                                                  a_maxFallbackCandidates,
                                const std::optional<int>                                   a_maxCellDistance) noexcept
  {
    CH_TIME("ParticleManagement::findNearestNeighborCandidates");

    a_edges.clear();

    const std::size_t nLocal = a_localValid.size();
    const int         K      = 1 + a_maxFallbackCandidates;

    const NNCellBuckets<Packed>& buckets = a_spatialIndex.buckets;

    // a_maxCellDistance's physical bound, computed once (patch/round-level, not per query) -- see
    // this function's own a_maxCellDistance docs for why the sqrt(SpaceDim) factor is required for
    // safety (a Euclidean-distance prune must never be tighter than the true per-axis Chebyshev
    // cutoff resolveTrivialTier()/judgeProposals() apply downstream, or it could wrongly exclude a
    // candidate they would still have accepted).
    std::optional<Real> maxDistBoundSquared;
    if (a_maxCellDistance.has_value()) {
      Real minDx = a_dx[0];
      for (int dir = 1; dir < SpaceDim; dir++) {
        minDx = std::min(minDx, a_dx[dir]);
      }
      const Real bound = static_cast<Real>(a_maxCellDistance.value()) * minDx * std::sqrt(static_cast<Real>(SpaceDim));
      maxDistBoundSquared = bound * bound;
    }

    struct RankedCandidate
    {
      Real        d2;
      std::size_t idx;
      bool        isLocal;
    };

    for (std::size_t qi = 0; qi < nLocal; qi++) {
      const NNMergeParticle<Packed>& q = a_localValid[qi];

      // Lazy deletion -- a_localValid is the FULL per-patch set (see this function's own
      // contract-change docs), so a query already consumed by an earlier pass/query this round
      // must be skipped explicitly here (the old design's pre-filtered aliveLocal made this
      // implicit).
      if (a_consumedIDs.count(q.globalID) > 0) {
        continue;
      }

      // a_dx is this patch's own level's dx -- correct here since every query in a_localValid is
      // necessarily on this same level (see this function's a_dx docs).
      const IntVect qCell = cellKeyOf(q.position, a_probLo, a_dx);
      const auto    cntIt = a_liveCellCount.find(NNCellKey{q.level, qCell});
      const int     count = (cntIt != a_liveCellCount.end()) ? cntIt->second : 0;

      if (count <= a_numParticlesPerCellThresh) {
        continue;
      }

      std::vector<RankedCandidate> best;
      best.reserve(static_cast<std::size_t>(K) + 1);

      auto considerCandidate = [&](const Real a_d2, const std::size_t a_idx, const bool a_isLocal) {
        if (static_cast<int>(best.size()) == K && a_d2 >= best.back().d2) {
          return;
        }
        const auto insertPos = std::upper_bound(best.begin(),
                                                best.end(),
                                                a_d2,
                                                [](const Real a_val, const RankedCandidate& a_c) {
                                                  return a_val < a_c.d2;
                                                });
        best.insert(insertPos, RankedCandidate{a_d2, a_idx, a_isLocal});
        if (static_cast<int>(best.size()) > K) {
          best.pop_back();
        }
      };

      g_nnDiagQueries++;

      // Single, whole-patch tree walk -- see NNSpatialIndex::tree/buildNNSpatialIndex()'s own
      // docs. Empty either because this patch had no particles at all, or because
      // buildNNSpatialIndex()'s average-density pre-filter decided this patch wasn't worth a tree
      // -- the latter is a heuristic, not a guarantee (see its own docs), so a query CAN still
      // legitimately exist here despite it.
      if (!a_spatialIndex.tree.empty()) {
        auto shouldPrune = [&](const RealVect& a_aabbLo, const RealVect& a_aabbHi) -> bool {
          const Real d2 = nnMergePointToAABBDistSquared(q.position, a_aabbLo, a_aabbHi);
          if (maxDistBoundSquared.has_value() && d2 > maxDistBoundSquared.value()) {
            return true;
          }
          if (static_cast<int>(best.size()) < K) {
            return false;
          }
          return d2 >= best.back().d2;
        };
        nnMergeWalkSubNodes(a_spatialIndex.tree,
                            0,
                            q.position,
                            qi,
                            a_localValid,
                            a_ghosts,
                            buckets.sorted,
                            a_consumedIDs,
                            considerCandidate,
                            shouldPrune);
      }
      else if (!buckets.sorted.empty()) {
        // Fallback: no tree was built for this patch (the average-density pre-filter said skip),
        // but this specific query turned out to be over-threshold anyway -- always correct, just
        // not fast; expected to almost never trigger (see buildNNSpatialIndex()'s own docs).
        for (const NNIndexedParticle& entry : buckets.sorted) {
          if (entry.isLocal && entry.idx == qi) {
            continue; // self
          }

          const NNMergeParticle<Packed>& cand = entry.isLocal ? a_localValid[entry.idx] : a_ghosts[entry.idx];

          if (a_consumedIDs.count(cand.globalID) > 0) {
            continue;
          }

          considerCandidate(nnMergeDistSquared(q.position, cand.position), entry.idx, entry.isLocal);
        }
      }

      if (best.empty()) {
        continue;
      }

      NNMergeEdge<Packed> edge;
      edge.distanceSquared  = best[0].d2;
      edge.queryID          = q.globalID;
      edge.candidateID      = best[0].isLocal ? a_localValid[best[0].idx].globalID : a_ghosts[best[0].idx].globalID;
      edge.candidateIsLocal = best[0].isLocal;

      for (std::size_t k = 1; k < best.size(); k++) {
        NNFallbackCandidate fc;
        fc.candidateID      = best[k].isLocal ? a_localValid[best[k].idx].globalID : a_ghosts[best[k].idx].globalID;
        fc.distanceSquared  = best[k].d2;
        fc.candidateIsLocal = best[k].isLocal;
        edge.fallbackCandidates.push_back(fc);
      }

      a_edges.push_back(std::move(edge));
    }
  }

  template <typename Packed, typename Combine, typename PosValid, typename Allocator, typename OnMerge>
  inline void
  resolveTrivialTier(const std::vector<NNMergeEdge<Packed>>&              a_edges,
                     const std::map<ParticleID, bool>&                    a_exposed,
                     const std::map<ParticleID, NNMergeParticle<Packed>>& a_particlesByID,
                     const Combine&                                       a_combine,
                     std::unordered_map<NNCellKey, int, NNCellKeyHasher>& a_liveCellCount,
                     std::set<ParticleID>&                                a_consumedIDs,
                     const int                                            a_numParticlesPerCellThresh,
                     std::set<ParticleID>&                                a_hasOutgoingCommitment,
                     std::vector<NNMergeResult<Packed>>&                  a_results,
                     const RealVect&                                      a_probLo,
                     const std::vector<RealVect>&                         a_dxByLevel,
                     const PosValid&                                      a_isPositionValid,
                     const Allocator&                                     a_allocateID,
                     const int                                            a_maxFallbackCandidates,
                     const std::optional<int>                             a_maxCellDistance,
                     unsigned long long*                                  a_crossLevelMergeCount,
                     const OnMerge&                                       a_onMergeCommitted) noexcept
  {
    CH_TIME("ParticleManagement::resolveTrivialTier");

    a_results.clear();

    // 1. Sort a copy of a_edges by (distanceSquared, queryID, candidateID) -- see this function's
    //    docs for why the deterministic tiebreak is required, not optional.
    std::vector<const NNMergeEdge<Packed>*> sorted;
    sorted.reserve(a_edges.size());
    for (const auto& e : a_edges) {
      sorted.push_back(&e);
    }
    std::sort(sorted.begin(), sorted.end(), [](const NNMergeEdge<Packed>* a_a, const NNMergeEdge<Packed>* a_b) {
      if (a_a->distanceSquared != a_b->distanceSquared) {
        return a_a->distanceSquared < a_b->distanceSquared;
      }
      if (a_a->queryID != a_b->queryID) {
        return a_a->queryID < a_b->queryID;
      }
      return a_a->candidateID < a_b->candidateID;
    });

    for (const NNMergeEdge<Packed>* edgePtr : sorted) {
      const NNMergeEdge<Packed>& edge = *edgePtr;

      // a. Stale query -- already consumed by an earlier, closer edge this same pass.
      if (a_consumedIDs.count(edge.queryID) > 0) {
        continue;
      }

      const auto queryExposedIt = a_exposed.find(edge.queryID);
      const bool queryExposed   = (queryExposedIt != a_exposed.end()) && queryExposedIt->second;

      const auto queryDataIt = a_particlesByID.find(edge.queryID);
      if (queryDataIt == a_particlesByID.end()) {
        continue;
      }
      const NNMergeParticle<Packed>& query = queryDataIt->second;

      bool resolved           = false; // an outcome (commit, or a legitimate "leave unmatched") was reached
      bool originalStale      = false;
      bool originalIneligible = false; // too far or invalid position, evaluated on the ORIGINAL candidate

      const int numAttempts = a_maxFallbackCandidates + 1;
      for (int attempt = 0; attempt < numAttempts; attempt++) {
        ParticleID candidateID;
        bool       candidateIsLocal;

        if (attempt == 0) {
          candidateID      = edge.candidateID;
          candidateIsLocal = edge.candidateIsLocal;
        }
        else {
          const std::size_t fbIdx = static_cast<std::size_t>(attempt - 1);
          if (fbIdx >= edge.fallbackCandidates.size()) {
            break;
          }
          candidateID      = edge.fallbackCandidates[fbIdx].candidateID;
          candidateIsLocal = edge.fallbackCandidates[fbIdx].candidateIsLocal;
        }

        const bool candidateStale = candidateIsLocal && (a_consumedIDs.count(candidateID) > 0);
        const bool candidateBusy  = candidateIsLocal && !candidateStale &&
                                   (a_hasOutgoingCommitment.count(candidateID) > 0);

        bool candidateExposed = false;
        if (candidateIsLocal) {
          const auto exposedIt = a_exposed.find(candidateID);
          candidateExposed     = (exposedIt != a_exposed.end()) && exposedIt->second;
        }

        const auto candDataIt      = a_particlesByID.find(candidateID);
        const bool candidateExists = candDataIt != a_particlesByID.end();

        bool     tooFar     = false;
        bool     invalidPos = false;
        RealVect mergedPos(D_DECL(0.0, 0.0, 0.0));
        Real     mergedWeight = 0.0;

        if (!candidateStale && candidateExists) {
          const NNMergeParticle<Packed>& cand = candDataIt->second;

          if (a_maxCellDistance.has_value()) {
            tooFar = nnMergeCrossLevelTooFar(query.position,
                                             query.level,
                                             cand.position,
                                             cand.level,
                                             a_probLo,
                                             a_dxByLevel,
                                             a_maxCellDistance.value());
          }

          if (!tooFar) {
            nnMergeWeightedCentroid(query.position, query.weight, cand.position, cand.weight, mergedPos, mergedWeight);
            invalidPos = !a_isPositionValid(mergedPos);
          }
        }

        const bool ineligible = tooFar || invalidPos;

        if (attempt == 0) {
          originalStale      = candidateStale;
          originalIneligible = ineligible;
        }

        const bool trivialOk = candidateIsLocal && candidateExists && !candidateStale && !candidateBusy &&
                               !candidateExposed && !queryExposed && !ineligible;

        if (trivialOk) {
          const NNMergeParticle<Packed>& cand = candDataIt->second;

          // Each participant's OWN level's dx -- see a_dxByLevel's docs. Using the wrong level's
          // dx here (e.g. a single shared dx) silently breaks this recheck for every level but
          // one, since the resulting cell key would not match the one a_liveCellCount was
          // actually populated under for that participant.
          const NNCellKey qKey(query.level, cellKeyOf(query.position, a_probLo, a_dxByLevel[query.level]));
          const NNCellKey cKey(cand.level, cellKeyOf(cand.position, a_probLo, a_dxByLevel[cand.level]));

          const auto qCountIt = a_liveCellCount.find(qKey);
          const auto cCountIt = a_liveCellCount.find(cKey);
          const int  qCount   = (qCountIt != a_liveCellCount.end()) ? qCountIt->second : 0;
          const int  cCount   = (cCountIt != a_liveCellCount.end()) ? cCountIt->second : 0;

          if (qCount > a_numParticlesPerCellThresh && cCount > a_numParticlesPerCellThresh) {
            // d. COMMIT.
            NNMergeResult<Packed> result;
            result.firstID  = query.globalID;
            result.secondID = cand.globalID;

            result.merged.position  = mergedPos;
            result.merged.weight    = mergedWeight;
            result.merged.globalID  = a_allocateID();
            result.merged.ownerRank = query.ownerRank;
            // Provisional -- placeMergedParticles() determines the TRUE destination level (and
            // patch/rank) from the merged position via findDestination(), never from this field.
            result.merged.level   = query.level;
            result.merged.payload = a_combine(query.payload, query.weight, cand.payload, cand.weight);

            a_results.push_back(result);

            a_consumedIDs.insert(query.globalID);
            a_consumedIDs.insert(cand.globalID);

            a_liveCellCount[qKey] -= 1;
            a_liveCellCount[cKey] -= 1;

            if (a_crossLevelMergeCount != nullptr && query.level != cand.level) {
              (*a_crossLevelMergeCount)++;
            }
            a_onMergeCommitted(query, cand, result.merged);

            resolved = true;
            break;
          }
          else {
            // e. Threshold recheck failed -- legitimately left unmatched, NOT fallback-eligible
            //    (a different candidate wouldn't change the query's own cell occupancy).
            resolved = true;
            break;
          }
        }
        // else: fall through to the next fallback attempt, if any remain.
      } // end attempt loop

      if (!resolved) {
        // f./g.: not resolved trivially -- see this function's docs for why only the ORIGINAL
        // candidate's staleness/ineligibility (not any fallback's) decides whether a proposal is
        // generated.
        if (!originalStale && !originalIneligible) {
          a_hasOutgoingCommitment.insert(query.globalID);
        }
      }
    }
  }

  template <typename Packed>
  inline void
  generateProposals(const std::vector<NNMergeEdge<Packed>>&              a_edges,
                    const std::set<ParticleID>&                          a_hasOutgoingCommitment,
                    const std::map<ParticleID, NNMergeParticle<Packed>>& a_particlesByID,
                    const int                                            a_numRanks,
                    std::vector<std::vector<NNMergeProposal<Packed>>>&   a_proposalsByDestRank) noexcept
  {
    CH_TIME("ParticleManagement::generateProposals");

    a_proposalsByDestRank.clear();
    a_proposalsByDestRank.resize(a_numRanks);

    for (const auto& edge : a_edges) {
      if (a_hasOutgoingCommitment.count(edge.queryID) == 0) {
        continue;
      }

      const auto sourceIt = a_particlesByID.find(edge.queryID);
      const auto targetIt = a_particlesByID.find(edge.candidateID);
      if (sourceIt == a_particlesByID.end() || targetIt == a_particlesByID.end()) {
        continue;
      }

      NNMergeProposal<Packed> proposal;
      proposal.source          = sourceIt->second;
      proposal.targetID        = edge.candidateID;
      proposal.distanceSquared = edge.distanceSquared;

      const int destRank = static_cast<int>(targetIt->second.ownerRank);
      CH_assert(destRank >= 0 && destRank < a_numRanks);
      a_proposalsByDestRank[destRank].push_back(proposal);
    }
  }

  template <typename Packed, typename Combine, typename PosValid, typename Allocator, typename OnMerge>
  inline void
  judgeProposals(const std::vector<NNMergeProposal<Packed>>&          a_incomingProposals,
                 const std::map<ParticleID, NNMergeParticle<Packed>>& a_myParticlesByID,
                 const std::set<ParticleID>&                          a_hasOutgoingCommitment,
                 const std::map<ParticleID, ParticleID>&              a_outgoingTargetOf,
                 const Combine&                                       a_combine,
                 std::unordered_map<NNCellKey, int, NNCellKeyHasher>& a_liveCellCount,
                 const int                                            a_numParticlesPerCellThresh,
                 const RealVect&                                      a_probLo,
                 const std::vector<RealVect>&                         a_dxByLevel,
                 const PosValid&                                      a_isPositionValid,
                 const Allocator&                                     a_allocateID,
                 std::vector<NNMergeResult<Packed>>&                  a_results,
                 std::vector<std::vector<NNMergeVerdict>>&            a_verdictsByDestRank,
                 const std::optional<int>                             a_maxCellDistance,
                 unsigned long long*                                  a_crossLevelMergeCount,
                 const OnMerge&                                       a_onMergeCommitted) noexcept
  {
    CH_TIME("ParticleManagement::judgeProposals");

    a_results.clear();
    // NOTE: a_verdictsByDestRank is APPENDED to, never resized here -- the caller must have
    // already sized it to the total rank count (same convention as placeMergedParticles's
    // a_scatterByDestRank).

    // 1. Group incoming proposals by target id.
    std::map<ParticleID, std::vector<const NNMergeProposal<Packed>*>> byTarget;
    for (const auto& p : a_incomingProposals) {
      byTarget[p.targetID].push_back(&p);
    }

    struct TargetDecision
    {
      ParticleID                     targetID;
      const NNMergeProposal<Packed>* winner;
    };

    std::vector<TargetDecision> decisions;
    decisions.reserve(byTarget.size());

    for (auto& kv : byTarget) {
      const ParticleID targetID  = kv.first;
      auto&            proposals = kv.second;

      const auto targetDataIt = a_myParticlesByID.find(targetID);
      if (targetDataIt == a_myParticlesByID.end()) {
        // Target no longer alive -- defensive skip, reject every incoming proposal.
        for (const NNMergeProposal<Packed>* p : proposals) {
          const int destRank = static_cast<int>(p->source.ownerRank);
          a_verdictsByDestRank[destRank].push_back(NNMergeVerdict{p->source.globalID, false});
        }
        continue;
      }

      const bool targetHasOutgoing = a_hasOutgoingCommitment.count(targetID) > 0;

      const NNMergeProposal<Packed>* winner = nullptr;

      if (!targetHasOutgoing) {
        // 2. Argmin over incoming proposals, deterministic id tiebreak.
        for (const NNMergeProposal<Packed>* p : proposals) {
          if (winner == nullptr || p->distanceSquared < winner->distanceSquared ||
              (p->distanceSquared == winner->distanceSquared && p->source.globalID < winner->source.globalID)) {
            winner = p;
          }
        }
      }
      else {
        // 3. Normally reject everything, UNLESS a mutual-match proposer is present.
        const auto ownTargetIt = a_outgoingTargetOf.find(targetID);
        if (ownTargetIt != a_outgoingTargetOf.end()) {
          const ParticleID ownTarget = ownTargetIt->second;
          for (const NNMergeProposal<Packed>* p : proposals) {
            if (p->source.globalID == ownTarget) {
              if (targetID < ownTarget) {
                winner = p;
              }
              break;
            }
          }
        }
      }

      if (winner != nullptr) {
        decisions.push_back(TargetDecision{targetID, winner});
      }

      // Reject every OTHER proposal to this target immediately -- a rejection is never
      // contingent on processing order, only an ACCEPT is (see rule 4).
      for (const NNMergeProposal<Packed>* p : proposals) {
        if (p != winner) {
          const int destRank = static_cast<int>(p->source.ownerRank);
          a_verdictsByDestRank[destRank].push_back(NNMergeVerdict{p->source.globalID, false});
        }
      }
    }

    // 4. Resolve all target decisions in ONE globally distance-sorted pass.
    std::sort(decisions.begin(), decisions.end(), [](const TargetDecision& a_a, const TargetDecision& a_b) {
      if (a_a.winner->distanceSquared != a_b.winner->distanceSquared) {
        return a_a.winner->distanceSquared < a_b.winner->distanceSquared;
      }
      return a_a.targetID < a_b.targetID;
    });

    for (const TargetDecision& decision : decisions) {
      const NNMergeProposal<Packed>& proposal = *decision.winner;
      const auto                     targetIt = a_myParticlesByID.find(decision.targetID);
      CH_assert(targetIt != a_myParticlesByID.end());
      const NNMergeParticle<Packed>& target = targetIt->second;
      const NNMergeParticle<Packed>& source = proposal.source;

      const int sourceDestRank = static_cast<int>(source.ownerRank);

      // 6. a_maxCellDistance / a_isPositionValid: a physical-validity rule, applies here exactly
      // as in resolveTrivialTier(), regardless of the mutual-match exception above.
      bool tooFar = false;
      if (a_maxCellDistance.has_value()) {
        tooFar = nnMergeCrossLevelTooFar(target.position,
                                         target.level,
                                         source.position,
                                         source.level,
                                         a_probLo,
                                         a_dxByLevel,
                                         a_maxCellDistance.value());
      }
      if (tooFar) {
        a_verdictsByDestRank[sourceDestRank].push_back(NNMergeVerdict{source.globalID, false});
        continue;
      }

      RealVect mergedPos(D_DECL(0.0, 0.0, 0.0));
      Real     mergedWeight = 0.0;
      nnMergeWeightedCentroid(target.position, target.weight, source.position, source.weight, mergedPos, mergedWeight);

      if (!a_isPositionValid(mergedPos)) {
        a_verdictsByDestRank[sourceDestRank].push_back(NNMergeVerdict{source.globalID, false});
        continue;
      }

      // Dynamic crowding recheck for BOTH participants (mirrors resolveTrivialTier()'s symmetric
      // "qCount > thresh && cCount > thresh" discipline exactly -- each participant's OWN level's
      // dx, see a_dxByLevel's docs). Checking only the target here was a real bug: a source
      // particle's own cell can drop below threshold between proposal generation and judging (an
      // earlier merge within this same round can drain it, same reasoning as resolveTrivialTier()'s
      // own recheck), and accepting anyway would silently violate the crowding trigger on the
      // source side. If this rank has no live visibility into the source's cell (never shipped
      // here as a ghost), sCount defaults to 0 -- rejected -- the same conservative default already
      // used when a target's own count is unknown.
      const NNCellKey tKey(target.level, cellKeyOf(target.position, a_probLo, a_dxByLevel[target.level]));
      const NNCellKey sKey(source.level, cellKeyOf(source.position, a_probLo, a_dxByLevel[source.level]));
      const auto      tCountIt = a_liveCellCount.find(tKey);
      const auto      sCountIt = a_liveCellCount.find(sKey);
      const int       tCount   = (tCountIt != a_liveCellCount.end()) ? tCountIt->second : 0;
      const int       sCount   = (sCountIt != a_liveCellCount.end()) ? sCountIt->second : 0;

      if (tCount <= a_numParticlesPerCellThresh || sCount <= a_numParticlesPerCellThresh) {
        a_verdictsByDestRank[sourceDestRank].push_back(NNMergeVerdict{source.globalID, false});
        continue;
      }

      // Accept.
      NNMergeResult<Packed> result;
      result.firstID  = target.globalID;
      result.secondID = source.globalID;

      result.merged.position  = mergedPos;
      result.merged.weight    = mergedWeight;
      result.merged.globalID  = a_allocateID();
      result.merged.ownerRank = target.ownerRank;
      // Provisional -- placeMergedParticles() determines the TRUE destination level (and
      // patch/rank) from the merged position via findDestination(), never from this field.
      result.merged.level   = target.level;
      result.merged.payload = a_combine(target.payload, target.weight, source.payload, source.weight);

      a_results.push_back(result);

      // Symmetric decrement (mirrors resolveTrivialTier()'s own commit) -- only meaningful when
      // this rank actually has a live entry for the source's cell (sCountIt found above); if it
      // doesn't, there is nothing here to decrement (the source's own rank's copy of its cell
      // count, if it has one, is this rank's problem to worry about only via the accepted verdict
      // it will receive, not via a_liveCellCount here).
      a_liveCellCount[tKey] -= 1;
      if (sCountIt != a_liveCellCount.end()) {
        a_liveCellCount[sKey] -= 1;
      }

      if (a_crossLevelMergeCount != nullptr && target.level != source.level) {
        (*a_crossLevelMergeCount)++;
      }
      a_onMergeCommitted(target, source, result.merged);

      a_verdictsByDestRank[sourceDestRank].push_back(NNMergeVerdict{source.globalID, true});
    }
  }

  inline void
  applyVerdicts(const std::vector<NNMergeVerdict>& a_incomingVerdicts, std::set<ParticleID>& a_consumedIDs) noexcept
  {
    CH_TIME("ParticleManagement::applyVerdicts");

    for (const auto& v : a_incomingVerdicts) {
      if (v.accepted) {
        a_consumedIDs.insert(v.proposerID);
      }
    }
  }

  template <typename Packed>
  inline void
  placeMergedParticles(const std::vector<NNMergeResult<Packed>>&                                     a_results,
                       const PositionLocator&                                                        a_locate,
                       const int                                                                     a_thisRank,
                       std::map<std::pair<int, unsigned int>, std::vector<NNMergeParticle<Packed>>>& a_localByPatch,
                       std::vector<std::vector<NNMergeParticle<Packed>>>& a_scatterByDestRank) noexcept
  {
    CH_TIME("ParticleManagement::placeMergedParticles");

    a_localByPatch.clear();
    // a_scatterByDestRank is APPENDED to, never resized here -- see its own docs.

    for (const auto& result : a_results) {
      // O(1) point->block query -- see PositionLocator/NNParticleLocation. NEVER a linear scan
      // over patches here; that is exactly the bottleneck this design avoids at large patch counts.
      const NNParticleLocation loc = a_locate(result.merged.position);

      // A merged particle landing outside every known patch (e.g. exactly on/outside the domain
      // boundary due to floating-point roundoff in the weighted centroid) is not expected in
      // normal operation; MayDay rather than silently dropping mass.
      if (!loc.valid) {
        MayDay::Error("ParticleManagement::placeMergedParticles -- merged particle position not "
                      "found in any patch");
      }

      if (loc.ownerRank == a_thisRank) {
        a_localByPatch[std::make_pair(loc.level, loc.gridIndex)].push_back(result.merged);
      }
      else {
        a_scatterByDestRank[loc.ownerRank].push_back(result.merged);
      }
    }
  }

  template <typename P,
            typename Packed,
            typename Traits,
            typename Gather,
            typename Combine,
            typename Scatter,
            typename Allocator,
            typename PosValid,
            typename OnMerge>
  inline void
  mergeNearestNeighborsRound(AmrMesh&                      a_amr,
                             ParticleContainer<P, Traits>& a_particles,
                             const int                     a_numParticlesPerCellThresh,
                             const Gather&                 a_gather,
                             const Combine&                a_combine,
                             const Scatter&                a_scatter,
                             const Allocator&              a_allocateID,
                             const bool                    a_iterateLocalTierToConvergence,
                             const int                     a_maxFallbackCandidates,
                             const std::optional<int>      a_maxCellDistance,
                             const PosValid&               a_isPositionValid,
                             unsigned long long*           a_crossLevelMergeCount,
                             const OnMerge&                a_onMergeCommitted,
                             const int                     a_subPartitionThresh)
  {
    CH_TIME("ParticleManagement::mergeNearestNeighborsRound");

    // TEMPORARY PERF INSTRUMENTATION -- do not commit. Per-phase wall-clock timing, requested to
    // locate the bottleneck behind mergeNearestNeighborsRound()'s cost relative to
    // ItoSolver::makeSuperparticles() (equal_weight_kd). Deliberately raw std::chrono, NOT
    // CH_TIME/CH_TIMER (Chombo's own TraceTimer profiler is reportedly slow enough itself to
    // distort a comparison like this) and not even this codebase's own Timer::wallClock()
    // wrapper -- std::chrono::steady_clock directly, per-rank, printed unconditionally (no
    // ParallelOps reduction -- this is a 1-rank investigation; would need one to be meaningful at
    // >1 rank).
    using NNPerfClock = std::chrono::steady_clock;
    auto nnPerfNow    = []() {
      return NNPerfClock::now();
    };
    auto nnPerfSeconds = [](const NNPerfClock::time_point& a_t0, const NNPerfClock::time_point& a_t1) -> double {
      return std::chrono::duration<double>(a_t1 - a_t0).count();
    };
    const NNPerfClock::time_point nnPerfRoundStart = nnPerfNow();

    const int myRank   = procID();
    const int numRanks = static_cast<int>(numProc());

    const std::string   realm       = a_particles.getRealm();
    const int           finestLevel = a_amr.getFinestLevel();
    const RealVect      probLo      = a_amr.getProbLo();
    const Vector<Real>& dxScalar    = a_amr.getDx();

    // 1. Ghost-fill -- fresh, exactly once, strictly before anything below. The search radius
    // this algorithm needs is one cell in every direction (the crowding trigger and the trivial-
    // tier candidate search only ever look at a particle's own cell and its immediate
    // neighborhood), so a width of 1 must already have been registered -- see this function's own
    // docs: registration is deliberately the CALLER's responsibility, not done here. A particle
    // ghost mask is only actually built during a regrid (AmrMesh::registerParticleGhostMask()
    // itself just records that a width should be tracked from then on); registering it here,
    // inside a function that runs long after grid setup, would be too late unless a regrid
    // happened to follow immediately afterward. The caller must register width 1 on this realm
    // BEFORE the relevant grids are (re)built (e.g. from a stepper's registerOperators()).
    constexpr int ghostWidth = 1;

    const NNPerfClock::time_point nnPerfGhostFillStart = nnPerfNow();
    a_particles.clearGhostParticles();
    a_particles.fillGhostParticles(a_amr.getParticleGhostMask(realm, ghostWidth),
                                   a_amr.getParticleGhostMaskCoarToFine(realm, ghostWidth),
                                   a_amr.getParticleGhostMaskFineToCoar(realm, ghostWidth));
    const double nnPerfGhostFillSeconds = nnPerfSeconds(nnPerfGhostFillStart, nnPerfNow());

    struct PatchWork
    {
      int       level;
      DataIndex din;
      RealVect  dx;
      int       blockSize;
      Box       haloBlockBox;
    };

    std::vector<PatchWork>                            patchWork;
    std::vector<std::vector<NNMergeParticle<Packed>>> patchLocalValid;   // parallel to patchWork
    std::vector<std::vector<NNMergeParticle<Packed>>> patchGhosts;       // parallel to patchWork
    std::vector<NNSpatialIndex<Packed>>               patchSpatialIndex; // parallel to patchWork

    std::map<ParticleID, NNMergeParticle<Packed>>       particlesByID; // pooled, this rank's own patches
    std::map<ParticleID, bool>                          exposed;       // pooled, this rank's own valid particles
    std::unordered_map<NNCellKey, int, NNCellKeyHasher> liveCellCount; // pooled, level-qualified cell key

    // Per-level dx, indexed by level -- every cell-key computation below and downstream (in
    // resolveTrivialTier()/judgeProposals()) MUST use the dx of the PARTICLE'S OWN level (see
    // NNMergeParticle::level's docs), never a single dx shared across levels.
    std::vector<RealVect> dxByLevel(finestLevel + 1);
    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      dxByLevel[lvl] = RealVect(D_DECL(dxScalar[lvl], dxScalar[lvl], dxScalar[lvl]));
    }

    const NNPerfClock::time_point nnPerfGatherStart = nnPerfNow();
    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      const RealVect& dx = dxByLevel[lvl];

      const DisjointBoxLayout& dbl = a_amr.getGrids(realm)[lvl];

      const AMRParticleGhostMask& maskSame = a_amr.getParticleGhostMask(realm, ghostWidth);
      const AMRParticleGhostMask& maskC2F  = a_amr.getParticleGhostMaskCoarToFine(realm, ghostWidth);
      const AMRParticleGhostMask& maskF2C  = a_amr.getParticleGhostMaskFineToCoar(realm, ghostWidth);

      const RefCountedPtr<MergeBlockGrid>& mergeBlockGrid = a_amr.getMergeBlockGrid(realm)[lvl];

      const DataIterator& dit  = dbl.dataIterator();
      const int           nbox = dit.size();

      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

        std::vector<NNMergeParticle<Packed>> localValid;
        std::vector<NNMergeParticle<Packed>> ghosts;

        for (std::size_t i = 0; i < leaf.size(); i++) {
          const bool isGhost = leaf.isGhost(i);

          // A ghost's TRUE origin level is NOT this receiving patch's level `lvl` -- it is
          // lvl-1/lvl+1 for a Coarse/Fine ghost (see ParticleSoA::GhostType). Getting this wrong
          // computes the ghost's cell key with the WRONG dx, exactly the bug that motivated
          // NNMergeParticle::level existing at all.
          int originLevel = lvl;
          if (isGhost) {
            switch (leaf.ghost(i)) {
            case GhostType::Coarse:
              originLevel = lvl - 1;
              break;
            case GhostType::Fine:
              originLevel = lvl + 1;
              break;
            default:
              originLevel = lvl; // GhostType::SameLevel
              break;
            }
          }

          NNMergeParticle<Packed> p;
          p.position  = leaf.position(i);
          p.weight    = leaf.weight(i);
          p.globalID  = leaf.particleID(i);
          p.ownerRank = leaf.rankID(i);
          p.level     = originLevel;
          p.payload   = a_gather(leaf, i);

          const IntVect cell = cellKeyOf(p.position, probLo, dxByLevel[originLevel]);

          if (isGhost) {
            ghosts.push_back(p);
            // REQUIRED -- see generateProposals()'s own documented contract for a_particlesByID:
            // "every currently-alive local valid particle and every ghost". Without this, a
            // genuinely cross-rank ghost candidate (one whose true owning patch is NOT among this
            // rank's own patches, so it's never seen as "native" anywhere in this loop) has no
            // entry here at all, generateProposals()'s lookup silently fails, and the whole
            // proposal is dropped -- confirmed empirically: this masks itself as "correct" at 1
            // rank purely because every patch (hence every ghost's true owner too) happens to be
            // this same rank, so the ghost's id is ALSO reachable via its owning patch's OWN
            // native entry processed elsewhere in this same loop; at higher rank counts, a ghost
            // whose owner is a genuinely different rank has no such fallback and the proposal is
            // silently lost.
            particlesByID[p.globalID] = p;
          }
          else {
            localValid.push_back(p);
            particlesByID[p.globalID] = p;

            bool isExposed = !maskSame[lvl].isNull() && isBoundaryExposed(cell, (*maskSame[lvl])[din]);
            if (!isExposed && lvl < finestLevel && !maskC2F[lvl].isNull()) {
              isExposed = isBoundaryExposed(cell, (*maskC2F[lvl])[din]);
            }
            if (!isExposed && lvl > 0 && !maskF2C[lvl].isNull()) {
              isExposed = isBoundaryExposed(cell, (*maskF2C[lvl])[din]);
            }
            exposed[p.globalID] = isExposed;
          }
        }

        // Built once here, per patch, per round -- over the FULL localValid/ghosts sets, before
        // any consumption happens this round -- and reused across every internal
        // a_iterateLocalTierToConvergence pass via lazy deletion (see
        // findNearestNeighborCandidates()'s own docs for why rebuilding per pass would reintroduce
        // real, avoidable cost on the up-to-1000x-per-round loop this redesign exists to speed up).
        const MergeBlockGrid::PatchBlocks& patchBlocks = mergeBlockGrid->getPatchBlocks(din);

        NNSpatialIndex<Packed> spatialIndex;
        buildNNSpatialIndex(localValid,
                            ghosts,
                            probLo,
                            dx,
                            mergeBlockGrid->getBlockSize(),
                            patchBlocks.haloBlockBox,
                            a_subPartitionThresh,
                            a_numParticlesPerCellThresh,
                            static_cast<long>(dbl[din].numPts()),
                            spatialIndex);

        patchWork.push_back(PatchWork{lvl, din, dx, mergeBlockGrid->getBlockSize(), patchBlocks.haloBlockBox});
        patchLocalValid.push_back(std::move(localValid));
        patchGhosts.push_back(std::move(ghosts));
        patchSpatialIndex.push_back(std::move(spatialIndex));
      }
    }

    // Populate liveCellCount from particlesByID -- deliberately NOT incremented inline during the
    // per-particle gather loop above. A single physical particle can appear as a ghost in MULTIPLE
    // of this rank's own patches at once (e.g. any particle near a corner/edge shared by several
    // patches this rank happens to own -- routine with tiled grid generation, which favors keeping
    // neighboring patches on the same rank), so incrementing a counter once per occurrence there
    // would double- (or triple-, quadruple-...) count that one particle's contribution to its own
    // cell's live count. particlesByID has no such problem: it is a map keyed by globalID, so
    // repeated occurrences of the same particle (native once, ghost N times) collapse to exactly
    // one entry -- deriving the count from it here, in one pass after the map is fully built,
    // counts every distinct particle exactly once regardless of how many patches it was seen in.
    for (const auto& kv : particlesByID) {
      const NNMergeParticle<Packed>& p    = kv.second;
      const IntVect                  cell = cellKeyOf(p.position, probLo, dxByLevel[p.level]);
      liveCellCount[NNCellKey(p.level, cell)] += 1;
    }
    const double nnPerfGatherSeconds = nnPerfSeconds(nnPerfGatherStart, nnPerfNow());

    std::set<ParticleID>             consumedIDs;
    std::set<ParticleID>             hasOutgoingCommitment;
    std::map<ParticleID, ParticleID> outgoingTargetOf;

    std::vector<NNMergeResult<Packed>> allResults;

    // 2. Per-patch candidate search + trivial-tier resolution, pooled rank-wide.
    //
    // lastPooledEdges deliberately survives past the loop below: it is reused (NOT
    // recomputed) for outgoingTargetOf/generateProposals afterward. findNearestNeighborCandidates()
    // is a brute-force O(n^2) scan per patch (a real spatial index is a separate, still-pending
    // hardening item -- see TODO.md), so re-running it a second time over the same alive set
    // purely to re-derive edges already computed here would double the single most expensive step
    // in this function for no new information: nothing consumes any more particles between the
    // last runLocalTierPass() and generateProposals(), so the edges are still exactly current.
    std::vector<NNMergeEdge<Packed>> lastPooledEdges;

    auto runLocalTierPass = [&]() -> bool {
      lastPooledEdges.clear();

      for (std::size_t pw = 0; pw < patchWork.size(); pw++) {
        // patchLocalValid[pw]/patchGhosts[pw] are this patch's FULL sets (not pre-filtered) and
        // patchSpatialIndex[pw] was built once, over those same full sets, before this round's
        // local-tier passes began -- findNearestNeighborCandidates() does its own lazy-deletion
        // filtering against consumedIDs internally now (see its own docs), so nothing needs to be
        // rebuilt or re-filtered here on every pass.
        std::vector<NNMergeEdge<Packed>> patchEdges;
        findNearestNeighborCandidates(patchLocalValid[pw],
                                      patchGhosts[pw],
                                      consumedIDs,
                                      patchSpatialIndex[pw],
                                      liveCellCount,
                                      probLo,
                                      patchWork[pw].dx,
                                      a_numParticlesPerCellThresh,
                                      patchEdges,
                                      a_maxFallbackCandidates,
                                      a_maxCellDistance);

        lastPooledEdges.insert(lastPooledEdges.end(), patchEdges.begin(), patchEdges.end());
      }

      std::vector<NNMergeResult<Packed>> passResults;
      resolveTrivialTier(lastPooledEdges,
                         exposed,
                         particlesByID,
                         a_combine,
                         liveCellCount,
                         consumedIDs,
                         a_numParticlesPerCellThresh,
                         hasOutgoingCommitment,
                         passResults,
                         probLo,
                         dxByLevel,
                         a_isPositionValid,
                         a_allocateID,
                         a_maxFallbackCandidates,
                         a_maxCellDistance,
                         a_crossLevelMergeCount,
                         a_onMergeCommitted);

      allResults.insert(allResults.end(), passResults.begin(), passResults.end());

      return !passResults.empty();
    };

    const NNPerfClock::time_point nnPerfLocalTierStart = nnPerfNow();
    if (a_iterateLocalTierToConvergence) {
      constexpr int safetyCap = 1000;
      for (int iter = 0; iter < safetyCap; iter++) {
        // hasOutgoingCommitment must NOT accumulate across internal iterations -- each iteration
        // re-derives it from scratch via a fresh candidate search over the currently-alive set.
        hasOutgoingCommitment.clear();
        if (!runLocalTierPass()) {
          break;
        }
      }
    }
    else {
      runLocalTierPass();
    }
    const double nnPerfLocalTierSeconds = nnPerfSeconds(nnPerfLocalTierStart, nnPerfNow());

    // outgoingTargetOf: for every query id with an outgoing commitment, remember its own true
    // nearest-neighbor candidate id (needed by judgeProposals()'s mutual-match rule). Reuses
    // lastPooledEdges from the local-tier pass above -- see the comment there for why this is
    // correct without a fresh search: resolveTrivialTier() only ever adds a queryID to
    // hasOutgoingCommitment when its edge's candidate was valid (not stale) at the moment that
    // edge was processed, and nothing gets consumed between then and here.
    for (const auto& edge : lastPooledEdges) {
      if (hasOutgoingCommitment.count(edge.queryID) > 0) {
        outgoingTargetOf[edge.queryID] = edge.candidateID;
      }
    }

    // 3. generateProposals() for everything the trivial tier could not resolve.
    const NNPerfClock::time_point                     nnPerfGenProposalsStart = nnPerfNow();
    std::vector<std::vector<NNMergeProposal<Packed>>> proposalsByDestRank;
    generateProposals(lastPooledEdges, hasOutgoingCommitment, particlesByID, numRanks, proposalsByDestRank);
    const double nnPerfGenProposalsSeconds = nnPerfSeconds(nnPerfGenProposalsStart, nnPerfNow());

    // 4. Exchange proposals.
    const NNPerfClock::time_point              nnPerfExchProposalsStart = nnPerfNow();
    const std::vector<NNMergeProposal<Packed>> incomingProposals        = nnMergeExchangeByRank(proposalsByDestRank);
    const double nnPerfExchProposalsSeconds = nnPerfSeconds(nnPerfExchProposalsStart, nnPerfNow());

    // 5. judgeProposals() over the pooled incoming proposals.
    const NNPerfClock::time_point                 nnPerfMyParticlesCopyStart = nnPerfNow();
    std::map<ParticleID, NNMergeParticle<Packed>> myParticlesByID;
    for (const auto& kv : particlesByID) {
      if (consumedIDs.count(kv.first) == 0) {
        myParticlesByID[kv.first] = kv.second;
      }
    }
    const double nnPerfMyParticlesCopySeconds = nnPerfSeconds(nnPerfMyParticlesCopyStart, nnPerfNow());

    std::vector<NNMergeResult<Packed>>       judgeResults;
    std::vector<std::vector<NNMergeVerdict>> verdictsByDestRank(numRanks);

    const NNPerfClock::time_point nnPerfJudgeStart = nnPerfNow();
    judgeProposals(incomingProposals,
                   myParticlesByID,
                   hasOutgoingCommitment,
                   outgoingTargetOf,
                   a_combine,
                   liveCellCount,
                   a_numParticlesPerCellThresh,
                   probLo,
                   dxByLevel,
                   a_isPositionValid,
                   a_allocateID,
                   judgeResults,
                   verdictsByDestRank,
                   a_maxCellDistance,
                   a_crossLevelMergeCount,
                   a_onMergeCommitted);
    const double nnPerfJudgeSeconds = nnPerfSeconds(nnPerfJudgeStart, nnPerfNow());

    allResults.insert(allResults.end(), judgeResults.begin(), judgeResults.end());
    for (const auto& r : judgeResults) {
      consumedIDs.insert(r.firstID);
    }

    // 6. Exchange verdicts and applyVerdicts().
    const NNPerfClock::time_point     nnPerfVerdictsStart = nnPerfNow();
    const std::vector<NNMergeVerdict> incomingVerdicts    = nnMergeExchangeByRank(verdictsByDestRank);
    applyVerdicts(incomingVerdicts, consumedIDs);
    const double nnPerfVerdictsSeconds = nnPerfSeconds(nnPerfVerdictsStart, nnPerfNow());

    // 7. placeMergedParticles(), then exchange and apply the scatter bucket.
    //
    // O(1) point->block query -- see NNParticleLocation/PositionLocator's docs. This is precisely
    // what avoids the "linear scan over every patch" bottleneck: ParticleContainer::findDestination()
    // is the SAME tile-hash lookup a_particles already uses internally to route its own particles
    // during remap() -- reusing it here (rather than adding a parallel AmrMesh/Realm-level
    // forwarding method) keeps this one mapping the single source of truth for "which box owns
    // this position", instead of a second implementation that could drift from it.
    const PositionLocator locate = [&a_particles](const RealVect& a_pos) -> NNParticleLocation {
      const auto dst = a_particles.findDestination(a_pos);
      return NNParticleLocation{dst.found, dst.level, dst.gridIndex, dst.rank};
    };

    const NNPerfClock::time_point                                                nnPerfPlaceStart = nnPerfNow();
    std::map<std::pair<int, unsigned int>, std::vector<NNMergeParticle<Packed>>> localByPatch;
    std::vector<std::vector<NNMergeParticle<Packed>>>                            scatterByDestRank(numRanks);

    placeMergedParticles(allResults, locate, myRank, localByPatch, scatterByDestRank);

    const std::vector<NNMergeParticle<Packed>> incomingScattered  = nnMergeExchangeByRank(scatterByDestRank);
    const double                               nnPerfPlaceSeconds = nnPerfSeconds(nnPerfPlaceStart, nnPerfNow());

    const NNPerfClock::time_point nnPerfApplyStart = nnPerfNow();
    // ---- Apply everything to the real ParticleSoA containers ----
    // Remove every consumed id (trivial-tier participants and accepted judge-side participants --
    // applyVerdicts() already folded accepted proposers into consumedIDs above).
    for (std::size_t pw = 0; pw < patchWork.size(); pw++) {
      ParticleSoA<P, Traits>& leaf = a_particles[patchWork[pw].level][patchWork[pw].din];

      for (std::size_t i = 0; i < leaf.size();) {
        if (!leaf.isGhost(i) && consumedIDs.count(leaf.particleID(i)) > 0) {
          leaf.remove(i);
        }
        else {
          i++;
        }
      }
    }
    a_particles.clearGhostParticles();

    // Scatter in every locally-placed merged particle. (level, gridIndex) -> DataIndex is a single
    // O(1) hash lookup via LevelTiles::getMyGrids() (global grid index -> local DataIndex, for
    // this rank's own boxes only), done once per bucket here -- not per particle, and never a scan.
    for (const auto& kv : localByPatch) {
      const int          level     = kv.first.first;
      const unsigned int gridIndex = kv.first.second;

      const DataIndex din = a_amr.getLevelTiles(realm)[level]->getMyGrids().at(gridIndex);

      ParticleSoA<P, Traits>& leaf = a_particles[level][din];
      for (const auto& merged : kv.second) {
        // Correct the provisional ownerRank (see resolveTrivialTier()/judgeProposals()'s own
        // "Provisional" comments -- it's inherited from whichever parent committed the merge, NOT
        // the true final destination) to the TRUE rank now that placeMergedParticles() has
        // resolved it via findDestination(): this IS that rank. Leaving the stale, provisional
        // value here is a real bug -- confirmed empirically as the root cause of a mass
        // conservation failure that only appeared with BOTH more than one MPI rank AND more than
        // one round: if this particle later becomes a proposer in a subsequent round and its
        // proposal is accepted, judgeProposals() routes the acceptance verdict to
        // NNMergeParticle::ownerRank, which -- left uncorrected -- could be a rank other than the
        // one that actually holds this particle. That true owner then never learns to remove it,
        // so it survives un-consumed while its weight has ALSO been folded into a new merged
        // particle elsewhere: mass is silently created.
        NNMergeParticle<Packed> corrected = merged;
        corrected.ownerRank               = myRank;
        a_scatter(leaf, corrected);
      }
    }

    // Ingest every merged particle that was shipped in from another rank -- same O(1) locate, no
    // scan over this rank's patches.
    for (const auto& merged : incomingScattered) {
      const NNParticleLocation loc = locate(merged.position);

      if (!loc.valid || loc.ownerRank != myRank) {
        // placeMergedParticles() already routed this to this rank via loc.ownerRank; a result
        // that no longer locates to this rank (or locates nowhere) here indicates an inconsistency
        // between the two locate() calls, which must never happen since the grid does not change
        // mid-round -- MayDay rather than silently dropping mass.
        MayDay::Error("ParticleManagement::mergeNearestNeighborsRound -- incoming scattered "
                      "particle not found in any of this rank's own patches");
      }

      const DataIndex din = a_amr.getLevelTiles(realm)[loc.level]->getMyGrids().at(loc.gridIndex);

      // Same ownerRank correction as the locally-placed loop above -- this rank is the TRUE owner
      // (already checked via loc.ownerRank == myRank above), not whatever provisional value the
      // committing rank happened to stamp on it.
      NNMergeParticle<Packed> corrected = merged;
      corrected.ownerRank               = myRank;

      ParticleSoA<P, Traits>& leaf = a_particles[loc.level][din];
      a_scatter(leaf, corrected);
    }
    const double nnPerfApplySeconds = nnPerfSeconds(nnPerfApplyStart, nnPerfNow());
    const double nnPerfTotalSeconds = nnPerfSeconds(nnPerfRoundStart, nnPerfNow());

    pout() << "NNPERF rank=" << myRank << " ghostFill=" << nnPerfGhostFillSeconds
           << " gatherAndIndex=" << nnPerfGatherSeconds << " localTier=" << nnPerfLocalTierSeconds
           << " generateProposals=" << nnPerfGenProposalsSeconds << " exchangeProposals=" << nnPerfExchProposalsSeconds
           << " myParticlesCopy=" << nnPerfMyParticlesCopySeconds << " judgeProposals=" << nnPerfJudgeSeconds
           << " exchangeVerdicts=" << nnPerfVerdictsSeconds << " placeAndExchange=" << nnPerfPlaceSeconds
           << " applyToContainers=" << nnPerfApplySeconds << " TOTAL=" << nnPerfTotalSeconds << endl;

    // TEMPORARY DIAGNOSTIC -- do not commit.
    pout() << "NNTREEDIAG queries=" << g_nnDiagQueries << " nodeVisits=" << g_nnDiagNodeVisits
           << " leafCandidateScans=" << g_nnDiagLeafCandidates << " leavesVisited=" << g_nnDiagLeavesVisited
           << " avgNodeVisitsPerQuery="
           << (g_nnDiagQueries > 0 ? static_cast<double>(g_nnDiagNodeVisits) / static_cast<double>(g_nnDiagQueries)
                                   : 0.0)
           << " avgLeafScansPerQuery="
           << (g_nnDiagQueries > 0 ? static_cast<double>(g_nnDiagLeafCandidates) / static_cast<double>(g_nnDiagQueries)
                                   : 0.0)
           << " avgLeavesVisitedPerQuery="
           << (g_nnDiagQueries > 0 ? static_cast<double>(g_nnDiagLeavesVisited) / static_cast<double>(g_nnDiagQueries)
                                   : 0.0)
           << " avgLeafSizeInTree="
           << (g_nnDiagLeafSizeCount > 0
                 ? static_cast<double>(g_nnDiagLeafSizeSum) / static_cast<double>(g_nnDiagLeafSizeCount)
                 : 0.0)
           << " maxLeafSizeInTree=" << g_nnDiagMaxLeafSize << " avgTreeSize="
           << (g_nnDiagTreeCount > 0 ? static_cast<double>(g_nnDiagTreeSizeSum) / static_cast<double>(g_nnDiagTreeCount)
                                     : 0.0)
           << " numTrees=" << g_nnDiagTreeCount << endl;
    g_nnDiagQueries        = 0;
    g_nnDiagNodeVisits     = 0;
    g_nnDiagLeafCandidates = 0;
    g_nnDiagLeavesVisited  = 0;
    g_nnDiagLeafSizeSum    = 0;
    g_nnDiagLeafSizeCount  = 0;
    g_nnDiagMaxLeafSize    = 0;
    g_nnDiagTreeSizeSum    = 0;
    g_nnDiagTreeCount      = 0;
  }

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#endif
