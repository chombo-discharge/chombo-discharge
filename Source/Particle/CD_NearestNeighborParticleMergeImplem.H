/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_NearestNeighborParticleMergeImplem.H
 * @brief  Implementation of CD_NearestNeighborParticleMerge.H
 * @author Robert Marskar
 */

#ifndef CD_NEARESTNEIGHBORPARTICLEMERGEIMPLEM_H
#define CD_NEARESTNEIGHBORPARTICLEMERGEIMPLEM_H

// Std includes
#include <algorithm>
#include <cmath>

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

namespace detail {

/**
 * @brief Squared Euclidean distance between two physical positions.
 * @param[in] a_a First position.
 * @param[in] a_b Second position.
 * @return Squared distance between a_a and a_b.
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
 * @brief Chebyshev (max-per-axis) distance between two cell-index vectors -- see
 * CD_NearestNeighborParticleMerge.H's a_maxCellDistance docs.
 * @param[in] a_a First cell index.
 * @param[in] a_b Second cell index.
 * @return Max-per-axis distance between a_a and a_b.
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
 * @brief a_maxCellDistance check for two positions that may be on DIFFERENT AMR levels.
 * @details Raw cell indices computed with two different dx values are not directly comparable, so
 * both positions are re-expressed in the finer of the two levels' cell-size units (higher AMR
 * level = finer) before taking the Chebyshev distance. Reduces to the same-level case when
 * a_levelA == a_levelB.
 * @param[in] a_posA         First participant's position.
 * @param[in] a_levelA       First participant's AMR level.
 * @param[in] a_posB         Second participant's position.
 * @param[in] a_levelB       Second participant's AMR level.
 * @param[in] a_probLo       Physical position of the domain's low corner.
 * @param[in] a_dxByLevel    Cell spacing per AMR level.
 * @param[in] a_maxCellDistance The cap (see CD_NearestNeighborParticleMerge.H's a_maxCellDistance
 * docs).
 * @return True iff the pair exceeds a_maxCellDistance.
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
 * @brief Weighted-centroid merge of two positions/weights -- shared by resolveTrivialTier() and
 * judgeProposals() so the formula is written exactly once.
 * @param[out] a_mergedPos   Weighted-centroid position.
 * @param[out] a_mergedWeight Summed weight.
 * @param[in]  a_posA        First particle's position.
 * @param[in]  a_weightA     First particle's weight.
 * @param[in]  a_posB        Second particle's position.
 * @param[in]  a_weightB     Second particle's weight.
 */
inline void
nnMergeWeightedCentroid(RealVect&       a_mergedPos,
                        Real&           a_mergedWeight,
                        const RealVect& a_posA,
                        const Real      a_weightA,
                        const RealVect& a_posB,
                        const Real      a_weightB) noexcept
{
  a_mergedWeight = a_weightA + a_weightB;

  for (int dir = 0; dir < SpaceDim; dir++) {
    a_mergedPos[dir] = (a_posA[dir] * a_weightA + a_posB[dir] * a_weightB) / a_mergedWeight;
  }
}

/**
 * @brief Squared distance from a physical point to a cell's axis-aligned bounding box (0 if the
 * point is inside).
 * @param[in] a_point  Physical position.
 * @param[in] a_cell   Cell index.
 * @param[in] a_probLo Physical position of the domain's low corner.
 * @param[in] a_dx     Cell spacing.
 * @return Squared distance from a_point to [a_probLo + a_cell*a_dx, a_probLo + (a_cell + 1)*a_dx].
 */
inline Real
nnMergePointToCellDistSquared(const RealVect& a_point,
                              const IntVect&  a_cell,
                              const RealVect& a_probLo,
                              const RealVect& a_dx) noexcept
{
  Real d2 = 0.0;

  for (int dir = 0; dir < SpaceDim; dir++) {
    const Real lo = a_probLo[dir] + a_cell[dir] * a_dx[dir];
    const Real hi = lo + a_dx[dir];

    Real delta = 0.0;

    if (a_point[dir] < lo) {
      delta = lo - a_point[dir];
    }
    else if (a_point[dir] > hi) {
      delta = a_point[dir] - hi;
    }

    d2 += delta * delta;
  }

  return d2;
}

/**
 * @brief Generic alltoallv-style exchange: send a per-destination-rank bucket of trivially-
 * copyable records, return the flat concatenation of everything received.
 * @details Packs one bucket per rank and exchanges via MPI_Alltoall (byte counts) + MPI_Alltoallv
 * (payload). In a serial (non-MPI) build, or with a single rank, this returns this rank's own
 * bucket unchanged.
 * @tparam T Trivially-copyable record type (NNMergeProposal<Packed>, NNMergeVerdict, or
 * MergeParticle<Packed> -- all satisfy this given Packed itself is required to).
 * @param[in] a_sendByRank One bucket per destination rank (size == number of ranks).
 * @return Flat vector of every record received from every rank (including this rank's own
 * bucket addressed to itself).
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

inline IntVect
cellKeyOf(const RealVect& a_position, const RealVect& a_probLo, const RealVect& a_dx) noexcept
{
  IntVect iv;

  for (int dir = 0; dir < SpaceDim; dir++) {
    iv[dir] = static_cast<int>(std::floor((a_position[dir] - a_probLo[dir]) / a_dx[dir]));
  }

  return iv;
}

/**
 * @brief Exact per-cell crowding pre-filter shared by buildNNSpatialIndex()/buildNNCellSpatialIndex(): is any
 * cell in this patch's local+ghost particle set over its a_cellBudget?
 * @details A count-only tally keyed as the crowding trigger (cellKeyOf(pos, probLo, dx), local and ghost
 * alike); checks the true per-cell maximum, not an average density (a low-average patch can still hold one
 * over-threshold cell).
 * @tparam Packed See MergeParticle.
 * @param[in] a_localValid This patch's own full valid-particle set.
 * @param[in] a_ghosts Every ghost shipped to this patch this round.
 * @param[in] a_probLo Physical position of the domain's low corner.
 * @param[in] a_dx Cell spacing of THIS PATCH's own level.
 * @param[in] a_cellBudget The count each cell is drained down to -- see NNCellBudget.
 * @return True iff some cell holds strictly more particles than its own budget allows.
 */
template <typename Packed>
inline bool
nnMergeAnyCellCrowded(const std::vector<MergeParticle<Packed>>& a_localValid,
                      const std::vector<MergeParticle<Packed>>& a_ghosts,
                      const RealVect&                           a_probLo,
                      const RealVect&                           a_dx,
                      const NNCellBudget&                       a_cellBudget) noexcept
{
  const std::size_t total = a_localValid.size() + a_ghosts.size();

  std::unordered_map<NNWalkCell, int, LevelTiles::TileHasher> cellCount;
  cellCount.reserve(total);

  for (const auto& p : a_localValid) {
    cellCount[cellKeyOf(p.position, a_probLo, a_dx)]++;
  }

  for (const auto& p : a_ghosts) {
    cellCount[cellKeyOf(p.position, a_probLo, a_dx)]++;
  }

  // Compared against the SMALLEST budget any cell has rather than this cell's own. The count above
  // is pooled per bare cell, in the query patch's single frame, so it has no level to look a budget
  // up with -- and this is only a pre-filter guarding whether the patch gets a spatial index at all.
  // Taking the most permissive budget keeps it strictly conservative: it can never skip a patch the
  // real per-cell tests downstream would have found work in. For a uniform budget this is exactly
  // the plain threshold comparison it replaced.
  const int minBudget = a_cellBudget.minimum();

  for (const auto& kv : cellCount) {
    if (kv.second > minBudget) {
      return true;
    }
  }

  return false;
}

template <typename Packed, typename Cloud>
inline void
buildNNSpatialIndex(NNSpatialIndex<Packed, Cloud>&            a_index,
                    const std::vector<MergeParticle<Packed>>& a_localValid,
                    const std::vector<MergeParticle<Packed>>& a_ghosts,
                    const RealVect&                           a_probLo,
                    const RealVect&                           a_dx,
                    const NNCellBudget&                       a_cellBudget,
                    const int                                 a_kNearest) noexcept
{
  CH_TIME("ParticleManagement::buildNNSpatialIndex");

  a_index.pointCloud = nullptr;
  a_index.pointCloudGraph.clear();
  a_index.pointCloudNLocal = 0;
  a_index.pointCloudK      = a_kNearest;

  const std::size_t total = a_localValid.size() + a_ghosts.size();
  if (total == 0) {
    return;
  }

  // Exact per-cell pre-filter (see this function's own docs): build nothing unless some cell is
  // over threshold.
  if (!nnMergeAnyCellCrowded(a_localValid, a_ghosts, a_probLo, a_dx, a_cellBudget)) {
    return;
  }

  // Assemble the whole-patch point cloud (local valid first, then ghosts). Metadata is each
  // point's globalID; a hit's cloud index decodes ownership (< nLocal => local valid, else
  // ghost). Positions are embedded into 3D (z = 0 in 2D) for EBGeometry's Vec3T.
  const std::size_t nLocal = a_localValid.size();

  std::vector<EBGeometry::Vec3T<Real>> positions;
  std::vector<ParticleID>              metadata;

  positions.reserve(total);
  metadata.reserve(total);

  auto pushPoint = [&](const MergeParticle<Packed>& a_p) {
    EBGeometry::Vec3T<Real> v(0.0, 0.0, 0.0);

    for (int dir = 0; dir < SpaceDim; dir++) {
      v[dir] = a_p.position[dir];
    }

    positions.push_back(v);
    metadata.push_back(a_p.globalID);
  };

  for (const auto& p : a_localValid) {
    pushPoint(p);
  }

  for (const auto& p : a_ghosts) {
    pushPoint(p);
  }

  a_index.pointCloud       = std::make_shared<Cloud>(positions, metadata);
  a_index.pointCloudNLocal = nLocal;

  // Precompute the whole kNN graph in one batched pass, reused across every internal
  // a_iterateLocalTierToConvergence pass (positions never change within a round; only consumed
  // status does, and that is filtered at lookup time). Ghost rows are computed but never read;
  // see NNSpatialIndex::pointCloudGraph.
  a_index.pointCloudGraph = a_index.pointCloud->allNearestNeighbors(static_cast<std::size_t>(a_kNearest));
}

template <typename Packed, typename Cloud>
inline void
findNearestNeighborCandidates(std::vector<NNMergeEdge<Packed>>&                          a_edges,
                              const std::vector<MergeParticle<Packed>>&                  a_localValid,
                              const std::vector<MergeParticle<Packed>>&                  a_ghosts,
                              const std::unordered_set<ParticleID>&                      a_consumedIDs,
                              const NNSpatialIndex<Packed, Cloud>&                       a_spatialIndex,
                              const std::unordered_map<NNCellKey, int, NNCellKeyHasher>& a_liveCellCount,
                              const RealVect&                                            a_probLo,
                              const RealVect&                                            a_dx,
                              const NNCellBudget&                                        a_cellBudget,
                              const int                                                  a_maxFallbackCandidates,
                              const std::optional<int>                                   a_maxCellDistance) noexcept
{
  CH_TIME("ParticleManagement::findNearestNeighborCandidates");

  a_edges.clear();

  const std::size_t nLocal = a_localValid.size();
  const int         K      = 1 + a_maxFallbackCandidates;

  // a_maxCellDistance's physical bound, computed once per call. The sqrt(SpaceDim) factor keeps
  // this Euclidean prune no tighter than the per-axis Chebyshev cutoff applied downstream by
  // resolveTrivialTier()/judgeProposals(), so it never excludes a candidate they would accept
  // (see this function's own a_maxCellDistance docs).
  std::optional<Real> maxDistBoundSquared;
  if (a_maxCellDistance.has_value()) {
    Real minDx = a_dx[0];

    for (int dir = 1; dir < SpaceDim; dir++) {
      minDx = std::min(minDx, a_dx[dir]);
    }

    const Real bound    = static_cast<Real>(a_maxCellDistance.value()) * minDx * std::sqrt(static_cast<Real>(SpaceDim));
    maxDistBoundSquared = bound * bound;
  }

  struct RankedCandidate
  {
    Real        d2;
    std::size_t idx;
    bool        isLocal;
  };

  for (std::size_t qi = 0; qi < nLocal; qi++) {
    const MergeParticle<Packed>& q = a_localValid[qi];

    // Lazy deletion -- a_localValid is the FULL per-patch set (see this function's own docs), so
    // a query already consumed by an earlier pass/query this round is skipped explicitly here.
    if (a_consumedIDs.count(q.globalID) > 0) {
      continue;
    }

    // a_dx is this patch's own level's dx -- correct here since every query in a_localValid is
    // necessarily on this same level (see this function's a_dx docs).
    const IntVect   qCell = cellKeyOf(q.position, a_probLo, a_dx);
    const NNCellKey qKey{q.level, qCell};
    const auto      cntIt = a_liveCellCount.find(qKey);
    const int       count = (cntIt != a_liveCellCount.end()) ? cntIt->second : 0;

    if (count <= a_cellBudget(qKey)) {
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

    // Read this query's precomputed kNN row out of the point-cloud graph (see
    // buildNNSpatialIndex()). The query's own cloud index equals qi (local valid points packed
    // first); allNearestNeighbors() already excluded it as self. Each hit's cloud index decodes
    // ownership (< pointCloudNLocal => local valid, else ghost); apply the a_maxCellDistance bound
    // and the a_consumedIDs skip before feeding survivors into the shared top-K logic. Filtering
    // the distance-sorted k-nearest by the bound yields exactly the within-bound k-nearest.
    // pointCloud is null only for a patch with no crowded cell, which generates no query.
    if (a_spatialIndex.pointCloud) {
      const std::size_t nLocalPC = a_spatialIndex.pointCloudNLocal;
      const int         Kgraph   = a_spatialIndex.pointCloudK;
      const std::size_t rowStart = qi * static_cast<std::size_t>(Kgraph);

      for (int h = 0; h < Kgraph; h++) {
        const typename Cloud::Hit& hit = a_spatialIndex.pointCloudGraph[rowStart + static_cast<std::size_t>(h)];

        // A hit with distanceSquared == max() is an empty slot (fewer than k neighbors existed).
        if (hit.distanceSquared == std::numeric_limits<Real>::max()) {
          break;
        }
        const std::size_t c       = hit.index;
        const bool        isLocal = c < nLocalPC;
        const std::size_t idx     = isLocal ? c : c - nLocalPC;
        const Real        d2      = hit.distanceSquared;

        if (maxDistBoundSquared.has_value() && d2 > maxDistBoundSquared.value()) {
          continue;
        }

        const ParticleID gid = isLocal ? a_localValid[idx].globalID : a_ghosts[idx].globalID;

        if (a_consumedIDs.count(gid) > 0) {
          continue;
        }

        considerCandidate(d2, idx, isLocal);
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

template <typename Packed>
inline void
buildNNCellSpatialIndex(NNCellSpatialIndex&                       a_index,
                        const std::vector<MergeParticle<Packed>>& a_localValid,
                        const std::vector<MergeParticle<Packed>>& a_ghosts,
                        const RealVect&                           a_probLo,
                        const RealVect&                           a_dx,
                        const NNCellBudget&                       a_cellBudget) noexcept
{
  CH_TIME("ParticleManagement::buildNNCellSpatialIndex");

  a_index.clear();

  const std::size_t total = a_localValid.size() + a_ghosts.size();
  if (total == 0) {
    return;
  }

  // Exact per-cell pre-filter, same contract as buildNNSpatialIndex()'s own: build nothing for
  // this patch unless some cell is over threshold.
  if (!nnMergeAnyCellCrowded(a_localValid, a_ghosts, a_probLo, a_dx, a_cellBudget)) {
    return;
  }

  // Bucket every particle (local and ghost alike) by cell, then build one cloud per occupied
  // cell -- not just crowded ones, since an under-threshold cell can still be a valid Moore-
  // neighbor candidate source (see findNearestNeighborCandidatesOneCell()).
  std::unordered_map<NNWalkCell,
                     std::pair<std::vector<EBGeometry::Vec3T<Real>>, std::vector<NNCellMeta>>,
                     LevelTiles::TileHasher>
    perCell;
  perCell.reserve(total);

  auto pushPoint = [&](const MergeParticle<Packed>& a_p, const bool a_isLocal) {
    const NNWalkCell cell = cellKeyOf(a_p.position, a_probLo, a_dx);

    EBGeometry::Vec3T<Real> v(0.0, 0.0, 0.0);

    for (int dir = 0; dir < SpaceDim; dir++) {
      v[dir] = a_p.position[dir];
    }

    auto& bucket = perCell[cell];

    bucket.first.push_back(v);
    bucket.second.push_back(NNCellMeta{a_p.globalID, a_isLocal, a_p.level});
  };

  for (const auto& p : a_localValid) {
    pushPoint(p, true);
  }

  for (const auto& p : a_ghosts) {
    pushPoint(p, false);
  }

  a_index.reserve(perCell.size());

  for (auto& kv : perCell) {
    a_index[kv.first] = std::make_shared<NNCellCloud>(kv.second.first, kv.second.second);
  }
}

template <typename Packed>
inline void
findNearestNeighborCandidatesOneCell(std::vector<NNMergeEdge<Packed>>&                          a_edges,
                                     const std::vector<MergeParticle<Packed>>&                  a_localValid,
                                     const std::unordered_set<ParticleID>&                      a_consumedIDs,
                                     const NNCellSpatialIndex&                                  a_cellIndex,
                                     const std::unordered_map<NNCellKey, int, NNCellKeyHasher>& a_liveCellCount,
                                     const RealVect&                                            a_probLo,
                                     const RealVect&                                            a_dx,
                                     const std::vector<RealVect>&                               a_dxByLevel,
                                     const NNCellBudget&                                        a_cellBudget,
                                     const int a_maxFallbackCandidates) noexcept
{
  CH_TIME("ParticleManagement::findNearestNeighborCandidatesOneCell");

  a_edges.clear();

  const int K = 1 + a_maxFallbackCandidates;

  struct RankedCandidate
  {
    Real       d2;
    ParticleID id;
    bool       isLocal;
  };

  // The 3^SpaceDim Moore-neighborhood offsets (own cell included), computed once and reused for
  // every query in this call.
  std::vector<IntVect> offsets;
  offsets.reserve(27);

  const Box neighborhood(-IntVect::Unit, IntVect::Unit);

  for (BoxIterator bit(neighborhood); bit.ok(); ++bit) {
    offsets.push_back(bit());
  }

  using PCBVH = NNCellCloud;

  for (const MergeParticle<Packed>& q : a_localValid) {
    // Lazy deletion -- see findNearestNeighborCandidates()'s own docs.
    if (a_consumedIDs.count(q.globalID) > 0) {
      continue;
    }

    const IntVect   qCell = cellKeyOf(q.position, a_probLo, a_dx);
    const NNCellKey qKey{q.level, qCell};
    const auto      cntIt = a_liveCellCount.find(qKey);
    const int       count = (cntIt != a_liveCellCount.end()) ? cntIt->second : 0;

    if (count <= a_cellBudget(qKey)) {
      continue;
    }

    EBGeometry::Vec3T<Real> queryVec(0.0, 0.0, 0.0);

    for (int dir = 0; dir < SpaceDim; dir++) {
      queryVec[dir] = q.position[dir];
    }

    // Rank the Moore-neighborhood cells by their EXACT minimum possible distance to the query,
    // ascending. This is a true lower bound (point-to-AABB distance), so once a full top-K set is
    // found, no cell beyond the current worst-of-K's distance can ever improve it -- see the break
    // below.
    std::vector<std::pair<Real, IntVect>> rankedCells;
    rankedCells.reserve(offsets.size());

    for (const IntVect& offset : offsets) {
      const IntVect cell = qCell + offset;
      const Real d2 = (offset == IntVect::Zero) ? 0.0 : nnMergePointToCellDistSquared(q.position, cell, a_probLo, a_dx);

      rankedCells.push_back(std::make_pair(d2, cell));
    }

    std::sort(rankedCells.begin(),
              rankedCells.end(),
              [](const std::pair<Real, IntVect>& a_a, const std::pair<Real, IntVect>& a_b) {
                if (a_a.first != a_b.first) {
                  return a_a.first < a_b.first;
                }

                // Deterministic tiebreak -- IntVect::operator< is a component-wise dominance
                // partial order, not a strict weak ordering (see NNCellKeyHasher's own docs), so a
                // plain lexicographic compare is used instead.
                for (int dir = 0; dir < SpaceDim; dir++) {
                  if (a_a.second[dir] != a_b.second[dir]) {
                    return a_a.second[dir] < a_b.second[dir];
                  }
                }

                return false;
              });

    std::vector<RankedCandidate> best;
    best.reserve(static_cast<std::size_t>(K) + 1);

    auto considerCandidate = [&](const Real a_d2, const ParticleID a_id, const bool a_isLocal) {
      if (static_cast<int>(best.size()) == K && a_d2 >= best.back().d2) {
        return;
      }

      const auto insertPos = std::upper_bound(best.begin(),
                                              best.end(),
                                              a_d2,
                                              [](const Real a_val, const RankedCandidate& a_c) {
                                                return a_val < a_c.d2;
                                              });

      best.insert(insertPos, RankedCandidate{a_d2, a_id, a_isLocal});

      if (static_cast<int>(best.size()) > K) {
        best.pop_back();
      }
    };

    std::vector<typename PCBVH::Hit> hitBuf;

    for (const auto& ranked : rankedCells) {
      const Real cellDistSquared = ranked.first;

      if (static_cast<int>(best.size()) == K && cellDistSquared > best.back().d2) {
        break;
      }

      const auto it = a_cellIndex.find(ranked.second);

      if (it == a_cellIndex.end() || it->second == nullptr) {
        continue;
      }

      const PCBVH& cloud = *(it->second);

      // The own cell also matches the query itself, so one extra neighbor is requested and the
      // self match filtered out below; a neighbor cell can never contain the query, so no extra
      // slot is needed there.
      const bool        isOwnCell = (ranked.second == qCell);
      const std::size_t kQuery    = isOwnCell ? static_cast<std::size_t>(K + 1) : static_cast<std::size_t>(K);

      hitBuf.resize(kQuery);

      const std::size_t found = cloud.closestPoints(queryVec, kQuery, hitBuf.data());

      for (std::size_t h = 0; h < found; h++) {
        const typename PCBVH::Hit& hit  = hitBuf[h];
        const NNCellMeta&          meta = cloud.metadata(hit.index);

        if (meta.globalID == q.globalID) {
          continue;
        }

        if (a_consumedIDs.count(meta.globalID) > 0) {
          continue;
        }

        // Own-cell/Moore-adjacency in the query patch's own-level frame does not by itself bound
        // the Chebyshev distance in the finer-of-the-two-levels frame that resolveTrivialTier()/
        // judgeProposals() actually enforce (see this function's own docs) -- a same-level
        // candidate needs no extra check (query and candidate already share one frame), but a
        // cross-level one is re-checked here, exactly, before ever being considered.
        if (meta.level != q.level) {
          const EBGeometry::Vec3T<Real>& candVec = cloud.position(hit.index);

          RealVect candPos(D_DECL(candVec[0], candVec[1], candVec[2]));

          if (nnMergeCrossLevelTooFar(q.position, q.level, candPos, meta.level, a_probLo, a_dxByLevel, 1)) {
            continue;
          }
        }

        considerCandidate(hit.distanceSquared, meta.globalID, meta.isLocal);
      }
    }

    if (best.empty()) {
      continue;
    }

    NNMergeEdge<Packed> edge;

    edge.distanceSquared  = best[0].d2;
    edge.queryID          = q.globalID;
    edge.candidateID      = best[0].id;
    edge.candidateIsLocal = best[0].isLocal;

    for (std::size_t k = 1; k < best.size(); k++) {
      NNFallbackCandidate fc;

      fc.candidateID      = best[k].id;
      fc.distanceSquared  = best[k].d2;
      fc.candidateIsLocal = best[k].isLocal;

      edge.fallbackCandidates.push_back(fc);
    }

    a_edges.push_back(std::move(edge));
  }
}

template <typename Packed, typename Combine, typename PosValid, typename Allocator, typename OnMerge>
inline void
resolveTrivialTier(std::unordered_map<NNCellKey, int, NNCellKeyHasher>&            a_liveCellCount,
                   std::unordered_set<ParticleID>&                                 a_consumedIDs,
                   std::unordered_set<ParticleID>&                                 a_hasOutgoingCommitment,
                   std::vector<NNMergeResult<Packed>>&                             a_results,
                   const std::vector<NNMergeEdge<Packed>>&                         a_edges,
                   const std::unordered_map<ParticleID, NNPooledParticle<Packed>>& a_particlesByID,
                   const Combine&                                                  a_combine,
                   const NNCellBudget&                                             a_cellBudget,
                   const RealVect&                                                 a_probLo,
                   const std::vector<RealVect>&                                    a_dxByLevel,
                   const PosValid&                                                 a_isPositionValid,
                   const Allocator&                                                a_allocateID,
                   const int                                                       a_maxFallbackCandidates,
                   const std::optional<int>                                        a_maxCellDistance,
                   unsigned long long*                                             a_crossLevelMergeCount,
                   const OnMerge&                                                  a_onMergeCommitted) noexcept
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

    const auto queryDataIt = a_particlesByID.find(edge.queryID);

    if (queryDataIt == a_particlesByID.end()) {
      continue;
    }

    const MergeParticle<Packed>& query        = queryDataIt->second.particle;
    const bool                   queryExposed = queryDataIt->second.exposed;

    bool resolved           = false; // an outcome (commit, or a legitimate "leave unmatched") was reached
    bool originalStale      = false;
    bool originalIneligible = false; // too far or invalid position, evaluated on the ORIGINAL candidate

    const int numAttempts = a_maxFallbackCandidates + 1;

    for (int attempt = 0; attempt < numAttempts; attempt++) {
      // A boundary-exposed query can never merge trivially -- trivialOk requires !queryExposed
      // unconditionally, so no fallback candidate could ever help. Run the ORIGINAL candidate
      // (attempt 0) to set originalStale/originalIneligible, then stop rather than walking doomed
      // fallbacks: each would cost consumed/busy lookups plus, for an alive candidate, a
      // cross-level distance calc and an a_isPositionValid() call (an EB point query for EB-aware
      // callers). Outcome is identical -- the query still proposes off its original candidate.
      if (queryExposed && attempt > 0) {
        break;
      }

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

      const auto candDataIt      = a_particlesByID.find(candidateID);
      const bool candidateExists = candDataIt != a_particlesByID.end();

      // exposed is meaningful only for local candidates (ghosts carry false); the single
      // particlesByID lookup above now also yields it, so there is no separate exposed map.
      const bool candidateExposed = candidateIsLocal && candidateExists && candDataIt->second.exposed;

      bool tooFar     = false;
      bool invalidPos = false;

      RealVect mergedPos(D_DECL(0.0, 0.0, 0.0));
      Real     mergedWeight = 0.0;

      if (!candidateStale && candidateExists) {
        const MergeParticle<Packed>& cand = candDataIt->second.particle;

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
          nnMergeWeightedCentroid(mergedPos, mergedWeight, query.position, query.weight, cand.position, cand.weight);
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
        const MergeParticle<Packed>& cand = candDataIt->second.particle;

        // Each participant's OWN level's dx -- see a_dxByLevel's docs. The cell key must match
        // the one a_liveCellCount was populated under for that participant, so a single shared
        // dx would break the recheck for every level but one.
        const NNCellKey qKey(query.level, cellKeyOf(query.position, a_probLo, a_dxByLevel[query.level]));
        const NNCellKey cKey(cand.level, cellKeyOf(cand.position, a_probLo, a_dxByLevel[cand.level]));

        const auto qCountIt = a_liveCellCount.find(qKey);
        const auto cCountIt = a_liveCellCount.find(cKey);

        const int qCount = (qCountIt != a_liveCellCount.end()) ? qCountIt->second : 0;
        const int cCount = (cCountIt != a_liveCellCount.end()) ? cCountIt->second : 0;

        if (qCount > a_cellBudget(qKey) && cCount > a_cellBudget(cKey)) {
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
generateProposals(std::vector<std::vector<NNMergeProposal<Packed>>>&              a_proposalsByDestRank,
                  const std::vector<NNMergeEdge<Packed>>&                         a_edges,
                  const std::unordered_set<ParticleID>&                           a_hasOutgoingCommitment,
                  const std::unordered_map<ParticleID, NNPooledParticle<Packed>>& a_particlesByID,
                  const int                                                       a_numRanks) noexcept
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

    proposal.source          = sourceIt->second.particle;
    proposal.targetID        = edge.candidateID;
    proposal.distanceSquared = edge.distanceSquared;

    const int destRank = static_cast<int>(targetIt->second.particle.ownerRank);

    CH_assert(destRank >= 0 && destRank < a_numRanks);

    a_proposalsByDestRank[destRank].push_back(proposal);
  }
}

template <typename Packed, typename Combine, typename PosValid, typename Allocator, typename OnMerge>
inline void
judgeProposals(std::unordered_map<NNCellKey, int, NNCellKeyHasher>&            a_liveCellCount,
               std::vector<NNMergeResult<Packed>>&                             a_results,
               std::vector<std::vector<NNMergeVerdict>>&                       a_verdictsByDestRank,
               const std::vector<NNMergeProposal<Packed>>&                     a_incomingProposals,
               const std::unordered_map<ParticleID, NNPooledParticle<Packed>>& a_particlesByID,
               const std::unordered_set<ParticleID>&                           a_consumedIDs,
               const std::unordered_set<ParticleID>&                           a_hasOutgoingCommitment,
               const std::map<ParticleID, ParticleID>&                         a_outgoingTargetOf,
               const Combine&                                                  a_combine,
               const NNCellBudget&                                             a_cellBudget,
               const RealVect&                                                 a_probLo,
               const std::vector<RealVect>&                                    a_dxByLevel,
               const PosValid&                                                 a_isPositionValid,
               const Allocator&                                                a_allocateID,
               const std::optional<int>                                        a_maxCellDistance,
               unsigned long long*                                             a_crossLevelMergeCount,
               const OnMerge&                                                  a_onMergeCommitted) noexcept
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

    const auto targetDataIt = a_particlesByID.find(targetID);

    if (targetDataIt == a_particlesByID.end() || a_consumedIDs.count(targetID) > 0) {
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
    const auto                     targetIt = a_particlesByID.find(decision.targetID);

    CH_assert(targetIt != a_particlesByID.end() && a_consumedIDs.count(decision.targetID) == 0);

    const MergeParticle<Packed>& target = targetIt->second.particle;
    const MergeParticle<Packed>& source = proposal.source;

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

    nnMergeWeightedCentroid(mergedPos, mergedWeight, target.position, target.weight, source.position, source.weight);

    if (!a_isPositionValid(mergedPos)) {
      a_verdictsByDestRank[sourceDestRank].push_back(NNMergeVerdict{source.globalID, false});

      continue;
    }

    // Dynamic crowding recheck for BOTH participants (mirrors resolveTrivialTier()'s symmetric
    // recheck -- each participant's OWN level's dx, see a_dxByLevel's docs). The source's cell
    // must be rechecked too: an earlier merge this round can drain it below threshold between
    // proposal generation and judging, and accepting anyway would violate the crowding trigger on
    // the source side. If this rank has no live entry for the source's cell (never shipped here as
    // a ghost), sCount defaults to 0 -- rejected -- the same conservative default used for an
    // unknown target count.
    const NNCellKey tKey(target.level, cellKeyOf(target.position, a_probLo, a_dxByLevel[target.level]));
    const NNCellKey sKey(source.level, cellKeyOf(source.position, a_probLo, a_dxByLevel[source.level]));

    const auto tCountIt = a_liveCellCount.find(tKey);
    const auto sCountIt = a_liveCellCount.find(sKey);

    const int tCount = (tCountIt != a_liveCellCount.end()) ? tCountIt->second : 0;
    const int sCount = (sCountIt != a_liveCellCount.end()) ? sCountIt->second : 0;

    if (tCount <= a_cellBudget(tKey) || sCount <= a_cellBudget(sKey)) {
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

    // Symmetric decrement (mirrors resolveTrivialTier()'s commit), but only when this rank has a
    // live entry for the source's cell (sCountIt found above); otherwise there is nothing local
    // to decrement -- the source's owning rank updates its own copy on the accepted verdict.
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
applyVerdicts(std::unordered_set<ParticleID>&    a_consumedIDs,
              const std::vector<NNMergeVerdict>& a_incomingVerdicts) noexcept
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
placeMergedParticles(std::map<std::pair<int, unsigned int>, std::vector<MergeParticle<Packed>>>& a_localByPatch,
                     std::vector<std::vector<MergeParticle<Packed>>>&                            a_scatterByDestRank,
                     const std::vector<NNMergeResult<Packed>>&                                   a_results,
                     const PositionLocator&                                                      a_locate,
                     const int                                                                   a_thisRank) noexcept
{
  CH_TIME("ParticleManagement::placeMergedParticles");

  a_localByPatch.clear();
  // a_scatterByDestRank is APPENDED to, never resized here -- see its own docs.

  for (const auto& result : a_results) {
    // O(1) point->block query -- see PositionLocator. Never a linear scan over
    // patches (the bottleneck this avoids at large patch counts).
    const LevelTiles::LevelAndBox loc = a_locate(result.merged.position);

    // A merged particle landing outside every known patch (e.g. exactly on/outside the domain
    // boundary due to floating-point roundoff in the weighted centroid) is not expected in
    // normal operation; MayDay rather than silently dropping mass.
    if (!loc.valid) {
      MayDay::Error("ParticleManagement::placeMergedParticles -- merged particle position not "
                    "found in any patch");
    }

    if (loc.rank == a_thisRank) {
      a_localByPatch[std::make_pair(loc.level, loc.gridIndex)].push_back(result.merged);
    }
    else {
      a_scatterByDestRank[loc.rank].push_back(result.merged);
    }
  }
}

template <typename P, typename Packed, typename Traits, typename Gather>
inline void
gatherMergeParticles(ParticleContainer<P, Traits>&                             a_particles,
                     std::vector<NNPatchWork>&                                 a_patchWork,
                     std::vector<std::vector<MergeParticle<Packed>>>&          a_patchLocalValid,
                     std::vector<std::vector<MergeParticle<Packed>>>&          a_patchGhosts,
                     std::unordered_map<ParticleID, NNPooledParticle<Packed>>& a_particlesByID,
                     std::unordered_map<NNCellKey, int, NNCellKeyHasher>&      a_liveCellCount,
                     const AmrMesh&                                            a_amr,
                     const std::string&                                        a_realm,
                     const int                                                 a_finestLevel,
                     const RealVect&                                           a_probLo,
                     const std::vector<RealVect>&                              a_dxByLevel,
                     const int                                                 a_ghostWidth,
                     const Gather&                                             a_gather) noexcept
{
  CH_TIME("ParticleManagement::gatherMergeParticles");

  a_patchWork.clear();
  a_patchLocalValid.clear();
  a_patchGhosts.clear();
  a_particlesByID.clear();
  a_liveCellCount.clear();

  // Which cells hold a particle that some other box can also see -- the trivial tier's eligibility test.
  // Read from the realm, which derives it from the three ghost masks when they are built (see
  // Realm::m_particleGhostExposure), rather than re-deriving it here per particle per timestep.
  const AMRMask& exposure = a_amr.getParticleGhostExposure(a_realm, a_ghostWidth);

  for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
    const RealVect& dx = a_dxByLevel[lvl];

    const DisjointBoxLayout& dbl = a_amr.getGrids(a_realm)[lvl];

    const DataIterator& dit  = dbl.dataIterator();
    const int           nbox = dit.size();

    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      std::vector<MergeParticle<Packed>> localValid;
      std::vector<MergeParticle<Packed>> ghosts;

      localValid.reserve(leaf.size());
      ghosts.reserve(leaf.size());

      a_particlesByID.reserve(a_particlesByID.size() + leaf.size());

      const BaseFab<bool>& exposureDin = (*exposure[lvl])[din];

      for (std::size_t i = 0; i < leaf.size(); i++) {
        const bool isGhost = leaf.isGhost(i);

        int originLevel = lvl;
        if (isGhost) {
          switch (leaf.ghost(i)) {
          case GhostType::Coarse: {
            originLevel = lvl - 1;

            break;
          }
          case GhostType::Fine: {
            originLevel = lvl + 1;

            break;
          }
          default: {
            originLevel = lvl; // GhostType::SameLevel

            break;
          }
          }
        }

        MergeParticle<Packed> p;

        p.position  = leaf.position(i);
        p.weight    = leaf.weight(i);
        p.globalID  = leaf.particleID(i);
        p.ownerRank = leaf.rankID(i);
        p.level     = originLevel;
        p.payload   = a_gather(leaf, i);

        const IntVect cell = cellKeyOf(p.position, a_probLo, a_dxByLevel[originLevel]);

        if (isGhost) {
          ghosts.push_back(p);
          a_particlesByID[p.globalID] = NNPooledParticle<Packed>{p, false};
        }
        else {
          // A non-ghost particle is resident in this patch, so its cell lies in this box -- which is
          // exactly the region the exposure mask is defined over (it carries no ghost cells).
          CH_assert(exposureDin.box().contains(cell));

          const bool isExposed = exposureDin(cell, 0);

          localValid.push_back(p);
          a_particlesByID[p.globalID] = NNPooledParticle<Packed>{p, isExposed};
        }
      }

      a_patchWork.push_back(NNPatchWork{lvl, din, dx});
      a_patchLocalValid.push_back(std::move(localValid));
      a_patchGhosts.push_back(std::move(ghosts));
    }
  }

  // Derived from the pooled a_particlesByID, NOT inline in the gather loop above -- a particle can
  // appear as a ghost in several of this rank's own patches at once, so counting per occurrence
  // would multiply-count it; particlesByID collapses those to one entry per id.
  for (const auto& kv : a_particlesByID) {
    const MergeParticle<Packed>& p    = kv.second.particle;
    const IntVect                cell = cellKeyOf(p.position, a_probLo, a_dxByLevel[p.level]);

    a_liveCellCount[NNCellKey(p.level, cell)] += 1;
  }
}

template <typename P,
          typename Packed,
          typename Traits,
          typename Combine,
          typename Scatter,
          typename PosValid,
          typename Allocator,
          typename OnMerge>
inline void
finishMergeRound(ParticleContainer<P, Traits>&                             a_particles,
                 std::unordered_map<ParticleID, NNPooledParticle<Packed>>& a_particlesByID,
                 std::unordered_set<ParticleID>&                           a_consumedIDs,
                 std::unordered_map<NNCellKey, int, NNCellKeyHasher>&      a_liveCellCount,
                 std::vector<NNMergeResult<Packed>>&                       a_allResults,
                 unsigned long long*                                       a_crossLevelMergeCount,
                 const AmrMesh&                                            a_amr,
                 const std::string&                                        a_realm,
                 const NNCellBudget&                                       a_cellBudget,
                 const Combine&                                            a_combine,
                 const Scatter&                                            a_scatter,
                 const Allocator&                                          a_allocateID,
                 const std::optional<int>                                  a_maxCellDistance,
                 const PosValid&                                           a_isPositionValid,
                 const RealVect&                                           a_probLo,
                 const std::vector<RealVect>&                              a_dxByLevel,
                 const std::vector<NNPatchWork>&                           a_patchWork,
                 const std::vector<NNMergeEdge<Packed>>&                   a_pooledEdges,
                 const std::unordered_set<ParticleID>&                     a_hasOutgoingCommitment,
                 const OnMerge&                                            a_onMergeCommitted) noexcept
{
  CH_TIME("ParticleManagement::finishMergeRound");

  const int myRank   = procID();
  const int numRanks = static_cast<int>(numProc());

  // outgoingTargetOf: for every query id with an outgoing commitment, remember its own nearest-
  // neighbor candidate id (needed by judgeProposals()'s mutual-match rule). Nothing is consumed
  // between the local tier's last pass and here, so a_pooledEdges is still exactly current.
  std::map<ParticleID, ParticleID> outgoingTargetOf;

  for (const auto& edge : a_pooledEdges) {
    if (a_hasOutgoingCommitment.count(edge.queryID) > 0) {
      outgoingTargetOf[edge.queryID] = edge.candidateID;
    }
  }

  std::vector<std::vector<NNMergeProposal<Packed>>> proposalsByDestRank;
  generateProposals(proposalsByDestRank, a_pooledEdges, a_hasOutgoingCommitment, a_particlesByID, numRanks);

  const std::vector<NNMergeProposal<Packed>> incomingProposals = nnMergeExchangeByRank(proposalsByDestRank);

  std::vector<NNMergeResult<Packed>>       judgeResults;
  std::vector<std::vector<NNMergeVerdict>> verdictsByDestRank(numRanks);

  judgeProposals(a_liveCellCount,
                 judgeResults,
                 verdictsByDestRank,
                 incomingProposals,
                 a_particlesByID,
                 a_consumedIDs,
                 a_hasOutgoingCommitment,
                 outgoingTargetOf,
                 a_combine,
                 a_cellBudget,
                 a_probLo,
                 a_dxByLevel,
                 a_isPositionValid,
                 a_allocateID,
                 a_maxCellDistance,
                 a_crossLevelMergeCount,
                 a_onMergeCommitted);

  a_allResults.insert(a_allResults.end(), judgeResults.begin(), judgeResults.end());

  for (const auto& r : judgeResults) {
    a_consumedIDs.insert(r.firstID);
  }

  const std::vector<NNMergeVerdict> incomingVerdicts = nnMergeExchangeByRank(verdictsByDestRank);

  applyVerdicts(a_consumedIDs, incomingVerdicts);

  // O(1) point->block query. ParticleContainer::findDestination() delegates to the shared
  // LevelTiles::findDestination core, so this is the same mapping remap() routes its own particles
  // through -- never a linear scan over patches.
  const PositionLocator locate = [&a_particles](const RealVect& a_pos) -> LevelTiles::LevelAndBox {
    return a_particles.findDestination(a_pos);
  };

  std::map<std::pair<int, unsigned int>, std::vector<MergeParticle<Packed>>> localByPatch;
  std::vector<std::vector<MergeParticle<Packed>>>                            scatterByDestRank(numRanks);

  placeMergedParticles(localByPatch, scatterByDestRank, a_allResults, locate, myRank);

  const std::vector<MergeParticle<Packed>> incomingScattered = nnMergeExchangeByRank(scatterByDestRank);

  // ---- Apply everything to the real ParticleSoA containers ----
  for (std::size_t pw = 0; pw < a_patchWork.size(); pw++) {
    ParticleSoA<P, Traits>& leaf = a_particles[a_patchWork[pw].level][a_patchWork[pw].din];

    for (std::size_t i = 0; i < leaf.size();) {
      if (!leaf.isGhost(i) && a_consumedIDs.count(leaf.particleID(i)) > 0) {
        leaf.remove(i);
      }
      else {
        i++;
      }
    }
  }

  a_particles.clearGhostParticles();

  for (const auto& kv : localByPatch) {
    const int          level     = kv.first.first;
    const unsigned int gridIndex = kv.first.second;

    const DataIndex din = a_amr.getLevelTiles(a_realm)[level]->getMyGrids().at(gridIndex);

    ParticleSoA<P, Traits>& leaf = a_particles[level][din];

    for (const auto& merged : kv.second) {
      // Correct the provisional ownerRank to this rank -- see mergeNearestNeighborsTree()'s own
      // docs on why this matters (a stale rank would misroute a later merge's verdict).
      MergeParticle<Packed> corrected = merged;

      corrected.ownerRank = myRank;
      a_scatter(leaf, corrected);
    }
  }

  for (const auto& merged : incomingScattered) {
    const LevelTiles::LevelAndBox loc = locate(merged.position);

    if (!loc.valid || loc.rank != myRank) {
      MayDay::Error("ParticleManagement::finishMergeRound -- incoming scattered particle not "
                    "found in any of this rank's own patches");
    }

    const DataIndex din = a_amr.getLevelTiles(a_realm)[loc.level]->getMyGrids().at(loc.gridIndex);

    MergeParticle<Packed> corrected = merged;
    corrected.ownerRank             = myRank;

    ParticleSoA<P, Traits>& leaf = a_particles[loc.level][din];
    a_scatter(leaf, corrected);
  }
}

template <typename P,
          typename Packed,
          typename Traits,
          typename Cloud,
          typename Gather,
          typename Combine,
          typename Scatter,
          typename Allocator,
          typename PosValid,
          typename OnMerge>
inline void
mergeNearestNeighborsRoundImpl(ParticleContainer<P, Traits>& a_particles,
                               unsigned long long*           a_crossLevelMergeCount,
                               const AmrMesh&                a_amr,
                               const NNCellBudget&           a_cellBudget,
                               const Gather&                 a_gather,
                               const Combine&                a_combine,
                               const Scatter&                a_scatter,
                               const Allocator&              a_allocateID,
                               const bool                    a_iterateLocalTierToConvergence,
                               const int                     a_maxFallbackCandidates,
                               const std::optional<int>      a_maxCellDistance,
                               const int                     a_ghostWidth,
                               const PosValid&               a_isPositionValid,
                               const OnMerge&                a_onMergeCommitted)
{
  CH_TIME("ParticleManagement::mergeNearestNeighborsRoundImpl");

  const std::string realm       = a_particles.getRealm();
  const int         finestLevel = a_amr.getFinestLevel();
  const RealVect    probLo      = a_amr.getProbLo();

  // Per-level dx, indexed by level -- every cell-key computation below and downstream (in
  // resolveTrivialTier()/judgeProposals()) MUST use the dx of the PARTICLE'S OWN level (see
  // MergeParticle::level's docs), never a single dx shared across levels.
  const std::vector<RealVect>& dxByLevel = a_amr.getDxAsRealVect().constStdVector();

  CH_assert(a_ghostWidth >= 1);
  const int ghostWidth = a_ghostWidth;

  // Ghost particles are assumed ALREADY FILLED by the caller (an a_ghostWidth-cell halo -- one cell
  // per unit of a_ghostWidth in every direction, which is what the crowding trigger and the candidate
  // search look at) and are used here as candidates. This function does not fill them: filling is the
  // caller's responsibility, on equal footing with the caller cell-sorting before a cell-based merge
  // (see ItoSolver / ParticleManagement::ParticleMergeKind). A particle ghost mask of this same width
  // must have been registered on this realm before the grids were (re)built, and every rank must use
  // the same a_ghostWidth (see the a_ghostWidth docs -- it defines boundary-exposure).
  std::vector<NNPatchWork>                                 patchWork;
  std::vector<std::vector<MergeParticle<Packed>>>          patchLocalValid; // parallel to patchWork
  std::vector<std::vector<MergeParticle<Packed>>>          patchGhosts;     // parallel to patchWork
  std::unordered_map<ParticleID, NNPooledParticle<Packed>> particlesByID;
  std::unordered_map<NNCellKey, int, NNCellKeyHasher>      liveCellCount;

  gatherMergeParticles<P, Packed, Traits>(a_particles,
                                          patchWork,
                                          patchLocalValid,
                                          patchGhosts,
                                          particlesByID,
                                          liveCellCount,
                                          a_amr,
                                          realm,
                                          finestLevel,
                                          probLo,
                                          dxByLevel,
                                          ghostWidth,
                                          a_gather);

  // Backend-specific: one whole-patch point-cloud index (see Cloud) per patch, built once here, per
  // round, over the FULL localValid/ghosts sets, before any consumption this round -- and reused
  // across every internal a_iterateLocalTierToConvergence pass via lazy deletion (see
  // findNearestNeighborCandidates()'s docs).
  std::vector<NNSpatialIndex<Packed, Cloud>> patchSpatialIndex;
  patchSpatialIndex.reserve(patchWork.size());

  for (std::size_t pw = 0; pw < patchWork.size(); pw++) {
    NNSpatialIndex<Packed, Cloud> spatialIndex;

    buildNNSpatialIndex(spatialIndex,
                        patchLocalValid[pw],
                        patchGhosts[pw],
                        probLo,
                        patchWork[pw].dx,
                        a_cellBudget,
                        1 + a_maxFallbackCandidates);

    patchSpatialIndex.push_back(std::move(spatialIndex));
  }

  std::unordered_set<ParticleID> consumedIDs;
  std::unordered_set<ParticleID> hasOutgoingCommitment;

  std::vector<NNMergeResult<Packed>> allResults;

  // 2. Per-patch candidate search + trivial-tier resolution, pooled rank-wide.
  //
  // lastPooledEdges deliberately survives past the loop below: it is reused (NOT recomputed) by
  // finishMergeRound() afterward. Nothing consumes any particles between the last runLocalTierPass()
  // and finishMergeRound(), so the edges are still exactly current; re-running the candidate search
  // purely to re-derive them would repeat the most expensive step for no new information.
  std::vector<NNMergeEdge<Packed>> lastPooledEdges;

  auto runLocalTierPass = [&]() -> bool {
    lastPooledEdges.clear();

    for (std::size_t pw = 0; pw < patchWork.size(); pw++) {
      // patchLocalValid[pw]/patchGhosts[pw] are this patch's FULL sets and patchSpatialIndex[pw]
      // was built once over them before the local-tier passes began. findNearestNeighborCandidates()
      // filters against consumedIDs itself (see its docs), so nothing is rebuilt or re-filtered here.
      std::vector<NNMergeEdge<Packed>> patchEdges;

      findNearestNeighborCandidates(patchEdges,
                                    patchLocalValid[pw],
                                    patchGhosts[pw],
                                    consumedIDs,
                                    patchSpatialIndex[pw],
                                    liveCellCount,
                                    probLo,
                                    patchWork[pw].dx,
                                    a_cellBudget,
                                    a_maxFallbackCandidates,
                                    a_maxCellDistance);

      lastPooledEdges.insert(lastPooledEdges.end(), patchEdges.begin(), patchEdges.end());
    }

    std::vector<NNMergeResult<Packed>> passResults;

    resolveTrivialTier(liveCellCount,
                       consumedIDs,
                       hasOutgoingCommitment,
                       passResults,
                       lastPooledEdges,
                       particlesByID,
                       a_combine,
                       a_cellBudget,
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

  // 3-7. Propose/judge/verdict, placement, and container application -- shared with
  // mergeNearestNeighborsOneCell() (see finishMergeRound()'s own docs).
  finishMergeRound<P, Packed, Traits>(a_particles,
                                      particlesByID,
                                      consumedIDs,
                                      liveCellCount,
                                      allResults,
                                      a_crossLevelMergeCount,
                                      a_amr,
                                      realm,
                                      a_cellBudget,
                                      a_combine,
                                      a_scatter,
                                      a_allocateID,
                                      a_maxCellDistance,
                                      a_isPositionValid,
                                      probLo,
                                      dxByLevel,
                                      patchWork,
                                      lastPooledEdges,
                                      hasOutgoingCommitment,
                                      a_onMergeCommitted);
}

} // namespace detail

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
mergeNearestNeighborsTree(ParticleContainer<P, Traits>& a_particles,
                          const AmrMesh&                a_amr,
                          const NNCellBudget&           a_cellBudget,
                          const Gather&                 a_gather,
                          const Combine&                a_combine,
                          const Scatter&                a_scatter,
                          const Allocator&              a_allocateID,
                          const bool                    a_iterateLocalTierToConvergence,
                          const int                     a_maxFallbackCandidates,
                          const std::optional<int>      a_maxCellDistance,
                          const int                     a_ghostWidth,
                          const PosValid&               a_isPositionValid,
                          unsigned long long*           a_crossLevelMergeCount,
                          const OnMerge&                a_onMergeCommitted)
{
  CH_TIME("ParticleManagement::mergeNearestNeighborsTree");

  detail::mergeNearestNeighborsRoundImpl<P, Packed, Traits, EBGeometry::PointCloudBVH<Real, ParticleID>>(
    a_particles,
    a_crossLevelMergeCount,
    a_amr,
    a_cellBudget,
    a_gather,
    a_combine,
    a_scatter,
    a_allocateID,
    a_iterateLocalTierToConvergence,
    a_maxFallbackCandidates,
    a_maxCellDistance,
    a_ghostWidth,
    a_isPositionValid,
    a_onMergeCommitted);
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
mergeNearestNeighborsHash(ParticleContainer<P, Traits>& a_particles,
                          const AmrMesh&                a_amr,
                          const NNCellBudget&           a_cellBudget,
                          const Gather&                 a_gather,
                          const Combine&                a_combine,
                          const Scatter&                a_scatter,
                          const Allocator&              a_allocateID,
                          const bool                    a_iterateLocalTierToConvergence,
                          const int                     a_maxFallbackCandidates,
                          const std::optional<int>      a_maxCellDistance,
                          const int                     a_ghostWidth,
                          const PosValid&               a_isPositionValid,
                          unsigned long long*           a_crossLevelMergeCount,
                          const OnMerge&                a_onMergeCommitted)
{
  CH_TIME("ParticleManagement::mergeNearestNeighborsHash");

  detail::mergeNearestNeighborsRoundImpl<P, Packed, Traits, EBGeometry::PointCloudHashGrid<Real, ParticleID>>(
    a_particles,
    a_crossLevelMergeCount,
    a_amr,
    a_cellBudget,
    a_gather,
    a_combine,
    a_scatter,
    a_allocateID,
    a_iterateLocalTierToConvergence,
    a_maxFallbackCandidates,
    a_maxCellDistance,
    a_ghostWidth,
    a_isPositionValid,
    a_onMergeCommitted);
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
mergeNearestNeighborsOneCell(ParticleContainer<P, Traits>& a_particles,
                             const AmrMesh&                a_amr,
                             const NNCellBudget&           a_cellBudget,
                             const Gather&                 a_gather,
                             const Combine&                a_combine,
                             const Scatter&                a_scatter,
                             const Allocator&              a_allocateID,
                             const bool                    a_iterateLocalTierToConvergence,
                             const int                     a_maxFallbackCandidates,
                             const PosValid&               a_isPositionValid,
                             unsigned long long*           a_crossLevelMergeCount,
                             const OnMerge&                a_onMergeCommitted)
{
  CH_TIME("ParticleManagement::mergeNearestNeighborsOneCell");

  using namespace detail;

  const std::string realm       = a_particles.getRealm();
  const int         finestLevel = a_amr.getFinestLevel();
  const RealVect    probLo      = a_amr.getProbLo();

  // Per-level dx, indexed by level -- see mergeNearestNeighborsRoundImpl()'s own note on why every cell
  // key must use the dx of the particle's OWN level.
  const std::vector<RealVect>& dxByLevel = a_amr.getDxAsRealVect().constStdVector();

  // Structurally fixed at 1 -- findNearestNeighborCandidatesOneCell() only ever visits the query's
  // own cell and its Moore neighbors, so a width-1 halo already contains every reachable candidate
  // (see this function's own docs).
  constexpr int            ghostWidth      = 1;
  const std::optional<int> maxCellDistance = 1;

  std::vector<NNPatchWork>                                 patchWork;
  std::vector<std::vector<MergeParticle<Packed>>>          patchLocalValid;
  std::vector<std::vector<MergeParticle<Packed>>>          patchGhosts;
  std::unordered_map<ParticleID, NNPooledParticle<Packed>> particlesByID;
  std::unordered_map<NNCellKey, int, NNCellKeyHasher>      liveCellCount;

  gatherMergeParticles<P, Packed, Traits>(a_particles,
                                          patchWork,
                                          patchLocalValid,
                                          patchGhosts,
                                          particlesByID,
                                          liveCellCount,
                                          a_amr,
                                          realm,
                                          finestLevel,
                                          probLo,
                                          dxByLevel,
                                          ghostWidth,
                                          a_gather);

  // Backend-specific: one PointCloudBVH per OCCUPIED cell (not per patch), built once here, per
  // round, over the FULL localValid/ghosts sets -- see buildNNCellSpatialIndex()'s own docs.
  std::vector<NNCellSpatialIndex> patchCellSpatialIndex;
  patchCellSpatialIndex.reserve(patchWork.size());

  for (std::size_t pw = 0; pw < patchWork.size(); pw++) {
    NNCellSpatialIndex cellIndex;

    buildNNCellSpatialIndex(cellIndex, patchLocalValid[pw], patchGhosts[pw], probLo, patchWork[pw].dx, a_cellBudget);

    patchCellSpatialIndex.push_back(std::move(cellIndex));
  }

  std::unordered_set<ParticleID> consumedIDs;
  std::unordered_set<ParticleID> hasOutgoingCommitment;

  std::vector<NNMergeResult<Packed>> allResults;
  std::vector<NNMergeEdge<Packed>>   lastPooledEdges;

  auto runLocalTierPass = [&]() -> bool {
    lastPooledEdges.clear();

    for (std::size_t pw = 0; pw < patchWork.size(); pw++) {
      std::vector<NNMergeEdge<Packed>> patchEdges;

      findNearestNeighborCandidatesOneCell(patchEdges,
                                           patchLocalValid[pw],
                                           consumedIDs,
                                           patchCellSpatialIndex[pw],
                                           liveCellCount,
                                           probLo,
                                           patchWork[pw].dx,
                                           dxByLevel,
                                           a_cellBudget,
                                           a_maxFallbackCandidates);

      lastPooledEdges.insert(lastPooledEdges.end(), patchEdges.begin(), patchEdges.end());
    }

    std::vector<NNMergeResult<Packed>> passResults;

    resolveTrivialTier(liveCellCount,
                       consumedIDs,
                       hasOutgoingCommitment,
                       passResults,
                       lastPooledEdges,
                       particlesByID,
                       a_combine,
                       a_cellBudget,
                       probLo,
                       dxByLevel,
                       a_isPositionValid,
                       a_allocateID,
                       a_maxFallbackCandidates,
                       maxCellDistance,
                       a_crossLevelMergeCount,
                       a_onMergeCommitted);

    allResults.insert(allResults.end(), passResults.begin(), passResults.end());

    return !passResults.empty();
  };

  if (a_iterateLocalTierToConvergence) {
    constexpr int safetyCap = 1000;

    for (int iter = 0; iter < safetyCap; iter++) {
      hasOutgoingCommitment.clear();

      if (!runLocalTierPass()) {
        break;
      }
    }
  }
  else {
    runLocalTierPass();
  }

  finishMergeRound<P, Packed, Traits>(a_particles,
                                      particlesByID,
                                      consumedIDs,
                                      liveCellCount,
                                      allResults,
                                      a_crossLevelMergeCount,
                                      a_amr,
                                      realm,
                                      a_cellBudget,
                                      a_combine,
                                      a_scatter,
                                      a_allocateID,
                                      maxCellDistance,
                                      a_isPositionValid,
                                      probLo,
                                      dxByLevel,
                                      patchWork,
                                      lastPooledEdges,
                                      hasOutgoingCommitment,
                                      a_onMergeCommitted);
}

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#endif
