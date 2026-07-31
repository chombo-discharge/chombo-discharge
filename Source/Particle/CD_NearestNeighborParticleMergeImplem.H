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
      @details Raw cell indices computed with two different dx values are not directly comparable, so
      both positions are re-expressed in the finer of the two levels' cell-size units (higher AMR
      level = finer) before taking the Chebyshev distance. Reduces to the same-level case when
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
      @details Packs one bucket per rank and exchanges via MPI_Alltoall (byte counts) + MPI_Alltoallv
      (payload). In a serial (non-MPI) build, or with a single rank, this returns this rank's own
      bucket unchanged.
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

template <typename Packed>
inline void
buildNNSpatialIndex(const std::vector<NNMergeParticle<Packed>>& a_localValid,
                    const std::vector<NNMergeParticle<Packed>>& a_ghosts,
                    const RealVect&                             a_probLo,
                    const RealVect&                             a_dx,
                    const int                                   a_numParticlesPerCellThresh,
                    const int                                   a_kNearest,
                    NNSpatialIndex<Packed>&                     a_index) noexcept
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

  // Exact per-cell pre-filter (see this function's own docs): a count-only tally keyed as the
  // crowding trigger (cellKeyOf(pos, probLo, dx), local and ghost alike). Build nothing unless
  // some cell is over threshold. Checks the true per-cell maximum, not an average density (a
  // low-average patch can still hold one over-threshold cell).
  std::unordered_map<NNWalkCell, int, LevelTiles::TileHasher> cellCount;
  cellCount.reserve(total);

  for (const auto& p : a_localValid) {
    cellCount[cellKeyOf(p.position, a_probLo, a_dx)]++;
  }

  for (const auto& p : a_ghosts) {
    cellCount[cellKeyOf(p.position, a_probLo, a_dx)]++;
  }

  bool anyCellCrowded = false;

  for (const auto& kv : cellCount) {
    if (kv.second > a_numParticlesPerCellThresh) {
      anyCellCrowded = true;
      break;
    }
  }

  if (!anyCellCrowded) {
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

  auto pushPoint = [&](const NNMergeParticle<Packed>& a_p) {
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

  a_index.pointCloud       = std::make_shared<EBGeometry::PointCloudBVH<Real, ParticleID>>(positions, metadata);
  a_index.pointCloudNLocal = nLocal;

  // Precompute the whole kNN graph in one batched pass, reused across every internal
  // a_iterateLocalTierToConvergence pass (positions never change within a round; only consumed
  // status does, and that is filtered at lookup time). Ghost rows are computed but never read;
  // see NNSpatialIndex::pointCloudGraph.
  a_index.pointCloudGraph = a_index.pointCloud->allNearestNeighbors(static_cast<std::size_t>(a_kNearest));
}

template <typename Packed>
inline void
findNearestNeighborCandidates(const std::vector<NNMergeParticle<Packed>>&                a_localValid,
                              const std::vector<NNMergeParticle<Packed>>&                a_ghosts,
                              const std::unordered_set<ParticleID>&                      a_consumedIDs,
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
    const NNMergeParticle<Packed>& q = a_localValid[qi];

    // Lazy deletion -- a_localValid is the FULL per-patch set (see this function's own docs), so
    // a query already consumed by an earlier pass/query this round is skipped explicitly here.
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

    // Read this query's precomputed kNN row out of the PointCloudBVH graph (see
    // buildNNSpatialIndex()). The query's own cloud index equals qi (local valid points packed
    // first); allNearestNeighbors() already excluded it as self. Each hit's cloud index decodes
    // ownership (< pointCloudNLocal => local valid, else ghost); apply the a_maxCellDistance bound
    // and the a_consumedIDs skip before feeding survivors into the shared top-K logic. Filtering
    // the distance-sorted k-nearest by the bound yields exactly the within-bound k-nearest.
    // pointCloud is null only for a patch with no crowded cell, which generates no query.
    if (a_spatialIndex.pointCloud) {
      using PCBVH = EBGeometry::PointCloudBVH<Real, ParticleID>;

      const std::size_t nLocalPC = a_spatialIndex.pointCloudNLocal;
      const int         Kgraph   = a_spatialIndex.pointCloudK;
      const std::size_t rowStart = qi * static_cast<std::size_t>(Kgraph);

      for (int h = 0; h < Kgraph; h++) {
        const typename PCBVH::Hit& hit = a_spatialIndex.pointCloudGraph[rowStart + static_cast<std::size_t>(h)];

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

template <typename Packed, typename Combine, typename PosValid, typename Allocator, typename OnMerge>
inline void
resolveTrivialTier(const std::vector<NNMergeEdge<Packed>>&                         a_edges,
                   const std::unordered_map<ParticleID, NNPooledParticle<Packed>>& a_particlesByID,
                   const Combine&                                                  a_combine,
                   std::unordered_map<NNCellKey, int, NNCellKeyHasher>&            a_liveCellCount,
                   std::unordered_set<ParticleID>&                                 a_consumedIDs,
                   const int                                                       a_numParticlesPerCellThresh,
                   std::unordered_set<ParticleID>&                                 a_hasOutgoingCommitment,
                   std::vector<NNMergeResult<Packed>>&                             a_results,
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

    const NNMergeParticle<Packed>& query        = queryDataIt->second.particle;
    const bool                     queryExposed = queryDataIt->second.exposed;

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
        const NNMergeParticle<Packed>& cand = candDataIt->second.particle;

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
        const NNMergeParticle<Packed>& cand = candDataIt->second.particle;

        // Each participant's OWN level's dx -- see a_dxByLevel's docs. The cell key must match
        // the one a_liveCellCount was populated under for that participant, so a single shared
        // dx would break the recheck for every level but one.
        const NNCellKey qKey(query.level, cellKeyOf(query.position, a_probLo, a_dxByLevel[query.level]));
        const NNCellKey cKey(cand.level, cellKeyOf(cand.position, a_probLo, a_dxByLevel[cand.level]));

        const auto qCountIt = a_liveCellCount.find(qKey);
        const auto cCountIt = a_liveCellCount.find(cKey);

        const int qCount = (qCountIt != a_liveCellCount.end()) ? qCountIt->second : 0;
        const int cCount = (cCountIt != a_liveCellCount.end()) ? cCountIt->second : 0;

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
generateProposals(const std::vector<NNMergeEdge<Packed>>&                         a_edges,
                  const std::unordered_set<ParticleID>&                           a_hasOutgoingCommitment,
                  const std::unordered_map<ParticleID, NNPooledParticle<Packed>>& a_particlesByID,
                  const int                                                       a_numRanks,
                  std::vector<std::vector<NNMergeProposal<Packed>>>&              a_proposalsByDestRank) noexcept
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
judgeProposals(const std::vector<NNMergeProposal<Packed>>&                     a_incomingProposals,
               const std::unordered_map<ParticleID, NNPooledParticle<Packed>>& a_particlesByID,
               const std::unordered_set<ParticleID>&                           a_consumedIDs,
               const std::unordered_set<ParticleID>&                           a_hasOutgoingCommitment,
               const std::map<ParticleID, ParticleID>&                         a_outgoingTargetOf,
               const Combine&                                                  a_combine,
               std::unordered_map<NNCellKey, int, NNCellKeyHasher>&            a_liveCellCount,
               const int                                                       a_numParticlesPerCellThresh,
               const RealVect&                                                 a_probLo,
               const std::vector<RealVect>&                                    a_dxByLevel,
               const PosValid&                                                 a_isPositionValid,
               const Allocator&                                                a_allocateID,
               std::vector<NNMergeResult<Packed>>&                             a_results,
               std::vector<std::vector<NNMergeVerdict>>&                       a_verdictsByDestRank,
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

    const NNMergeParticle<Packed>& target = targetIt->second.particle;
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
applyVerdicts(const std::vector<NNMergeVerdict>& a_incomingVerdicts,
              std::unordered_set<ParticleID>&    a_consumedIDs) noexcept
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
    // O(1) point->block query -- see PositionLocator/NNParticleLocation. Never a linear scan over
    // patches (the bottleneck this avoids at large patch counts).
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
                           const int                     a_ghostWidth,
                           const PosValid&               a_isPositionValid,
                           unsigned long long*           a_crossLevelMergeCount,
                           const OnMerge&                a_onMergeCommitted)
{
  CH_TIME("ParticleManagement::mergeNearestNeighborsRound");

  // This function is the public orchestrator; the step helpers and their types live in detail::.
  using namespace detail;

  const int myRank   = procID();
  const int numRanks = static_cast<int>(numProc());

  const std::string   realm       = a_particles.getRealm();
  const int           finestLevel = a_amr.getFinestLevel();
  const RealVect      probLo      = a_amr.getProbLo();
  const Vector<Real>& dxScalar    = a_amr.getDx();

  // Ghost particles are assumed ALREADY FILLED by the caller (an a_ghostWidth-cell halo -- one cell
  // per unit of a_ghostWidth in every direction, which is what the crowding trigger and the candidate
  // search look at) and are used here as candidates. This function does not fill them: filling is the
  // caller's responsibility, on equal footing with the caller cell-sorting before a cell-based merge
  // (see ItoSolver / ParticleManagement::ParticleMergeKind). A particle ghost mask of this same width
  // must have been registered on this realm before the grids were (re)built, and every rank must use
  // the same a_ghostWidth (see the a_ghostWidth docs -- it defines boundary-exposure).

  struct PatchWork
  {
    int       level;
    DataIndex din;
    RealVect  dx;
  };

  std::vector<PatchWork>                            patchWork;
  std::vector<std::vector<NNMergeParticle<Packed>>> patchLocalValid;   // parallel to patchWork
  std::vector<std::vector<NNMergeParticle<Packed>>> patchGhosts;       // parallel to patchWork
  std::vector<NNSpatialIndex<Packed>>               patchSpatialIndex; // parallel to patchWork

  // pooled over this rank's own patches; each value carries the particle plus its boundary-exposure
  // flag (folded in from what used to be a separate exposed map, halving the gather's per-particle
  // inserts and resolveTrivialTier()'s per-query/candidate lookups)
  std::unordered_map<ParticleID, NNPooledParticle<Packed>> particlesByID;
  std::unordered_map<NNCellKey, int, NNCellKeyHasher>      liveCellCount; // pooled, level-qualified cell key

  // Per-level dx, indexed by level -- every cell-key computation below and downstream (in
  // resolveTrivialTier()/judgeProposals()) MUST use the dx of the PARTICLE'S OWN level (see
  // NNMergeParticle::level's docs), never a single dx shared across levels.
  std::vector<RealVect> dxByLevel(finestLevel + 1);

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    dxByLevel[lvl] = RealVect(D_DECL(dxScalar[lvl], dxScalar[lvl], dxScalar[lvl]));
  }

  CH_assert(a_ghostWidth >= 1);
  const int ghostWidth = a_ghostWidth;

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    const RealVect& dx = dxByLevel[lvl];

    const DisjointBoxLayout& dbl = a_amr.getGrids(realm)[lvl];

    const AMRParticleGhostMask& maskSame = a_amr.getParticleGhostMask(realm, ghostWidth);
    const AMRParticleGhostMask& maskC2F  = a_amr.getParticleGhostMaskCoarToFine(realm, ghostWidth);
    const AMRParticleGhostMask& maskF2C  = a_amr.getParticleGhostMaskFineToCoar(realm, ghostWidth);

    const DataIterator& dit  = dbl.dataIterator();
    const int           nbox = dit.size();

    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      std::vector<NNMergeParticle<Packed>> localValid;
      std::vector<NNMergeParticle<Packed>> ghosts;

      localValid.reserve(leaf.size());
      ghosts.reserve(leaf.size());

      // Reserve to the running total -- particlesByID is pooled across every patch this rank owns,
      // so reserving up front avoids repeated rehashes as it grows. leaf.size() is O(1).
      particlesByID.reserve(particlesByID.size() + leaf.size());

      // Boundary-exposure masks for THIS patch (din), resolved ONCE here rather than per particle:
      // din/lvl are fixed within this box, so the RefCountedPtr null-checks and the LayoutData[din]
      // lookups are loop-invariant. A null pointer means the mask does not apply at this level (no
      // coarser/finer neighbour, or none defined).
      const ParticleGhostMask* maskSameDin = maskSame[lvl].isNull() ? nullptr : &(*maskSame[lvl])[din];
      const ParticleGhostMask* maskC2FDin  = (lvl < finestLevel && !maskC2F[lvl].isNull()) ? &(*maskC2F[lvl])[din]
                                                                                           : nullptr;
      const ParticleGhostMask* maskF2CDin  = (lvl > 0 && !maskF2C[lvl].isNull()) ? &(*maskF2C[lvl])[din] : nullptr;

      for (std::size_t i = 0; i < leaf.size(); i++) {
        const bool isGhost = leaf.isGhost(i);

        // A ghost's TRUE origin level is NOT this receiving patch's level `lvl` -- it is
        // lvl-1/lvl+1 for a Coarse/Fine ghost (see ParticleSoA::GhostType). Its cell key must be
        // computed with that origin level's dx (see NNMergeParticle::level).
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
          // REQUIRED -- generateProposals()'s a_particlesByID contract is "every currently-alive
          // local valid particle AND every ghost". A ghost whose true owning patch is not among
          // this rank's own patches has no native entry anywhere in this loop, so without this
          // insert generateProposals()'s lookup fails and the proposal is dropped. Ghosts are never
          // exposed (a ghost never initiates a query).
          particlesByID[p.globalID] = NNPooledParticle<Packed>{p, false};
        }
        else {
          bool isExposed = (maskSameDin != nullptr) && isBoundaryExposed(cell, *maskSameDin);

          if (!isExposed && maskC2FDin != nullptr) {
            isExposed = isBoundaryExposed(cell, *maskC2FDin);
          }

          if (!isExposed && maskF2CDin != nullptr) {
            isExposed = isBoundaryExposed(cell, *maskF2CDin);
          }

          localValid.push_back(p);
          particlesByID[p.globalID] = NNPooledParticle<Packed>{p, isExposed};
        }
      }

      // Built once here, per patch, per round -- over the FULL localValid/ghosts sets, before any
      // consumption this round -- and reused across every internal a_iterateLocalTierToConvergence
      // pass via lazy deletion (see findNearestNeighborCandidates()'s docs).
      NNSpatialIndex<Packed> spatialIndex;

      buildNNSpatialIndex(localValid,
                          ghosts,
                          probLo,
                          dx,
                          a_numParticlesPerCellThresh,
                          1 + a_maxFallbackCandidates,
                          spatialIndex);

      patchWork.push_back(PatchWork{lvl, din, dx});
      patchLocalValid.push_back(std::move(localValid));
      patchGhosts.push_back(std::move(ghosts));
      patchSpatialIndex.push_back(std::move(spatialIndex));
    }
  }

  // Populate liveCellCount from particlesByID, NOT inline in the gather loop above. A single
  // particle can appear as a ghost in several of this rank's own patches at once, so counting per
  // occurrence would multiply-count it. particlesByID is keyed by globalID, so those repeated
  // occurrences collapse to one entry -- deriving the count from it counts each particle once.
  for (const auto& kv : particlesByID) {
    const NNMergeParticle<Packed>& p    = kv.second.particle;
    const IntVect                  cell = cellKeyOf(p.position, probLo, dxByLevel[p.level]);

    liveCellCount[NNCellKey(p.level, cell)] += 1;
  }

  std::unordered_set<ParticleID>   consumedIDs;
  std::unordered_set<ParticleID>   hasOutgoingCommitment;
  std::map<ParticleID, ParticleID> outgoingTargetOf;

  std::vector<NNMergeResult<Packed>> allResults;

  // 2. Per-patch candidate search + trivial-tier resolution, pooled rank-wide.
  //
  // lastPooledEdges deliberately survives past the loop below: it is reused (NOT recomputed) for
  // outgoingTargetOf/generateProposals afterward. Nothing consumes any particles between the last
  // runLocalTierPass() and generateProposals(), so the edges are still exactly current; re-running
  // the candidate search purely to re-derive them would repeat the most expensive step for no new
  // information.
  std::vector<NNMergeEdge<Packed>> lastPooledEdges;

  auto runLocalTierPass = [&]() -> bool {
    lastPooledEdges.clear();

    for (std::size_t pw = 0; pw < patchWork.size(); pw++) {
      // patchLocalValid[pw]/patchGhosts[pw] are this patch's FULL sets and patchSpatialIndex[pw]
      // was built once over them before the local-tier passes began. findNearestNeighborCandidates()
      // filters against consumedIDs itself (see its docs), so nothing is rebuilt or re-filtered here.
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

  // outgoingTargetOf: for every query id with an outgoing commitment, remember its own nearest-
  // neighbor candidate id (needed by judgeProposals()'s mutual-match rule). Reuses lastPooledEdges
  // (see above) -- correct without a fresh search since nothing is consumed between the last pass
  // and here.
  for (const auto& edge : lastPooledEdges) {
    if (hasOutgoingCommitment.count(edge.queryID) > 0) {
      outgoingTargetOf[edge.queryID] = edge.candidateID;
    }
  }

  // 3. generateProposals() for everything the trivial tier could not resolve.
  std::vector<std::vector<NNMergeProposal<Packed>>> proposalsByDestRank;
  generateProposals(lastPooledEdges, hasOutgoingCommitment, particlesByID, numRanks, proposalsByDestRank);

  // 4. Exchange proposals.
  const std::vector<NNMergeProposal<Packed>> incomingProposals = nnMergeExchangeByRank(proposalsByDestRank);

  // 5. judgeProposals() over the pooled incoming proposals. Judge targets are looked up directly
  // in particlesByID and filtered against consumedIDs at the point of use inside judgeProposals(),
  // rather than from a separate pre-filtered copy.
  std::vector<NNMergeResult<Packed>>       judgeResults;
  std::vector<std::vector<NNMergeVerdict>> verdictsByDestRank(numRanks);

  judgeProposals(incomingProposals,
                 particlesByID,
                 consumedIDs,
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

  allResults.insert(allResults.end(), judgeResults.begin(), judgeResults.end());

  for (const auto& r : judgeResults) {
    consumedIDs.insert(r.firstID);
  }

  // 6. Exchange verdicts and applyVerdicts().
  const std::vector<NNMergeVerdict> incomingVerdicts = nnMergeExchangeByRank(verdictsByDestRank);

  applyVerdicts(incomingVerdicts, consumedIDs);

  // 7. placeMergedParticles(), then exchange and apply the scatter bucket.
  //
  // O(1) point->block query (see NNParticleLocation/PositionLocator). ParticleContainer::
  // findDestination() is the same tile-hash lookup a_particles uses to route its own particles
  // during remap(); reusing it keeps one source of truth for "which box owns this position"
  // instead of a linear scan or a parallel implementation that could drift.
  const PositionLocator locate = [&a_particles](const RealVect& a_pos) -> NNParticleLocation {
    const auto dst = a_particles.findDestination(a_pos);

    return NNParticleLocation{dst.found, dst.level, dst.gridIndex, dst.rank};
  };

  std::map<std::pair<int, unsigned int>, std::vector<NNMergeParticle<Packed>>> localByPatch;
  std::vector<std::vector<NNMergeParticle<Packed>>>                            scatterByDestRank(numRanks);

  placeMergedParticles(allResults, locate, myRank, localByPatch, scatterByDestRank);

  const std::vector<NNMergeParticle<Packed>> incomingScattered = nnMergeExchangeByRank(scatterByDestRank);

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
      // Correct the provisional ownerRank to this rank -- placeMergedParticles() has resolved the
      // TRUE owner via findDestination() and it is this rank. The provisional value (inherited
      // from a parent, see resolveTrivialTier()/judgeProposals()) must not survive: if this
      // particle later proposes and is accepted, judgeProposals() routes the verdict to
      // NNMergeParticle::ownerRank, so a stale rank would send the removal to the wrong rank and
      // the particle would survive while its weight is also merged elsewhere -- created mass.
      NNMergeParticle<Packed> corrected = merged;

      corrected.ownerRank = myRank;
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
}

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#endif
