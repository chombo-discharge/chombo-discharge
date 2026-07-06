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
#include <cstdlib>

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>
#include <SPMD.H>
#include <DataIterator.H>

// Our includes
#include <CD_NearestNeighborParticleMerge.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace ParticleManagement {

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

  template <typename Packed>
  inline void
  findNearestNeighborCandidates(const std::vector<NNMergeParticle<Packed>>& a_localValid,
                                const std::vector<NNMergeParticle<Packed>>& a_ghosts,
                                const std::map<IntVect, int>&               a_liveCellCount,
                                const RealVect&                             a_probLo,
                                const RealVect&                             a_dx,
                                const int                                   a_numParticlesPerCellThresh,
                                std::vector<NNMergeEdge<Packed>>&           a_edges,
                                const int                                   a_maxFallbackCandidates) noexcept
  {
    CH_TIME("ParticleManagement::findNearestNeighborCandidates");

    a_edges.clear();

    const std::size_t nLocal = a_localValid.size();
    const std::size_t nGhost = a_ghosts.size();
    const int         K      = 1 + a_maxFallbackCandidates;

    struct RankedCandidate
    {
      Real        d2;
      std::size_t idx;
      bool        isLocal;
    };

    for (std::size_t qi = 0; qi < nLocal; qi++) {
      const NNMergeParticle<Packed>& q = a_localValid[qi];

      const IntVect qCell = cellKeyOf(q.position, a_probLo, a_dx);
      const auto    cntIt = a_liveCellCount.find(qCell);
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

      for (std::size_t li = 0; li < nLocal; li++) {
        if (li == qi) {
          continue;
        }
        considerCandidate(nnMergeDistSquared(q.position, a_localValid[li].position), li, true);
      }
      for (std::size_t gi = 0; gi < nGhost; gi++) {
        considerCandidate(nnMergeDistSquared(q.position, a_ghosts[gi].position), gi, false);
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

  template <typename Packed, typename Combine, typename PosValid, typename Allocator>
  inline void
  resolveTrivialTier(const std::vector<NNMergeEdge<Packed>>&              a_edges,
                     const std::map<ParticleID, bool>&                    a_exposed,
                     const std::map<ParticleID, NNMergeParticle<Packed>>& a_particlesByID,
                     const Combine&                                       a_combine,
                     std::map<IntVect, int>&                              a_liveCellCount,
                     std::set<ParticleID>&                                a_consumedIDs,
                     const int                                            a_numParticlesPerCellThresh,
                     std::set<ParticleID>&                                a_hasOutgoingCommitment,
                     std::vector<NNMergeResult<Packed>>&                  a_results,
                     const RealVect&                                      a_probLo,
                     const RealVect&                                      a_dx,
                     const PosValid&                                      a_isPositionValid,
                     const Allocator&                                     a_allocateID,
                     const int                                            a_maxFallbackCandidates,
                     const std::optional<int>                             a_maxCellDistance) noexcept
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
            const IntVect qCell = cellKeyOf(query.position, a_probLo, a_dx);
            const IntVect cCell = cellKeyOf(cand.position, a_probLo, a_dx);
            tooFar              = nnMergeChebyshevCellDistance(qCell, cCell) > a_maxCellDistance.value();
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

          const IntVect qCell = cellKeyOf(query.position, a_probLo, a_dx);
          const IntVect cCell = cellKeyOf(cand.position, a_probLo, a_dx);

          const auto qCountIt = a_liveCellCount.find(qCell);
          const auto cCountIt = a_liveCellCount.find(cCell);
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
            result.merged.payload   = a_combine(query.payload, query.weight, cand.payload, cand.weight);

            a_results.push_back(result);

            a_consumedIDs.insert(query.globalID);
            a_consumedIDs.insert(cand.globalID);

            a_liveCellCount[qCell] -= 1;
            a_liveCellCount[cCell] -= 1;

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

  template <typename Packed, typename Combine, typename PosValid, typename Allocator>
  inline void
  judgeProposals(const std::vector<NNMergeProposal<Packed>>&          a_incomingProposals,
                 const std::map<ParticleID, NNMergeParticle<Packed>>& a_myParticlesByID,
                 const std::set<ParticleID>&                          a_hasOutgoingCommitment,
                 const std::map<ParticleID, ParticleID>&              a_outgoingTargetOf,
                 const Combine&                                       a_combine,
                 std::map<IntVect, int>&                              a_liveCellCount,
                 const int                                            a_numParticlesPerCellThresh,
                 const RealVect&                                      a_probLo,
                 const RealVect&                                      a_dx,
                 const PosValid&                                      a_isPositionValid,
                 const Allocator&                                     a_allocateID,
                 std::vector<NNMergeResult<Packed>>&                  a_results,
                 std::vector<std::vector<NNMergeVerdict>>&            a_verdictsByDestRank,
                 const std::optional<int>                             a_maxCellDistance) noexcept
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
        const IntVect tCell = cellKeyOf(target.position, a_probLo, a_dx);
        const IntVect sCell = cellKeyOf(source.position, a_probLo, a_dx);
        tooFar              = nnMergeChebyshevCellDistance(tCell, sCell) > a_maxCellDistance.value();
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

      // Dynamic per-target threshold recheck (mirrors resolveTrivialTier()'s discipline).
      const IntVect tCell    = cellKeyOf(target.position, a_probLo, a_dx);
      const auto    tCountIt = a_liveCellCount.find(tCell);
      const int     tCount   = (tCountIt != a_liveCellCount.end()) ? tCountIt->second : 0;

      if (tCount <= a_numParticlesPerCellThresh) {
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
      result.merged.payload   = a_combine(target.payload, target.weight, source.payload, source.weight);

      a_results.push_back(result);

      a_liveCellCount[tCell] -= 1;

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
            typename PosValid>
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
                             const PosValid&               a_isPositionValid)
  {
    CH_TIME("ParticleManagement::mergeNearestNeighborsRound");

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

    a_particles.clearGhostParticles();
    a_particles.fillGhostParticles(a_amr.getParticleGhostMask(realm, ghostWidth),
                                   a_amr.getParticleGhostMaskCoarToFine(realm, ghostWidth),
                                   a_amr.getParticleGhostMaskFineToCoar(realm, ghostWidth));

    struct PatchWork
    {
      int       level;
      DataIndex din;
      RealVect  dx;
    };

    std::vector<PatchWork>                            patchWork;
    std::vector<std::vector<NNMergeParticle<Packed>>> patchLocalValid; // parallel to patchWork
    std::vector<std::vector<NNMergeParticle<Packed>>> patchGhosts;     // parallel to patchWork

    std::map<ParticleID, NNMergeParticle<Packed>> particlesByID; // pooled, this rank's own patches
    std::map<ParticleID, bool>                    exposed;       // pooled, this rank's own valid particles
    std::map<IntVect, int>                        liveCellCount; // pooled, GLOBAL cell key

    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      const RealVect dx(D_DECL(dxScalar[lvl], dxScalar[lvl], dxScalar[lvl]));

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

        for (std::size_t i = 0; i < leaf.size(); i++) {
          NNMergeParticle<Packed> p;
          p.position  = leaf.position(i);
          p.weight    = leaf.weight(i);
          p.globalID  = leaf.particleID(i);
          p.ownerRank = leaf.rankID(i);
          p.payload   = a_gather(leaf, i);

          const IntVect cell = cellKeyOf(p.position, probLo, dx);

          if (leaf.isGhost(i)) {
            ghosts.push_back(p);
            liveCellCount[cell] += 1;
          }
          else {
            localValid.push_back(p);
            liveCellCount[cell] += 1;
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

        patchWork.push_back(PatchWork{lvl, din, dx});
        patchLocalValid.push_back(std::move(localValid));
        patchGhosts.push_back(std::move(ghosts));
      }
    }

    const RealVect dxLevel0(D_DECL(dxScalar[0], dxScalar[0], dxScalar[0]));

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
        std::vector<NNMergeParticle<Packed>> aliveLocal;
        aliveLocal.reserve(patchLocalValid[pw].size());
        for (const auto& p : patchLocalValid[pw]) {
          if (consumedIDs.count(p.globalID) == 0) {
            aliveLocal.push_back(p);
          }
        }
        std::vector<NNMergeParticle<Packed>> aliveGhosts;
        aliveGhosts.reserve(patchGhosts[pw].size());
        for (const auto& p : patchGhosts[pw]) {
          if (consumedIDs.count(p.globalID) == 0) {
            aliveGhosts.push_back(p);
          }
        }

        std::vector<NNMergeEdge<Packed>> patchEdges;
        findNearestNeighborCandidates(aliveLocal,
                                      aliveGhosts,
                                      liveCellCount,
                                      probLo,
                                      patchWork[pw].dx,
                                      a_numParticlesPerCellThresh,
                                      patchEdges,
                                      a_maxFallbackCandidates);

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
                         dxLevel0,
                         a_isPositionValid,
                         a_allocateID,
                         a_maxFallbackCandidates,
                         a_maxCellDistance);

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
    std::vector<std::vector<NNMergeProposal<Packed>>> proposalsByDestRank;
    generateProposals(lastPooledEdges, hasOutgoingCommitment, particlesByID, numRanks, proposalsByDestRank);

    // 4. Exchange proposals.
    const std::vector<NNMergeProposal<Packed>> incomingProposals = nnMergeExchangeByRank(proposalsByDestRank);

    // 5. judgeProposals() over the pooled incoming proposals.
    std::map<ParticleID, NNMergeParticle<Packed>> myParticlesByID;
    for (const auto& kv : particlesByID) {
      if (consumedIDs.count(kv.first) == 0) {
        myParticlesByID[kv.first] = kv.second;
      }
    }

    std::vector<NNMergeResult<Packed>>       judgeResults;
    std::vector<std::vector<NNMergeVerdict>> verdictsByDestRank(numRanks);

    judgeProposals(incomingProposals,
                   myParticlesByID,
                   hasOutgoingCommitment,
                   outgoingTargetOf,
                   a_combine,
                   liveCellCount,
                   a_numParticlesPerCellThresh,
                   probLo,
                   dxLevel0,
                   a_isPositionValid,
                   a_allocateID,
                   judgeResults,
                   verdictsByDestRank,
                   a_maxCellDistance);

    allResults.insert(allResults.end(), judgeResults.begin(), judgeResults.end());
    for (const auto& r : judgeResults) {
      consumedIDs.insert(r.firstID);
    }

    // 6. Exchange verdicts and applyVerdicts().
    const std::vector<NNMergeVerdict> incomingVerdicts = nnMergeExchangeByRank(verdictsByDestRank);
    applyVerdicts(incomingVerdicts, consumedIDs);

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
        a_scatter(leaf, merged);
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

      ParticleSoA<P, Traits>& leaf = a_particles[loc.level][din];
      a_scatter(leaf, merged);
    }
  }

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#endif
