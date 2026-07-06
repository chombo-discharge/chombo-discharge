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

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_NearestNeighborParticleMerge.H>
#include <CD_NamespaceHeader.H>

namespace ParticleManagement {

  inline bool
  isBoundaryExposed(const IntVect& a_cellIndex, const ParticleGhostMask& a_ghostMask) noexcept
  {
    // A particle is boundary-exposed iff its own cell has at least one outgoing ghost-shipment
    // target registered on the mask that drives THIS patch's own ghost-fill. See RULES.md,
    // "Boundary-exposed classification" -- this is deliberately a lookup, not a geometric
    // recomputation:
    //
    //   return a_ghostMask.numTargets(a_cellIndex) > 0;
    //
    // Note for mergeNearestNeighborsRound(): a particle can be exposed toward a same-level
    // neighbor, a coarser level, or a finer level -- three DIFFERENT masks
    // (AmrMesh::getParticleGhostMask()/getParticleGhostMaskCoarToFine()/
    // getParticleGhostMaskFineToCoar()). This function only tests ONE mask; the orchestrator is
    // responsible for OR-ing all three together per particle (exposed iff ANY of the three say
    // so), not this function.
    MayDay::Error("ParticleManagement::isBoundaryExposed -- not yet implemented, see RULES.md");
    return false;
  }

  template <typename Packed>
  inline void
  findNearestNeighborCandidates(const std::vector<NNMergeParticle<Packed>>& a_localValid,
                                const std::vector<NNMergeParticle<Packed>>& a_ghosts,
                                const std::map<IntVect, int>&               a_liveCellCount,
                                const int                                   a_numParticlesPerCellThresh,
                                std::vector<NNMergeEdge<Packed>>&           a_edges) noexcept
  {
    CH_TIME("ParticleManagement::findNearestNeighborCandidates");

    a_edges.clear();

    // 1. Build a single spatial index (e.g. a flat, scratch-buffer KD-tree in the style of
    //    ParticleManagement::buildEqualWeightKDLeaves in CD_ParticleManagementImplem.H -- no
    //    per-node heap allocation, reused thread_local scratch across calls) over the COMBINED
    //    position set: every particle in a_localValid, followed by every particle in a_ghosts.
    //    Remember which half of the combined index each entry came from (local vs ghost) and its
    //    original NNMergeParticle, so results can be reported back by id.
    //
    // 2. For each particle p in a_localValid (query sources -- ghosts never query):
    //      a. cellKey = cellKeyOf(p.position)   -- unclamped, position-derived (see RULES.md,
    //         "Crowding trigger": must NOT be clamped to any one patch's box, so a ghost bucket
    //         into its true origin cell).
    //      b. if a_liveCellCount.at(cellKey) <= a_numParticlesPerCellThresh: skip -- p's own
    //         cell is not (yet) crowded enough to be a candidate query source.
    //      c. Query the spatial index for p's single nearest OTHER point (excluding itself),
    //         among the combined local+ghost set. This is the ONE nearest-neighbor query per
    //         particle this whole algorithm performs -- see the file-level docs for why it is
    //         always run against the FULL neighborhood, never a restricted one.
    //      d. If a candidate was found, append one NNMergeEdge: distanceSquared, p's id, the
    //         candidate's id, and whether the candidate came from the local half of the index
    //         (NNMergeEdge::candidateIsLocal).
    //
    // 3. Do NOT filter on the CANDIDATE's own cell threshold here -- that is deliberately
    //    deferred to resolveTrivialTier()/judgeProposals(), which re-check it dynamically at
    //    processing time against live counts that may already reflect earlier merges this round
    //    (see RULES.md, "Crowding trigger").
    MayDay::Error("ParticleManagement::findNearestNeighborCandidates -- not yet implemented, see RULES.md");
  }

  template <typename Packed, typename Combine>
  inline void
  resolveTrivialTier(const std::vector<NNMergeEdge<Packed>>& a_edges,
                     const std::map<ParticleID, bool>&       a_exposed,
                     const Combine&                          a_combine,
                     std::map<IntVect, int>&                 a_liveCellCount,
                     std::set<ParticleID>&                   a_consumedIDs,
                     std::set<ParticleID>&                   a_hasOutgoingCommitment,
                     std::vector<NNMergeResult<Packed>>&     a_results,
                     const int                               a_maxFallbackCandidates,
                     const std::optional<int>                a_maxCellDistance) noexcept
  {
    CH_TIME("ParticleManagement::resolveTrivialTier");

    a_results.clear();

    // 1. Sort a COPY of a_edges by (distanceSquared, queryID, candidateID) -- the deterministic
    //    tiebreak on ids is required, not optional: an unbroken exact-distance tie (e.g. two
    //    true mutual nearest neighbors) left to whatever order the edges happened to be
    //    generated/inserted in is a real, observed source of nondeterminism that cascades into
    //    unrelated pairs several hops away (see RULES.md, "Trivial tier").
    //
    // 2. For each edge (queryID, candidateID, d2, candidateIsLocal), in that sorted order:
    //      a. if queryID is already in a_consumedIDs: stale (already consumed by an earlier,
    //         closer edge in this same pass) -- skip.
    //      b. Evaluate the CURRENT candidate (initially candidateID, possibly replaced by a
    //         fallback below) with the following per-candidate check, "evalCandidate":
    //           candidateIsStale = candidateIsLocal && candidate already in a_consumedIDs.
    //           candidateIsBusy  = candidateIsLocal && !candidateIsStale &&
    //                              candidate is in a_hasOutgoingCommitment.
    //           candidateTooFar  = a_maxCellDistance.has_value() &&
    //                              cellChebyshevDistance(cellKeyOf(query.position),
    //                                                     cellKeyOf(candidate.position))
    //                                > *a_maxCellDistance
    //             -- cellKeyOf() here is the SAME unclamped, position-derived cell key used for
    //             the crowding trigger (floor(pos/dx)), NOT either particle's own patch-local
    //             IntVect cell index -- using a shared, position-derived key is what makes this
    //             comparison meaningful across a patch boundary at all. cellChebyshevDistance is
    //             the max-per-axis difference between the two integer cell-index vectors.
    //           trivialOk = candidateIsLocal && !candidateIsStale && !candidateIsBusy &&
    //                       !candidateTooFar &&
    //                       !a_exposed.at(queryID) && !a_exposed.at(candidate).
    //         (a ghost candidate, i.e. !candidateIsLocal, is never trivial-tier eligible --
    //         trivialOk stays false; this is unaffected by any of the above.)
    //      c. If trivialOk: dynamic per-cell threshold recheck, using a_liveCellCount AS IT
    //         CURRENTLY STANDS (i.e. already reflecting decrements from edges processed earlier
    //         in this SAME sorted pass, not the generation-time snapshot) -- both queryID's and
    //         the candidate's current cell counts must still exceed the threshold.
    //      d. if trivialOk and the threshold recheck passed: COMMIT --
    //           - merged = a_combine(query, candidate) folded with the position/weight
    //             weighted-centroid formula (position/weight are handled directly by this
    //             function, not by the Combine callback -- see PayloadCombine's docs).
    //           - append {queryID, candidateID, merged} to a_results.
    //           - insert queryID and candidateID into a_consumedIDs.
    //           - decrement a_liveCellCount at both of their cells.
    //           - done with this edge, move to the next.
    //      e. else if the threshold recheck failed specifically (trivialOk was true but crowding
    //         no longer holds): queryID is simply left unmatched this round -- neither a trivial
    //         merge nor a proposal; NOT a fallback-eligible case (a different candidate wouldn't
    //         change the query's OWN cell occupancy). See RULES.md, "Crowding trigger" for why
    //         this is a legitimate, expected outcome, and "What's proven vs. what's empirically
    //         validated" for the caveat that a persistently isolated crowded cell is possible and
    //         correct when nothing moves between rounds. Done with this edge.
    //      f. else if candidateIsStale OR candidateTooFar (evaluated against the ORIGINAL,
    //         un-retried candidate -- see below): do NOT insert queryID into
    //         a_hasOutgoingCommitment, do NOT generate a proposal for it. For candidateIsStale,
    //         there is no meaningful id to propose to. For candidateTooFar, the pair is exactly
    //         as physically invalid via judgeProposals() as it would be trivially -- there is
    //         nothing valid to hand off to a judge either, so a proposal would be meaningless.
    //         Either way, marking queryID "busy" here would incorrectly block some OTHER
    //         particle from using it as an ordinary trivial candidate later, even though queryID
    //         hasn't committed to anything. If a_maxFallbackCandidates > 0 and queryID is not
    //         itself exposed (if it is, no local candidate could ever be trivially eligible, so
    //         don't bother): immediately re-query for the next-nearest still-available candidate,
    //         excluding every candidate already tried for this edge, and go back to step (b) with
    //         that new candidate (consuming one of a_maxFallbackCandidates attempts). If no
    //         further candidate is found, or attempts are exhausted: queryID is left unmatched
    //         this round, same as (e) -- no proposal. Note that whether a PROPOSAL is generated at
    //         the end depends only on the ORIGINAL candidate (the one from
    //         findNearestNeighborCandidates(), before any fallback retries) being stale or too
    //         far -- a fallback candidate that is itself too far only ends that ONE retry attempt
    //         (try the next fallback instead), it does not by itself suppress the eventual
    //         proposal if the original candidate was otherwise propose-eligible.
    //      g. else (candidateIsBusy, or candidate is exposed, or candidate is a ghost -- i.e.
    //         trivialOk was false for a reason OTHER than staleness or distance): if
    //         a_maxFallbackCandidates > 0 and queryID is not itself exposed: same fallback retry
    //         as (f) -- try the next candidate instead of giving up immediately. Once fallback
    //         attempts are exhausted (or a_maxFallbackCandidates is 0): insert queryID into
    //         a_hasOutgoingCommitment, so generateProposals() picks it up -- using queryID's
    //         ORIGINAL edge (the one from findNearestNeighborCandidates(), NOT any fallback
    //         candidate that was also tried and rejected) as the proposal target. The candidate
    //         genuinely still exists (and is within range) in every one of these cases, so
    //         judgeProposals() -- which already has a rule for "target has its own outgoing
    //         commitment" (accept only on a mutual match) -- should get a chance to decide, not
    //         be silently bypassed.
    //
    // Empirical note (Python prototype, Scripts/ParticleMergePrototype/): instrumenting exactly
    // why edges failed to commit in a deliberately extreme single-round test (a cell seeded at
    // 64 particles against a threshold of 16) found "candidate already busy" (case g) to be the
    // SINGLE LARGEST loss category -- larger than genuine boundary exposure. An earlier revision
    // of this design treated case (g)'s busy sub-case the same as case (f) (silent skip, no
    // proposal) -- this was a real completeness gap relative to what's documented above, found
    // and fixed empirically. Enabling a_maxFallbackCandidates (even just 1) recovered the large
    // majority of the resulting loss; see RULES.md for the full writeup.
    //
    // a_maxCellDistance verified directly with a hand-placed pair 2 cells apart (Chebyshev):
    // uncapped, they merge trivially at the exact weighted centroid; capped at 1 cell, both
    // survive completely untouched (no merge, no proposal, no side effects); capped at exactly 2
    // cells, they merge again -- confirming the cap is inclusive (<=, not <).
    MayDay::Error("ParticleManagement::resolveTrivialTier -- not yet implemented, see RULES.md");
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

    // For each edge whose queryID is in a_hasOutgoingCommitment (i.e. resolveTrivialTier()
    // could not resolve it locally):
    //   - look up the full NNMergeParticle for queryID in a_particlesByID -- this becomes
    //     NNMergeProposal::source.
    //   - look up the candidate's owning rank via a_particlesByID.at(candidateID).ownerRank --
    //     THIS, not any id comparison, determines the destination bucket (see NNMergeProposal's
    //     docs for why an id-based routing rule was tried and found unsound: it can send the two
    //     directions of a genuinely mutual pair to two different, uncoordinated judges).
    //   - append {source, candidateID, distanceSquared} to
    //     a_proposalsByDestRank[candidate's owning rank] -- including when that rank equals this
    //     rank's own (see this function's header docs: the caller hands that specific bucket to
    //     judgeProposals() directly, in-process, rather than shipping it over MPI).
    //
    // Each queryID appears in at most ONE outgoing proposal (exactly the single edge
    // findNearestNeighborCandidates() produced for it) -- this function does not aggregate or
    // deduplicate across queries, since each query already has at most one edge by construction.
    MayDay::Error("ParticleManagement::generateProposals -- not yet implemented, see RULES.md");
  }

  template <typename Packed, typename Combine>
  inline void
  judgeProposals(const std::vector<NNMergeProposal<Packed>>&          a_incomingProposals,
                 const std::map<ParticleID, NNMergeParticle<Packed>>& a_myParticlesByID,
                 const std::set<ParticleID>&                          a_hasOutgoingCommitment,
                 const std::map<ParticleID, ParticleID>&              a_outgoingTargetOf,
                 const Combine&                                       a_combine,
                 std::map<IntVect, int>&                              a_liveCellCount,
                 const int                                            a_numParticlesPerCellThresh,
                 std::vector<NNMergeResult<Packed>>&                  a_results,
                 std::vector<std::vector<NNMergeVerdict>>&            a_verdictsByDestRank) noexcept
  {
    CH_TIMERS("ParticleManagement::judgeProposals");
    CH_TIMER("ParticleManagement::judgeProposals::phase1_candidates", t1);
    CH_TIMER("ParticleManagement::judgeProposals::phase2_sortedAccept", t2);
    CH_TIMER("ParticleManagement::judgeProposals::phase3_verdicts", t3);

    a_results.clear();

    CH_START(t1);
    // Phase 1: group a_incomingProposals by targetID, then for EACH target determine its
    // candidate winner and sort key WITHOUT committing anything yet:
    //   - if targetID is not in a_myParticlesByID (shouldn't normally happen -- a proposal was
    //     addressed here by construction -- but guard it anyway): every incoming proposal for it
    //     is rejected immediately (append {source.id=proposer, accepted=false} to
    //     a_verdictsByDestRank[proposer's rank] for each), do not add a phase-2 candidate.
    //   - else if targetID is in a_hasOutgoingCommitment (this target is itself proposing to
    //     someone else this round):
    //       myTarget = a_outgoingTargetOf.at(targetID)
    //       mutual = the subset of this target's incoming proposals whose source id == myTarget
    //       if mutual is non-empty AND targetID < myTarget (deterministic lower-id-performs-the-
    //         accept tiebreak -- see this function's header docs for why the mutual-match
    //         exception exists at all, and why exactly one side must act on it):
    //           record a phase-2 candidate (distanceSquared = mutual[0]'s, targetID, winner =
    //           mutual[0], target's data).
    //       else:
    //           reject every incoming proposal for this target immediately (same as the
    //           not-in-a_myParticlesByID branch above) -- either it is genuinely spoken for
    //           elsewhere with no reciprocity, or the reciprocal match will be accepted by the
    //           OTHER side (myTarget's own rank), which independently reaches the same
    //           conclusion from the mirrored data it has.
    //   - else (no outgoing commitment): record a phase-2 candidate using the plain
    //     minimum-(distanceSquared, source.id) proposal among this target's incoming proposals.
    CH_STOP(t1);

    CH_START(t2);
    // Phase 2: sort ALL phase-1 candidates (mutual-match and plain alike, across every target
    // this rank is judging) by (distanceSquared, targetID) and process them in that order:
    //   - re-check targetID has not already been consumed (should not happen within this single
    //     phase-2 pass since each target appears at most once, but mirrors the discipline in
    //     resolveTrivialTier() for symmetry/robustness).
    //   - dynamic per-cell threshold recheck on targetID's CURRENT cell count in a_liveCellCount
    //     (already reflecting decrements from candidates accepted earlier in this same sorted
    //     pass -- this is the fix for the bug where accepting one target could silently push a
    //     DIFFERENT target sharing its cell below threshold, and the outcome depended on
    //     whatever arbitrary order targets were iterated in).
    //   - if the recheck passes: commit -- combine target + winner's source into a merged
    //     particle (weighted centroid/weight directly, a_combine for the opaque payload), append
    //     {targetID, winner.source.id, merged} to a_results, decrement a_liveCellCount at
    //     targetID's cell, and remember that this target's winning proposal is `winner` (for
    //     phase 3).
    //   - if the recheck fails: this target is not accepted after all -- every one of its
    //     incoming proposals gets a rejected verdict in phase 3, including the one that looked
    //     like the winner in phase 1.
    CH_STOP(t2);

    CH_START(t3);
    // Phase 3: for every target that had at least one phase-1 candidate (i.e. was NOT rejected
    // immediately in phase 1 -- those already got their verdicts emitted directly), emit exactly
    // one verdict per ORIGINAL incoming proposal for that target: accepted=true for the specific
    // proposal recorded as the winner in phase 2 (only if that target was actually committed),
    // accepted=false for every other proposal against the same target. Route each verdict to
    // a_verdictsByDestRank[that proposal's source.ownerRank].
    CH_STOP(t3);

    MayDay::Error("ParticleManagement::judgeProposals -- not yet implemented, see RULES.md");
  }

  inline void
  applyVerdicts(const std::vector<NNMergeVerdict>& a_incomingVerdicts, std::set<ParticleID>& a_consumedIDs) noexcept
  {
    CH_TIME("ParticleManagement::applyVerdicts");

    // For each verdict in a_incomingVerdicts:
    //   if verdict.accepted: insert verdict.proposerID into a_consumedIDs.
    //   else: no-op -- the proposer's particle is simply left alone (see RULES.md, "Single pass,
    //   caller loops for convergence": reconsidered next round's fresh candidate search, never
    //   retried against a second-best candidate within this round).
    //
    // The caller (mergeNearestNeighborsRound()) is responsible for actually removing every id in
    // a_consumedIDs from its own valid particle lists after this call returns.
    MayDay::Error("ParticleManagement::applyVerdicts -- not yet implemented, see RULES.md");
  }

  template <typename Packed>
  inline void
  placeMergedParticles(const std::vector<NNMergeResult<Packed>>&            a_results,
                       const std::vector<Box>&                              a_myPatchBoxes,
                       std::map<int, std::vector<NNMergeParticle<Packed>>>& a_localByPatch,
                       std::vector<std::vector<NNMergeParticle<Packed>>>&   a_scatterByDestRank) noexcept
  {
    CH_TIME("ParticleManagement::placeMergedParticles");

    a_localByPatch.clear();

    // For each result in a_results (pooled from BOTH resolveTrivialTier() and judgeProposals()
    // -- placement makes no distinction between the two, see this function's header docs):
    //   - point-locate result.merged.position against a_myPatchBoxes (a linear scan is fine at
    //     small patch counts; the real AmrMesh/DisjointBoxLayout integration may expose a
    //     cheaper spatial index -- left to the implementation).
    //   - if found in patch index i: append result.merged to a_localByPatch[i].
    //   - if not found in ANY of a_myPatchBoxes: this rank does not own the destination patch --
    //     the caller (mergeNearestNeighborsRound()) is expected to have a way to resolve "which
    //     rank owns the patch containing this position" globally (mirroring how a
    //     DisjointBoxLayout is already globally known/replicated) and bucket the particle into
    //     a_scatterByDestRank accordingly. This function signature intentionally leaves that
    //     global patch-to-rank lookup abstract (only local point-location against a_myPatchBoxes
    //     is done here) since it depends on whatever the real AmrMesh integration exposes for
    //     it.
    //
    // Note this runs even for trivial-tier (same-rank) results: the union of two patches is not
    // generally convex (T-junctions, corners, refinement boundaries), so a locally-decided
    // merge's centroid can still legitimately fall outside every one of a_myPatchBoxes.
    MayDay::Error("ParticleManagement::placeMergedParticles -- not yet implemented, see RULES.md");
  }

  template <typename P, typename Packed, typename Traits, typename Gather, typename Combine, typename Scatter>
  inline void
  mergeNearestNeighborsRound(AmrMesh&                      a_amr,
                             ParticleContainer<P, Traits>& a_particles,
                             const int                     a_numParticlesPerCellThresh,
                             const Gather&                 a_gather,
                             const Combine&                a_combine,
                             const Scatter&                a_scatter,
                             const bool                    a_iterateLocalTierToConvergence,
                             const int                     a_maxFallbackCandidates,
                             const std::optional<int>      a_maxCellDistance)
  {
    CH_TIMERS("ParticleManagement::mergeNearestNeighborsRound");
    CH_TIMER("ParticleManagement::mergeNearestNeighborsRound::ghostFill", t1);
    CH_TIMER("ParticleManagement::mergeNearestNeighborsRound::localPhase", t2);
    CH_TIMER("ParticleManagement::mergeNearestNeighborsRound::exchangeProposals", t3);
    CH_TIMER("ParticleManagement::mergeNearestNeighborsRound::judge", t4);
    CH_TIMER("ParticleManagement::mergeNearestNeighborsRound::exchangeVerdicts", t5);
    CH_TIMER("ParticleManagement::mergeNearestNeighborsRound::placement", t6);

    CH_START(t1);
    // 1. Ghost-fill, fresh, exactly once, strictly before anything below (see RULES.md,
    //    "Sequencing" -- this is a hard invariant, not an optimization detail):
    //      - fetch the three ghost masks this merge needs from a_amr: same-level
    //        (getParticleGhostMask), coarse-to-fine, and fine-to-coarse, all at whatever ghost
    //        width was registered via a_amr.registerParticleGhostMask() for this particle
    //        container's realm (the CALLER's responsibility to have registered a width at least
    //        as large as this merge's own search radius -- this function does not register
    //        anything itself).
    //      - a_particles.fillGhostParticles(maskSame, maskCoarToFine, maskFineToCoar).
    //      - a_particles.clearGhostParticles() must NOT be called again until after THIS round's
    //        placement step has fully committed -- ghosts from this fill are read throughout
    //        steps 2-3 below.
    CH_STOP(t1);

    CH_START(t2);
    // 2. Per patch this rank owns:
    //      a. Gather a_particles' valid particles (via a_gather, packed into NNMergeParticle<
    //         Packed>) and its ghosts (same gather, applied to the ghost-tagged slots
    //         ParticleSoA already carries -- see ParticleSoA::isGhost()/ghost()).
    //      b. Compute per-cell live counts (valid + ghost) using the SAME unclamped cellKeyOf()
    //         used everywhere else this round.
    //      c. Compute isBoundaryExposed() for every locally-valid particle, OR-ing across all
    //         three ghost masks from step 1 (a particle exposed toward EITHER a same-level,
    //         coarser, or finer neighbor counts as exposed).
    //      d. findNearestNeighborCandidates() for this patch.
    //      e. Pool this patch's edges into one rank-wide edge list (do NOT resolve per-patch in
    //         isolation -- resolveTrivialTier() needs the full, rank-wide, globally
    //         distance-sorted pass; see that function's docs).
    //    After all of this rank's own patches have been processed: resolveTrivialTier() over the
    //    pooled edge list, applying its results directly to the corresponding ParticleSoA
    //    containers via a_scatter for merged particles and ordinary removal for consumed ids
    //    (removal mechanics -- e.g. swap-and-pop via ParticleSoA::remove() -- are the caller's
    //    responsibility, not resolveTrivialTier()'s, which only reports which ids were consumed).
    //
    //    resolveTrivialTier() is always called with a_maxFallbackCandidates and a_maxCellDistance
    //    forwarded unchanged (see their own docs -- both are orthogonal to the iteration option
    //    below and combine with it rather than competing with it).
    //
    //    If a_iterateLocalTierToConvergence is false (the default): steps (a)-(e) plus the one
    //    resolveTrivialTier() call run exactly once, as above.
    //
    //    If a_iterateLocalTierToConvergence is true: after that first resolveTrivialTier() call,
    //    check whether it committed at least one merge. If so, repeat steps (a)-(e) plus
    //    resolveTrivialTier() again -- re-gathering valid particles (now reflecting the previous
    //    local iteration's removals/insertions), reusing the SAME ghosts and SAME boundary-exposed
    //    flags from step 1/this function's own first pass (neither changes between local
    //    iterations: ghost data depends only on the ghost-fill from step 1, which is not redone,
    //    and a particle's exposure is a fixed geometric property of its cell, not of who else is
    //    still alive) -- until a full local iteration commits zero further trivial merges.
    //    Terminate the loop as soon as one full iteration yields zero commits (not a fixed
    //    iteration count), since further iterations after that point cannot change anything
    //    without new ghost data. generateProposals()/judgeProposals() (steps 3-5 below) still run
    //    exactly once regardless of this option -- only the LOCAL portion loops.
    //
    //    NOTE on the ACTUAL measured relationship between these two options (see RULES.md for the
    //    full writeup): a_iterateLocalTierToConvergence was originally hypothesized to close most
    //    of the single-round completion gap by itself, on the reasoning that a particle whose
    //    candidate gets consumed by a closer commit elsewhere is left unmatched for the rest of
    //    that pass even if a different partner exists, and re-running the search should let it
    //    find one. Measured effect was much smaller than that reasoning predicted (~5% more
    //    merges in the validation scenario) -- re-running the search only helps when a
    //    particle's BEST choice actually changes between iterations, which turned out to be rare.
    //    The dominant loss mechanism was instead a STABLE condition (a candidate marked "busy" by
    //    its own separate, unrelated merge attempt) that persists identically across iterations,
    //    which a_maxFallbackCandidates (via resolveTrivialTier(), see its docs) is what actually
    //    addresses, and recovers far more of the gap. Both options remain independently useful
    //    and compose (this function applies a_maxFallbackCandidates on every local iteration,
    //    however many there are), but a_maxFallbackCandidates is the higher-value one.
    CH_STOP(t2);

    CH_START(t3);
    // 3. generateProposals() over whatever resolveTrivialTier() left as
    //    a_hasOutgoingCommitment. Exchange: hand a_proposalsByDestRank[this rank] to
    //    judgeProposals() directly, in-process; ship every other bucket to its destination rank
    //    and collect the incoming proposals addressed to this rank's own particles. See this
    //    file's header-level @details for the open question this step depends on: no existing
    //    general-purpose "ship arbitrary payloads to a specific rank" primitive is currently
    //    exposed outside ParticleContainer's own protected, ParticleSoA-specific
    //    gatherToPool()/gatherGhostsFromMasks()/distributeFromPool() machinery -- this exchange
    //    needs either a new general utility, or NNMergeProposal records represented as
    //    specially-tagged particles riding the existing transport. Resolve this before writing
    //    real code here.
    CH_STOP(t3);

    CH_START(t4);
    // 4. judgeProposals() over the pooled incoming proposals from step 3, using this rank's own
    //    currently-alive valid particles (post-step-2 removals already applied) as
    //    a_myParticlesByID.
    CH_STOP(t4);

    CH_START(t5);
    // 5. Exchange verdicts (same open transport question as step 3, mirrored: this rank's own
    //    in-process bucket plus real MPI for everyone else) and applyVerdicts() against this
    //    rank's own outgoing proposals from step 3. Remove every consumed id from the real
    //    ParticleSoA containers.
    CH_STOP(t5);

    CH_START(t6);
    // 6. Pool every NNMergeResult from steps 2 and 4 and call placeMergedParticles(). Append
    //    every locally-placed merged particle into the correct ParticleSoA via a_scatter.
    //    Exchange the scatter bucket (same open transport question again) and apply incoming
    //    scattered particles the same way. Only after this step completes may
    //    a_particles.clearGhostParticles() be called for the next round's fresh fill.
    CH_STOP(t6);

    MayDay::Error("ParticleManagement::mergeNearestNeighborsRound -- not yet implemented, see RULES.md");
  }

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#endif
