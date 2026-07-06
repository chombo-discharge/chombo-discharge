/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_NearestNeighborParticleMerge.H
  @brief  Declaration of a distributed, MPI-safe nearest-neighbor particle-merge algorithm.
  @author Robert Marskar
*/

#ifndef CD_NEARESTNEIGHBORPARTICLEMERGE_H
#define CD_NEARESTNEIGHBORPARTICLEMERGE_H

// Std includes
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>

// Chombo includes
#include <RealVect.H>
#include <IntVect.H>

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_ParticleContainer.H>
#include <CD_ParticleGhostMask.H>
#include <CD_NamespaceHeader.H>

class AmrMesh;

/**
  @brief Namespace for various particle management tools.
  @details This file adds a nearest-neighbor super-particle merge algorithm that is safe across
  patch and rank boundaries -- see RULES.md at the repository root for the full design rationale
  and the rules this implementation must follow. In short: candidates are always found against a
  particle's TRUE full neighborhood (own patch plus every ghost, filled once per round, never
  restricted), so boundary-adjacent particles are never biased toward worse matches than interior
  ones; pairs where neither particle is reachable from outside the patch resolve immediately
  ("trivial tier"); everything else goes through a propose/judge/verdict protocol that is safe for
  any number of simultaneous cross-patch suitors (not just pairwise boundaries) and free of the
  ordering-dependent bugs a naive implementation of the same idea runs into (see RULES.md,
  "Judge", for the specific failure modes this design closes).

  This algorithm was designed and validated (mass conservation, no double-consumption, and
  determinism under adversarially randomized processing order, across six 2D/3D domain topologies)
  with a Python/mpi4py prototype before being translated into this header -- see
  Scripts/ParticleMergePrototype/. The prototype's `merge_lib.py` is the direct analogue of the
  "pure logic" functions declared below (dimension-generic, transport-agnostic); its
  `run_scenario.py` is the analogue of mergeNearestNeighborsRound()'s orchestration, minus the
  real MPI exchange, which the prototype stood in for with mpi4py's `alltoall`.
*/
namespace ParticleManagement {

  /**
    @brief Minimal, payload-agnostic description of one particle as input to the nearest-neighbor
    merge algorithm.
    @details Carries exactly what the matching/decision logic needs -- position, weight, a
    globally unique id, and the particle's CURRENT owning rank -- plus an OPAQUE, caller-defined
    reconciliation payload. The matching/decision logic (candidate search, trivial-tier
    resolution, proposal routing, judging) never inspects @c payload; it is only unpacked by the
    caller-supplied combine callback when a pair is actually finalized. This mirrors how the
    existing mergers in CD_ParticleManagement.H gather SoA columns into a small packed
    intermediate (e.g. NonCommParticle<M,N>) rather than operating on the full particle type --
    see RULES.md, "Data model", for why this is deliberately NOT a static per-type column
    declaration (i.e. NOT reusing the H5Part column-descriptor pattern): different merge
    algorithms may legitimately want different reconciliation data, which only a caller-supplied
    combine callback can express.
    @tparam Packed Caller-defined payload type (e.g. a single Real energy for ItoParticle, or an
    empty struct if no type-specific reconciliation data is required). Must be copyable and
    trivially linearizable (see NNMergeProposal, which ships two of these as one atomic pair).
  */
  template <typename Packed>
  struct NNMergeParticle
  {
    /**
      @brief Particle position.
    */
    RealVect position;

    /**
      @brief Particle weight (used for the crowding count and the merged-particle centroid; not
      required to be an integer).
    */
    Real weight;

    /**
      @brief Globally unique particle id. Used to break cross-patch merge ties deterministically
      (see judgeProposals()) and, together with @c payload, as the audit key for a lineage-based
      mass-conservation check (see the prototype's per-particle ledger audit in
      Scripts/ParticleMergePrototype/run_scenario.py).
    */
    ParticleID globalID;

    /**
      @brief The rank that currently owns this particle. ALWAYS read from the particle's own
      data -- see RULES.md, "Data model": this is the single source of truth and must never be
      inferred from MPI transport/source-rank metadata (e.g. the sender of an MPI message is not
      necessarily authoritative about who logically owns the particle it carries).
    */
    RankID ownerRank;

    /**
      @brief The AMR level this particle currently lives on.
      @details REQUIRED for correct multi-level behavior: cellKeyOf() needs the cell spacing (dx)
      of the LEVEL A PARTICLE ACTUALLY LIVES ON to produce a meaningful cell index, and a ghost
      candidate can legitimately be on a different level than the querying patch (a
      GhostType::Coarse or GhostType::Fine ghost, not just GhostType::SameLevel -- see
      ParticleSoA::GhostType). Every function that computes a cell key for crowding/threshold
      purposes must look up the per-level dx for THIS field, never assume a single global dx --
      see resolveTrivialTier()/judgeProposals()'s a_dxByLevel parameter.
    */
    int level;

    /**
      @brief Opaque, caller-defined reconciliation payload. Never inspected by this algorithm's
      matching/decision logic -- only carried through to the caller-supplied combine callback.
    */
    Packed payload;
  };

  /**
    @brief One additional ranked candidate beyond a query's primary nearest neighbor.
    @details Only ever populated when a_maxFallbackCandidates > 0 is requested from
    findNearestNeighborCandidates() -- see NNMergeEdge::fallbackCandidates and
    resolveTrivialTier()'s docs for how this is consumed. Fetched via a single K-nearest spatial
    query at candidate-generation time (K = 1 + a_maxFallbackCandidates), NOT via a fresh query
    per retry attempt -- resolveTrivialTier() simply walks further down this pre-fetched,
    already-distance-sorted list when a closer candidate turns out to be blocked, which avoids
    re-traversing the spatial index on every retry.
  */
  struct NNFallbackCandidate
  {
    /**
      @brief Global id of this candidate.
    */
    ParticleID candidateID;

    /**
      @brief Squared distance between the query particle and this candidate.
    */
    Real distanceSquared;

    /**
      @brief True iff this candidate is a locally-valid particle in the SAME patch as the query
      (i.e. not a ghost) -- same meaning as NNMergeEdge::candidateIsLocal.
    */
    bool candidateIsLocal;
  };

  /**
    @brief One candidate merge edge: a single particle's own nearest-neighbor query result.
    @details Produced once per over-threshold, still-unconsumed valid particle by
    findNearestNeighborCandidates() -- exactly one edge per query particle, its single TRUE
    nearest neighbor among {this patch's own other valid particles} union {every ghost}. Distance
    is carried explicitly (never recomputed downstream) so that both the trivial-tier resolution
    and the judge can apply the SAME "process globally closest-first, with a deterministic id
    tiebreak on exact ties" discipline -- see RULES.md, "Trivial tier" and "Judge", for why both
    of those needed this (an unbroken exact-distance tie left to arbitrary
    iteration/insertion order was a real, observed source of nondeterminism that cascaded into
    unrelated downstream decisions).
    @tparam Packed See NNMergeParticle.
  */
  template <typename Packed>
  struct NNMergeEdge
  {
    /**
      @brief Squared distance between the query particle and its candidate.
    */
    Real distanceSquared;

    /**
      @brief Global id of the query particle (always locally valid/owned -- ghosts never
      initiate a query, see RULES.md, "Propose").
    */
    ParticleID queryID;

    /**
      @brief Global id of the candidate (query's nearest neighbor). May be a local valid
      particle (same patch) or a ghost.
    */
    ParticleID candidateID;

    /**
      @brief True iff the candidate is a locally-valid particle in the SAME patch as the query
      (i.e. not a ghost). Note this is NOT sufficient on its own to decide trivial-tier
      eligibility -- see isBoundaryExposed() and RULES.md, "Boundary-exposed classification": a
      same-patch candidate can still require the propose/judge/verdict protocol if either
      participant is boundary-exposed.
    */
    bool candidateIsLocal;

    /**
      @brief Next-nearest candidates beyond candidateID, in ascending distanceSquared order.
      @details Empty unless a_maxFallbackCandidates > 0 was passed to
      findNearestNeighborCandidates(); size at most a_maxFallbackCandidates. See
      resolveTrivialTier() for how these are consumed: if candidateID turns out to be blocked
      (stale, busy, exposed, a ghost, or too far -- see a_maxCellDistance), the query particle
      tries fallbackCandidates[0], then [1], and so on, before giving up -- never re-querying the
      spatial index for this edge. This member's own content is never re-derived or refreshed once
      populated; it reflects the position/ownership snapshot at generation time, same as
      candidateID itself -- validity (still alive, not busy) is always rechecked dynamically
      against a_consumedIDs/a_hasOutgoingCommitment when a candidate is actually tried, never
      trusted from anything computed here.
    */
    std::vector<NNFallbackCandidate> fallbackCandidates;
  };

  /**
    @brief A cross-patch merge proposal: one particle's own nearest-neighbor query result,
    bundled as a self-contained, atomic pair for shipment to whichever rank owns the candidate.
    @details Carries FULL data for both participants (not just a bare id reference) so that (a)
    the judge's decision is entirely self-contained -- no need to cross-reference separately
    tracked local state -- and (b) the judge is structurally prevented from ever trying to
    recombine pieces across DIFFERENT incoming proposals; a proposal is always an opaque,
    pre-formed pair with a distance already attached. See RULES.md, "Propose".

    Routing rule (see generateProposals()): a proposal is ALWAYS sent to the owner of whatever
    the querying particle found as its own best candidate -- never by any id comparison. An
    earlier "route to whoever owns the lower-id particle" framing was found to be unsound: it can
    send the two directions of a genuinely mutual best-match pair to two different,
    uncoordinated judges with no way to reconcile them.
    @tparam Packed See NNMergeParticle.
  */
  template <typename Packed>
  struct NNMergeProposal
  {
    /**
      @brief The proposer's own particle (full data) -- always locally valid/owned by the
      sending rank.
    */
    NNMergeParticle<Packed> source;

    /**
      @brief Global id of the particle being proposed to. Owned by the receiving (judge) rank.
    */
    ParticleID targetID;

    /**
      @brief Squared distance between @c source and the target. Duplicated from the originating
      NNMergeEdge so the judge never needs to look anything up to sort/compare proposals.
    */
    Real distanceSquared;
  };

  /**
    @brief The judge's reply to one proposal.
    @details See RULES.md, "Verdict and placement": accepted implies the judge has already
    created the merged particle and removed its own (target) particle; the proposer, on
    receiving an accepted verdict, removes its own (source) particle. Rejected implies no action
    on either side -- the proposer's particle is simply reconsidered next round.
  */
  struct NNMergeVerdict
  {
    /**
      @brief Global id of the proposer this verdict is a reply to (i.e. NNMergeProposal::source's
      id). Used by the proposer to match an incoming verdict back to one of its own outgoing
      proposals.
    */
    ParticleID proposerID;

    /**
      @brief True iff the judge accepted this specific proposal (i.e. it was the winning
      candidate for its target -- see judgeProposals()).
    */
    bool accepted;
  };

  /**
    @brief One resolved merge: the two consumed particles' ids and the newly created merged
    particle (not yet placed into any patch -- see placeMergedParticles()).
    @tparam Packed See NNMergeParticle.
  */
  template <typename Packed>
  struct NNMergeResult
  {
    /**
      @brief Global id of the first consumed particle. For a trivial-tier result, this is the
      query particle. For a judge-side result, this is the target.
    */
    ParticleID firstID;

    /**
      @brief Global id of the second consumed particle. For a trivial-tier result, this is the
      local candidate. For a judge-side result, this is the winning proposer.
    */
    ParticleID secondID;

    /**
      @brief The newly created merged particle: weighted-centroid position, weighted-average
      opaque payload (computed by the caller-supplied combine callback), summed weight, and a
      FRESH globally-unique id (never colliding with any existing particle's id -- see
      mergeNearestNeighborsRound() for the id-allocation requirement this implies).
    */
    NNMergeParticle<Packed> merged;
  };

  /**
    @brief Combine callback signature: produce a merged particle's opaque payload from its two
    parents' payloads and weights.
    @details Deliberately a template parameter wherever it is consumed (see
    resolveTrivialTier(), judgeProposals()), NOT a std::function -- this sits on a hot path (one
    call per resolved merge) and must actually inline, unlike the std::function-based
    a_gather/a_combine/a_scatter callbacks the EXISTING mergers in CD_ParticleManagement.H
    accept. Position and weight are handled by this algorithm directly (weighted-centroid,
    summed weight); this callback is responsible ONLY for the opaque @c Packed payload (e.g.
    weighted-average energy for ItoParticle).
    @param[in] a_first  First parent's payload.
    @param[in] a_firstWeight First parent's weight.
    @param[in] a_second Second parent's payload.
    @param[in] a_secondWeight Second parent's weight.
    @return The merged particle's payload.
  */
  template <typename Packed>
  using PayloadCombine = std::function<
    Packed(const Packed& a_first, const Real a_firstWeight, const Packed& a_second, const Real a_secondWeight)>;

  /**
    @brief Position-validity predicate for a would-be merged particle.
    @details Deliberately a template parameter wherever consumed (see resolveTrivialTier(),
    judgeProposals()), NOT a std::function, for the same inlining reason as PayloadCombine.
    Called with the weighted-centroid position ALREADY computed, exactly once per candidate pair
    that is otherwise eligible to merge -- never on a pair that would already be rejected for an
    unrelated reason (staleness, busy, exposure, a_maxCellDistance), so a caller doing real
    geometric work here (e.g. an embedded-boundary point query) never pays for it needlessly.
    Default is an always-true predicate (matches the original, no-op behavior). A pair that fails
    this check is treated exactly like a_maxCellDistance's "too far" case: never proposed either
    (there is nothing valid to hand a judge if the position itself is invalid), a fallback
    candidate is tried instead if available, and otherwise the particle is simply left unmatched
    this round with no side effects on anyone else's processing.
    @param[in] a_mergedPosition The candidate merge's weighted-centroid position.
    @return True iff a particle is allowed to exist at this position (e.g. false if it would fall
    inside an embedded boundary).
  */
  using PositionValid = std::function<bool(const RealVect& a_mergedPosition)>;

  /**
    @brief Fresh-id allocator callback: produce one globally-unique id for a newly created merged
    particle.
    @details Deliberately a template parameter wherever consumed (see resolveTrivialTier(),
    judgeProposals()), NOT a std::function, for the same inlining reason as PayloadCombine.
    Called exactly once per committed merge (trivial-tier or judge-side). The caller (ultimately
    mergeNearestNeighborsRound()) is responsible for ensuring ids never collide -- not just within
    one call, but across every rank and every past/future call for the life of the simulation; a
    simple rank-namespaced counter (this rank's index times a stride, plus a monotonically
    increasing local count) is sufficient and is what the Python prototype used (see
    Scripts/ParticleMergePrototype/merge_lib.py's make_particle()) -- there is no existing
    general-purpose global-id allocation scheme elsewhere in the particle infrastructure to reuse
    instead.
    @return A fresh globally-unique particle id.
  */
  using IDAllocator = std::function<ParticleID()>;

  /**
    @brief The unclamped, position-derived cell index a given physical position falls in.
    @details This is the SAME cell-key convention used throughout this algorithm for the crowding
    trigger and a_maxCellDistance (see RULES.md, "Crowding trigger") -- deliberately NOT clamped to
    any one patch's own box, so that a ghost particle buckets into its true origin cell and two
    particles in different patches' own local frames remain directly comparable.
    @warning a_dx MUST be the cell spacing of the AMR level the position's particle actually lives
    on (see NNMergeParticle::level) -- NEVER a single fixed dx shared across every level. An
    IntVect returned by this function is only meaningful together with the level it was computed
    for; see NNCellKey, which pairs the two so they can never be silently conflated (a bug found
    empirically: pooling bare IntVect cell keys from different levels into one crowding-count map,
    or computing a dynamic recheck with the wrong level's dx, causes coarse- and fine-level cells
    to alias each other's counts, silently suppressing merges on every level but one).
    @param[in] a_position Physical position.
    @param[in] a_probLo   Physical position of the domain's low corner.
    @param[in] a_dx       Cell spacing (per direction) of the position's OWN level.
    @return floor((a_position - a_probLo) / a_dx), component-wise.
  */
  inline IntVect
  cellKeyOf(const RealVect& a_position, const RealVect& a_probLo, const RealVect& a_dx) noexcept;

  /**
    @brief A cell key together with the AMR level it was computed on -- see cellKeyOf()'s warning.
    @details Two NNCellKeys are equal only if BOTH the level and the cell index match; a coarse-
    level cell and an unrelated fine-level cell that happen to share the same bare IntVect are
    correctly kept distinct because their level differs. This is the key type for every live
    per-cell count pooled across multiple levels (see a_liveCellCount in resolveTrivialTier()/
    judgeProposals()/findNearestNeighborCandidates()).
  */
  using NNCellKey = std::pair<int, IntVect>;

  /**
    @brief Whether a particle is boundary-exposed -- i.e. would be shipped as a ghost to at
    least one neighboring patch.
    @details This is the eligibility test for the trivial tier (see resolveTrivialTier()) -- NOT
    "same rank" and NOT "same ParticleSoA" (both were tried during design and found unsound: a
    particle can be exposed to a different patch owned by the very same rank, and a particle can
    be exposed while its candidate happens to sit in the same ParticleSoA -- see RULES.md,
    "Boundary-exposed classification", for the failure modes this rule closes). This is a direct
    lookup against the patch's own OUTGOING ghost mask -- i.e. whether this particle's cell has
    any registered destination -- not a new geometric computation: @c ParticleGhostMask::numTargets()
    already answers exactly this question for the cell a particle lives in.
    @param[in] a_cellIndex Cell index (IntVect) of the particle, in this patch's own local frame
                            -- i.e. whatever a_ghostMask itself is keyed by (see
                            ParticleGhostMask::define()/addTarget()).
    @param[in] a_ghostMask This patch's own OUTGOING particle ghost mask, as returned by
                            AmrMesh::getParticleGhostMask() (or getParticleGhostMaskCoarToFine()/
                            getParticleGhostMaskFineToCoar() for cross-level exposure -- a
                            particle can be boundary-exposed toward a coarser or finer level
                            too, not just the same level; see mergeNearestNeighborsRound() for
                            how the three mask kinds are combined into one boundary-exposed
                            test).
    @return True iff the particle's cell has at least one outgoing ghost-shipment target.
  */
  inline bool
  isBoundaryExposed(const IntVect& a_cellIndex, const ParticleGhostMask& a_ghostMask) noexcept;

  /**
    @brief For every over-threshold, currently-alive valid particle in a patch, find its single
    TRUE nearest neighbor among {this patch's own other valid particles} union {every ghost
    shipped to this patch}.
    @details This is the ONLY place a nearest-neighbor search happens, and it is ALWAYS run
    against the full (local + ghost) view -- never a restricted, local-only view -- see RULES.md,
    "Sequencing": restricting this search for some particles (e.g. "same-patch candidates only,
    resolve those first") is exactly the bug this design's predecessor had ("compression toward
    patch interior" -- boundary-adjacent particles got locked into worse matches than interior
    ones purely because of where the domain happened to be cut). Crowding eligibility (is the
    QUERY particle's own cell currently over threshold) is checked here, at generation time,
    using the live per-cell counts computed from BOTH valid and ghost occupancy (a ghost is a
    faithful stand-in for its true origin cell's population); the SYMMETRIC check on the
    candidate's cell is deliberately NOT done here -- it is re-checked dynamically, at
    processing time, by resolveTrivialTier() and judgeProposals(), since earlier merges within
    the same round can change it (see RULES.md, "Crowding trigger").

    Implementation is expected to build a spatial index (e.g. a KD-tree, matching the "kd-tree
    over the particles in the patch" framing this whole feature started from) over the combined
    local-valid-plus-ghost position set once per patch per round, then query it once per
    over-threshold valid particle -- mirroring the allocation-light, thread-local-scratch style
    already used by the existing mergers in CD_ParticleManagementImplem.H (e.g.
    buildEqualWeightKDLeaves(), mergeAdjacentNearest()). That existing infrastructure builds a
    DIFFERENT kind of tree, though (equal-weight spatial partitioning -- its "leaves" are groups
    of particles, not point-query results), so it isn't directly reusable here; only its
    underlying philosophy carries over (no per-node heap allocation, no pointer-chasing, reused
    thread-local scratch across calls).

    @b Leaf contents: each point entered into the spatial index needs only its position
    (duplicated inline, since that's what the traversal actually compares -- the entire point of
    building an index rather than following pointers) plus a single index tagged local-vs-ghost,
    referring back into a_localValid or a_ghosts. Nothing else needs to be duplicated into the
    tree -- id, owner rank, and the opaque payload are all already available in O(1) from that
    index once a query returns, since a_localValid/a_ghosts are already the small, purpose-built
    NNMergeParticle records this decision layer needs (not the full ParticleSoA column set).
    Deliberately NOT tagged with anything about a candidate's current busy/consumed status: that
    status doesn't exist yet when the tree is built (nothing has been resolved at this point in
    the round), and it changes continuously as resolveTrivialTier() processes its sorted edge
    list afterward -- baking it into the tree would mean mutating tree data on every single
    commit, which defeats the purpose of building a static index once. Busy/consumed status is
    tracked externally (a_consumedIDs, a_hasOutgoingCommitment) and rechecked dynamically wherever
    a candidate is actually about to be used, never trusted from anything attached to the tree.

    When a_maxFallbackCandidates (forwarded from resolveTrivialTier(), see its own docs) is
    greater than zero: query the tree for the K=1+a_maxFallbackCandidates nearest points instead
    of just 1, and populate NNMergeEdge::fallbackCandidates with candidates 2..K (already
    distance-sorted by construction of a K-nearest query). This is a single spatial-index
    traversal per query particle regardless of K -- resolveTrivialTier() walks the pre-fetched
    fallbackCandidates list when the primary candidate is blocked, rather than this function (or
    the tree) being asked again per retry attempt.

    @tparam Packed See NNMergeParticle.
    @param[in]  a_localValid       This patch's own valid particles (candidates AND query
                                    sources).
    @param[in]  a_ghosts           Every ghost shipped to this patch this round (candidates
                                    only; ghosts never initiate a query).
    @param[in]  a_liveCellCount    Live particle count (valid + ghost) per level-qualified cell
                                    key (see NNCellKey), already reflecting this round's frozen
                                    ghost-fill snapshot.
    @param[in]  a_probLo           Physical position of the domain's low corner (see cellKeyOf()).
    @param[in]  a_dx               Cell spacing of THIS PATCH's own level (see cellKeyOf()'s
                                    warning) -- valid here because every query in a_localValid is
                                    necessarily on this same level; candidates never have their
                                    crowding checked by this function, so a cross-level ghost
                                    candidate needs no dx here.
    @param[in]  a_numParticlesPerCellThresh Crowding threshold N -- a cell is eligible while it
                                    holds strictly more than this many particles.
    @param[in]  a_maxFallbackCandidates Forwarded from resolveTrivialTier() (see its docs); 0
                                    (default) means an ordinary single-nearest query and every
                                    edge's fallbackCandidates stays empty.
    @param[out] a_edges            One entry per over-threshold, currently-alive query particle
                                    that has at least one candidate (empty if none found).
  */
  template <typename Packed>
  inline void
  findNearestNeighborCandidates(const std::vector<NNMergeParticle<Packed>>& a_localValid,
                                const std::vector<NNMergeParticle<Packed>>& a_ghosts,
                                const std::map<NNCellKey, int>&             a_liveCellCount,
                                const RealVect&                             a_probLo,
                                const RealVect&                             a_dx,
                                const int                                   a_numParticlesPerCellThresh,
                                std::vector<NNMergeEdge<Packed>>&           a_edges,
                                const int                                   a_maxFallbackCandidates = 0) noexcept;

  /**
    @brief Resolve every candidate edge that is safely mergeable with zero cross-rank
    coordination, and mark which remaining particles have committed to a non-trivial outgoing
    proposal.
    @details An edge is trivial-tier eligible iff NEITHER the query nor the candidate is
    boundary-exposed (see isBoundaryExposed()) -- evaluated per RULES.md, "Trivial tier": this
    check happens AFTER finding the true nearest neighbor, never as a restriction on the search
    itself.

    All eligible edges (across every patch this rank owns) are processed in ONE globally
    distance-sorted order, with a deterministic secondary tiebreak on participant ids -- exact
    distance ties (e.g. true mutual nearest neighbors) are common enough that leaving them to
    arbitrary iteration/insertion order produces real, observed nondeterminism that cascades into
    unrelated pairs several hops away.

    At the moment each edge is actually processed (not at generation time): re-verify both
    participants are still alive (lazy invalidation -- skip if either was already consumed by a
    closer edge processed earlier in this same pass); if the candidate is local, re-verify it has
    NOT, in the meantime, committed to its OWN outgoing (non-trivial) proposal -- a purely
    geometric "exposed" check on the candidate alone is insufficient, since a non-exposed
    particle's true nearest neighbor can still be a foreign ghost even though the particle itself
    is unreachable externally; consuming such a candidate trivially anyway is a real
    double-consumption bug this dynamic recheck closes. Finally, re-verify both participants'
    CURRENT cell occupancy (reflecting decrements from edges already processed this pass, not the
    generation-time snapshot) still exceeds the threshold.

    On a successful trivial-tier commit: remove both participants from the caller's own valid
    particle list (conceptually -- see a_consumedIDs), decrement both cells' live counts, and record the merge via
    a_combine/into a_results. On an edge that cannot be resolved trivially, exactly one of two
    things happens, and the distinction matters:
      - If the candidate is a ghost, OR either participant is boundary-exposed, OR the candidate
        turned out to already have its own outgoing commitment: mark the query particle as
        having an outgoing commitment (see a_hasOutgoingCommitment) so generateProposals() can
        pick it up. In all three cases the candidate genuinely still exists (it just isn't
        trivially available right now), so a proposal is meaningful -- judgeProposals() already
        has a rule for "target has its own outgoing commitment" (accept only on a mutual match,
        reject otherwise), and that rule is exactly what should get a chance to run here, not a
        silent no-op. (An earlier revision of this design silently dropped the busy-candidate
        case instead of proposing -- confirmed empirically, in the Python prototype, to be a
        genuine completeness gap: it left the query particle with no outcome at all for the
        round, not even a clean reject, and measurably hurt overall completion once
        a_maxFallbackCandidates below was introduced and made the omission's effect visible.)
      - If the candidate is stale (already consumed -- found dead by the lazy-invalidation check
        above): do NOT mark the query particle as having an outgoing commitment, and do NOT
        generate a proposal. There is no meaningful id to propose to. Critically, this also means
        the query particle is NOT marked "busy" -- it remains fully available for some OTHER
        particle to claim as an ordinary trivial candidate later in this same sorted pass, or via
        a_maxFallbackCandidates below.

    @tparam Packed  See NNMergeParticle.
    @tparam Combine See PayloadCombine.
    @param[in]     a_edges           Candidate edges from findNearestNeighborCandidates(), for
                                      ALL patches this rank owns (pooled together so the global
                                      distance sort spans the whole rank, not just one patch).
    @param[in]     a_exposed         Boundary-exposed flag (isBoundaryExposed()) per particle id,
                                      for every locally-valid particle across all of this rank's
                                      patches.
    @param[in]     a_particlesByID   Full particle data (position/weight/id/rank/payload) for
                                      every currently-alive local valid particle and every ghost
                                      across all of this rank's patches, keyed by id -- needed to
                                      actually compute a merge's weighted-centroid position and
                                      call a_combine (an edge only carries ids/distance, not full
                                      particle data).
    @param[in]     a_combine         Opaque-payload combine callback (see PayloadCombine).
    @param[in,out] a_liveCellCount   Live per-cell counts; decremented in place as trivial merges
                                      commit, so later edges in the same pass see the up-to-date
                                      occupancy.
    @param[in,out] a_consumedIDs     Ids already consumed this round (by an earlier trivial merge,
                                      or found dead for any other reason); an edge referencing a
                                      consumed id is stale and skipped. Updated in place.
    @param[in]     a_numParticlesPerCellThresh Crowding threshold N, for the dynamic per-cell
                                      recheck (same value used everywhere else this round).
    @param[out]    a_hasOutgoingCommitment Set of query-particle ids that were NOT resolved
                                      trivially (excluding the stale-candidate case, see above)
                                      and must go on to generateProposals().
    @param[out]    a_results         One entry per trivial-tier merge actually committed, in the
                                      order they were processed.
    @param[in]     a_maxFallbackCandidates User-selectable option, default 0 (disabled, matches
                                      the originally documented single-candidate-per-particle
                                      behavior). When a query particle's chosen candidate is
                                      blocked -- already consumed, OR already has its own outgoing
                                      commitment, OR exposed/a ghost -- and the query particle
                                      itself is NOT exposed (if it is, no local candidate could
                                      ever be trivially eligible, so retrying is pointless): try
                                      the edge's own NNMergeEdge::fallbackCandidates in order
                                      (already fetched by findNearestNeighborCandidates() via one
                                      K-nearest spatial-index query, K = 1+a_maxFallbackCandidates
                                      -- see its docs), up to this many extra attempts, before
                                      giving up. There is no fresh spatial-index query here, and
                                      no "excluding every candidate already tried" bookkeeping to
                                      maintain -- the pre-fetched list is already in ascending
                                      distance order with duplicates impossible, so trying the
                                      next one is just advancing an index. Retries are resolved
                                      immediately at the point of failure, not deferred to
                                      preserve strict global distance-sort order across every
                                      particle's retries -- a deliberate simplification, justified
                                      because a retry only ever fires for an edge whose first
                                      choice already failed, so it was never competing on equal
                                      footing with anyone's first attempt anyway. The proposal
                                      generated for a query particle
                                      that exhausts every fallback (see a_hasOutgoingCommitment
                                      above) still uses its TRUE, original nearest neighbor --
                                      fallbacks only affect local trivial-tier eligibility, never
                                      what gets proposed cross-patch. Measured in the Python
                                      prototype (Scripts/ParticleMergePrototype/) to be the single
                                      largest lever found so far for closing the single-round
                                      completion gap: in the 64-particles-per-cell/threshold-16
                                      scenario also used to validate a_iterateLocalTierToConvergence,
                                      instrumentation showed "candidate already has an outgoing
                                      commitment" was the single largest cause of lost local
                                      merges -- larger than genuine boundary exposure. A single
                                      fallback candidate recovered the great majority of that loss
                                      in one pass; additional fallback attempts beyond that gave
                                      steeply diminishing returns. See RULES.md, "Optional: retry
                                      with a fallback candidate when the first choice is blocked"
                                      for the full empirical writeup.
    @param[in]     a_maxCellDistance User-selectable physical merge-eligibility cap, in whole
                                      grid cells, using Chebyshev (max-per-axis) distance between
                                      the SAME unclamped, position-derived cell key already used
                                      for the crowding trigger (`floor(pos/dx)`, see
                                      findNearestNeighborCandidates()'s docs) -- same cell = 0,
                                      any of the 8 (2D) / 26 (3D) Moore-neighborhood-adjacent
                                      cells including diagonals = 1, and so on. Using this
                                      position-derived key rather than each patch's own local
                                      IntVect cell index is what makes the comparison meaningful
                                      across a patch boundary in the first place -- two particles
                                      in different patches' own local frames aren't directly
                                      comparable otherwise. `std::nullopt` (default) disables
                                      this and matches the original, unrestricted-search behavior
                                      (any candidate, however far, is eligible). A query particle
                                      whose TRUE nearest neighbor exceeds this cap is NOT proposed
                                      to -- unlike the busy/exposed/ghost cases, this is a
                                      physical-validity rule, not a coordination one: a pair this
                                      far apart is exactly as invalid via judgeProposals() as it
                                      would be trivially, so there is nothing valid to hand off to
                                      a judge either. If a_maxFallbackCandidates above is nonzero,
                                      the query particle gets a fallback retry for a CLOSER
                                      candidate instead of giving up outright; otherwise it is
                                      simply left unmatched this round, with no side effects on
                                      any other particle's processing (same treatment as a stale
                                      candidate, see a_hasOutgoingCommitment above). Verified
                                      directly in the Python prototype with a hand-placed pair 2
                                      cells apart: uncapped, they merge at the exact weighted
                                      centroid; capped at 1 cell, both survive completely
                                      untouched; capped at exactly 2 cells (the boundary), they
                                      merge again -- confirming the cap is inclusive.
    @param[in]     a_probLo          Physical position of the domain's low corner (see
                                      cellKeyOf()); needed for the dynamic crowding recheck and
                                      a_maxCellDistance.
    @param[in]     a_dxByLevel       Cell spacing PER AMR LEVEL (index 0 = level 0, etc.) -- see
                                      cellKeyOf()'s warning. Each participant's OWN cell key is
                                      computed with a_dxByLevel[participant.level], never a single
                                      shared dx: a candidate can legitimately be a ghost from a
                                      coarser or finer level than the query's own patch, and
                                      getting this wrong silently breaks the dynamic crowding
                                      recheck for every level but one (confirmed empirically).
    @param[in]     a_isPositionValid Position-validity predicate (see PositionValid), called once
                                      per otherwise-eligible pair with the weighted-centroid
                                      position already computed. Treated exactly like
                                      a_maxCellDistance's "too far" case on failure: no proposal
                                      for the original candidate, a fallback retry if available,
                                      otherwise left unmatched with no side effects.
    @param[in]     a_allocateID      Fresh-id allocator (see IDAllocator), called once per
                                      committed trivial-tier merge.
    @param[in,out] a_crossLevelMergeCount Optional diagnostic counter (nullptr, the default,
                                      disables it entirely -- zero cost when not needed):
                                      incremented once per committed merge whose two participants
                                      were on DIFFERENT AMR levels (query.level != cand.level).
                                      NOT reset by this function -- the caller owns its lifetime,
                                      so counts can be accumulated across a_iterateLocalTierToConvergence's
                                      internal iterations, or across this rank's judgeProposals()
                                      acceptances too (see that function's own parameter).
  */
  template <typename Packed, typename Combine, typename PosValid, typename Allocator>
  inline void
  resolveTrivialTier(const std::vector<NNMergeEdge<Packed>>&              a_edges,
                     const std::map<ParticleID, bool>&                    a_exposed,
                     const std::map<ParticleID, NNMergeParticle<Packed>>& a_particlesByID,
                     const Combine&                                       a_combine,
                     std::map<NNCellKey, int>&                            a_liveCellCount,
                     std::set<ParticleID>&                                a_consumedIDs,
                     const int                                            a_numParticlesPerCellThresh,
                     std::set<ParticleID>&                                a_hasOutgoingCommitment,
                     std::vector<NNMergeResult<Packed>>&                  a_results,
                     const RealVect&                                      a_probLo,
                     const std::vector<RealVect>&                         a_dxByLevel,
                     const PosValid&                                      a_isPositionValid,
                     const Allocator&                                     a_allocateID,
                     const int                                            a_maxFallbackCandidates = 0,
                     const std::optional<int>                             a_maxCellDistance       = std::nullopt,
                     unsigned long long*                                  a_crossLevelMergeCount  = nullptr) noexcept;

  /**
    @brief Package one outgoing proposal for every edge that resolveTrivialTier() could not
    resolve locally, bucketed by destination rank.
    @details Routing rule: ALWAYS to the owner of the candidate (NNMergeEdge::candidateID's
    owning rank) -- never by any id comparison; see NNMergeProposal's documentation for why an
    id-based routing rule was tried and found unsound. If the destination rank equals this
    rank's own (a boundary-exposed candidate that happens to live in a DIFFERENT patch owned by
    the same rank), the proposal still goes into the normal per-destination-rank bucket -- the
    caller (mergeNearestNeighborsRound()) is expected to hand the "bucket addressed to myself"
    directly to judgeProposals() in-process, without an actual MPI round-trip, while every other
    bucket goes out over MPI. This keeps the decision logic in judgeProposals() completely
    unaware of whether a given proposal arrived locally or over the network -- see RULES.md,
    "Propose": "this still goes through the exact same pooled judge logic ... never a shortcut."
    @tparam Packed See NNMergeParticle.
    @param[in]  a_edges                 Candidate edges from findNearestNeighborCandidates(), for
                                         all of this rank's patches.
    @param[in]  a_hasOutgoingCommitment Query-particle ids that resolveTrivialTier() could not
                                         resolve (i.e. exactly the edges this function packages).
    @param[in]  a_particlesByID         Full particle data (position/weight/id/rank/payload) for
                                         every currently-alive local valid particle and every
                                         ghost, keyed by id -- needed to build the full
                                         NNMergeParticle bundle for each proposal.
    @param[in]  a_numRanks              Total MPI rank count (size of the returned bucket vector).
    @param[out] a_proposalsByDestRank   One bucket per rank; a_proposalsByDestRank[r] holds every
                                         proposal whose candidate is owned by rank r (including
                                         r == this rank's own).
  */
  template <typename Packed>
  inline void
  generateProposals(const std::vector<NNMergeEdge<Packed>>&              a_edges,
                    const std::set<ParticleID>&                          a_hasOutgoingCommitment,
                    const std::map<ParticleID, NNMergeParticle<Packed>>& a_particlesByID,
                    const int                                            a_numRanks,
                    std::vector<std::vector<NNMergeProposal<Packed>>>&   a_proposalsByDestRank) noexcept;

  /**
    @brief Judge every incoming proposal that targets a particle this rank owns, and produce the
    verdicts to send back.
    @details See RULES.md, "Judge", for the full rule set this implements; summarized:
      1. Group incoming proposals by target id.
      2. For a target with NO outgoing commitment of its own: the eligible candidate is the
         minimum-distance incoming proposal (deterministic id tiebreak on exact ties).
      3. For a target that DOES have its own outgoing commitment (i.e. it is simultaneously
         proposing to something else): normally reject every incoming proposal unconditionally
         (it is "spoken for") -- EXCEPT when its own outgoing target is itself one of its current
         incoming proposers, a genuinely mutual best-match pair. Blindly rejecting that case
         causes a PERMANENT deadlock (a stable mutual pair proposes to each other identically
         every round and is rejected every round, forever, since nothing changes to break the
         symmetry) -- so it must be accepted instead. Because BOTH sides of a mutual pair detect
         the match independently, a deterministic tiebreak (lower-id side performs the accept,
         higher-id side falls back to the normal reject) is required, or both sides would accept
         and double-merge the same two particles.
      4. ALL eligible target decisions (plain-argmin and mutual-match alike) are resolved in ONE
         globally distance-sorted pass across every target this rank is judging (with an id
         tiebreak) -- NOT by iterating targets in whatever order the incoming-proposal list
         happens to be in. Accepting one target decrements its cell's live count, which can push
         a DIFFERENT target sharing that cell below threshold; judging targets in arbitrary order
         makes that outcome depend on arbitrary processing order. Sorting closest-first with a
         dynamic per-target threshold recheck (mirroring resolveTrivialTier()'s discipline
         exactly) removes that dependency.
      5. This handles any number of simultaneous suitors for one particle, not just pairwise
         boundaries -- a corner where 3+ patches meet is resolved by this exact same logic, with
         no special-casing, since every proposal naming a given target funnels into the same
         "group by target id" step regardless of how many different ranks it came from.
      6. A_maxCellDistance and a_isPositionValid apply here exactly as they do in
         resolveTrivialTier() -- a physical-validity rule, not a coordination one. An incoming
         proposal that is too far, or whose weighted-centroid position fails a_isPositionValid, is
         never eligible to be accepted (skipped as though it did not exist for argmin purposes),
         REGARDLESS of the mutual-match exception in rule 3 -- an invalid position is invalid on
         both sides of a pair, mutual match or not.

    Proposals addressed to this rank's own particles arrive from two sources, pooled identically
    here: proposals received over MPI from other ranks, and this rank's own
    generateProposals()-produced bucket addressed to itself (an in-process handoff, no network
    round-trip -- see generateProposals()).

    @tparam Packed  See NNMergeParticle.
    @tparam Combine See PayloadCombine.
    @tparam PosValid See PositionValid.
    @param[in]     a_incomingProposals All proposals targeting a particle this rank owns, pooled
                                       from every source rank (including this rank's own
                                       in-process bucket).
    @param[in]     a_myParticlesByID   This rank's own currently-alive valid particles, keyed by
                                       id (i.e. potential judge targets).
    @param[in]     a_hasOutgoingCommitment Ids (among a_myParticlesByID) that have their own
                                       outgoing commitment this round (needed for rule 3 above).
    @param[in]     a_outgoingTargetOf  For an id in a_hasOutgoingCommitment, the id of the
                                       particle it is itself proposing to (needed to detect the
                                       mutual-match exception in rule 3).
    @param[in]     a_combine           Opaque-payload combine callback (see PayloadCombine).
    @param[in,out] a_liveCellCount     Live per-cell counts for this rank's own patches;
                                       decremented in place as targets are accepted, in the same
                                       distance-sorted processing order (see rule 4 above).
    @param[in]     a_numParticlesPerCellThresh Crowding threshold N (same value used everywhere
                                       else this round).
    @param[in]     a_probLo            Physical position of the domain's low corner (see
                                       cellKeyOf()).
    @param[in]     a_dxByLevel         Cell spacing PER AMR LEVEL -- see resolveTrivialTier()'s
                                       a_dxByLevel docs for why this must never be a single shared
                                       dx (a proposal's source and this rank's target can be on
                                       different levels).
    @param[in]     a_isPositionValid   Position-validity predicate (see PositionValid and rule 6
                                       above).
    @param[in]     a_allocateID        Fresh-id allocator (see IDAllocator), called once per
                                       accepted target.
    @param[in]     a_maxCellDistance   Same physical merge-eligibility cap as
                                       resolveTrivialTier()'s (see rule 6 above); std::nullopt
                                       (default) disables it.
    @param[out]    a_results           One entry per target actually accepted, in the order
                                       processed.
    @param[out]    a_verdictsByDestRank One bucket per proposer rank (pre-sized to the total rank
                                       count by the caller); every proposal (accepted or not) gets
                                       exactly one verdict, addressed back to
                                       NNMergeProposal::source's owning rank.
    @param[in,out] a_crossLevelMergeCount Optional diagnostic counter -- see
                                       resolveTrivialTier()'s docs; incremented once per accepted
                                       target whose target/source were on different AMR levels.
  */
  template <typename Packed, typename Combine, typename PosValid, typename Allocator>
  inline void
  judgeProposals(const std::vector<NNMergeProposal<Packed>>&          a_incomingProposals,
                 const std::map<ParticleID, NNMergeParticle<Packed>>& a_myParticlesByID,
                 const std::set<ParticleID>&                          a_hasOutgoingCommitment,
                 const std::map<ParticleID, ParticleID>&              a_outgoingTargetOf,
                 const Combine&                                       a_combine,
                 std::map<NNCellKey, int>&                            a_liveCellCount,
                 const int                                            a_numParticlesPerCellThresh,
                 const RealVect&                                      a_probLo,
                 const std::vector<RealVect>&                         a_dxByLevel,
                 const PosValid&                                      a_isPositionValid,
                 const Allocator&                                     a_allocateID,
                 std::vector<NNMergeResult<Packed>>&                  a_results,
                 std::vector<std::vector<NNMergeVerdict>>&            a_verdictsByDestRank,
                 const std::optional<int>                             a_maxCellDistance      = std::nullopt,
                 unsigned long long*                                  a_crossLevelMergeCount = nullptr) noexcept;

  /**
    @brief Apply every verdict addressed to one of this rank's own outgoing proposals: remove
    the source particle for every accepted proposal, do nothing for rejected ones.
    @details Rejected proposals are deliberately a no-op -- the source particle is simply left
    alone and reconsidered next round's fresh candidate search (see RULES.md, "Single pass,
    caller loops for convergence": this algorithm never internally retries a rejected proposal
    against a second-best candidate within the same round).
    @param[in]  a_incomingVerdicts All verdicts addressed to this rank, pooled from every judge
                                   rank (including this rank's own in-process verdict for
                                   proposals it judged for itself).
    @param[out] a_consumedIDs      Ids of every source particle whose proposal was accepted --
                                   the caller removes these from its own valid particle lists.
  */
  inline void
  applyVerdicts(const std::vector<NNMergeVerdict>& a_incomingVerdicts, std::set<ParticleID>& a_consumedIDs) noexcept;

  /**
    @brief Result of locating a physical position against the AMR grid -- see placeMergedParticles().
    @details Deliberately independent of any AmrMesh/Realm/ParticleContainer type (this header
    stays Amr-agnostic, see RULES.md) -- mirrors ParticleContainer::Destination's fields exactly,
    since that IS the intended real implementation: an O(1) point->block query on the realm's
    LevelTiles tile->box hash (the SAME lookup a ParticleContainer already uses internally to
    route its own particles during remap()), never a linear scan over patches. Reusing
    ParticleContainer::findDestination() rather than adding a second, parallel point-location
    method keeps this one mapping the single source of truth for "which box owns this position".
  */
  struct NNParticleLocation
  {
    /**
      @brief True iff a covering box was found (false near/outside the domain boundary).
    */
    bool valid;

    /**
      @brief AMR level containing the point.
    */
    int level;

    /**
      @brief Global grid index of the covering box within that level (NOT a local DataIndex --
      the caller is expected to convert this to whatever local indexing it needs, e.g. via
      LevelTiles::getMyGrids(), only when ownerRank is this rank).
    */
    unsigned int gridIndex;

    /**
      @brief The rank that owns the covering box.
    */
    int ownerRank;
  };

  /**
    @brief Locate-by-position callback -- see NNParticleLocation and placeMergedParticles().
    @details Deliberately a template parameter wherever consumed, NOT required to be exactly this
    std::function type -- this alias exists purely to document the expected signature (mirrors
    PayloadCombine's role). Unlike PayloadCombine, call frequency here is bounded by the number of
    MERGE RESULTS this round (not by the number of candidate pairs squared), so a non-inlined
    indirect call is an acceptable cost -- this is not a hot inner-loop path.
    @param[in] a_position Physical position to locate.
    @return Which box (level, global grid index) and rank covers it.
  */
  using PositionLocator = std::function<NNParticleLocation(const RealVect& a_position)>;

  /**
    @brief Point-locate a newly created merged particle against the AMR grid, and bucket it for a
    scatter exchange if it belongs to someone else.
    @details Applied UNIFORMLY to every merge result, trivial-tier and judge-decided alike -- see
    RULES.md, "Placement": the union of two patches is not generally convex (T-junctions,
    corners, refinement boundaries), so even a purely local (same-rank) merge's result can
    legitimately belong to a third patch, possibly owned by a different rank. Point-location must
    be O(1) per result (see PositionLocator/NNParticleLocation) -- NEVER a linear scan over every
    patch, which would turn this function into a real bottleneck at large patch counts.
    @tparam Packed See NNMergeParticle.
    @param[in]  a_results        Every merge result produced this round by resolveTrivialTier()
                                  and judgeProposals(), pooled together.
    @param[in]  a_locate         Locate-by-position callback (see PositionLocator).
    @param[in]  a_thisRank       This rank's own number -- an NNParticleLocation whose ownerRank
                                  equals this is what makes a result "local" rather than needing a
                                  scatter.
    @param[out] a_localByPatch   Results whose merged particle lands in one of this rank's own
                                  patches, bucketed by (level, NNParticleLocation::gridIndex).
    @param[out] a_scatterByDestRank Results whose merged particle belongs to a different rank,
                                  bucketed by destination rank (pre-sized to the total rank count
                                  by the caller), for one final exchange.
  */
  template <typename Packed>
  inline void
  placeMergedParticles(const std::vector<NNMergeResult<Packed>>&                                     a_results,
                       const PositionLocator&                                                        a_locate,
                       const int                                                                     a_thisRank,
                       std::map<std::pair<int, unsigned int>, std::vector<NNMergeParticle<Packed>>>& a_localByPatch,
                       std::vector<std::vector<NNMergeParticle<Packed>>>& a_scatterByDestRank) noexcept;

  /**
    @brief Run one full round (one timestep's worth) of the nearest-neighbor merge algorithm.
    @details Orchestrates, in order (see RULES.md for the full rationale behind this exact
    sequencing):
      1. Ghost-fill -- fresh, exactly once, strictly before any of the below. Uses
         AmrMesh::getParticleGhostMask()/getParticleGhostMaskCoarToFine()/
         getParticleGhostMaskFineToCoar() and ParticleContainer::fillGhostParticles(), matching
         the ALREADY-EXISTING ghost-fill machinery -- this function does not reimplement ghost
         transport, it only requires that a ghost width of (at least) 1 has ALREADY been
         registered on this realm via AmrMesh::registerParticleGhostMask(), and does NOT call
         registerParticleGhostMask() itself. @b Timing matters here: a particle ghost mask is
         only actually built during a regrid (registerParticleGhostMask() itself just records
         that a width should be tracked from then on) -- registering it from inside this
         function, which typically runs long after the grids were last (re)built, would silently
         leave the mask undefined/empty until the next regrid happens to occur. The caller must
         register the width EARLY, before the relevant grids are (re)built (e.g. from a
         TimeStepper's registerOperators(), which runs before initial grid setup) -- see the
         ItoKMC-specific merger config for a concrete example.
      2. Per patch: findNearestNeighborCandidates(), then resolveTrivialTier(). If
         a_iterateLocalTierToConvergence is true, this step (candidate search + trivial-tier
         resolution, using the SAME ghost-fill from step 1 -- no new ghost data is needed since
         nothing about local re-matching requires it) repeats until no further local commits
         occur, rather than running exactly once. See a_iterateLocalTierToConvergence's own
         documentation below for why this is worth having as a caller-selected option rather
         than baked in one way or the other.
      3. generateProposals() for everything the trivial tier could not resolve.
      4. Exchange proposals: this rank's own bucket (addressed to itself) is handed directly to
         judgeProposals() in-process; every other bucket must be shipped to its destination rank
         and the corresponding incoming proposals collected. @b Implementation note: no existing
         general-purpose "send arbitrary payloads to a specific destination rank and collect what
         arrives" primitive is currently exposed by this codebase outside
         ParticleContainer's OWN protected gatherToPool()/gatherGhostsFromMasks()/
         distributeFromPool() machinery (which operates on ParticleSoA/MoverPool specifically,
         not on arbitrary NNMergeProposal records). Implementing this exchange will require
         either a new general-purpose alltoallv-style utility, or representing a proposal as a
         specially-tagged particle so it can ride the existing ParticleContainer transport --
         this decision is left open for the implementation phase.
      5. judgeProposals() over the pooled incoming proposals.
      6. Exchange verdicts (same open transport question as step 4, mirrored in the other
         direction) and applyVerdicts().
      7. placeMergedParticles(), then exchange and apply the scatter bucket (same open transport
         question again).

    Fresh, globally-unique ids for newly created particles must be allocated by whatever calls
    this function such that they can NEVER collide with an existing particle's id or another
    rank's freshly allocated id within the same round -- the prototype used a rank-namespaced
    counter (see Scripts/ParticleMergePrototype/merge_lib.py's make_particle()); the real
    implementation should use whatever global id allocation scheme the rest of the particle
    infrastructure already relies on, if one exists, rather than inventing a second one.

    This function does a SINGLE cross-patch exchange (propose/judge/verdict happens exactly
    once) regardless of a_iterateLocalTierToConvergence -- see RULES.md, "Single pass, caller
    loops for convergence": a cell more than roughly 2x over threshold is not necessarily fully
    drained by the cross-patch protocol in one call, and the caller is still expected to invoke
    this once per timestep, relying on repeated calls (and, in a real simulation, particle motion
    between them) for further convergence of anything the cross-patch protocol itself left
    unresolved. What a_iterateLocalTierToConvergence changes is only the LOCAL (trivial-tier)
    portion of a single call -- see its own documentation below.

    @tparam P       SoA payload type (e.g. ItoParticle).
    @tparam Packed  Opaque reconciliation payload type for the merge decision layer (see
                    NNMergeParticle) -- typically much smaller than P's full column set (e.g.
                    just energy for ItoParticle; velocity/mobility/diffusion are already
                    established as safe to drop across a merge, since they are recomputed from
                    field interpolation immediately afterward).
    @tparam Traits  Column descriptor for P (defaults to ParticleTraits<P>).
    @tparam Gather  Callable `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed` -- packs one
                    SoA slot's payload columns into the opaque Packed representation.
    @tparam Combine See PayloadCombine.
    @tparam Scatter Callable `(ParticleSoA<P,Traits>&, const NNMergeParticle<Packed>&) -> void` --
                    appends one surviving/merged particle back into the real SoA container.
    @param[in,out] a_amr                       AMR mesh -- source of ghost masks and (eventually)
                                                the MPI communicator/patch layout used for the
                                                open transport questions noted above.
    @param[in,out] a_particles                 The real particle container to merge, for every
                                                patch this rank owns.
    @param[in]     a_numParticlesPerCellThresh Crowding threshold N.
    @param[in]     a_gather                    Per-slot gather callback (see Gather).
    @param[in]     a_combine                   Opaque-payload combine callback (see
                                                PayloadCombine).
    @param[in]     a_scatter                   Per-surviving-particle scatter callback (see
                                                Scatter).
    @param[in]     a_iterateLocalTierToConvergence User-selectable option, default false (matches
                                                the originally documented, simplest behavior).
                                                When false: findNearestNeighborCandidates() and
                                                resolveTrivialTier() each run exactly once, like
                                                every other step in this function. When true: that
                                                candidate-search-plus-trivial-tier pair repeats,
                                                using the SAME ghost-fill from step 1 (no new
                                                ghost data is needed -- only which local particles
                                                are still alive has changed), until a full local
                                                pass commits zero further trivial merges, and only
                                                THEN does the function proceed to
                                                generateProposals()/judgeProposals() (still exactly
                                                once, unaffected by this option).

                                                Why this option exists: a single
                                                findNearestNeighborCandidates() pass fixes each
                                                query particle's nearest-neighbor edge once, at
                                                the start of the pass -- if that candidate is
                                                consumed by a different, closer edge before this
                                                edge is reached, the query particle is left
                                                unmatched for the rest of that pass, even though a
                                                different, still-alive partner might exist nearby.
                                                Measured in the Python prototype
                                                (Scripts/ParticleMergePrototype/): a single
                                                isolated, non-exposed cell seeded at 64 particles
                                                against a threshold of 16 (so trivial tier alone
                                                has everything it needs to drain the cell fully,
                                                with no cross-patch coordination required at all)
                                                completed only ~12 merges in one pass -- far below
                                                the 24 that draining 64 down to 16 requires, and
                                                also below the 32-merge hard ceiling one
                                                non-iterated pass can ever reach (each participant
                                                gets only one query per pass, capping any single
                                                pass at n/2 merges regardless of match quality).
                                                @b Measured effect turned out smaller than the
                                                reasoning above predicts: re-running the search
                                                only helps a particle whose BEST available choice
                                                actually changes between iterations, and empirically
                                                that's rare -- in the same 64-vs-16 test, enabling
                                                this option recovered only roughly 5% more merges,
                                                converging in just 3-4 internal iterations. The
                                                dominant loss turned out to be a different, more
                                                stable effect that iteration alone cannot fix -- see
                                                a_maxFallbackCandidates below, which was
                                                specifically introduced to address it and recovers
                                                far more. This option is still worth having (it is
                                                a real, if modest, free improvement, and combines
                                                with a_maxFallbackCandidates rather than competing
                                                with it), just not for the reason originally
                                                hypothesized. See RULES.md, "Single pass, caller
                                                loops for convergence" for the full empirical
                                                writeup and the tradeoff (extra local re-search
                                                compute per call, in exchange for needing fewer
                                                external timesteps to clear a severe local
                                                backlog).
    @param[in]     a_maxFallbackCandidates     Forwarded directly to resolveTrivialTier() -- see
                                                its own documentation for the full mechanism and
                                                rationale. Default 0 (disabled, matches the
                                                originally documented behavior). This is the
                                                option that actually closes most of the gap
                                                a_iterateLocalTierToConvergence was originally
                                                (and only partially correctly) hypothesized to
                                                close: instrumentation traced the dominant loss to
                                                a particle being blocked from a perfectly good
                                                local partner simply because that partner's OWN,
                                                separate merge attempt had already failed and
                                                marked it "busy" -- a stable condition that does
                                                NOT change between iterations of
                                                a_iterateLocalTierToConvergence, since the busy
                                                particle's own true nearest neighbor doesn't
                                                change either. Giving the BLOCKED particle a
                                                fallback candidate instead directly addresses this;
                                                a single fallback candidate recovered the large
                                                majority of the achievable improvement in testing,
                                                with steeply diminishing returns beyond that.
    @param[in]     a_maxCellDistance           Forwarded directly to resolveTrivialTier() -- see
                                                its own documentation for the full mechanism.
                                                `std::nullopt` (default) disables this and matches
                                                the original, unrestricted-search behavior. A
                                                physical-validity cap, in whole grid cells, on how
                                                far apart two particles may be and still merge --
                                                e.g. a value of 1 permits merging within the same
                                                cell or any immediately (Moore-neighborhood)
                                                adjacent cell, but never farther. Independent of,
                                                and composes with, both options above.
    @param[in]     a_isPositionValid           Forwarded directly to resolveTrivialTier() and
                                                judgeProposals() -- see PositionValid. Default is
                                                an always-true predicate (matches the originally
                                                documented, unrestricted behavior). Independent of,
                                                and composes with, every option above -- e.g. an
                                                embedded-boundary point query for an EB-aware
                                                caller (see the ItoKMC-specific merger config for a
                                                concrete example).
    @param[in]     a_allocateID                Fresh-id allocator (see IDAllocator). Mandatory --
                                                there is no universally safe default; the caller
                                                must supply a scheme that stays collision-free
                                                across every rank and every past/future call.
    @param[in,out] a_crossLevelMergeCount      Optional diagnostic counter (nullptr, the default,
                                                disables it entirely): incremented once for every
                                                merge this call commits (trivial-tier or judge-side)
                                                whose two participants were on DIFFERENT AMR levels.
                                                NOT reset by this function -- accumulates across
                                                a_iterateLocalTierToConvergence's internal
                                                iterations and across repeated calls if the caller
                                                keeps reusing the same counter, so it is safe to use
                                                to track cross-level activity over several rounds.
  */
  template <typename P,
            typename Packed,
            typename Traits = ParticleTraits<P>,
            typename Gather,
            typename Combine,
            typename Scatter,
            typename Allocator,
            typename PosValid = std::function<bool(const RealVect&)>>
  inline void
  mergeNearestNeighborsRound(
    AmrMesh&                      a_amr,
    ParticleContainer<P, Traits>& a_particles,
    const int                     a_numParticlesPerCellThresh,
    const Gather&                 a_gather,
    const Combine&                a_combine,
    const Scatter&                a_scatter,
    const Allocator&              a_allocateID,
    const bool                    a_iterateLocalTierToConvergence = false,
    const int                     a_maxFallbackCandidates         = 0,
    const std::optional<int>      a_maxCellDistance               = std::nullopt,
    const PosValid&               a_isPositionValid =
      [](const RealVect&) {
        return true;
      },
    unsigned long long* a_crossLevelMergeCount = nullptr);

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#include <CD_NearestNeighborParticleMergeImplem.H>

#endif
