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
#include <limits>
#include <map>
#include <optional>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

// Chombo includes
#include <RealVect.H>
#include <IntVect.H>

// EBGeometry includes
#include <EBGeometry.hpp>

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_ParticleContainer.H>
#include <CD_ParticleGhostMask.H>
#include <CD_NamespaceHeader.H>

class AmrMesh;

/**
  @brief Namespace for various particle management tools.
  @details This file adds a nearest-neighbor super-particle merge that is safe across patch and rank
  boundaries. Candidates are found against a particle's full neighborhood (its own patch plus every
  ghost, filled once per round); a pair where neither participant is reachable from outside the patch
  resolves immediately (the "trivial tier"), and everything else goes through a propose/judge/verdict
  protocol that is safe for any number of simultaneous cross-patch suitors.
*/
namespace ParticleManagement {

  /**
    @brief Minimal, payload-agnostic description of one particle as input to the nearest-neighbor
    merge.
    @details Carries what the matching/decision logic needs -- position, weight, a globally unique
    id, the particle's current owning rank, and its AMR level -- plus an opaque, caller-defined
    reconciliation payload. The matching/decision logic never inspects @c payload; it is unpacked
    only by the caller-supplied combine callback when a pair is finalized.
    @tparam Packed Caller-defined payload type (e.g. a single Real energy for ItoParticle, or an
    empty struct if no type-specific reconciliation data is required). Must be copyable and
    trivially linearizable (NNMergeProposal ships two of these as one atomic pair).
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
      (see judgeProposals()).
    */
    ParticleID globalID;

    /**
      @brief The rank that currently owns this particle. Always read from the particle's own data --
      the single source of truth, never inferred from MPI transport/source-rank metadata.
    */
    RankID ownerRank;

    /**
      @brief The AMR level this particle currently lives on.
      @details Required for multi-level correctness: cellKeyOf() needs the dx of the level the
      particle lives on, and a ghost candidate can be on a different level than the querying patch
      (GhostType::Coarse/Fine, not just SameLevel). Every cell-key computation uses the per-level dx
      for this field, never a single global dx (see a_dxByLevel in resolveTrivialTier()/
      judgeProposals()).
    */
    int level;

    /**
      @brief Opaque, caller-defined reconciliation payload. Never inspected by this algorithm's
      matching/decision logic -- only carried through to the caller-supplied combine callback.
    */
    Packed payload;
  };

  namespace detail {

    /**
      @brief One additional ranked candidate beyond a query's primary nearest neighbor.
      @details Populated only when a_maxFallbackCandidates > 0 is requested from
      findNearestNeighborCandidates() (see NNMergeEdge::fallbackCandidates). Fetched via a single
      K-nearest query at generation time (K = 1 + a_maxFallbackCandidates); resolveTrivialTier() walks
      this pre-fetched, distance-sorted list when a closer candidate is blocked, without re-querying
      the spatial index.
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
      @details Produced once per over-threshold, unconsumed valid particle by
      findNearestNeighborCandidates() -- exactly one edge per query particle, its single nearest
      neighbor among {this patch's own other valid particles} union {every ghost}. Distance is carried
      explicitly (never recomputed downstream) so trivial-tier resolution and the judge apply the same
      "closest-first, with a deterministic id tiebreak on exact ties" ordering.
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
        initiate a query).
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
        eligibility -- see isBoundaryExposed(): a same-patch candidate can still require the
        propose/judge/verdict protocol if either participant is boundary-exposed.
      */
      bool candidateIsLocal;

      /**
        @brief Next-nearest candidates beyond candidateID, in ascending distanceSquared order.
        @details Empty unless a_maxFallbackCandidates > 0 was passed to
        findNearestNeighborCandidates(); size at most a_maxFallbackCandidates. resolveTrivialTier()
        tries these in order when candidateID is blocked (stale, busy, exposed, a ghost, or too far),
        without re-querying the spatial index. Content is a generation-time snapshot; a candidate's
        validity is rechecked against a_consumedIDs/a_hasOutgoingCommitment when it is actually tried.
      */
      std::vector<NNFallbackCandidate> fallbackCandidates;
    };

    /**
      @brief A cross-patch merge proposal: one particle's own nearest-neighbor query result,
      bundled as a self-contained, atomic pair for shipment to whichever rank owns the candidate.
      @details Carries full data for both participants (not a bare id reference) so the judge's
      decision is self-contained and it can never recombine pieces across different incoming
      proposals. Routing rule (see generateProposals()): a proposal is always sent to the owner of the
      querying particle's best candidate, never by any id comparison.
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
      @details Accepted implies the judge has already created the merged particle and removed its own
      (target) particle; the proposer, on receiving an accepted verdict, removes its own (source)
      particle. Rejected implies no action on either side -- the proposer's particle is simply
      reconsidered next round.
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

  } // namespace detail

  /**
    @brief Illustrative signature of the combine callback: produce a merged particle's opaque
    payload from its two parents' payloads and weights.
    @details This alias documents the expected callback shape only; it is NOT the type used at the
    call sites. resolveTrivialTier()/judgeProposals() take the combine callback as a template
    parameter (not this std::function) so it inlines on the per-merge hot path. Position and weight
    are handled by the algorithm directly (weighted centroid, summed weight); this callback produces
    only the opaque @c Packed payload (e.g. weighted-average energy for ItoParticle).
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
    @brief Illustrative signature of the position-validity predicate for a would-be merged particle.
    @details Documents the expected callback shape; consumed as a template parameter (not this
    std::function) at the call sites, for the same inlining reason as PayloadCombine. Called with the
    weighted-centroid position already computed, once per candidate pair that is otherwise eligible
    to merge (never on a pair already rejected for another reason), so a caller doing real geometric
    work here pays only when it matters. Default is always-true. A pair that fails is treated like
    a_maxCellDistance's "too far" case: not proposed, a fallback tried if available, otherwise left
    unmatched with no side effects.
    @param[in] a_mergedPosition The candidate merge's weighted-centroid position.
    @return True iff a particle may exist at this position (e.g. false inside an embedded boundary).
  */
  using PositionValid = std::function<bool(const RealVect& a_mergedPosition)>;

  /**
    @brief Illustrative signature of the fresh-id allocator: produce one globally-unique id for a
    newly created merged particle.
    @details Documents the expected callback shape; consumed as a template parameter (not this
    std::function). Called once per committed merge (trivial-tier or judge-side). The caller must
    ensure ids never collide -- not just within one call, but across every rank and every
    past/future call for the life of the simulation; a rank-namespaced counter (this rank's index
    times a stride plus a monotonically increasing local count) suffices.
    @return A fresh globally-unique particle id.
  */
  using IDAllocator = std::function<ParticleID()>;

  /**
    @brief Optional per-merge diagnostic callback: observe every individual merge as it commits.
    @details Consumed as a template parameter (not required to be exactly this std::function).
    Called once per committed merge (trivial-tier or judge-side) with both parents' full pre-merge
    data and the resulting merged particle (finalized but not yet placed; see
    placeMergedParticles()). Default is a no-op. Intended for invariant checks, e.g. recomputing the
    weighted centroid from the parents and comparing against a_merged.position.
    @param[in] a_first  First parent (pre-merge, full data).
    @param[in] a_second Second parent (pre-merge, full data).
    @param[in] a_merged The resulting merged particle.
  */
  template <typename Packed>
  using MergeCommittedCallback = std::function<void(const NNMergeParticle<Packed>& a_first,
                                                    const NNMergeParticle<Packed>& a_second,
                                                    const NNMergeParticle<Packed>& a_merged)>;

  namespace detail {

    /**
      @brief The unclamped, position-derived cell index a given physical position falls in.
      @details The cell-key convention used throughout for the crowding trigger and a_maxCellDistance.
      Not clamped to any patch's box, so a ghost buckets into its true origin cell and particles in
      different patches' local frames stay directly comparable.
      @warning a_dx MUST be the cell spacing of the AMR level the position's particle lives on (see
      NNMergeParticle::level), never a single fixed dx shared across levels. The returned IntVect is
      only meaningful together with that level; NNCellKey pairs the two so they cannot be conflated
      (pooling bare cell keys from different levels aliases coarse- and fine-level cell counts).
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
      @brief Hash functor for NNCellKey -- REQUIRED whenever NNCellKey is used as a hash-map key;
      NNCellKey must never be used as an ordered (%std::map or %std::set) key.
      @details IntVect::operator< (which std::pair<int,IntVect>::operator< calls for the second
      component) is a component-wise "strictly dominates in every direction" partial order, not a
      strict weak ordering. Two unrelated cells where neither dominates the other (e.g. (76,104) vs
      (80,101)) compare as neither-less-than, which %std::map/%std::set treat as equal keys -- pooling
      unrelated cells' counts into one entry. Hashing depends only on IntVect::operator== (full
      component-wise equality), which is correct.
    */
    struct NNCellKeyHasher
    {
      /**
        @brief Hash operator.
        @param[in] a_key Cell key to hash.
        @return Hash value.
      */
      std::size_t
      operator()(const NNCellKey& a_key) const noexcept
      {
        std::size_t hash = LevelTiles::TileHasher()(a_key.second);
        hash ^= static_cast<std::size_t>(static_cast<unsigned int>(a_key.first)) + 0x9e3779b9 + (hash << 6) +
                (hash >> 2);
        return hash;
      }
    };

    /**
      @brief A non-level-qualified cell index, distinct from NNCellKey, used where a cell is expressed
      in a single consistent frame (the query patch's own level).
      @details A separate alias from NNCellKey (which is level-qualified for the cross-level
      crowding-count map) so the two cannot be passed interchangeably. NNWalkCell is always computed
      with the query patch's own level's dx, even for a cross-level ghost, so candidates are
      geometrically comparable to the query in one frame (see cellKeyOf()'s warning).
    */
    using NNWalkCell = IntVect;

    /**
      @brief A patch's spatial index for one round: an EBGeometry SIMD PointCloudBVH over the whole
      patch (every local valid particle followed by every ghost) plus its precomputed k-nearest-
      neighbor graph.
      @details Built once per patch per mergeNearestNeighborsRound() call via buildNNSpatialIndex()
      and reused across every internal a_iterateLocalTierToConvergence pass via lazy deletion (see
      findNearestNeighborCandidates()'s own docs). A patch with no crowded cell gets no index at all
      (pointCloud stays null -- see buildNNSpatialIndex()'s per-cell pre-filter). The point cloud's
      metadata is each point's globalID; a hit's cloud index c decodes ownership directly:
      c < pointCloudNLocal => local valid a_localValid[c], else ghost a_ghosts[c - pointCloudNLocal].
    */
    template <typename Packed>
    struct NNSpatialIndex
    {
      /**
        @brief The SIMD PointCloudBVH over this patch's whole point set (local valid first, then
        ghosts). Null iff there were no particles at all, or no cell in this patch could possibly be
        crowded (see buildNNSpatialIndex()'s per-cell pre-filter).
      */
      std::shared_ptr<EBGeometry::PointCloudBVH<Real, ParticleID>> pointCloud;

      /**
        @brief Number of local-valid points at the front of the PointCloudBVH's cloud order (used to
        decode a hit's local-vs-ghost ownership from its cloud index). Only meaningful when
        pointCloud is non-null.
      */
      std::size_t pointCloudNLocal = 0;

      /**
        @brief The whole k-nearest-neighbor graph, precomputed once via PointCloudBVH::
        allNearestNeighbors(pointCloudK) at build time (batched, leaf-order / cache-warm, seeded from
        each point's own leaf -- strictly cheaper than per-query calls). Row-major: entry
        [c * pointCloudK + j] is the j-th nearest OTHER point of cloud index c. Only the local-valid
        rows (c < pointCloudNLocal) are ever read by findNearestNeighborCandidates(); ghost rows are
        computed but unused. Empty iff pointCloud is null.
      */
      std::vector<typename EBGeometry::PointCloudBVH<Real, ParticleID>::Hit> pointCloudGraph;

      /**
        @brief k used to build pointCloudGraph (= 1 + a_maxFallbackCandidates); the row stride.
      */
      int pointCloudK = 0;
    };

    /**
      @brief Build a patch's spatial index for one round: an EBGeometry SIMD PointCloudBVH over the
      whole patch (every local valid particle then every ghost) plus its precomputed k-nearest-
      neighbor graph -- UNLESS no cell in this patch could possibly be crowded, in which case nothing
      is built.
      @details The point cloud carries each point's globalID as metadata; a hit's cloud index decodes
      ownership (c < nLocal => local valid, else ghost). After the cloud is built, its whole
      a_kNearest-nearest-neighbor graph is precomputed once via PointCloudBVH::allNearestNeighbors()
      (batched, leaf-order/cache-warm, own-leaf-seeded) into a_index.pointCloudGraph and reused across
      every internal a_iterateLocalTierToConvergence pass (positions never change within a round; only
      consumed status does, and that is filtered at lookup time -- see
      findNearestNeighborCandidates()).

      @b Per-cell pre-filter: if NO cell in this patch has more than a_numParticlesPerCellThresh
      particles, nothing is built (a_index.pointCloud stays null, exactly as if the patch had zero
      particles) -- the SAME cheap "is this even worth it" check equal_weight_kd effectively gets for
      free by dispatching per cell. This is an EXACT check (not an average-density heuristic, which
      would be a real bug: a low-average patch can still hold one badly-clustered over-threshold cell):
      a lightweight per-cell tally of the patch's own particles, keyed exactly as the crowding trigger
      (cellKeyOf(pos, probLo, dx), local and ghost alike), decides it directly.
      @tparam Packed See NNMergeParticle.
      @param[in]  a_localValid                This patch's own full valid-particle set.
      @param[in]  a_ghosts                    Every ghost shipped to this patch this round.
      @param[in]  a_probLo                    Physical position of the domain's low corner.
      @param[in]  a_dx                        Cell spacing of THIS PATCH's own level.
      @param[in]  a_numParticlesPerCellThresh The physics-facing crowding threshold used elsewhere
                                              this round (see mergeNearestNeighborsRound()'s own
                                              a_numParticlesPerCellThresh) -- used ONLY for the
                                              per-cell pre-filter above.
      @param[in]  a_kNearest                  k = 1 + a_maxFallbackCandidates -- the neighbor count to
                                              precompute per point into NNSpatialIndex::pointCloudGraph.
      @param[out] a_index                     The resulting spatial index.
    */
    template <typename Packed>
    inline void
    buildNNSpatialIndex(const std::vector<NNMergeParticle<Packed>>& a_localValid,
                        const std::vector<NNMergeParticle<Packed>>& a_ghosts,
                        const RealVect&                             a_probLo,
                        const RealVect&                             a_dx,
                        const int                                   a_numParticlesPerCellThresh,
                        const int                                   a_kNearest,
                        NNSpatialIndex<Packed>&                     a_index) noexcept;

    /**
      @brief Whether a particle is boundary-exposed -- i.e. would be shipped as a ghost to at
      least one neighboring patch.
      @details The eligibility test for the trivial tier (see resolveTrivialTier()). A direct lookup
      against the patch's own outgoing ghost mask -- whether this particle's cell has any registered
      destination (@c ParticleGhostMask::numTargets()) -- not "same rank" and not "same ParticleSoA": a
      particle can be exposed to a different patch owned by the same rank, and can be exposed while its
      candidate sits in the same ParticleSoA.
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
      @details This is the only place a nearest-neighbor search happens, and it always runs against
      the full (local + ghost) view, never a local-only view: restricting it biases boundary-adjacent
      particles toward worse matches than interior ones. Crowding eligibility (is the query's own cell
      over threshold) is checked here at generation time from the live per-cell counts over valid +
      ghost occupancy; the symmetric check on the candidate's cell is done later, dynamically, by
      resolveTrivialTier()/judgeProposals(), since earlier merges this round can change it.

      a_localValid/a_ghosts are this patch's full, unfiltered sets for the round; this function does the
      alive-filtering itself via a_consumedIDs (checked against both the query and every candidate).
      a_spatialIndex is built once per patch per mergeNearestNeighborsRound() call over the full set and
      reused across every a_iterateLocalTierToConvergence pass via this same lazy-deletion check, so it
      is never rebuilt per pass.

      Each query's precomputed row is read from a_spatialIndex's k-nearest-neighbor graph (see
      NNSpatialIndex::pointCloudGraph / buildNNSpatialIndex()), computed once for the whole patch via
      the EBGeometry SIMD PointCloudBVH's batched allNearestNeighbors(). Per-query work is a row lookup
      plus filtering: a hit is skipped if its id is in a_consumedIDs, or (when a_maxCellDistance is set)
      if it lies beyond the bound; survivors feed the shared top-K logic that builds the edge. Self is
      excluded by the graph (each row is the nearest OTHER points).

      a_maxCellDistance, when set, filters the graph row: rows are distance-sorted, so an out-of-bound
      candidate is farther than the bound and can never displace a within-bound one, and filtering
      yields exactly the within-bound k-nearest (see its param docs for the exact, safe-but-not-tight
      bound).

      When a_maxFallbackCandidates > 0, the graph is built with k = 1 + a_maxFallbackCandidates, so each
      row holds the K nearest and NNMergeEdge::fallbackCandidates is populated with candidates 2..K
      (distance-sorted); resolveTrivialTier() walks that pre-fetched list when the primary is blocked,
      without re-querying the index.

      @tparam Packed See NNMergeParticle.
      @param[in]  a_localValid       This patch's own FULL valid-particle set for this round (query
                                      sources AND candidates) -- see the contract-change note above;
                                      not pre-filtered by a_consumedIDs.
      @param[in]  a_ghosts           Every ghost shipped to this patch this round (candidates
                                      only; ghosts never initiate a query) -- also not pre-filtered.
      @param[in]  a_consumedIDs      Lazy-deletion skip set: ids already consumed earlier this
                                      round (by an earlier a_iterateLocalTierToConvergence pass, or
                                      an earlier query within this same pass). Checked against both
                                      the query and every candidate encountered during the walk.
      @param[in]  a_spatialIndex     This patch's spatial index for this round (see NNSpatialIndex),
                                      built once via buildNNSpatialIndex() over the FULL a_localValid/
                                      a_ghosts sets before any consumption happens.
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
      @param[out] a_edges            One entry per over-threshold, currently-alive query particle
                                      that has at least one candidate (empty if none found).
      @param[in]  a_maxFallbackCandidates Forwarded from resolveTrivialTier() (see its docs); 0
                                      (default) means an ordinary single-nearest query and every
                                      edge's fallbackCandidates stays empty.
      @param[in]  a_maxCellDistance  Same hard cutoff as resolveTrivialTier()'s own a_maxCellDistance
                                      (see its docs); std::nullopt (default) disables it. Applied here
                                      purely as an additional, always-safe filter on the graph row --
                                      resolveTrivialTier()/judgeProposals() still apply the exact,
                                      level-aware, per-axis Chebyshev check (nnMergeCrossLevelTooFar())
                                      downstream regardless, unchanged. The bound used here is
                                      deliberately looser: a candidate is dropped only once its
                                      EUCLIDEAN distance exceeds a_maxCellDistance * min(a_dx) *
                                      sqrt(SpaceDim) -- the sqrt(SpaceDim) factor is required because a
                                      per-axis Chebyshev cell distance of N does not bound Euclidean
                                      physical distance by N * dx alone (a diagonal candidate exactly N
                                      cells away in both directions is N * dx * sqrt(SpaceDim) away in a
                                      straight line), and this bound must never be tighter than the
                                      true Chebyshev cutoff or it could wrongly drop a candidate
                                      resolveTrivialTier()/judgeProposals() would still have accepted.
                                      Safe but not maximally tight.
    */
    template <typename Packed>
    inline void
    findNearestNeighborCandidates(const std::vector<NNMergeParticle<Packed>>&                a_localValid,
                                  const std::vector<NNMergeParticle<Packed>>&                a_ghosts,
                                  const std::unordered_set<ParticleID>&                      a_consumedIDs,
                                  const NNSpatialIndex<Packed>&                              a_spatialIndex,
                                  const std::unordered_map<NNCellKey, int, NNCellKeyHasher>& a_liveCellCount,
                                  const RealVect&                                            a_probLo,
                                  const RealVect&                                            a_dx,
                                  const int                         a_numParticlesPerCellThresh,
                                  std::vector<NNMergeEdge<Packed>>& a_edges,
                                  const int                         a_maxFallbackCandidates = 0,
                                  const std::optional<int>          a_maxCellDistance       = std::nullopt) noexcept;

    /**
      @brief Resolve every candidate edge that is safely mergeable with zero cross-rank
      coordination, and mark which remaining particles have committed to a non-trivial outgoing
      proposal.
      @details An edge is trivial-tier eligible iff NEITHER the query nor the candidate is
      boundary-exposed (see isBoundaryExposed()) -- this check happens AFTER finding the true nearest
      neighbor, never as a restriction on the search itself.

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
          silent no-op.
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
      @param[in]     a_maxFallbackCandidates Default 0 (disabled). When a query particle's chosen
                                        candidate is blocked (already consumed, or has its own outgoing
                                        commitment, or exposed/a ghost) and the query particle itself is
                                        NOT exposed, try the edge's NNMergeEdge::fallbackCandidates in
                                        order (pre-fetched by findNearestNeighborCandidates() via one
                                        K-nearest query, K = 1 + a_maxFallbackCandidates), up to this
                                        many extra attempts, before giving up. No fresh query and no
                                        already-tried bookkeeping -- the list is ascending-distance with
                                        no duplicates, so trying the next is advancing an index. Retries
                                        resolve immediately at the point of failure rather than being
                                        re-sorted into the global order (a retry only fires for an edge
                                        whose first choice already failed). A query particle that
                                        exhausts every fallback still proposes with its true, original
                                        nearest neighbor -- fallbacks affect only local trivial-tier
                                        eligibility, never what is proposed cross-patch.
      @param[in]     a_maxCellDistance User-selectable physical merge-eligibility cap, in whole grid
                                        cells, using Chebyshev (max-per-axis) distance between the same
                                        unclamped, position-derived cell key used for the crowding
                                        trigger (`floor(pos/dx)`, see findNearestNeighborCandidates())
                                        -- same cell = 0, any Moore-adjacent cell (incl. diagonals) = 1,
                                        and so on. The position-derived key (not each patch's local
                                        IntVect) is what makes the comparison meaningful across a patch
                                        boundary. `std::nullopt` (default) disables it. A query whose
                                        nearest neighbor exceeds the cap is NOT proposed to -- this is a
                                        physical-validity rule, not a coordination one: the pair is
                                        exactly as invalid via judgeProposals() as it would be
                                        trivially. If a_maxFallbackCandidates is nonzero the query gets a
                                        fallback retry for a closer candidate; otherwise it is left
                                        unmatched with no side effects. The cap is inclusive.
      @param[in]     a_probLo          Physical position of the domain's low corner (see
                                        cellKeyOf()); needed for the dynamic crowding recheck and
                                        a_maxCellDistance.
      @param[in]     a_dxByLevel       Cell spacing PER AMR LEVEL (index 0 = level 0, etc.) -- see
                                        cellKeyOf()'s warning. Each participant's OWN cell key is
                                        computed with a_dxByLevel[participant.level], never a single
                                        shared dx: a candidate can legitimately be a ghost from a
                                        coarser or finer level than the query's own patch, and
                                        getting this wrong silently breaks the dynamic crowding
                                        recheck for every level but one.
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
      @param[in]     a_onMergeCommitted Optional per-merge diagnostic callback (see
                                        MergeCommittedCallback); default is a no-op.
    */
    template <typename Packed, typename Combine, typename PosValid, typename Allocator, typename OnMerge>
    inline void
    resolveTrivialTier(
      const std::vector<NNMergeEdge<Packed>>&                        a_edges,
      const std::unordered_map<ParticleID, bool>&                    a_exposed,
      const std::unordered_map<ParticleID, NNMergeParticle<Packed>>& a_particlesByID,
      const Combine&                                                 a_combine,
      std::unordered_map<NNCellKey, int, NNCellKeyHasher>&           a_liveCellCount,
      std::unordered_set<ParticleID>&                                a_consumedIDs,
      const int                                                      a_numParticlesPerCellThresh,
      std::unordered_set<ParticleID>&                                a_hasOutgoingCommitment,
      std::vector<NNMergeResult<Packed>>&                            a_results,
      const RealVect&                                                a_probLo,
      const std::vector<RealVect>&                                   a_dxByLevel,
      const PosValid&                                                a_isPositionValid,
      const Allocator&                                               a_allocateID,
      const int                                                      a_maxFallbackCandidates = 0,
      const std::optional<int>                                       a_maxCellDistance       = std::nullopt,
      unsigned long long*                                            a_crossLevelMergeCount  = nullptr,
      const OnMerge&                                                 a_onMergeCommitted =
        [](const NNMergeParticle<Packed>&, const NNMergeParticle<Packed>&, const NNMergeParticle<Packed>&) {
        }) noexcept;

    /**
      @brief Package one outgoing proposal for every edge that resolveTrivialTier() could not
      resolve locally, bucketed by destination rank.
      @details Routing rule: ALWAYS to the owner of the candidate (NNMergeEdge::candidateID's
      owning rank) -- never by any id comparison (see NNMergeProposal). If the destination rank equals
      this rank's own (a boundary-exposed candidate that happens to live in a DIFFERENT patch owned by
      the same rank), the proposal still goes into the normal per-destination-rank bucket -- the
      caller (mergeNearestNeighborsRound()) is expected to hand the "bucket addressed to myself"
      directly to judgeProposals() in-process, without an actual MPI round-trip, while every other
      bucket goes out over MPI. This keeps the decision logic in judgeProposals() completely
      unaware of whether a given proposal arrived locally or over the network -- it still goes through
      the exact same pooled judge logic, never a shortcut.
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
    generateProposals(const std::vector<NNMergeEdge<Packed>>&                        a_edges,
                      const std::unordered_set<ParticleID>&                          a_hasOutgoingCommitment,
                      const std::unordered_map<ParticleID, NNMergeParticle<Packed>>& a_particlesByID,
                      const int                                                      a_numRanks,
                      std::vector<std::vector<NNMergeProposal<Packed>>>&             a_proposalsByDestRank) noexcept;

    /**
      @brief Judge every incoming proposal that targets a particle this rank owns, and produce the
      verdicts to send back.
      @details The full rule set this function implements, summarized:
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
      @param[in]     a_particlesByID     Full particle data (position/weight/id/rank/payload) for
                                         every currently-alive local valid particle and every ghost
                                         across all of this rank's patches, keyed by id -- a target
                                         is a potential judge target only if its id is present here
                                         AND absent from a_consumedIDs (see that param).
      @param[in]     a_consumedIDs       Ids already consumed this round (by an earlier trivial or
                                         judged merge). A target id present in a_particlesByID but
                                         ALSO in a_consumedIDs is treated as not-alive -- this
                                         replaces what used to be a separate, pre-filtered
                                         a_myParticlesByID copy of a_particlesByID; checking both
                                         containers at the point of use instead avoids rebuilding
                                         that copy every round.
      @param[in]     a_hasOutgoingCommitment Ids (among a_particlesByID, not in a_consumedIDs) that
                                         have their own outgoing commitment this round (needed for
                                         rule 3 above).
      @param[in]     a_outgoingTargetOf  For an id in a_hasOutgoingCommitment, the id of the
                                         particle it is itself proposing to (needed to detect the
                                         mutual-match exception in rule 3).
      @param[in]     a_combine           Opaque-payload combine callback (see PayloadCombine).
      @param[in,out] a_liveCellCount     Live per-cell counts for this rank's own patches;
                                         decremented in place as targets are accepted, in the same
                                         distance-sorted processing order (see rule 4 above). Used to
                                         dynamically recheck BOTH the target's AND the source's own
                                         cell -- mirrors resolveTrivialTier()'s symmetric recheck
                                         exactly; accepting a proposal whose source's own cell is no
                                         longer over threshold (e.g. drained by an earlier merge
                                         within this same round) would silently violate the crowding
                                         trigger on the source side. If this rank has no live
                                         visibility into the source's cell (it was never shipped here
                                         as a ghost), the source count defaults to 0 -- i.e. rejected
                                         -- the same conservative default already used when a
                                         target's own count is unknown.
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
      @param[in]     a_onMergeCommitted  Optional per-merge diagnostic callback -- see
                                         resolveTrivialTier()'s docs; called with (target, source,
                                         merged) for every accepted target.
    */
    template <typename Packed, typename Combine, typename PosValid, typename Allocator, typename OnMerge>
    inline void
    judgeProposals(
      const std::vector<NNMergeProposal<Packed>>&                    a_incomingProposals,
      const std::unordered_map<ParticleID, NNMergeParticle<Packed>>& a_particlesByID,
      const std::unordered_set<ParticleID>&                          a_consumedIDs,
      const std::unordered_set<ParticleID>&                          a_hasOutgoingCommitment,
      const std::map<ParticleID, ParticleID>&                        a_outgoingTargetOf,
      const Combine&                                                 a_combine,
      std::unordered_map<NNCellKey, int, NNCellKeyHasher>&           a_liveCellCount,
      const int                                                      a_numParticlesPerCellThresh,
      const RealVect&                                                a_probLo,
      const std::vector<RealVect>&                                   a_dxByLevel,
      const PosValid&                                                a_isPositionValid,
      const Allocator&                                               a_allocateID,
      std::vector<NNMergeResult<Packed>>&                            a_results,
      std::vector<std::vector<NNMergeVerdict>>&                      a_verdictsByDestRank,
      const std::optional<int>                                       a_maxCellDistance      = std::nullopt,
      unsigned long long*                                            a_crossLevelMergeCount = nullptr,
      const OnMerge&                                                 a_onMergeCommitted =
        [](const NNMergeParticle<Packed>&, const NNMergeParticle<Packed>&, const NNMergeParticle<Packed>&) {
        }) noexcept;

    /**
      @brief Apply every verdict addressed to one of this rank's own outgoing proposals: remove
      the source particle for every accepted proposal, do nothing for rejected ones.
      @details Rejected proposals are deliberately a no-op -- the source particle is simply left
      alone and reconsidered in next round's fresh candidate search. This algorithm never internally
      retries a rejected proposal against a second-best candidate within the same round.
      @param[in]  a_incomingVerdicts All verdicts addressed to this rank, pooled from every judge
                                     rank (including this rank's own in-process verdict for
                                     proposals it judged for itself).
      @param[out] a_consumedIDs      Ids of every source particle whose proposal was accepted --
                                     the caller removes these from its own valid particle lists.
    */
    inline void
    applyVerdicts(const std::vector<NNMergeVerdict>& a_incomingVerdicts,
                  std::unordered_set<ParticleID>&    a_consumedIDs) noexcept;

    /**
      @brief Result of locating a physical position against the AMR grid -- see placeMergedParticles().
      @details Deliberately independent of any AmrMesh/Realm/ParticleContainer type (this header
      stays Amr-agnostic) -- mirrors ParticleContainer::Destination's fields exactly,
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
      @details Applied UNIFORMLY to every merge result, trivial-tier and judge-decided alike: the
      union of two patches is not generally convex (T-junctions, corners, refinement boundaries), so
      even a purely local (same-rank) merge's result can legitimately belong to a third patch,
      possibly owned by a different rank. Point-location must be O(1) per result (see
      PositionLocator/NNParticleLocation) -- NEVER a linear scan over every patch, which would turn
      this function into a real bottleneck at large patch counts.
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

  } // namespace detail

  /**
    @brief Run one full round (one timestep's worth) of the nearest-neighbor merge algorithm.
    @details Orchestrates, in order:
      1. Ghost-fill -- fresh, exactly once, strictly before any of the below, via
         ParticleContainer::fillGhostParticles() with the AmrMesh particle ghost masks. This
         function does not fill ghosts itself: the caller fills them first, and must have registered
         a width-1 particle ghost mask before the grids were (re)built (a mask is only built during a
         regrid, so registering it late leaves it empty). See
         ItoSolver::mergeSuperparticlesNearestNeighbor for a concrete caller.
      2. Per patch: findNearestNeighborCandidates(), then resolveTrivialTier(). If
         a_iterateLocalTierToConvergence is true, this pair repeats (reusing the step-1 ghost-fill)
         until no further local commits occur. See a_iterateLocalTierToConvergence below.
      3. generateProposals() for everything the trivial tier could not resolve.
      4. Exchange proposals via nnMergeExchangeByRank() (an MPI_Alltoall/Alltoallv helper): this
         rank's own bucket (addressed to itself) is handed directly to judgeProposals() in-process,
         every other bucket is shipped to its destination rank and the incoming proposals collected.
      5. judgeProposals() over the pooled incoming proposals.
      6. Exchange verdicts (nnMergeExchangeByRank(), the other direction) and applyVerdicts().
      7. placeMergedParticles(), then exchange and apply the scatter bucket (nnMergeExchangeByRank()).

    Fresh, globally-unique ids for newly created particles must be allocated by the caller such that
    they can never collide with an existing particle's id or another rank's freshly allocated id
    within the same round; a rank-namespaced counter suffices (see a_allocateID / IDAllocator).

    This function does a SINGLE cross-patch exchange (propose/judge/verdict happens exactly
    once) regardless of a_iterateLocalTierToConvergence: a cell more than roughly 2x over threshold
    is not necessarily fully drained by the cross-patch protocol in one call, and the caller is
    still expected to invoke this once per timestep, relying on repeated calls (and, in a real
    simulation, particle motion between them) for further convergence of anything the cross-patch
    protocol itself left unresolved. What a_iterateLocalTierToConvergence changes is only the LOCAL
    (trivial-tier) portion of a single call -- see its own documentation below.

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
    @param[in,out] a_amr                       AMR mesh -- source of ghost masks, per-level dx, and
                                                the patch layout used to place merged particles.
    @param[in,out] a_particles                 The real particle container to merge, for every
                                                patch this rank owns.
    @param[in]     a_numParticlesPerCellThresh Crowding threshold N.
    @param[in]     a_gather                    Per-slot gather callback (see Gather).
    @param[in]     a_combine                   Opaque-payload combine callback (see
                                                PayloadCombine).
    @param[in]     a_scatter                   Per-surviving-particle scatter callback (see
                                                Scatter).
    @param[in]     a_iterateLocalTierToConvergence Default false. When false,
                                                findNearestNeighborCandidates() and
                                                resolveTrivialTier() each run once. When true, that
                                                candidate-search-plus-trivial-tier pair repeats
                                                (reusing the same ghost-fill from step 1; only which
                                                local particles are still alive changes) until a full
                                                local pass commits zero further trivial merges, and
                                                only then does generateProposals()/judgeProposals()
                                                run (still exactly once, unaffected by this option).
                                                It helps because a single pass fixes each query's edge
                                                once: if that candidate is consumed by a closer edge
                                                first, the query is left unmatched for the rest of the
                                                pass even though another partner may exist nearby. The
                                                gain is modest -- re-running only helps when a
                                                particle's best available choice actually changes
                                                between iterations; a_maxFallbackCandidates below
                                                addresses the more common, stable loss.
    @param[in]     a_maxFallbackCandidates     Forwarded to resolveTrivialTier() (see its docs).
                                                Default 0 (disabled). Lets a query blocked from a good
                                                local partner (because that partner is already "busy"
                                                with its own outgoing commitment -- a condition
                                                iteration alone does not clear) fall back to its next
                                                nearest candidates. A single fallback recovers most of
                                                the achievable improvement.
    @param[in]     a_maxCellDistance           Forwarded to resolveTrivialTier() (see its docs).
                                                `std::nullopt` (default) disables it. A physical cap,
                                                in whole grid cells, on how far apart two particles
                                                may be and still merge -- e.g. 1 permits merging within
                                                the same cell or any Moore-adjacent cell, never
                                                farther. Composes with the options above.
    @param[in]     a_isPositionValid           Forwarded to resolveTrivialTier()/judgeProposals()
                                                (see PositionValid). Default always-true. Composes
                                                with the options above -- e.g. an embedded-boundary
                                                point query for an EB-aware caller (see
                                                ItoSolver::mergeSuperparticlesNearestNeighbor).
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
    @param[in]     a_onMergeCommitted          Optional per-merge diagnostic callback (see
                                                MergeCommittedCallback); default is a no-op.
  */
  template <typename P,
            typename Packed,
            typename Traits = ParticleTraits<P>,
            typename Gather,
            typename Combine,
            typename Scatter,
            typename Allocator,
            typename PosValid = std::function<bool(const RealVect&)>,
            typename OnMerge  = std::function<
               void(const NNMergeParticle<Packed>&, const NNMergeParticle<Packed>&, const NNMergeParticle<Packed>&)>>
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
    unsigned long long* a_crossLevelMergeCount = nullptr,
    const OnMerge&      a_onMergeCommitted =
      [](const NNMergeParticle<Packed>&, const NNMergeParticle<Packed>&, const NNMergeParticle<Packed>&) {
      });

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#include <CD_NearestNeighborParticleMergeImplem.H>

#endif
