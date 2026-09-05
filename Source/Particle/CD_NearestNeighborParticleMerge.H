/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_NearestNeighborParticleMerge.H
 * @brief  Declaration of a distributed, MPI-safe nearest-neighbor particle-merge algorithm.
 * @details This file adds a nearest-neighbor super-particle merge that is safe across patch and rank
 * boundaries. Candidates are found against a particle's full neighborhood (its own patch plus every
 * ghost, filled once per round); a pair where neither participant is reachable from outside the patch
 * resolves immediately (the "trivial tier"), and everything else goes through a propose/judge/verdict
 * protocol that is safe for any number of simultaneous cross-patch suitors.
 * @author Robert Marskar
 */

#ifndef CD_NEARESTNEIGHBORPARTICLEMERGE_H
#define CD_NEARESTNEIGHBORPARTICLEMERGE_H

// Std includes
#include <algorithm>
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
#include <CD_ParticleManagement.H>
#include <CD_NamespaceHeader.H>

class AmrMesh;

namespace ParticleManagement {

namespace detail {

/**
 * @brief Pooled per-round record for one particle: the particle plus its cached boundary-exposure
 * flag.
 * @details particlesByID maps every local valid particle and every ghost to one of these. Folding
 * the exposure flag in (it used to live in a separate id->bool map filled in the same gather loop)
 * lets the gather do one hash insert per particle instead of two, and resolveTrivialTier() one
 * lookup per query/candidate instead of two. exposed is meaningful only for local valid particles;
 * ghosts always carry false (a ghost never initiates a query, and candidate-side exposure is gated
 * on locality).
 * @tparam Packed Caller-defined payload type (see MergeParticle).
 */
template <typename Packed>
struct NNPooledParticle
{
  /**
   * @brief The pooled particle (position/weight/ids/level/payload).
   */
  MergeParticle<Packed> particle;

  /**
   * @brief Whether this local valid particle's cell touches a patch boundary (has ghost targets).
   * Always false for ghosts.
   */
  bool exposed = false;
};

/**
 * @brief One additional ranked candidate beyond a query's primary nearest neighbor.
 * @details Populated only when a_maxFallbackCandidates > 0 is requested from
 * findNearestNeighborCandidates() (see NNMergeEdge::fallbackCandidates). Fetched via a single
 * K-nearest query at generation time (K = 1 + a_maxFallbackCandidates); resolveTrivialTier() walks
 * this pre-fetched, distance-sorted list when a closer candidate is blocked, without re-querying
 * the spatial index.
 */
struct NNFallbackCandidate
{
  /**
   * @brief Global id of this candidate.
   */
  ParticleID candidateID;

  /**
   * @brief Squared distance between the query particle and this candidate.
   */
  Real distanceSquared;

  /**
   * @brief True iff this candidate is a locally-valid particle in the SAME patch as the query
   * (i.e. not a ghost) -- same meaning as NNMergeEdge::candidateIsLocal.
   */
  bool candidateIsLocal;
};

/**
 * @brief One candidate merge edge: a single particle's own nearest-neighbor query result.
 * @details Produced once per over-threshold, unconsumed valid particle by
 * findNearestNeighborCandidates() -- exactly one edge per query particle, its single nearest
 * neighbor among {this patch's own other valid particles} union {every ghost}. Distance is carried
 * explicitly (never recomputed downstream) so trivial-tier resolution and the judge apply the same
 * "closest-first, with a deterministic id tiebreak on exact ties" ordering.
 * @tparam Packed See MergeParticle.
 */
template <typename Packed>
struct NNMergeEdge
{
  /**
   * @brief Squared distance between the query particle and its candidate.
   */
  Real distanceSquared;

  /**
   * @brief Global id of the query particle (always locally valid/owned -- ghosts never
   * initiate a query).
   */
  ParticleID queryID;

  /**
   * @brief Global id of the candidate (query's nearest neighbor). May be a local valid
   * particle (same patch) or a ghost.
   */
  ParticleID candidateID;

  /**
   * @brief True iff the candidate is a locally-valid particle in the SAME patch as the query
   * (i.e. not a ghost). Note this is NOT sufficient on its own to decide trivial-tier
   * eligibility -- see Realm::m_particleGhostExposure: a same-patch candidate can still require
   * the propose/judge/verdict protocol if either participant is boundary-exposed.
   */
  bool candidateIsLocal;

  /**
   * @brief Next-nearest candidates beyond candidateID, in ascending distanceSquared order.
   * @details Empty unless a_maxFallbackCandidates > 0 was passed to
   * findNearestNeighborCandidates(); size at most a_maxFallbackCandidates. resolveTrivialTier()
   * tries these in order when candidateID is blocked (stale, busy, exposed, a ghost, or too far),
   * without re-querying the spatial index. Content is a generation-time snapshot; a candidate's
   * validity is rechecked against a_consumedIDs/a_hasOutgoingCommitment when it is actually tried.
   */
  std::vector<NNFallbackCandidate> fallbackCandidates;
};

/**
 * @brief A cross-patch merge proposal: one particle's own nearest-neighbor query result,
 * bundled as a self-contained, atomic pair for shipment to whichever rank owns the candidate.
 * @details Carries full data for both participants (not a bare id reference) so the judge's
 * decision is self-contained and it can never recombine pieces across different incoming
 * proposals. Routing rule (see generateProposals()): a proposal is always sent to the owner of the
 * querying particle's best candidate, never by any id comparison.
 * @tparam Packed See MergeParticle.
 */
template <typename Packed>
struct NNMergeProposal
{
  /**
   * @brief The proposer's own particle (full data) -- always locally valid/owned by the
   * sending rank.
   */
  MergeParticle<Packed> source;

  /**
   * @brief Global id of the particle being proposed to. Owned by the receiving (judge) rank.
   */
  ParticleID targetID;

  /**
   * @brief Squared distance between @c source and the target. Duplicated from the originating
   * NNMergeEdge so the judge never needs to look anything up to sort/compare proposals.
   */
  Real distanceSquared;
};

/**
 * @brief The judge's reply to one proposal.
 * @details Accepted implies the judge has already created the merged particle and removed its own
 * (target) particle; the proposer, on receiving an accepted verdict, removes its own (source)
 * particle. Rejected implies no action on either side -- the proposer's particle is simply
 * reconsidered next round.
 */
struct NNMergeVerdict
{
  /**
   * @brief Global id of the proposer this verdict is a reply to (i.e. NNMergeProposal::source's
   * id). Used by the proposer to match an incoming verdict back to one of its own outgoing
   * proposals.
   */
  ParticleID proposerID;

  /**
   * @brief True iff the judge accepted this specific proposal (i.e. it was the winning
   * candidate for its target -- see judgeProposals()).
   */
  bool accepted;
};

/**
 * @brief One resolved merge: the two consumed particles' ids and the newly created merged
 * particle (not yet placed into any patch -- see placeMergedParticles()).
 * @tparam Packed See MergeParticle.
 */
template <typename Packed>
struct NNMergeResult
{
  /**
   * @brief Global id of the first consumed particle. For a trivial-tier result, this is the
   * query particle. For a judge-side result, this is the target.
   */
  ParticleID firstID;

  /**
   * @brief Global id of the second consumed particle. For a trivial-tier result, this is the
   * local candidate. For a judge-side result, this is the winning proposer.
   */
  ParticleID secondID;

  /**
   * @brief The newly created merged particle: weighted-centroid position, weighted-average
   * opaque payload (computed by the caller-supplied combine callback), summed weight, and a
   * FRESH globally-unique id (never colliding with any existing particle's id -- see
   * mergeNearestNeighborsTree() for the id-allocation requirement this implies).
   */
  MergeParticle<Packed> merged;
};

} // namespace detail

/**
 * @brief Illustrative signature of the combine callback: produce a merged particle's opaque
 * payload from its two parents' payloads and weights.
 * @details This alias documents the expected callback shape only; it is NOT the type used at the
 * call sites. resolveTrivialTier()/judgeProposals() take the combine callback as a template
 * parameter (not this std::function) so it inlines on the per-merge hot path. Position and weight
 * are handled by the algorithm directly (weighted centroid, summed weight); this callback produces
 * only the opaque @c Packed payload (e.g. weighted-average energy for ItoParticle).
 * @param[in] a_first First parent's payload.
 * @param[in] a_firstWeight First parent's weight.
 * @param[in] a_second Second parent's payload.
 * @param[in] a_secondWeight Second parent's weight.
 * @return The merged particle's payload.
 */
template <typename Packed>
using PayloadCombine = std::function<
  Packed(const Packed& a_first, const Real a_firstWeight, const Packed& a_second, const Real a_secondWeight)>;

/**
 * @brief Illustrative signature of the position-validity predicate for a would-be merged particle.
 * @details Documents the expected callback shape; consumed as a template parameter (not this
 * std::function) at the call sites, for the same inlining reason as PayloadCombine. Called with the
 * weighted-centroid position already computed, once per candidate pair that is otherwise eligible
 * to merge (never on a pair already rejected for another reason), so a caller doing real geometric
 * work here pays only when it matters. Default is always-true. A pair that fails is treated like
 * a_maxCellDistance's "too far" case: not proposed, a fallback tried if available, otherwise left
 * unmatched with no side effects.
 * @param[in] a_mergedPosition The candidate merge's weighted-centroid position.
 * @return True iff a particle may exist at this position (e.g. false inside an embedded boundary).
 */
using PositionValid = std::function<bool(const RealVect& a_mergedPosition)>;

/**
 * @brief Illustrative signature of the fresh-id allocator: produce one globally-unique id for a
 * newly created merged particle.
 * @details Documents the expected callback shape; consumed as a template parameter (not this
 * std::function). Called once per committed merge (trivial-tier or judge-side). The caller must
 * ensure ids never collide -- not just within one call, but across every rank and every
 * past/future call for the life of the simulation; a rank-namespaced counter (this rank's index
 * times a stride plus a monotonically increasing local count) suffices.
 * @return A fresh globally-unique particle id.
 */
using IDAllocator = std::function<ParticleID()>;

/**
 * @brief Optional per-merge diagnostic callback: observe every individual merge as it commits.
 * @details Consumed as a template parameter (not required to be exactly this std::function).
 * Called once per committed merge (trivial-tier or judge-side) with both parents' full pre-merge
 * data and the resulting merged particle (finalized but not yet placed; see
 * placeMergedParticles()). Default is a no-op. Intended for invariant checks, e.g. recomputing the
 * weighted centroid from the parents and comparing against a_merged.position.
 * @param[in] a_first First parent (pre-merge, full data).
 * @param[in] a_second Second parent (pre-merge, full data).
 * @param[in] a_merged The resulting merged particle.
 */
template <typename Packed>
using MergeCommittedCallback = std::function<void(const MergeParticle<Packed>& a_first,
                                                  const MergeParticle<Packed>& a_second,
                                                  const MergeParticle<Packed>& a_merged)>;

namespace detail {

/**
 * @brief The unclamped, position-derived cell index a given physical position falls in.
 * @details The cell-key convention used throughout for the crowding trigger and a_maxCellDistance.
 * Not clamped to any patch's box, so a ghost buckets into its true origin cell and particles in
 * different patches' local frames stay directly comparable.
 * @warning a_dx MUST be the cell spacing of the AMR level the position's particle lives on (see
 * MergeParticle::level), never a single fixed dx shared across levels. The returned IntVect is
 * only meaningful together with that level; NNCellKey pairs the two so they cannot be conflated
 * (pooling bare cell keys from different levels aliases coarse- and fine-level cell counts).
 * @param[in] a_position Physical position.
 * @param[in] a_probLo Physical position of the domain's low corner.
 * @param[in] a_dx Cell spacing (per direction) of the position's OWN level.
 * @return floor((a_position - a_probLo) / a_dx), component-wise.
 */
inline IntVect
cellKeyOf(const RealVect& a_position, const RealVect& a_probLo, const RealVect& a_dx) noexcept;

/**
 * @brief A cell key together with the AMR level it was computed on -- see cellKeyOf()'s warning.
 * @details Two NNCellKeys are equal only if BOTH the level and the cell index match; a coarse-
 * level cell and an unrelated fine-level cell that happen to share the same bare IntVect are
 * correctly kept distinct because their level differs. This is the key type for every live
 * per-cell count pooled across multiple levels (see a_liveCellCount in resolveTrivialTier()/
 * judgeProposals()/findNearestNeighborCandidates()).
 */
using NNCellKey = std::pair<int, IntVect>;

/**
 * @brief Hash functor for NNCellKey -- REQUIRED whenever NNCellKey is used as a hash-map key;
 * NNCellKey must never be used as an ordered (%std::map or %std::set) key.
 * @details IntVect::operator< (which std::pair<int,IntVect>::operator< calls for the second
 * component) is a component-wise "strictly dominates in every direction" partial order, not a
 * strict weak ordering. Two unrelated cells where neither dominates the other (e.g. (76,104) vs
 * (80,101)) compare as neither-less-than, which %std::map/%std::set treat as equal keys -- pooling
 * unrelated cells' counts into one entry. Hashing depends only on IntVect::operator== (full
 * component-wise equality), which is correct.
 */
struct NNCellKeyHasher
{
  /**
   * @brief Hash operator.
   * @param[in] a_key Cell key to hash.
   * @return Hash value.
   */
  std::size_t
  operator()(const NNCellKey& a_key) const noexcept
  {
    std::size_t hash = LevelTiles::TileHasher()(a_key.second);
    hash ^= static_cast<std::size_t>(static_cast<unsigned int>(a_key.first)) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    return hash;
  }
};

} // namespace detail

/**
 * @brief The per-cell particle count this merge drains each cell down to.
 * @details Every crowding decision in this algorithm asks the same question -- "does this cell hold
 * more particles than it is allowed to?" -- so this type is the single place that answers it. The
 * plain case is one uniform count for every cell, and an implicit constructor from @c int provides
 * it, so a caller with nothing special to say simply passes the count.
 *
 * A caller may additionally RESERVE slots in individual cells. A reserved slot is one already spoken
 * for by particles this call cannot see, which must nonetheless count against the cell's total. The
 * budget such a cell is drained to is the uniform count minus its reservations, floored at zero. The
 * motivating case is a two-container merge in which one container has already been reduced and its
 * results must still occupy room in the shared per-cell target.
 *
 * Reservations are keyed by (level, cell), not by cell alone, for the reason given in
 * detail::NNCellKey: a coarse and a fine cell can share a bare IntVect.
 * @note Lookup is a single branch when nothing is reserved, so the uniform case costs no more than
 * the bare integer it replaces.
 */
class NNCellBudget
{
public:
  /**
   * @brief Construct a uniform budget. Deliberately implicit -- a bare count is the common case.
   * @param[in] a_uniform Particle count every cell is drained down to.
   */
  NNCellBudget(const int a_uniform) noexcept : m_uniform(a_uniform)
  {}

  /**
   * @brief Reserve one slot in a cell, lowering that cell's budget by one.
   * @details Call once per particle that is invisible to this merge but still occupies the cell.
   * @param[in] a_level AMR level the cell lives on.
   * @param[in] a_cell Cell index, in the convention of detail::cellKeyOf().
   */
  inline void
  reserveOne(const int a_level, const IntVect& a_cell)
  {
    const int count = ++m_reserved[detail::NNCellKey(a_level, a_cell)];

    m_minBudget = std::min(m_minBudget, std::max(0, m_uniform - count));
  }

  /**
   * @brief Reserve one slot in the cell a position falls in.
   * @details Preferred over the IntVect overload for a caller holding positions: the cell index is
   * derived here, so it cannot disagree with the convention this algorithm buckets particles under.
   * @param[in] a_level AMR level the position's particle lives on.
   * @param[in] a_position Physical position of the particle occupying the slot.
   * @param[in] a_probLo Physical position of the domain's low corner.
   * @param[in] a_dx Cell spacing (per direction) of that particle's OWN level.
   */
  inline void
  reserveOne(const int a_level, const RealVect& a_position, const RealVect& a_probLo, const RealVect& a_dx)
  {
    this->reserveOne(a_level, detail::cellKeyOf(a_position, a_probLo, a_dx));
  }

  /**
   * @brief The budget for one cell: the uniform count less that cell's reservations, floored at 0.
   * @param[in] a_key Cell to look up.
   * @return Particle count this cell is drained down to.
   */
  inline int
  operator()(const detail::NNCellKey& a_key) const noexcept
  {
    if (m_reserved.empty()) {
      return m_uniform;
    }

    const auto it = m_reserved.find(a_key);

    return (it == m_reserved.end()) ? m_uniform : std::max(0, m_uniform - it->second);
  }

  /**
   * @brief The smallest budget any cell has.
   * @details For use where a cell's own budget cannot be looked up -- a conservative stand-in, since
   * no cell is drained to fewer particles than this. Equals the uniform count when nothing is
   * reserved.
   * @return The minimum over every cell's budget.
   */
  inline int
  minimum() const noexcept
  {
    return m_reserved.empty() ? m_uniform : std::min(m_uniform, m_minBudget);
  }

private:
  /**
   * @brief Particle count every cell is drained to before its own reservations are subtracted.
   */
  int m_uniform;

  /**
   * @brief Running minimum of every reserved cell's budget, maintained by reserveOne().
   */
  int m_minBudget = std::numeric_limits<int>::max();

  /**
   * @brief Slots already spoken for, per (level, cell). Empty for a uniform budget.
   */
  std::unordered_map<detail::NNCellKey, int, detail::NNCellKeyHasher> m_reserved;
};

namespace detail {

/**
 * @brief A non-level-qualified cell index, distinct from NNCellKey, used where a cell is expressed
 * in a single consistent frame (the query patch's own level).
 * @details A separate alias from NNCellKey (which is level-qualified for the cross-level
 * crowding-count map) so the two cannot be passed interchangeably. NNWalkCell is always computed
 * with the query patch's own level's dx, even for a cross-level ghost, so candidates are
 * geometrically comparable to the query in one frame (see cellKeyOf()'s warning).
 */
using NNWalkCell = IntVect;

/**
 * @brief A patch's spatial index for one round: an EBGeometry point-cloud index (PointCloudBVH or
 * PointCloudHashGrid, selected via @c Cloud) over the whole patch (every local valid particle
 * followed by every ghost) plus its precomputed k-nearest-neighbor graph.
 * @details Built once per patch per mergeNearestNeighborsTree()/mergeNearestNeighborsHash() call
 * via buildNNSpatialIndex() and reused across every internal a_iterateLocalTierToConvergence pass
 * via lazy deletion (see findNearestNeighborCandidates()'s own docs). A patch with no crowded cell
 * gets no index at all (pointCloud stays null -- see buildNNSpatialIndex()'s per-cell pre-filter).
 * The point cloud's metadata is each point's globalID; a hit's cloud index c decodes ownership
 * directly: c < pointCloudNLocal => local valid a_localValid[c], else ghost
 * a_ghosts[c - pointCloudNLocal].
 * @tparam Packed See MergeParticle.
 * @tparam Cloud EBGeometry point-cloud index type -- @c EBGeometry::PointCloudBVH<Real,ParticleID> (default, used
 * by mergeNearestNeighborsTree()) or @c EBGeometry::PointCloudHashGrid<Real,ParticleID> (used by
 * mergeNearestNeighborsHash()). The two expose an identical query interface (same @c Hit shape, same method names),
 * so every consumer of this struct is written once and works for either.
 */
template <typename Packed, typename Cloud = EBGeometry::PointCloudBVH<Real, ParticleID>>
struct NNSpatialIndex
{
  /**
   * @brief The point-cloud index over this patch's whole point set (local valid first, then
   * ghosts). Null iff there were no particles at all, or no cell in this patch could possibly be
   * crowded (see buildNNSpatialIndex()'s per-cell pre-filter).
   */
  std::shared_ptr<Cloud> pointCloud;

  /**
   * @brief Number of local-valid points at the front of the point cloud's order (used to decode a
   * hit's local-vs-ghost ownership from its cloud index). Only meaningful when pointCloud is
   * non-null.
   */
  std::size_t pointCloudNLocal = 0;

  /**
   * @brief The whole k-nearest-neighbor graph, precomputed once via @c Cloud::allNearestNeighbors(pointCloudK) at
   * build time (batched, cache-warm, seeded from each point's own spatial locality -- strictly cheaper than
   * per-query calls). Row-major: entry [c * pointCloudK + j] is the j-th nearest OTHER point of cloud index c. Only
   * the local-valid rows (c < pointCloudNLocal) are ever read by findNearestNeighborCandidates(); ghost rows are
   * computed but unused. Empty iff pointCloud is null.
   */
  std::vector<typename Cloud::Hit> pointCloudGraph;

  /**
   * @brief k used to build pointCloudGraph (= 1 + a_maxFallbackCandidates); the row stride.
   */
  int pointCloudK = 0;
};

/**
 * @brief Build a patch's spatial index for one round: an EBGeometry point-cloud index (see @c Cloud) over the whole
 * patch (every local valid particle then every ghost) plus its precomputed k-nearest-neighbor graph -- UNLESS no
 * cell in this patch could possibly be crowded, in which case nothing is built.
 * @details The point cloud carries each point's globalID as metadata; a hit's cloud index decodes
 * ownership (c < nLocal => local valid, else ghost). After the cloud is built, its whole
 * a_kNearest-nearest-neighbor graph is precomputed once via @c Cloud::allNearestNeighbors()
 * (batched, own-locality-seeded) into a_index.pointCloudGraph and reused across
 * every internal a_iterateLocalTierToConvergence pass (positions never change within a round; only
 * consumed status does, and that is filtered at lookup time -- see
 * findNearestNeighborCandidates()).
 *
 * @b Per-cell pre-filter: if NO cell in this patch is over its a_cellBudget
 * particles, nothing is built (a_index.pointCloud stays null, exactly as if the patch had zero
 * particles) -- the SAME cheap "is this even worth it" check the kd merges effectively get for
 * free by dispatching per cell. This is an EXACT check (not an average-density heuristic, which
 * would be a real bug: a low-average patch can still hold one badly-clustered over-threshold cell):
 * a lightweight per-cell tally of the patch's own particles, keyed exactly as the crowding trigger
 * (cellKeyOf(pos, probLo, dx), local and ghost alike), decides it directly.
 * @tparam Packed See MergeParticle.
 * @tparam Cloud See NNSpatialIndex. Constructed with its two-argument constructor (positions, metadata), so this
 * function relies on each @c Cloud's own default for whatever third construction parameter it takes (target leaf
 * size for PointCloudBVH, target points per cell for PointCloudHashGrid) -- both defaults are reasonable, general
 * choices, so neither backend needs bespoke construction logic here.
 * @param[out] a_index The resulting spatial index.
 * @param[in] a_localValid This patch's own full valid-particle set.
 * @param[in] a_ghosts Every ghost shipped to this patch this round.
 * @param[in] a_probLo Physical position of the domain's low corner.
 * @param[in] a_dx Cell spacing of THIS PATCH's own level.
 * @param[in] a_cellBudget The physics-facing per-cell budget used elsewhere this round (see
 * mergeNearestNeighborsTree()'s own a_cellBudget) -- used ONLY for the per-cell pre-filter above.
 * @param[in] a_kNearest k = 1 + a_maxFallbackCandidates -- the neighbor count to precompute per point into
 * NNSpatialIndex::pointCloudGraph.
 */
template <typename Packed, typename Cloud = EBGeometry::PointCloudBVH<Real, ParticleID>>
inline void
buildNNSpatialIndex(NNSpatialIndex<Packed, Cloud>&            a_index,
                    const std::vector<MergeParticle<Packed>>& a_localValid,
                    const std::vector<MergeParticle<Packed>>& a_ghosts,
                    const RealVect&                           a_probLo,
                    const RealVect&                           a_dx,
                    const NNCellBudget&                       a_cellBudget,
                    const int                                 a_kNearest) noexcept;

/**
 * @brief For every over-threshold, currently-alive valid particle in a patch, find its single
 * TRUE nearest neighbor among {this patch's own other valid particles} union {every ghost
 * shipped to this patch}.
 * @details This is the only place a nearest-neighbor search happens, and it always runs against
 * the full (local + ghost) view, never a local-only view: restricting it biases boundary-adjacent
 * particles toward worse matches than interior ones. Crowding eligibility (is the query's own cell
 * over threshold) is checked here at generation time from the live per-cell counts over valid +
 * ghost occupancy; the symmetric check on the candidate's cell is done later, dynamically, by
 * resolveTrivialTier()/judgeProposals(), since earlier merges this round can change it.
 *
 * a_localValid/a_ghosts are this patch's full, unfiltered sets for the round; this function does the
 * alive-filtering itself via a_consumedIDs (checked against both the query and every candidate).
 * a_spatialIndex is built once per patch per mergeNearestNeighborsTree() call over the full set and
 * reused across every a_iterateLocalTierToConvergence pass via this same lazy-deletion check, so it
 * is never rebuilt per pass.
 *
 * Each query's precomputed row is read from a_spatialIndex's k-nearest-neighbor graph (see
 * NNSpatialIndex::pointCloudGraph / buildNNSpatialIndex()), computed once for the whole patch via
 * the EBGeometry SIMD PointCloudBVH's batched allNearestNeighbors(). Per-query work is a row lookup
 * plus filtering: a hit is skipped if its id is in a_consumedIDs, or (when a_maxCellDistance is set)
 * if it lies beyond the bound; survivors feed the shared top-K logic that builds the edge. Self is
 * excluded by the graph (each row is the nearest OTHER points).
 *
 * a_maxCellDistance, when set, filters the graph row: rows are distance-sorted, so an out-of-bound
 * candidate is farther than the bound and can never displace a within-bound one, and filtering
 * yields exactly the within-bound k-nearest (see its param docs for the exact, safe-but-not-tight
 * bound).
 *
 * When a_maxFallbackCandidates > 0, the graph is built with k = 1 + a_maxFallbackCandidates, so each
 * row holds the K nearest and NNMergeEdge::fallbackCandidates is populated with candidates 2..K
 * (distance-sorted); resolveTrivialTier() walks that pre-fetched list when the primary is blocked,
 * without re-querying the index.
 *
 * @tparam Packed See MergeParticle.
 * @tparam Cloud See NNSpatialIndex.
 * @param[out] a_edges One entry per over-threshold, currently-alive query particle that has at least one candidate
 * (empty if none found).
 * @param[in] a_localValid This patch's own FULL valid-particle set for this round (query sources AND candidates) --
 * see the contract-change note above; not pre-filtered by a_consumedIDs.
 * @param[in] a_ghosts Every ghost shipped to this patch this round (candidates only; ghosts never initiate a query)
 * -- also not pre-filtered.
 * @param[in] a_consumedIDs Lazy-deletion skip set: ids already consumed earlier this round (by an earlier
 * a_iterateLocalTierToConvergence pass, or an earlier query within this same pass). Checked against both the query
 * and every candidate encountered during the walk.
 * @param[in] a_spatialIndex This patch's spatial index for this round (see NNSpatialIndex), built once via
 * buildNNSpatialIndex() over the FULL a_localValid/ a_ghosts sets before any consumption happens.
 * @param[in] a_liveCellCount Live particle count (valid + ghost) per level-qualified cell key (see NNCellKey),
 * already reflecting this round's frozen ghost-fill snapshot.
 * @param[in] a_probLo Physical position of the domain's low corner (see cellKeyOf()).
 * @param[in] a_dx Cell spacing of THIS PATCH's own level (see cellKeyOf()'s warning) -- valid here because every
 * query in a_localValid is necessarily on this same level; candidates never have their crowding checked by this
 * function, so a cross-level ghost candidate needs no dx here.
 * @param[in] a_cellBudget Per-cell budget -- a cell is eligible while it holds strictly more
 * than this many particles.
 * @param[in] a_maxFallbackCandidates Forwarded from resolveTrivialTier() (see its docs); 0 (default) means an
 * ordinary single-nearest query and every edge's fallbackCandidates stays empty.
 * @param[in] a_maxCellDistance Same hard cutoff as resolveTrivialTier()'s own a_maxCellDistance (see its docs);
 * std::nullopt (default) disables it. Applied here purely as an additional, always-safe filter on the graph row --
 * resolveTrivialTier()/judgeProposals() still apply the exact, level-aware, per-axis Chebyshev check
 * (nnMergeCrossLevelTooFar()) downstream regardless, unchanged. The bound used here is deliberately looser: a
 * candidate is dropped only once its EUCLIDEAN distance exceeds a_maxCellDistance * min(a_dx) * sqrt(SpaceDim) --
 * the sqrt(SpaceDim) factor is required because a per-axis Chebyshev cell distance of N does not bound Euclidean
 * physical distance by N * dx alone (a diagonal candidate exactly N cells away in both directions is N * dx *
 * sqrt(SpaceDim) away in a straight line), and this bound must never be tighter than the true Chebyshev cutoff or it
 * could wrongly drop a candidate resolveTrivialTier()/judgeProposals() would still have accepted. Safe but not
 * maximally tight.
 */
template <typename Packed, typename Cloud = EBGeometry::PointCloudBVH<Real, ParticleID>>
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
                              const int                                                  a_maxFallbackCandidates = 0,
                              const std::optional<int> a_maxCellDistance = std::nullopt) noexcept;

/**
 * @brief Per-point metadata carried by a one-cell spatial index's PointCloudBVH: a particle's
 * global id and whether it is a local-valid particle of the patch the index was built for.
 * @details Unlike NNSpatialIndex (one whole-patch cloud, ownership decoded from cloud-index
 * position), each one-cell cloud is built independently per cell, so ownership is carried directly
 * as metadata instead.
 */
struct NNCellMeta
{
  /**
   * @brief Global id of this particle.
   */
  ParticleID globalID;

  /**
   * @brief True iff this particle is a locally-valid particle of the patch the enclosing
   * NNCellSpatialIndex was built for (i.e. not a ghost) -- same meaning as
   * NNMergeEdge::candidateIsLocal.
   */
  bool isLocal;

  /**
   * @brief This particle's TRUE AMR level (see MergeParticle::level) -- needed by
   * findNearestNeighborCandidatesOneCell() to re-check a cross-level candidate's Chebyshev
   * distance in the correct (finer-of-the-two-levels) frame, since the cell buckets themselves are
   * all keyed in the query patch's own level's frame.
   */
  int level;
};

/**
 * @brief One cell's PointCloudBVH for the one-cell nearest-neighbor merge, over every local valid
 * and ghost particle whose cellKeyOf() falls in that cell.
 */
using NNCellCloud = EBGeometry::PointCloudBVH<Real, NNCellMeta>;

/**
 * @brief A patch's one-cell spatial index: one NNCellCloud per occupied cell, keyed by NNWalkCell
 * (see its docs for the single-frame convention this uses).
 * @details Empty (no entries) iff no cell in the patch could possibly be crowded -- the same
 * per-cell pre-filter as buildNNSpatialIndex(); see buildNNCellSpatialIndex().
 */
using NNCellSpatialIndex = std::unordered_map<NNWalkCell, std::shared_ptr<NNCellCloud>, LevelTiles::TileHasher>;

/**
 * @brief Build a patch's one-cell spatial index for the one-cell nearest-neighbor merge: one
 * PointCloudBVH per occupied cell -- UNLESS no cell in this patch could possibly be crowded, in
 * which case nothing is built.
 * @details Every OCCUPIED cell gets a cloud, not just crowded ones -- an under-threshold cell can
 * still hold a valid Moore-neighbor candidate for a query in an adjacent crowded cell (see
 * findNearestNeighborCandidatesOneCell()). Uses the same exact per-cell tally as
 * buildNNSpatialIndex() to decide whether to build anything at all for this patch.
 * @tparam Packed See MergeParticle.
 * @param[out] a_index The resulting one-cell spatial index.
 * @param[in] a_localValid This patch's own full valid-particle set.
 * @param[in] a_ghosts Every ghost shipped to this patch this round.
 * @param[in] a_probLo Physical position of the domain's low corner.
 * @param[in] a_dx Cell spacing of THIS PATCH's own level.
 * @param[in] a_cellBudget The physics-facing per-cell budget used elsewhere this round; used
 * ONLY for the per-cell pre-filter above.
 */
template <typename Packed>
inline void
buildNNCellSpatialIndex(NNCellSpatialIndex&                       a_index,
                        const std::vector<MergeParticle<Packed>>& a_localValid,
                        const std::vector<MergeParticle<Packed>>& a_ghosts,
                        const RealVect&                           a_probLo,
                        const RealVect&                           a_dx,
                        const NNCellBudget&                       a_cellBudget) noexcept;

/**
 * @brief For every over-threshold, currently-alive valid particle in a patch, find its single TRUE
 * nearest neighbor restricted to its own cell and its Moore-adjacent neighbor cells (Chebyshev cell
 * distance <= 1) -- the one-cell counterpart of findNearestNeighborCandidates().
 * @details Visits the query's own cell plus every occupied Moore-neighbor cell (up to
 * 3^SpaceDim - 1 of them), ordered by ascending squared point-to-cell distance, querying each
 * cell's PointCloudBVH (the own cell via an external query that also matches the query particle
 * itself, filtered out by id; neighbor cells via an ordinary external query, no self match
 * possible). Stops visiting further cells once a full top-K candidate set has been found and the
 * next cell's minimum possible distance cannot beat or tie the current worst of K -- an exact, not
 * approximate, prune, since point-to-cell-AABB distance is a true lower bound on any point inside
 * that cell.
 *
 * The cell buckets themselves are all keyed in the query patch's OWN level's frame (see
 * buildNNCellSpatialIndex()), so Moore-adjacency in that frame does not, by itself, bound the
 * Chebyshev distance between the query and a candidate on a DIFFERENT AMR level in the
 * finer-of-the-two-levels frame that nnMergeCrossLevelTooFar() (and so resolveTrivialTier()/
 * judgeProposals()) actually enforce. A same-level candidate is unaffected (query and candidate
 * already share one frame). For a cross-level candidate, this function re-checks it exactly, in that
 * same finer-of-the-two-levels frame, before it is ever considered -- so every candidate this
 * function returns is guaranteed at Chebyshev distance <= 1 in the frame the downstream exact check
 * uses, not merely in the query patch's own frame. Callers MUST still forward that same cap (not
 * std::nullopt) as a_maxCellDistance to resolveTrivialTier()/judgeProposals().
 * @tparam Packed See MergeParticle.
 * @param[out] a_edges One entry per over-threshold, currently-alive query particle that has at least one candidate.
 * @param[in] a_localValid This patch's own FULL valid-particle set for this round (query sources).
 * @param[in] a_consumedIDs Lazy-deletion skip set: ids already consumed earlier this round.
 * @param[in] a_cellIndex This patch's one-cell spatial index (see NNCellSpatialIndex), built once via
 * buildNNCellSpatialIndex() over the FULL a_localValid/a_ghosts sets before any consumption happens.
 * @param[in] a_liveCellCount Live particle count (valid + ghost) per level-qualified cell key (see NNCellKey).
 * @param[in] a_probLo Physical position of the domain's low corner (see cellKeyOf()).
 * @param[in] a_dx Cell spacing of THIS PATCH's own level (see cellKeyOf()'s warning).
 * @param[in] a_dxByLevel Cell spacing per AMR level (see resolveTrivialTier()'s own docs) -- needed for the exact
 * cross-level re-check described above.
 * @param[in] a_cellBudget Per-cell budget -- a cell is eligible while it holds strictly more
 * than this many particles.
 * @param[in] a_maxFallbackCandidates Same meaning as findNearestNeighborCandidates()'s own parameter; 0 (default)
 * means an ordinary single-nearest query and every edge's fallbackCandidates stays empty.
 */
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
                                     const int a_maxFallbackCandidates = 0) noexcept;

/**
 * @brief Resolve every candidate edge that is safely mergeable with zero cross-rank
 * coordination, and mark which remaining particles have committed to a non-trivial outgoing
 * proposal.
 * @details An edge is trivial-tier eligible iff NEITHER the query nor the candidate is
 * boundary-exposed (see Realm::m_particleGhostExposure) -- this check happens AFTER finding the true nearest
 * neighbor, never as a restriction on the search itself.
 *
 * All eligible edges (across every patch this rank owns) are processed in ONE globally
 * distance-sorted order, with a deterministic secondary tiebreak on participant ids -- exact
 * distance ties (e.g. true mutual nearest neighbors) are common enough that leaving them to
 * arbitrary iteration/insertion order produces real, observed nondeterminism that cascades into
 * unrelated pairs several hops away.
 *
 * At the moment each edge is actually processed (not at generation time): re-verify both
 * participants are still alive (lazy invalidation -- skip if either was already consumed by a
 * closer edge processed earlier in this same pass); if the candidate is local, re-verify it has
 * NOT, in the meantime, committed to its OWN outgoing (non-trivial) proposal -- a purely
 * geometric "exposed" check on the candidate alone is insufficient, since a non-exposed
 * particle's true nearest neighbor can still be a foreign ghost even though the particle itself
 * is unreachable externally; consuming such a candidate trivially anyway is a real
 * double-consumption bug this dynamic recheck closes. Finally, re-verify both participants'
 * CURRENT cell occupancy (reflecting decrements from edges already processed this pass, not the
 * generation-time snapshot) still exceeds the threshold.
 *
 * On a successful trivial-tier commit: remove both participants from the caller's own valid
 * particle list (conceptually -- see a_consumedIDs), decrement both cells' live counts, and record the merge via
 * a_combine/into a_results. On an edge that cannot be resolved trivially, exactly one of two
 * things happens, and the distinction matters:
 *   - If the candidate is a ghost, OR either participant is boundary-exposed, OR the candidate
 *     turned out to already have its own outgoing commitment: mark the query particle as
 *     having an outgoing commitment (see a_hasOutgoingCommitment) so generateProposals() can
 *     pick it up. In all three cases the candidate genuinely still exists (it just isn't
 *     trivially available right now), so a proposal is meaningful -- judgeProposals() already
 *     has a rule for "target has its own outgoing commitment" (accept only on a mutual match,
 *     reject otherwise), and that rule is exactly what should get a chance to run here, not a
 *     silent no-op.
 *   - If the candidate is stale (already consumed -- found dead by the lazy-invalidation check
 *     above): do NOT mark the query particle as having an outgoing commitment, and do NOT
 *     generate a proposal. There is no meaningful id to propose to. Critically, this also means
 *     the query particle is NOT marked "busy" -- it remains fully available for some OTHER
 *     particle to claim as an ordinary trivial candidate later in this same sorted pass, or via
 *     a_maxFallbackCandidates below.
 *
 * @tparam Packed See MergeParticle.
 * @tparam Combine See PayloadCombine.
 * @param[in,out] a_liveCellCount Live per-cell counts; decremented in place as trivial merges commit, so later edges
 * in the same pass see the up-to-date occupancy.
 * @param[in,out] a_consumedIDs Ids already consumed this round (by an earlier trivial merge, or found dead for any
 * other reason); an edge referencing a consumed id is stale and skipped. Updated in place.
 * @param[out] a_hasOutgoingCommitment Set of query-particle ids that were NOT resolved trivially (excluding the
 * stale-candidate case, see above) and must go on to generateProposals().
 * @param[out] a_results One entry per trivial-tier merge actually committed, in the order they were processed.
 * @param[in] a_edges Candidate edges from findNearestNeighborCandidates(), for ALL patches this rank owns (pooled
 * together so the global distance sort spans the whole rank, not just one patch).
 * @param[in] a_particlesByID Full particle data for every currently-alive local valid particle and every ghost across
 * all of this rank's patches, keyed by id. Each value is an NNPooledParticle: the particle (position/weight/id/rank/
 * payload, needed to compute a merge's weighted-centroid position and call a_combine, since an edge only carries
 * ids/distance) plus its cached boundary-exposure flag (Realm::m_particleGhostExposure; false for ghosts).
 * @param[in] a_combine Opaque-payload combine callback (see PayloadCombine).
 * @param[in] a_cellBudget Per-cell budget, for the dynamic per-cell recheck (same budget used
 * everywhere else this round).
 * @param[in] a_probLo Physical position of the domain's low corner (see cellKeyOf()); needed for the dynamic
 * crowding recheck and a_maxCellDistance.
 * @param[in] a_dxByLevel Cell spacing PER AMR LEVEL (index 0 = level 0, etc.) -- see cellKeyOf()'s warning. Each
 * participant's OWN cell key is computed with a_dxByLevel[participant.level], never a single shared dx: a candidate
 * can legitimately be a ghost from a coarser or finer level than the query's own patch, and getting this wrong
 * silently breaks the dynamic crowding recheck for every level but one.
 * @param[in] a_isPositionValid Position-validity predicate (see PositionValid), called once per otherwise-eligible
 * pair with the weighted-centroid position already computed. Treated exactly like a_maxCellDistance's "too far" case
 * on failure: no proposal for the original candidate, a fallback retry if available, otherwise left unmatched with
 * no side effects.
 * @param[in] a_allocateID Fresh-id allocator (see IDAllocator), called once per committed trivial-tier merge.
 * @param[in] a_maxFallbackCandidates Default 0 (disabled). When a query particle's chosen candidate is blocked
 * (already consumed, or has its own outgoing commitment, or exposed/a ghost) and the query particle itself is NOT
 * exposed, try the edge's NNMergeEdge::fallbackCandidates in order (pre-fetched by findNearestNeighborCandidates()
 * via one K-nearest query, K = 1 + a_maxFallbackCandidates), up to this many extra attempts, before giving up. No
 * fresh query and no already-tried bookkeeping -- the list is ascending-distance with no duplicates, so trying the
 * next is advancing an index. Retries resolve immediately at the point of failure rather than being re-sorted into
 * the global order (a retry only fires for an edge whose first choice already failed). A query particle that
 * exhausts every fallback still proposes with its true, original nearest neighbor -- fallbacks affect only local
 * trivial-tier eligibility, never what is proposed cross-patch.
 * @param[in] a_maxCellDistance User-selectable physical merge-eligibility cap, in whole grid cells, using Chebyshev
 * (max-per-axis) distance between the same unclamped, position-derived cell key used for the crowding trigger
 * (`floor(pos/dx)`, see findNearestNeighborCandidates()) -- same cell = 0, any Moore-adjacent cell (incl. diagonals)
 * = 1, and so on. The position-derived key (not each patch's local IntVect) is what makes the comparison meaningful
 * across a patch boundary. `std::nullopt` (default) disables it. A query whose nearest neighbor exceeds the cap is
 * NOT proposed to -- this is a physical-validity rule, not a coordination one: the pair is exactly as invalid via
 * judgeProposals() as it would be trivially. If a_maxFallbackCandidates is nonzero the query gets a fallback retry
 * for a closer candidate; otherwise it is left unmatched with no side effects. The cap is inclusive.
 * @param[in,out] a_crossLevelMergeCount Optional diagnostic counter (nullptr, the default, disables it entirely --
 * zero cost when not needed): incremented once per committed merge whose two participants were on DIFFERENT AMR
 * levels (query.level != cand.level). NOT reset by this function -- the caller owns its lifetime, so counts can be
 * accumulated across a_iterateLocalTierToConvergence's internal iterations, or across this rank's judgeProposals()
 * acceptances too (see that function's own parameter).
 * @param[in] a_onMergeCommitted Optional per-merge diagnostic callback (see MergeCommittedCallback); default is a
 * no-op.
 */
template <typename Packed, typename Combine, typename PosValid, typename Allocator, typename OnMerge>
inline void
resolveTrivialTier(
  std::unordered_map<NNCellKey, int, NNCellKeyHasher>&            a_liveCellCount,
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
  const int                                                       a_maxFallbackCandidates = 0,
  const std::optional<int>                                        a_maxCellDistance       = std::nullopt,
  unsigned long long*                                             a_crossLevelMergeCount  = nullptr,
  const OnMerge&                                                  a_onMergeCommitted =
    [](const MergeParticle<Packed>&, const MergeParticle<Packed>&, const MergeParticle<Packed>&) {
    }) noexcept;

/**
 * @brief Package one outgoing proposal for every edge that resolveTrivialTier() could not
 * resolve locally, bucketed by destination rank.
 * @details Routing rule: ALWAYS to the owner of the candidate (NNMergeEdge::candidateID's
 * owning rank) -- never by any id comparison (see NNMergeProposal). If the destination rank equals
 * this rank's own (a boundary-exposed candidate that happens to live in a DIFFERENT patch owned by
 * the same rank), the proposal still goes into the normal per-destination-rank bucket -- the
 * caller (mergeNearestNeighborsTree()) is expected to hand the "bucket addressed to myself"
 * directly to judgeProposals() in-process, without an actual MPI round-trip, while every other
 * bucket goes out over MPI. This keeps the decision logic in judgeProposals() completely
 * unaware of whether a given proposal arrived locally or over the network -- it still goes through
 * the exact same pooled judge logic, never a shortcut.
 * @tparam Packed See MergeParticle.
 * @param[out] a_proposalsByDestRank One bucket per rank; a_proposalsByDestRank[r] holds every proposal whose
 * candidate is owned by rank r (including r == this rank's own).
 * @param[in] a_edges Candidate edges from findNearestNeighborCandidates(), for all of this rank's patches.
 * @param[in] a_hasOutgoingCommitment Query-particle ids that resolveTrivialTier() could not resolve (i.e. exactly
 * the edges this function packages).
 * @param[in] a_particlesByID Full particle data for every currently-alive local valid particle and every ghost,
 * keyed by id (each value an NNPooledParticle wrapping the particle) -- needed to build the full MergeParticle
 * bundle for each proposal.
 * @param[in] a_numRanks Total MPI rank count (size of the returned bucket vector).
 */
template <typename Packed>
inline void
generateProposals(std::vector<std::vector<NNMergeProposal<Packed>>>&              a_proposalsByDestRank,
                  const std::vector<NNMergeEdge<Packed>>&                         a_edges,
                  const std::unordered_set<ParticleID>&                           a_hasOutgoingCommitment,
                  const std::unordered_map<ParticleID, NNPooledParticle<Packed>>& a_particlesByID,
                  const int                                                       a_numRanks) noexcept;

/**
 * @brief Judge every incoming proposal that targets a particle this rank owns, and produce the
 * verdicts to send back.
 * @details The full rule set this function implements, summarized:
 *   1. Group incoming proposals by target id.
 *   2. For a target with NO outgoing commitment of its own: the eligible candidate is the
 *      minimum-distance incoming proposal (deterministic id tiebreak on exact ties).
 *   3. For a target that DOES have its own outgoing commitment (i.e. it is simultaneously
 *      proposing to something else): normally reject every incoming proposal unconditionally
 *      (it is "spoken for") -- EXCEPT when its own outgoing target is itself one of its current
 *      incoming proposers, a genuinely mutual best-match pair. Blindly rejecting that case
 *      causes a PERMANENT deadlock (a stable mutual pair proposes to each other identically
 *      every round and is rejected every round, forever, since nothing changes to break the
 *      symmetry) -- so it must be accepted instead. Because BOTH sides of a mutual pair detect
 *      the match independently, a deterministic tiebreak (lower-id side performs the accept,
 *      higher-id side falls back to the normal reject) is required, or both sides would accept
 *      and double-merge the same two particles.
 *   4. ALL eligible target decisions (plain-argmin and mutual-match alike) are resolved in ONE
 *      globally distance-sorted pass across every target this rank is judging (with an id
 *      tiebreak) -- NOT by iterating targets in whatever order the incoming-proposal list
 *      happens to be in. Accepting one target decrements its cell's live count, which can push
 *      a DIFFERENT target sharing that cell below threshold; judging targets in arbitrary order
 *      makes that outcome depend on arbitrary processing order. Sorting closest-first with a
 *      dynamic per-target threshold recheck (mirroring resolveTrivialTier()'s discipline
 *      exactly) removes that dependency.
 *   5. This handles any number of simultaneous suitors for one particle, not just pairwise
 *      boundaries -- a corner where 3+ patches meet is resolved by this exact same logic, with
 *      no special-casing, since every proposal naming a given target funnels into the same
 *      "group by target id" step regardless of how many different ranks it came from.
 *   6. A_maxCellDistance and a_isPositionValid apply here exactly as they do in
 *      resolveTrivialTier() -- a physical-validity rule, not a coordination one. An incoming
 *      proposal that is too far, or whose weighted-centroid position fails a_isPositionValid, is
 *      never eligible to be accepted (skipped as though it did not exist for argmin purposes),
 *      REGARDLESS of the mutual-match exception in rule 3 -- an invalid position is invalid on
 *      both sides of a pair, mutual match or not.
 *
 * Proposals addressed to this rank's own particles arrive from two sources, pooled identically
 * here: proposals received over MPI from other ranks, and this rank's own
 * generateProposals()-produced bucket addressed to itself (an in-process handoff, no network
 * round-trip -- see generateProposals()).
 *
 * @tparam Packed See MergeParticle.
 * @tparam Combine See PayloadCombine.
 * @tparam PosValid See PositionValid.
 * @param[in,out] a_liveCellCount Live per-cell counts for this rank's own patches; decremented in place as targets
 * are accepted, in the same distance-sorted processing order (see rule 4 above). Used to dynamically recheck BOTH
 * the target's AND the source's own cell -- mirrors resolveTrivialTier()'s symmetric recheck exactly; accepting a
 * proposal whose source's own cell is no longer over threshold (e.g. drained by an earlier merge within this same
 * round) would silently violate the crowding trigger on the source side. If this rank has no live visibility into
 * the source's cell (it was never shipped here as a ghost), the source count defaults to 0 -- i.e. rejected -- the
 * same conservative default already used when a target's own count is unknown.
 * @param[out] a_results One entry per target actually accepted, in the order processed.
 * @param[out] a_verdictsByDestRank One bucket per proposer rank (pre-sized to the total rank count by the caller);
 * every proposal (accepted or not) gets exactly one verdict, addressed back to NNMergeProposal::source's owning
 * rank.
 * @param[in] a_incomingProposals All proposals targeting a particle this rank owns, pooled from every source rank
 * (including this rank's own in-process bucket).
 * @param[in] a_particlesByID Full particle data for every currently-alive local valid particle and every ghost across
 * all of this rank's patches, keyed by id (each value an NNPooledParticle wrapping the particle) -- a target is a
 * potential judge target only if its id is present here AND absent from a_consumedIDs (see that param).
 * @param[in] a_consumedIDs Ids already consumed this round (by an earlier trivial or judged merge). A target id
 * present in a_particlesByID but ALSO in a_consumedIDs is treated as not-alive -- this replaces what used to be a
 * separate, pre-filtered a_myParticlesByID copy of a_particlesByID; checking both containers at the point of use
 * instead avoids rebuilding that copy every round.
 * @param[in] a_hasOutgoingCommitment Ids (among a_particlesByID, not in a_consumedIDs) that have their own outgoing
 * commitment this round (needed for rule 3 above).
 * @param[in] a_outgoingTargetOf For an id in a_hasOutgoingCommitment, the id of the particle it is itself proposing
 * to (needed to detect the mutual-match exception in rule 3).
 * @param[in] a_combine Opaque-payload combine callback (see PayloadCombine).
 * @param[in] a_cellBudget Per-cell budget (the same budget used everywhere else this round).
 * @param[in] a_probLo Physical position of the domain's low corner (see cellKeyOf()).
 * @param[in] a_dxByLevel Cell spacing PER AMR LEVEL -- see resolveTrivialTier()'s a_dxByLevel docs for why this must
 * never be a single shared dx (a proposal's source and this rank's target can be on different levels).
 * @param[in] a_isPositionValid Position-validity predicate (see PositionValid and rule 6 above).
 * @param[in] a_allocateID Fresh-id allocator (see IDAllocator), called once per accepted target.
 * @param[in] a_maxCellDistance Same physical merge-eligibility cap as resolveTrivialTier()'s (see rule 6 above);
 * std::nullopt (default) disables it.
 * @param[in,out] a_crossLevelMergeCount Optional diagnostic counter -- see resolveTrivialTier()'s docs; incremented
 * once per accepted target whose target/source were on different AMR levels.
 * @param[in] a_onMergeCommitted Optional per-merge diagnostic callback -- see resolveTrivialTier()'s docs; called
 * with (target, source, merged) for every accepted target.
 */
template <typename Packed, typename Combine, typename PosValid, typename Allocator, typename OnMerge>
inline void
judgeProposals(
  std::unordered_map<NNCellKey, int, NNCellKeyHasher>&            a_liveCellCount,
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
  const std::optional<int>                                        a_maxCellDistance      = std::nullopt,
  unsigned long long*                                             a_crossLevelMergeCount = nullptr,
  const OnMerge&                                                  a_onMergeCommitted =
    [](const MergeParticle<Packed>&, const MergeParticle<Packed>&, const MergeParticle<Packed>&) {
    }) noexcept;

/**
 * @brief Apply every verdict addressed to one of this rank's own outgoing proposals: remove
 * the source particle for every accepted proposal, do nothing for rejected ones.
 * @details Rejected proposals are deliberately a no-op -- the source particle is simply left
 * alone and reconsidered in next round's fresh candidate search. This algorithm never internally
 * retries a rejected proposal against a second-best candidate within the same round.
 * @param[out] a_consumedIDs Ids of every source particle whose proposal was accepted -- the caller removes these
 * from its own valid particle lists.
 * @param[in] a_incomingVerdicts All verdicts addressed to this rank, pooled from every judge rank (including this
 * rank's own in-process verdict for proposals it judged for itself).
 */
inline void
applyVerdicts(std::unordered_set<ParticleID>&    a_consumedIDs,
              const std::vector<NNMergeVerdict>& a_incomingVerdicts) noexcept;

/**
 * @brief Locate-by-position callback -- see placeMergedParticles().
 * @details Returns LevelTiles::LevelAndBox, which is what
 * ParticleContainer::findDestination() already produces: an O(1) point->block query on the realm's
 * LevelTiles tile->box hash, the SAME lookup a ParticleContainer uses internally to route its own
 * particles during remap(), never a linear scan over patches. Going through that one mapping keeps
 * it the single source of truth for "which box owns this position".
 *
 * The callback itself remains, so this header does not depend on ParticleContainer; only the result
 * type is shared. Call frequency is bounded by the number of MERGE RESULTS this round (not by the
 * number of candidate pairs), so a non-inlined indirect call is an acceptable cost -- this is not a
 * hot inner-loop path.
 * @param[in] a_position Physical position to locate.
 * @return Which box (level, global grid index) and rank covers it; @c valid is false off-domain.
 */
using PositionLocator = std::function<LevelTiles::LevelAndBox(const RealVect& a_position)>;

/**
 * @brief Point-locate a newly created merged particle against the AMR grid, and bucket it for a
 * scatter exchange if it belongs to someone else.
 * @details Applied UNIFORMLY to every merge result, trivial-tier and judge-decided alike: the
 * union of two patches is not generally convex (T-junctions, corners, refinement boundaries), so
 * even a purely local (same-rank) merge's result can legitimately belong to a third patch,
 * possibly owned by a different rank. Point-location must be O(1) per result (see
 * PositionLocator) -- NEVER a linear scan over every patch, which would turn
 * this function into a real bottleneck at large patch counts.
 * @tparam Packed See MergeParticle.
 * @param[out] a_localByPatch Results whose merged particle lands in one of this rank's own patches, bucketed by
 * (level, LevelTiles::LevelAndBox::gridIndex).
 * @param[out] a_scatterByDestRank Results whose merged particle belongs to a different rank, bucketed by destination
 * rank (pre-sized to the total rank count by the caller), for one final exchange.
 * @param[in] a_results Every merge result produced this round by resolveTrivialTier() and judgeProposals(), pooled
 * together.
 * @param[in] a_locate Locate-by-position callback (see PositionLocator).
 * @param[in] a_thisRank This rank's own number -- a located box whose rank equals this is what makes a
 * result "local" rather than needing a scatter.
 */
template <typename Packed>
inline void
placeMergedParticles(std::map<std::pair<int, unsigned int>, std::vector<MergeParticle<Packed>>>& a_localByPatch,
                     std::vector<std::vector<MergeParticle<Packed>>>&                            a_scatterByDestRank,
                     const std::vector<NNMergeResult<Packed>>&                                   a_results,
                     const PositionLocator&                                                      a_locate,
                     const int                                                                   a_thisRank) noexcept;

/**
 * @brief One patch's identity for the duration of a merge round: its AMR level, local DataIndex,
 * and its own level's cell spacing.
 * @details Shared by every backend (mergeNearestNeighborsTree()/mergeNearestNeighborsOneCell())
 * as the parallel index into their own per-patch vectors (patchLocalValid/patchGhosts/...).
 */
struct NNPatchWork
{
  /**
   * @brief AMR level this patch lives on.
   */
  int level;

  /**
   * @brief This patch's local DataIndex on its level.
   */
  DataIndex din;

  /**
   * @brief Cell spacing of this patch's own level.
   */
  RealVect dx;
};

/**
 * @brief Gather every patch this rank owns into the backend-agnostic per-patch particle sets and
 * pooled bookkeeping a merge round needs, shared by every search backend.
 * @details For every patch: splits its leaf into local-valid/ghost MergeParticle vectors (ghosts
 * carrying their TRUE origin AMR level, see MergeParticle::level), computes each local-valid
 * particle's boundary-exposure flag from the realm's same/coarse-to-fine/fine-to-coarse particle
 * ghost masks at the given a_ghostWidth, and pools everything rank-wide into a_particlesByID.
 * a_liveCellCount is derived from the pooled a_particlesByID afterward (once per particle, not
 * once per occurrence, since a ghost can appear in several of this rank's own patches at once).
 * Builds no spatial index -- that is backend-specific and is the caller's next step.
 * @tparam P SoA payload type (see mergeNearestNeighborsTree()).
 * @tparam Packed Opaque reconciliation payload type (see MergeParticle).
 * @tparam Traits Column descriptor for P.
 * @tparam Gather Callable `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed`.
 * @param[in,out] a_particles The real particle container being merged.
 * @param[out] a_patchWork One entry per patch this rank owns.
 * @param[out] a_patchLocalValid This patch's own local-valid particles, parallel to a_patchWork.
 * @param[out] a_patchGhosts Every ghost shipped to this patch, parallel to a_patchWork.
 * @param[out] a_particlesByID Every currently-alive local-valid particle and every ghost, pooled rank-wide by id.
 * @param[out] a_liveCellCount Live particle count (valid + ghost) per level-qualified cell key, derived from
 * a_particlesByID.
 * @param[in] a_amr AMR mesh -- source of ghost masks and the patch layout.
 * @param[in] a_realm Realm name (a_particles.getRealm()).
 * @param[in] a_finestLevel Finest AMR level (a_amr.getFinestLevel()).
 * @param[in] a_probLo Physical position of the domain's low corner.
 * @param[in] a_dxByLevel Cell spacing per AMR level.
 * @param[in] a_ghostWidth Particle ghost halo width already filled by the caller (see
 * mergeNearestNeighborsTree()'s own a_ghostWidth docs).
 * @param[in] a_gather Per-slot gather callback.
 */
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
                     const Gather&                                             a_gather) noexcept;

/**
 * @brief Finish a merge round: propose/judge/verdict the edges the local tier could not resolve,
 * place every merge result, and apply everything to the real particle containers -- the
 * backend-agnostic tail shared by mergeNearestNeighborsTree() and mergeNearestNeighborsOneCell().
 * @details Runs, in order: derive outgoingTargetOf from a_pooledEdges/a_hasOutgoingCommitment;
 * generateProposals(); exchange proposals; judgeProposals(); exchange verdicts; applyVerdicts();
 * placeMergedParticles(); exchange the scatter bucket; remove every consumed local-valid particle
 * from the real containers; scatter in every locally-placed and every incoming-scattered merged
 * particle. Identical to steps 3-7 of mergeNearestNeighborsTree()'s own docs.
 * @tparam P SoA payload type.
 * @tparam Packed Opaque reconciliation payload type.
 * @tparam Traits Column descriptor for P.
 * @tparam Combine See PayloadCombine.
 * @tparam Scatter Callable `(ParticleSoA<P,Traits>&, const MergeParticle<Packed>&) -> void`.
 * @tparam PosValid See PositionValid.
 * @tparam Allocator See IDAllocator.
 * @tparam OnMerge See MergeCommittedCallback.
 * @param[in,out] a_particles The real particle container being merged.
 * @param[in,out] a_particlesByID Every currently-alive local-valid particle and every ghost, pooled rank-wide by id.
 * @param[in,out] a_consumedIDs Ids already consumed this round; updated in place as judge-side merges commit.
 * @param[in,out] a_liveCellCount Live per-cell counts; decremented in place as judge-side merges commit.
 * @param[in,out] a_allResults Every merge result committed this round (trivial-tier results already inside on entry);
 * judge-side results are appended.
 * @param[in,out] a_crossLevelMergeCount Optional diagnostic counter (see MergeCommittedCallback's sibling docs).
 * @param[in] a_amr AMR mesh -- source of the LevelTiles used to place merged particles.
 * @param[in] a_realm Realm name.
 * @param[in] a_cellBudget The count each cell is drained down to -- see NNCellBudget.
 * @param[in] a_combine Opaque-payload combine callback.
 * @param[in] a_scatter Per-surviving-particle scatter callback.
 * @param[in] a_allocateID Fresh-id allocator.
 * @param[in] a_maxCellDistance Physical merge-eligibility cap forwarded to judgeProposals() (see its own docs);
 * the caller's own cap for a tree-backed round, or a hardcoded 1 for a one-cell round (see
 * findNearestNeighborCandidatesOneCell()'s docs on why that value is mandatory there).
 * @param[in] a_isPositionValid Position-validity predicate.
 * @param[in] a_probLo Physical position of the domain's low corner.
 * @param[in] a_dxByLevel Cell spacing per AMR level.
 * @param[in] a_patchWork One entry per patch this rank owns (see NNPatchWork), from gatherMergeParticles().
 * @param[in] a_pooledEdges Every candidate edge produced by the local tier's LAST candidate-search pass, pooled
 * across every patch this rank owns (used to derive outgoingTargetOf; nothing is consumed between the last pass and
 * here, so these are still current).
 * @param[in] a_hasOutgoingCommitment Query ids the local tier could not resolve trivially.
 * @param[in] a_onMergeCommitted Optional per-merge diagnostic callback.
 */
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
                 const OnMerge&                                            a_onMergeCommitted) noexcept;

/**
 * @brief Shared implementation behind mergeNearestNeighborsTree() and mergeNearestNeighborsHash(): one full round of
 * the whole-patch-tree nearest-neighbor merge, generic over which EBGeometry point-cloud index backs the per-patch
 * spatial search.
 * @details Identical algorithm and public contract to mergeNearestNeighborsTree() (see its own docs for the full
 * step-by-step description) -- the only difference between the two public entry points is which @c Cloud this
 * function is instantiated with: @c EBGeometry::PointCloudBVH<Real,ParticleID> for
 * mergeNearestNeighborsTree(), @c EBGeometry::PointCloudHashGrid<Real,ParticleID> for
 * mergeNearestNeighborsHash(). Both point-cloud types expose an identical query interface, so this one
 * implementation serves both without any backend-specific branching.
 * @tparam Cloud EBGeometry point-cloud index type (see NNSpatialIndex).
 * @tparam P,Packed,Traits,Gather,Combine,Scatter,Allocator,PosValid,OnMerge See mergeNearestNeighborsTree().
 * @param a_amr,a_particles,a_cellBudget,a_gather,a_combine,a_scatter,a_allocateID,
 * a_iterateLocalTierToConvergence,a_maxFallbackCandidates,a_maxCellDistance,a_ghostWidth,a_isPositionValid,
 * a_crossLevelMergeCount,a_onMergeCommitted See mergeNearestNeighborsTree().
 */
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
                               const OnMerge&                a_onMergeCommitted);

} // namespace detail

/**
 * @brief Run one full round (one timestep's worth) of the nearest-neighbor merge algorithm.
 * @details Orchestrates, in order:
 *   1. Ghost-fill -- fresh, exactly once, strictly before any of the below, via
 *      ParticleContainer::fillGhostParticles() with the AmrMesh particle ghost masks. This
 *      function does not fill ghosts itself: the caller fills them first, and must have registered
 *      a width-1 particle ghost mask before the grids were (re)built (a mask is only built during a
 *      regrid, so registering it late leaves it empty). See
 *      ItoSolver::mergeSuperparticlesNearestNeighbor for a concrete caller.
 *   2. Per patch: findNearestNeighborCandidates(), then resolveTrivialTier(). If
 *      a_iterateLocalTierToConvergence is true, this pair repeats (reusing the step-1 ghost-fill)
 *      until no further local commits occur. See a_iterateLocalTierToConvergence below.
 *   3. generateProposals() for everything the trivial tier could not resolve.
 *   4. Exchange proposals via nnMergeExchangeByRank() (an MPI_Alltoall/Alltoallv helper): this
 *      rank's own bucket (addressed to itself) is handed directly to judgeProposals() in-process,
 *      every other bucket is shipped to its destination rank and the incoming proposals collected.
 *   5. judgeProposals() over the pooled incoming proposals.
 *   6. Exchange verdicts (nnMergeExchangeByRank(), the other direction) and applyVerdicts().
 *   7. placeMergedParticles(), then exchange and apply the scatter bucket (nnMergeExchangeByRank()).
 *
 * Fresh, globally-unique ids for newly created particles must be allocated by the caller such that
 * they can never collide with an existing particle's id or another rank's freshly allocated id
 * within the same round; a rank-namespaced counter suffices (see a_allocateID / IDAllocator).
 *
 * This function does a SINGLE cross-patch exchange (propose/judge/verdict happens exactly
 * once) regardless of a_iterateLocalTierToConvergence: a cell more than roughly 2x over threshold
 * is not necessarily fully drained by the cross-patch protocol in one call, and the caller is
 * still expected to invoke this once per timestep, relying on repeated calls (and, in a real
 * simulation, particle motion between them) for further convergence of anything the cross-patch
 * protocol itself left unresolved. What a_iterateLocalTierToConvergence changes is only the LOCAL
 * (trivial-tier) portion of a single call -- see its own documentation below.
 *
 * @tparam P SoA payload type (e.g. ItoParticle).
 * @tparam Packed Opaque reconciliation payload type for the merge decision layer (see MergeParticle) -- typically
 * much smaller than P's full column set (e.g. just energy for ItoParticle; velocity/mobility/diffusion are already
 * established as safe to drop across a merge, since they are recomputed from field interpolation immediately
 * afterward).
 * @tparam Traits Column descriptor for P (defaults to ParticleTraits<P>).
 * @tparam Gather Callable `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed` -- packs one SoA slot's payload
 * columns into the opaque Packed representation.
 * @tparam Combine See PayloadCombine.
 * @tparam Scatter Callable `(ParticleSoA<P,Traits>&, const MergeParticle<Packed>&) -> void` -- appends one
 * surviving/merged particle back into the real SoA container.
 * @param[in,out] a_particles The real particle container to merge, for every patch this rank owns.
 * @param[in,out] a_crossLevelMergeCount Optional diagnostic counter (nullptr, the default, disables it entirely):
 * incremented once for every merge this call commits (trivial-tier or judge-side) whose two participants were on
 * DIFFERENT AMR levels. NOT reset by this function -- accumulates across a_iterateLocalTierToConvergence's internal
 * iterations and across repeated calls if the caller keeps reusing the same counter, so it is safe to use to track
 * cross-level activity over several rounds.
 * @param[in] a_amr AMR mesh -- source of ghost masks, per-level dx, and the patch layout used to place merged
 * particles.
 * @param[in] a_cellBudget The count each cell is drained down to -- see NNCellBudget.
 * @param[in] a_gather Per-slot gather callback (see Gather).
 * @param[in] a_combine Opaque-payload combine callback (see PayloadCombine).
 * @param[in] a_scatter Per-surviving-particle scatter callback (see Scatter).
 * @param[in] a_allocateID Fresh-id allocator (see IDAllocator). Mandatory -- there is no universally safe default; the
 * caller must supply a scheme that stays collision-free across every rank and every past/future call.
 * @param[in] a_iterateLocalTierToConvergence Default false. When false, findNearestNeighborCandidates() and
 * resolveTrivialTier() each run once. When true, that candidate-search-plus-trivial-tier pair repeats (reusing the
 * same ghost-fill from step 1; only which local particles are still alive changes) until a full local pass commits
 * zero further trivial merges, and only then does generateProposals()/judgeProposals() run (still exactly once,
 * unaffected by this option). It helps because a single pass fixes each query's edge once: if that candidate is
 * consumed by a closer edge first, the query is left unmatched for the rest of the pass even though another partner
 * may exist nearby. The gain is modest -- re-running only helps when a particle's best available choice actually
 * changes between iterations; a_maxFallbackCandidates below addresses the more common, stable loss.
 * @param[in] a_maxFallbackCandidates Forwarded to resolveTrivialTier() (see its docs). Default 0 (disabled). Lets a
 * query blocked from a good local partner (because that partner is already "busy" with its own outgoing commitment --
 * a condition iteration alone does not clear) fall back to its next nearest candidates. A single fallback recovers
 * most of the achievable improvement.
 * @param[in] a_maxCellDistance Forwarded to resolveTrivialTier() (see its docs). `std::nullopt` (default) disables it.
 * A physical cap, in whole grid cells, on how far apart two particles may be and still merge -- e.g. 1 permits merging
 * within the same cell or any Moore-adjacent cell, never farther. Composes with the options above.
 * @param[in] a_ghostWidth Width, in whole cells, of the particle ghost halo the CALLER has already filled and of the
 * registered particle ghost mask this function fetches to compute boundary-exposure. Must equal the filled/registered
 * width and be identical on every rank -- that is what makes a trivially-merged (non-exposed) pair provably invisible
 * to every other patch, so it is the load-bearing invariant against double-merging. Double-merge safety does NOT
 * depend on a_maxCellDistance; a_ghostWidth only needs to be >= a_maxCellDistance for cross-patch merges to actually
 * reach that far (a candidate that many cells across a boundary is a merge candidate only if it lies inside the halo).
 * Default 1, which -- with the default unbounded a_maxCellDistance -- caps cross-patch merges at one cell.
 * @param[in] a_isPositionValid Forwarded to resolveTrivialTier()/judgeProposals() (see PositionValid). Default
 * always-true. Composes with the options above -- e.g. an embedded-boundary point query for an EB-aware caller (see
 * ItoSolver::mergeSuperparticlesNearestNeighbor).
 * @param[in] a_onMergeCommitted Optional per-merge diagnostic callback (see MergeCommittedCallback); default is a
 * no-op.
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
             void(const MergeParticle<Packed>&, const MergeParticle<Packed>&, const MergeParticle<Packed>&)>>
inline void
mergeNearestNeighborsTree(
  ParticleContainer<P, Traits>& a_particles,
  const AmrMesh&                a_amr,
  const NNCellBudget&           a_cellBudget,
  const Gather&                 a_gather,
  const Combine&                a_combine,
  const Scatter&                a_scatter,
  const Allocator&              a_allocateID,
  const bool                    a_iterateLocalTierToConvergence = false,
  const int                     a_maxFallbackCandidates         = 0,
  const std::optional<int>      a_maxCellDistance               = std::nullopt,
  const int                     a_ghostWidth                    = 1,
  const PosValid&               a_isPositionValid =
    [](const RealVect&) {
      return true;
    },
  unsigned long long* a_crossLevelMergeCount = nullptr,
  const OnMerge&      a_onMergeCommitted =
    [](const MergeParticle<Packed>&, const MergeParticle<Packed>&, const MergeParticle<Packed>&) {
    });

/**
 * @brief Run one full round of the nearest-neighbor merge algorithm, identical in every respect to
 * mergeNearestNeighborsTree() except that the per-patch spatial search is backed by an
 * EBGeometry::PointCloudHashGrid (a uniform-grid index) instead of an EBGeometry::PointCloudBVH.
 * @details See mergeNearestNeighborsTree() for the full algorithm description, parameter contracts, and
 * defaults -- every one of them applies here unchanged; only the underlying point-cloud index differs (see
 * detail::mergeNearestNeighborsRoundImpl(), which both this function and mergeNearestNeighborsTree() forward to,
 * selecting the point-cloud type via its @c Cloud template parameter). PointCloudHashGrid trades PointCloudBVH's
 * adaptivity to non-uniform density for an O(1) expected build/query cost at near-uniform density; prefer it when a
 * crowded cell's neighborhood is not also badly clustered in density, and measure before relying on it being faster
 * in a given workload (see PointCloudHashGrid's own docs in EBGeometry for when it is and is not a good fit).
 * @tparam P,Packed,Traits,Gather,Combine,Scatter,Allocator,PosValid,OnMerge See mergeNearestNeighborsTree().
 * @param a_amr,a_particles,a_cellBudget,a_gather,a_combine,a_scatter,a_allocateID,
 * a_iterateLocalTierToConvergence,a_maxFallbackCandidates,a_maxCellDistance,a_ghostWidth,a_isPositionValid,
 * a_crossLevelMergeCount,a_onMergeCommitted See mergeNearestNeighborsTree().
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
             void(const MergeParticle<Packed>&, const MergeParticle<Packed>&, const MergeParticle<Packed>&)>>
inline void
mergeNearestNeighborsHash(
  ParticleContainer<P, Traits>& a_particles,
  const AmrMesh&                a_amr,
  const NNCellBudget&           a_cellBudget,
  const Gather&                 a_gather,
  const Combine&                a_combine,
  const Scatter&                a_scatter,
  const Allocator&              a_allocateID,
  const bool                    a_iterateLocalTierToConvergence = false,
  const int                     a_maxFallbackCandidates         = 0,
  const std::optional<int>      a_maxCellDistance               = std::nullopt,
  const int                     a_ghostWidth                    = 1,
  const PosValid&               a_isPositionValid =
    [](const RealVect&) {
      return true;
    },
  unsigned long long* a_crossLevelMergeCount = nullptr,
  const OnMerge&      a_onMergeCommitted =
    [](const MergeParticle<Packed>&, const MergeParticle<Packed>&, const MergeParticle<Packed>&) {
    });

/**
 * @brief Run one full round of the one-cell nearest-neighbor merge: the same distributed
 * propose/judge/verdict protocol as mergeNearestNeighborsTree(), but with candidates found via a
 * per-cell PointCloudBVH spatial index (see NNCellSpatialIndex) instead of one whole-patch tree.
 * @details Orchestrates, in order: 1. gatherMergeParticles() (ghosts assumed already filled by the
 * caller, exactly as mergeNearestNeighborsTree() requires -- see its own docs). 2. Per patch:
 * buildNNCellSpatialIndex(), then findNearestNeighborCandidatesOneCell() + resolveTrivialTier(),
 * optionally repeated to convergence (see a_iterateLocalTierToConvergence). 3. finishMergeRound()
 * for everything the local tier could not resolve.
 *
 * The merge distance is structurally fixed at Chebyshev cell distance 1 -- only the query's own
 * cell and its Moore-adjacent neighbors are ever visited (see
 * findNearestNeighborCandidatesOneCell()) -- so, unlike mergeNearestNeighborsTree(), this function
 * takes no a_maxCellDistance parameter, and requires exactly a width-1 particle ghost halo (already
 * filled by the caller) rather than a configurable a_ghostWidth.
 * @tparam P SoA payload type (e.g. ItoParticle).
 * @tparam Packed Opaque reconciliation payload type (see MergeParticle).
 * @tparam Traits Column descriptor for P (defaults to ParticleTraits<P>).
 * @tparam Gather Callable `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed`.
 * @tparam Combine See PayloadCombine.
 * @tparam Scatter Callable `(ParticleSoA<P,Traits>&, const MergeParticle<Packed>&) -> void`.
 * @param[in,out] a_particles The real particle container to merge, for every patch this rank owns.
 * @param[in] a_amr AMR mesh -- source of ghost masks, per-level dx, and the patch layout.
 * @param[in] a_cellBudget The count each cell is drained down to -- see NNCellBudget.
 * @param[in] a_gather Per-slot gather callback (see Gather).
 * @param[in] a_combine Opaque-payload combine callback (see PayloadCombine).
 * @param[in] a_scatter Per-surviving-particle scatter callback (see Scatter).
 * @param[in] a_allocateID Fresh-id allocator (see IDAllocator). Mandatory, same requirement as
 * mergeNearestNeighborsTree().
 * @param[in] a_iterateLocalTierToConvergence Default false. See mergeNearestNeighborsTree()'s own docs.
 * @param[in] a_maxFallbackCandidates Forwarded to resolveTrivialTier(). Default 0 (disabled). See
 * mergeNearestNeighborsTree()'s own docs.
 * @param[in] a_isPositionValid Forwarded to resolveTrivialTier()/judgeProposals() (see PositionValid). Default
 * always-true.
 * @param[in,out] a_crossLevelMergeCount Optional diagnostic counter (nullptr, the default, disables it). See
 * mergeNearestNeighborsTree()'s own docs.
 * @param[in] a_onMergeCommitted Optional per-merge diagnostic callback (see MergeCommittedCallback); default is a
 * no-op.
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
             void(const MergeParticle<Packed>&, const MergeParticle<Packed>&, const MergeParticle<Packed>&)>>
inline void
mergeNearestNeighborsOneCell(
  ParticleContainer<P, Traits>& a_particles,
  const AmrMesh&                a_amr,
  const NNCellBudget&           a_cellBudget,
  const Gather&                 a_gather,
  const Combine&                a_combine,
  const Scatter&                a_scatter,
  const Allocator&              a_allocateID,
  const bool                    a_iterateLocalTierToConvergence = false,
  const int                     a_maxFallbackCandidates         = 0,
  const PosValid&               a_isPositionValid =
    [](const RealVect&) {
      return true;
    },
  unsigned long long* a_crossLevelMergeCount = nullptr,
  const OnMerge&      a_onMergeCommitted =
    [](const MergeParticle<Packed>&, const MergeParticle<Packed>&, const MergeParticle<Packed>&) {
    });

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#include <CD_NearestNeighborParticleMergeImplem.H>

#endif
