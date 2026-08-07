/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_KDParticleMerge.H
 * @brief  Declaration of distributed, MPI-safe whole-patch kd-tree super-particle merge algorithms.
 * @author Robert Marskar
 */

#ifndef CD_KDPARTICLEMERGE_H
#define CD_KDPARTICLEMERGE_H

// Std includes
#include <cstdint>
#include <functional>
#include <vector>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_ParticleContainer.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace ParticleManagement {

/**
 * @brief Minimal, payload-agnostic description of one particle as input to the kd-tree carve
 * merge.
 * @details Carries what the tree build/carve protocol needs -- position, weight, a globally
 * unique id, and the particle's current owning rank -- plus an opaque, caller-defined payload
 * (e.g. energy). The build/carve logic never inspects @c payload; it is unpacked only by the
 * caller-supplied combine callback when a leaf/box is finalized. No AMR-level field is carried
 * (unlike the nearest-neighbor merge's NNMergeParticle): the tree build only ever needs physical
 * position (level-independent) and this patch's own cell spacing, never a per-particle level.
 * @tparam Packed Caller-defined payload type (e.g. a single Real energy for ItoParticle).
 */
template <typename Packed>
struct KDParticle
{
  /**
   * @brief Particle position.
   */
  RealVect position;

  /**
   * @brief Particle weight (used for the merged-particle centroid; not required to be an integer).
   */
  Real weight;

  /**
   * @brief Globally unique particle id -- also the source of a box's key anchor (see
   * detail::KDBoxKey) and the deterministic tiebreak needed since AABB volume alone is not a
   * strict total order.
   */
  ParticleID globalID;

  /**
   * @brief The rank that currently owns this particle. Always read from the particle's own data,
   * never inferred from MPI transport metadata. NOT the same question as isGhost below -- a
   * particle can be owned by this rank (via one of this rank's OTHER patches) while still being a
   * ghost in the patch whose build currently holds it.
   */
  RankID ownerRank;

  /**
   * @brief True iff this particle is a ghost copy in the patch whose build currently holds it
   * (i.e. not physically stored there). Build-internal bookkeeping, not part of the payload the
   * caller-supplied callbacks see meaningfully -- always false on a freshly merged result. This is
   * deliberately a per-particle field (not inferred from ownerRank) because ownership and physical
   * residency are different questions: a particle owned by this rank via one patch can be a ghost
   * in a NEIGHBORING patch this same rank also owns, and treating "owner == this rank" as "no
   * further contest possible" there is exactly the bug this field exists to prevent -- see
   * mergeKD()'s own docs.
   */
  bool isGhost = false;

  /**
   * @brief Opaque, caller-defined payload. Never inspected by the build/carve logic -- only
   * carried through to the caller-supplied combine callback.
   */
  Packed payload;
};

/**
 * @brief Illustrative signature of the N-ary combine callback: produce a merged particle's opaque
 * payload from an arbitrary-size group of parents' payloads and weights.
 * @details Deliberately N-ary, not pairwise like the nearest-neighbor merge's PayloadCombine --
 * leaves/boxes here are groups by construction (kd-tree leaves, carve boxes), not pairs off a
 * k-NN graph, so the natural callback shape takes a whole group at once rather than needing a
 * caller-chosen fold order for reproducibility. This alias documents the expected shape only; it
 * is NOT the type used at the call sites (consumed as a template parameter so it inlines on the
 * per-merge hot path). Position and total weight are handled by the algorithm directly (weighted
 * centroid, summed weight); this callback produces only the opaque @c Packed payload.
 * @param[in] a_payloads Pointer to the first of @c a_count parents' payloads.
 * @param[in] a_weights  Pointer to the first of @c a_count parents' weights, parallel to
 * @c a_payloads.
 * @param[in] a_count    Number of parents being combined (always >= 2 -- see the merge threshold
 * in mergeKD()'s own docs).
 * @return The merged particle's payload.
 */
template <typename Packed>
using KDPayloadCombine = std::function<
  Packed(const Packed* a_payloads, const Real* a_weights, const std::size_t a_count)>;

namespace detail {

/**
 * @brief One leaf of a kd tree: a contiguous index range into the (in-place reordered)
 * particle buffer, plus its axis-aligned bounding box.
 * @details Leaves are never grid-aligned by construction (see buildKDQuotaLeaves()) -- a
 * leaf's box can legitimately straddle a cell face. Whether a leaf is interior, a boundary
 * candidate, or unmergeable is NOT stored here; it is derived at classification time from the
 * leaf's own members/box (see mergeKD()'s implementation) so there is exactly one
 * source of truth for that decision.
 */
struct KDLeaf
{
  /**
   * @brief Start of this leaf's member range (inclusive) in the reordered particle buffer.
   */
  std::size_t lo;

  /**
   * @brief End of this leaf's member range (exclusive) in the reordered particle buffer.
   */
  std::size_t hi;

  /**
   * @brief Low corner of this leaf's axis-aligned bounding box.
   */
  RealVect boxLo;

  /**
   * @brief High corner of this leaf's axis-aligned bounding box.
   */
  RealVect boxHi;
};

/**
 * @brief Build one whole-patch kd tree: partition @a a_particles by position into leaves, each
 * of which becomes exactly one super-particle, subject to a live per-cell leaf quota of @a a_ppc.
 * @details @a a_particles is reordered IN PLACE, so each returned leaf's [lo,hi) range is contiguous
 * in the reordered buffer; no particle is ever copied or split synthetically.
 *
 * A leaf collapses to exactly ONE super-particle at its own weighted centroid, so the number of
 * leaves whose centroid falls in a cell is that cell's post-merge population. That count is tracked
 * live in a per-cell map while splitting, and a split that would push a cell past @a a_ppc is
 * refused. Splitting takes the heaviest candidate leaf first and stops when nothing may legally be
 * split.
 *
 * Splits are along the longest real-space axis and are deliberately UNSNAPPED to the grid, so a leaf
 * may straddle a cell face and draw members from either side. The split plane is the count median
 * while the node is wider than @a a_splitWeightLeafDx cell widths, and the weight median once it is
 * narrower: above that scale a split apportions leaves between cells, which is a question of counts;
 * below it the weight median drives the resulting super-particles toward equal weight.
 *
 * The quota constrains only leaves that can actually merge. mergeKD() vetoes any leaf
 * wider than @c s_kdMaxLeafExtent and leaves its members behind individually, so such a leaf
 * costs its cell its full member count rather than one particle; the quota therefore does not block
 * splitting it further. As a result a cell can finish slightly above @a a_ppc when a leaf is both
 * too wide to merge and in a cell already at quota, since both available choices overshoot.
 *
 * Particles in cells at or below @a a_ppc are excluded from the partition and emitted as singleton
 * leaves. A singleton never merges (the merge threshold is two members), so an uncrowded cell cannot
 * be drained by a neighbour's merge.
 * @tparam Packed See KDParticle.
 * @tparam CellTrueCount Callable `(const IntVect&) -> int` giving that cell's total particle count,
 * or 0 for an empty cell. A cell is crowded iff this exceeds @a a_ppc.
 * @param[in,out] a_particles Particles to partition; reordered in place.
 * @param[in] a_ppc Target particles per cell, enforced here as a per-cell ceiling on leaf count.
 * @param[in] a_splitWeightLeafDx Leaf size, in cell widths, at or below which the split plane
 * switches from the count median to the weight median. 0 disables the weight median.
 * @param[in] a_dx This patch's own cell spacing.
 * @param[in] a_probLo Physical position of the domain's low corner, for computing cell indices.
 * @param[in] a_cellTrueCount Per-cell particle-count lookup -- see @tparam CellTrueCount.
 * @param[out] a_leaves Every resulting leaf, in no particular order.
 */
template <typename Packed, typename CellTrueCount>
inline void
buildKDQuotaLeaves(std::vector<KDParticle<Packed>>& a_particles,
                   const int                        a_ppc,
                   const Real                       a_splitWeightLeafDx,
                   const RealVect&                  a_dx,
                   const RealVect&                  a_probLo,
                   const CellTrueCount&             a_cellTrueCount,
                   std::vector<KDLeaf>&             a_leaves) noexcept;

/**
 * @brief A boundary box's key: the total order the carve protocol's argmin resolves ownership
 * conflicts with.
 * @details Ascending lexicographic comparison on (volume, anchor): the tightest (smallest-volume)
 * box wins a contested particle; @c anchor (the lowest global particle id among the box's own
 * members) is the deterministic tiebreak volume alone cannot provide (every single-particle box
 * has zero volume, so volume alone is not a strict total order). Computed once by a box's
 * proposer and copied verbatim into every claim message, so comparing two keys for equality is
 * safe bit-for-bit -- no epsilon needed.
 */
struct KDBoxKey
{
  /**
   * @brief This box's AABB volume.
   */
  Real volume;

  /**
   * @brief Lowest global particle id among this box's own members.
   */
  ParticleID anchor;
};

/**
 * @brief Strict weak ordering for KDBoxKey -- see its own docs.
 * @param[in] a_lhs Left-hand key.
 * @param[in] a_rhs Right-hand key.
 * @return True iff @a a_lhs sorts strictly before @a a_rhs.
 */
inline bool
operator<(const KDBoxKey& a_lhs, const KDBoxKey& a_rhs) noexcept;

/**
 * @brief One member of a boundary box, as listed by its proposer: a global id plus that
 * particle's current owning rank.
 * @details Eligibility to join a box's final (post-carve) membership is by this explicit listing
 * only, never by raw geometric AABB containment -- a box's AABB can enclose a particle its
 * proposer never actually listed.
 */
struct KDMember
{
  /**
   * @brief This member's global particle id.
   */
  ParticleID id;

  /**
   * @brief This member's current owning rank.
   */
  RankID owner;
};

/**
 * @brief One outgoing claim: "my box, with this key, lists you as a member" -- the wire record
 * for the carve protocol's Phase 1 (claims routed to owners).
 * @details Trivially copyable (required for the MPI_Alltoallv-based exchange, see
 * kdExchangeByRank()).
 */
struct KDClaim
{
  /**
   * @brief Global id of the claimed member (the recipient rank must own this particle).
   */
  ParticleID memberID;

  /**
   * @brief The claiming box's key.
   */
  KDBoxKey key;

  /**
   * @brief The claiming box's proposer rank -- where the verdict should be sent back to.
   */
  RankID proposerRank;

  /**
   * @brief The claiming box's own index into the proposer's local box list, echoed back verbatim
   * in the reply. NOT redundant with proposerRank: two DIFFERENT boxes on the SAME rank (from two
   * different patches that same rank owns) can both claim the same foreign particle at once, and
   * the reply must let each of them individually learn whether IT specifically won -- a rank-level
   * "your rank won" answer cannot distinguish which of that rank's own competing boxes it was. See
   * mergeKD()'s own docs for the double-count this caused before this field existed.
   */
  int proposerBoxIdx;
};

/**
 * @brief One reply to a single claim: "your specific box nominally won or lost this member" --
 * the wire record for the carve protocol's Phase 2 (verdicts routed back to proposers, foreign
 * members only).
 * @details Trivially copyable (see KDClaim). Answers exactly the claim identified by
 * (memberID, claimantBoxIdx) -- one KDVerdict per incoming KDClaim, not one per rank, so
 * that two of the SAME rank's own boxes claiming the same particle each get their own, individually
 * correct answer (see KDClaim::proposerBoxIdx). NOMINAL, not final: winning the argmin for one
 * member does not by itself guarantee the winning box will actually commit -- it can independently
 * lose OTHER members to OTHER, unrelated lower-keyed competitors and end up under the merge
 * threshold itself. An owner must NOT delete a member on the strength of this verdict alone -- it
 * must wait for the corresponding KDCommit in Phase 3. See mergeKD()'s own docs
 * for a worked counterexample of the conservation bug that results from deleting here instead.
 */
struct KDVerdict
{
  /**
   * @brief Global id of the particle this verdict concerns.
   */
  ParticleID memberID;

  /**
   * @brief Echo of the claim's own KDClaim::proposerBoxIdx -- identifies which of the
   * recipient's own boxes this verdict answers for.
   */
  int claimantBoxIdx;

  /**
   * @brief True iff THIS specific box (memberID, claimantBoxIdx) won the argmin for this particle.
   */
  bool won;
};

/**
 * @brief One commit/release instruction: "the box that nominally won you did/didn't actually
 * commit" -- the wire record for the carve protocol's Phase 3 (sent by a proposer to each foreign
 * member's owner, once the proposer's own box's final outcome is known).
 * @details Trivially copyable (see KDClaim). This, not KDVerdict, is what actually
 * authorizes an owner to delete a particle -- see KDVerdict's own docs for why a nominal
 * win/loss result by itself is insufficient.
 */
struct KDCommit
{
  /**
   * @brief Global id of the particle this instruction concerns.
   */
  ParticleID memberID;

  /**
   * @brief True iff the box that nominally won this member actually committed (the owner should
   * delete its copy); false iff that box did not commit (>= 2 threshold or PosValid failed) and
   * the owner should leave this particle untouched, exactly as if it had never been claimed.
   */
  bool committed;
};

/**
 * @brief Generic Alltoallv-style exchange for a trivially-copyable record type: send a
 * per-destination-rank bucket, return the flat concatenation of everything received.
 * @details A deliberate, self-contained duplicate of the nearest-neighbor merge's
 * nnMergeExchangeByRank() (same MPI_Alltoall/Alltoallv shape) rather than a shared dependency --
 * that function lives in a header that transitively pulls in EBGeometry, which this file must
 * not depend on (see this file's own @details). In a serial (non-MPI) build, or with a single
 * rank, returns this rank's own bucket unchanged.
 * @tparam T Trivially-copyable record type (KDClaim, KDVerdict, KDCommit, or
 * KDParticle<Packed> with Packed itself trivially copyable).
 * @param[in] a_sendByRank One bucket per destination rank (size == number of ranks).
 * @return Flat vector of every record received from every rank (including this rank's own bucket
 * addressed to itself).
 */
template <typename T>
inline std::vector<T>
kdExchangeByRank(const std::vector<std::vector<T>>& a_sendByRank);

} // namespace detail

/**
 * @brief Run one non-iterative pass of the kd-tree carve merge over every patch this rank owns.
 * @details Orchestrates, in order:
 *   1. Build one kd tree per patch over its local particles plus every ghost. The caller must have
 *      filled ghosts at width 1; this function does not fill them itself.
 *   2. Classify each leaf as unmergeable (too wide to merge), interior (no boundary exposure and no
 *      ghost member), or a boundary box.
 *   3. Commit interior leaves immediately, with no communication.
 *   4. Resolve boundary boxes against competing patches in three fixed exchange rounds -- claim,
 *      nominal verdict, commit-or-release -- keyed by (AABB volume, anchor id).
 *   5. Collapse a leaf or box only if at least two members survive and @a a_isPositionValid accepts
 *      the centroid, then place the result. Anything else is left untouched.
 *
 * The round count is fixed, not a convergence loop; further reduction comes from calling again.
 * @tparam P SoA payload type.
 * @tparam Packed Opaque reconciliation payload type (see KDParticle).
 * @tparam Traits Column descriptor for P (defaults to ParticleTraits<P>).
 * @tparam Gather Callable `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed`.
 * @tparam Combine See KDPayloadCombine -- N-ary, not pairwise.
 * @tparam Scatter Callable `(ParticleSoA<P,Traits>&, const KDParticle<Packed>&) -> void`.
 * @tparam Allocator Callable `() -> ParticleID`, a fresh globally-unique id per committed merge.
 * @tparam PosValid Callable `(const RealVect&) -> bool`; false leaves the group's members untouched.
 * @param[in,out] a_amr AMR mesh -- source of the particle ghost masks and grid layout.
 * @param[in,out] a_particles The real particle container being merged.
 * @param[in] a_ppc Target particles per cell.
 * @param[in] a_splitWeightLeafDx Leaf size, in cell widths, at or below which the tree build
 * switches from count-median to weight-median splitting -- forwarded verbatim to
 * buildKDQuotaLeaves(). 0 disables the weight median entirely.
 * @param[in] a_gather Per-slot gather callback.
 * @param[in] a_combine N-ary opaque-payload combine callback.
 * @param[in] a_scatter Per-merged-particle scatter callback.
 * @param[in] a_allocateID Fresh-id allocator.
 * @param[in] a_isPositionValid Position-validity predicate.
 */
template <typename P,
          typename Packed,
          typename Traits = ParticleTraits<P>,
          typename Gather,
          typename Combine,
          typename Scatter,
          typename Allocator,
          typename PosValid>
inline void
mergeKDCarve(AmrMesh&                      a_amr,
             ParticleContainer<P, Traits>& a_particles,
             const int                     a_ppc,
             const Real                    a_splitWeightLeafDx,
             const Gather&                 a_gather,
             const Combine&                a_combine,
             const Scatter&                a_scatter,
             const Allocator&              a_allocateID,
             const PosValid&               a_isPositionValid);

/**
 * @brief Run one patch-local kd-tree merge over every patch this rank owns.
 * @details Same tree build and per-cell quota as mergeKDCarve(), with no boundary tier:
 *   1. Build one kd tree per patch over that patch's own particles. The caller must NOT have filled
 *      ghosts; a ghost here would be merged locally while its true owner merges it too.
 *   2. Every leaf whose largest per-axis extent is within detail::s_kdMaxLeafExtent commits
 *      immediately, regardless of boundary exposure. No leaf is ever contested, so there is no
 *      claim/verdict/commit exchange and no communication at all.
 *   3. A leaf collapses to one merged particle only if it has at least two members and
 *      @a a_isPositionValid accepts the centroid; otherwise its members are left untouched.
 *   4. Merged results are placed by a cheap point-in-the-owning-patch test, falling back to
 *      ParticleContainer::findDestination().
 *
 * Each patch therefore reduces only the particles it owns. Cheaper than mergeKDCarve() and free of
 * communication, at the cost of no coordination across a patch boundary: a cell straddled by two
 * patches is reduced by each of them independently.
 * @tparam P SoA payload type.
 * @tparam Packed Opaque reconciliation payload type (see KDParticle).
 * @tparam Traits Column descriptor for P (defaults to ParticleTraits<P>).
 * @tparam Gather Callable `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed`.
 * @tparam Combine See KDPayloadCombine -- N-ary, not pairwise.
 * @tparam Scatter Callable `(ParticleSoA<P,Traits>&, const KDParticle<Packed>&) -> void`.
 * @tparam Allocator Callable `() -> ParticleID`, a fresh globally-unique id per committed merge.
 * @tparam PosValid Callable `(const RealVect&) -> bool`; false leaves the leaf's members untouched.
 * @param[in,out] a_amr AMR mesh -- source of the grid layout.
 * @param[in,out] a_particles The real particle container being merged.
 * @param[in] a_ppc Target particles per cell.
 * @param[in] a_splitWeightLeafDx Leaf size, in cell widths, at or below which the tree build
 * switches from count-median to weight-median splitting -- forwarded verbatim to
 * buildKDQuotaLeaves(). 0 disables the weight median entirely.
 * @param[in] a_gather Per-slot gather callback.
 * @param[in] a_combine N-ary opaque-payload combine callback.
 * @param[in] a_scatter Per-merged-particle scatter callback.
 * @param[in] a_allocateID Fresh-id allocator.
 * @param[in] a_isPositionValid Position-validity predicate.
 */
template <typename P,
          typename Packed,
          typename Traits = ParticleTraits<P>,
          typename Gather,
          typename Combine,
          typename Scatter,
          typename Allocator,
          typename PosValid>
inline void
mergeKDPatch(AmrMesh&                      a_amr,
             ParticleContainer<P, Traits>& a_particles,
             const int                     a_ppc,
             const Real                    a_splitWeightLeafDx,
             const Gather&                 a_gather,
             const Combine&                a_combine,
             const Scatter&                a_scatter,
             const Allocator&              a_allocateID,
             const PosValid&               a_isPositionValid);

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#include <CD_KDParticleMergeImplem.H>

#endif
