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
#include <Box.H>
#include <FArrayBox.H>
#include <RealVect.H>

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_ParticleContainer.H>
#include <CD_ParticleManagement.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace ParticleManagement {

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
 * in mergeKDCarve()/mergeKDPatch()'s own docs).
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
 * leaf's own members/box (see mergeKDCarve()'s implementation) so there is exactly one
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
 * The quota constrains only leaves that can actually merge. mergeKDCarve() and mergeKDPatch()
 * veto any leaf
 * wider than @c s_kdMaxLeafExtent and leaves its members behind individually, so such a leaf
 * costs its cell its full member count rather than one particle; the quota therefore does not block
 * splitting it further. As a result a cell can finish slightly above @a a_ppc when a leaf is both
 * too wide to merge and in a cell already at quota, since both available choices overshoot.
 *
 * Particles in cells at or below @a a_ppc are excluded from the partition and emitted as singleton
 * leaves. A singleton never merges (the merge threshold is two members), so an uncrowded cell cannot
 * be drained by a neighbour's merge.
 * @tparam Packed See MergeParticle.
 * @param[in,out] a_particles Particles to partition; reordered in place.
 * @param[in,out] a_used Live per-cell leaf quota for this patch, zeroed on entry. Its box must cover
 * every cell a leaf centroid can fall in, i.e. it needs at least one ghost cell.
 * @param[out] a_leaves Every resulting leaf, in no particular order.
 * @param[in] a_ppc Target particles per cell, enforced here as a per-cell ceiling on leaf count.
 * @param[in] a_splitWeightLeafDx Leaf size, in cell widths, at or below which the split plane
 * switches from the count median to the weight median. 0 disables the weight median.
 * @param[in] a_dx This patch's own cell spacing.
 * @param[in] a_probLo Physical position of the domain's low corner, for computing cell indices.
 * @param[in] a_cellCounts Per-cell TRUE particle-count histogram (see kdFillCellHistogram()). A cell
 * is crowded, and its particles eligible to merge, iff its count exceeds @a a_ppc. A cell outside
 * this box reads as 0 (empty): unlike @a a_used, whose box coverage is a hard precondition, this
 * lookup guards rather than asserts, since a particle's own cell is always in-box but the guard
 * costs nothing and documents the intent.
 */
template <typename Packed>
inline void
buildKDQuotaLeaves(std::vector<MergeParticle<Packed>>& a_particles,
                   FArrayBox&                          a_used,
                   std::vector<KDLeaf>&                a_leaves,
                   const int                           a_ppc,
                   const Real                          a_splitWeightLeafDx,
                   const RealVect&                     a_dx,
                   const RealVect&                     a_probLo,
                   const FArrayBox&                    a_cellCounts) noexcept;

/**
 * @brief Generic Alltoallv-style exchange for a trivially-copyable record type: send a
 * per-destination-rank bucket, return the flat concatenation of everything received.
 * @details A deliberate, self-contained duplicate of the nearest-neighbor merge's
 * nnMergeExchangeByRank() (same MPI_Alltoall/Alltoallv shape) rather than a shared dependency --
 * that function lives in a header that transitively pulls in EBGeometry, which this file must
 * not depend on (see this file's own @details). In a serial (non-MPI) build, or with a single
 * rank, returns this rank's own bucket unchanged.
 * @tparam T Trivially-copyable record type (KDClaim, KDVerdict, KDCommit, or
 * MergeParticle<Packed> with Packed itself trivially copyable).
 * @param[in] a_sendByRank One bucket per destination rank (size == number of ranks).
 * @return Flat vector of every record received from every rank (including this rank's own bucket
 * addressed to itself).
 */
template <typename T>
inline std::vector<T>
kdExchangeByRank(const std::vector<std::vector<T>>& a_sendByRank);

/**
 * @brief Neighbor-restricted counterpart of kdExchangeByRank().
 * @details Two Isend/Irecv rounds -- record counts, then payloads -- over a_neighborRanks only, in place
 * of the MPI_Alltoall/Alltoallv pair. A rank outside the set is known to send nothing, so its count is
 * never asked for. Valid ONLY when every destination is bounded by that set; a destination outside it
 * would never be received, which for these buckets means a silently dropped merge and lost weight.
 * Asserted on the sending side in debug builds. The result is concatenated in rank order, this rank's
 * own bucket included at its own position, so it is indistinguishable from what kdExchangeByRank()
 * returns -- which matters, because the argmin arbitration downstream tie-breaks on record order.
 * @tparam T Trivially-copyable record type, as for kdExchangeByRank().
 * @param[in] a_sendByRank    One bucket per destination rank (size == number of ranks).
 * @param[in] a_neighborRanks Ranks this rank may exchange with.
 * @return Flat vector of every record received, plus this rank's own bucket.
 */
template <typename T>
inline std::vector<T>
kdExchangeByRankNeighbors(const std::vector<std::vector<T>>& a_sendByRank, const std::vector<int>& a_neighborRanks);

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
 * @tparam Packed Opaque reconciliation payload type (see MergeParticle).
 * @tparam Traits Column descriptor for P (defaults to ParticleTraits<P>).
 * @tparam Gather Callable `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed`.
 * @tparam Combine See KDPayloadCombine -- N-ary, not pairwise.
 * @tparam Scatter Callable `(ParticleSoA<P,Traits>&, const MergeParticle<Packed>&) -> void`.
 * @tparam Allocator Callable `() -> ParticleID`, a fresh globally-unique id per committed merge.
 * @tparam PosValid Callable `(const RealVect&) -> bool`; false leaves the group's members untouched.
 * @param[in] a_amr AMR mesh -- source of the particle ghost masks and grid layout.
 * @param[in,out] a_particles The real particle container being merged.
 * @param[in,out] a_cellHistogram Caller-owned per-cell scratch (1 component, >= 1 ghost cell, this
 * container's realm) holding the cell occupancy. Contents on entry are irrelevant -- it is overwritten
 * per patch. Owned by the caller so the allocation happens once per regrid rather than once per call.
 * @param[in,out] a_leafQuota Caller-owned per-cell scratch, same shape and lifetime as
 * @a a_cellHistogram, holding the live per-cell merge quota.
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
 * @tparam Packed Opaque reconciliation payload type (see MergeParticle).
 * @tparam Traits Column descriptor for P (defaults to ParticleTraits<P>).
 * @tparam Gather Callable `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed`.
 * @tparam Combine See KDPayloadCombine -- N-ary, not pairwise.
 * @tparam Scatter Callable `(ParticleSoA<P,Traits>&, const MergeParticle<Packed>&) -> void`.
 * @tparam Allocator Callable `() -> ParticleID`, a fresh globally-unique id per committed merge.
 * @tparam PosValid Callable `(const RealVect&) -> bool`; false leaves the leaf's members untouched.
 * @param[in,out] a_particles The real particle container being merged.
 * @param[in,out] a_cellHistogram Caller-owned per-cell scratch; see mergeKDCarve()'s own docs.
 * @param[in,out] a_leafQuota Caller-owned per-cell scratch; see mergeKDCarve()'s own docs.
 * @param[in] a_amr AMR mesh -- source of the grid layout.
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
             const PosValid&               a_isPositionValid);

/**
 * @brief Run the uncontested tier of the kd merge, splitting the input into merged and leftover.
 * @details Same tree build and per-cell quota as mergeKDCarve(), but with no arbitration tier at
 * all. In order:
 *   1. Build one kd tree per patch over that patch's own particles plus every ghost, which the
 *      caller must already have filled at width 1. Ghosts are never merged, but they take part in
 *      the partition: a ghost marks the leaf it lands in as contested, and the quota needs the
 *      cell's true occupancy.
 *   2. Commit every leaf that holds at least two members, is within detail::s_kdMaxLeafExtent, holds
 *      no ghost, and whose centroid @a a_isPositionValid accepts. Each such leaf becomes one
 *      super-particle in @a a_interior and its members are removed from @a a_particles.
 *   3. Everything else is left untouched in @a a_particles -- leaves holding a ghost, leaves too
 *      wide to collapse, and leaves with a single member.
 *
 * "Holds no ghost" is the whole safety condition, and it needs no communication to evaluate. Every
 * member of such a leaf is resident in this patch, and any other patch that draws one of those
 * members into a leaf of its own sees it as a ghost, so that leaf fails this same test there. No two
 * patches can commit the same particle and no weight is counted twice. Boundary exposure is
 * deliberately not consulted: an exposed but uncontested leaf merges here.
 *
 * The result is a split, not a finished merge -- @a a_particles is left holding exactly the
 * contested particles, for a caller that means to reduce them by some other means.
 * @tparam P SoA payload type.
 * @tparam Packed Opaque reconciliation payload type (see MergeParticle).
 * @tparam Traits Column descriptor for P (defaults to ParticleTraits<P>).
 * @tparam Gather Callable `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed`.
 * @tparam Combine See KDPayloadCombine -- N-ary, not pairwise.
 * @tparam Scatter Callable `(ParticleSoA<P,Traits>&, const MergeParticle<Packed>&) -> void`.
 * @tparam Allocator Callable `() -> ParticleID`, a fresh globally-unique id per committed merge.
 * @tparam PosValid Callable `(const RealVect&) -> bool`; false leaves the leaf's members untouched.
 * @param[in,out] a_particles Particles to reduce. On return holds only what was not committed.
 * @param[in,out] a_interior Receives one super-particle per committed leaf. Must be allocated on the
 * same realm as @a a_particles, since a leaf is written to the patch it was built on.
 * @param[in,out] a_cellHistogram Caller-owned per-cell scratch; see mergeKDCarve()'s own docs.
 * @param[in,out] a_leafQuota Caller-owned per-cell scratch; see mergeKDCarve()'s own docs.
 * @param[in] a_amr AMR mesh -- source of the grid layout.
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
                const PosValid&               a_isPositionValid);

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#include <CD_KDParticleMergeImplem.H>

#endif
