/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleManagement.H
 * @brief  Namespace containing various particle management utilities.
 * @author Robert Marskar
 */

#ifndef CD_PARTICLEMANAGEMENT_H
#define CD_PARTICLEMANAGEMENT_H

// Std includes
#include <functional>
#include <string>
#include <utility>
#include <vector>

// Chombo includes
#include <RealVect.H>
#include <MayDay.H>

// Our includes
#include <CD_CellInfo.H>
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Namespace for various particle management tools.
 */
namespace ParticleManagement {

/**
 * @brief How a kd merge divides a node into two children.
 * @details The two rules answer different questions. The weight median divides the node so the halves
 * carry as nearly equal a weight as possible, which equalizes the merged weights. The count median
 * divides so the halves hold as nearly equal a particle count as possible, at the cost of leaving the
 * merged weights uneven. Which one is right depends on the node's size: while a node still spans several
 * cells a split apportions leaves BETWEEN cells, which is a question of counts; once it is narrow the
 * weight median is what drives the super-particles toward equal weight. Hybrid is that rule, with the
 * crossover given as a length; Weight and Count are its two limits (crossover at infinity and at zero),
 * so all three are the same knob and every kd method reads it the same way.
 */
enum class KDPartition
{
  Weight, ///< Split at the weight median at every node; halves carry near-equal weight.
  Count,  ///< Split at the count median at every node; halves hold near-equal particle counts.
  Hybrid  ///< Count median while the node is wider than the crossover length, weight median below it.
};

/**
 * @brief Where a kd merge puts the particle a leaf reduces to.
 * @details The centroid keeps the leaf's mean position and discards its spread, so repeated merges
 * contract the sub-cell distribution onto a lattice of leaf centres. Sampling draws one of the leaf's own
 * particles with probability proportional to weight, which is distribution-preserving and inherits a
 * position some particle actually occupied. Random draws uniformly from the leaf's bounding box, which
 * also avoids the lattice but from a distribution that is not quite the leaf's own.
 */
enum class KDPlacement
{
  Centroid, ///< The leaf's weighted centroid.
  Sample,   ///< One of the leaf's own particles, drawn with probability proportional to weight.
  Random    ///< A uniformly random point in the leaf's bounding box.
};

/**
 * @brief How the AMR-scope kd merge resolves the leaves that touch a patch or rank boundary.
 * @details Both tiers see the same tree: the difference is only what happens to a leaf whose members are
 * not all owned by one patch. Carve arbitrates the contested particles between the competing patches in a
 * single non-iterative exchange, so a boundary leaf merges as one group. Nn leaves the boundary skin to a
 * nearest-neighbour pair merge instead, which needs no arbitration but drains the skin in rounds.
 */
enum class KDAmrBoundary
{
  Carve, ///< Arbitrate contested particles between patches (z-buffer carve).
  Nn     ///< Merge the interior with the kd tree and the boundary skin with nearest-neighbour pairs.
};

/**
 * @brief How the AMR-scope nearest-neighbour merge finds a particle's merge candidates.
 * @details The merge protocol is identical for all three; only the per-patch spatial index differs, and
 * with it the reachable merge distance.
 */
enum class NNSearch
{
  Tree,   ///< One whole-patch PointCloudBVH per patch.
  Hash,   ///< One whole-patch PointCloudHashGrid per patch.
  OneCell ///< One PointCloudBVH per occupied cell; merge distance structurally fixed at 1 cell.
};

/**
 * @brief The super-particle merge methods, named by the scope they group particles over.
 * @details The scope is the part that changes what the merge can do; how a method groups within that
 * scope is a separate axis (KDPartition/KDPlacement for the kd methods, NNSearch for the
 * nearest-neighbour one), so it is not spelled into the selector.
 */
enum class ParticleMergeMethod
{
  None,         ///< No merging.
  KdCell,       ///< Cell: makeKDCellMerger, one kd tree per cell.
  KdPatch,      ///< Patch: mergeKDPatch, one kd tree per patch, patch-local (no ghosts).
  KdAmr,        ///< AMR: one kd tree per patch, resolved across patches -- see KDAmrBoundary.
  NnAmr,        ///< AMR: nearest-neighbour pair merge over the hierarchy -- see NNSearch.
  Reinitialize, ///< Cell: makeReinitializeMerger; positions discarded and redrawn.
  External      ///< Cell: caller-supplied per-cell merger.
};

/**
 * @brief Map a merge-method selector string to a ParticleMergeMethod.
 * @details Valid input options are:
 * 'none' -> None
 * 'kd_cell' -> KdCell
 * 'kd_patch' -> KdPatch
 * 'kd_amr' -> KdAmr
 * 'nn_amr' -> NnAmr
 * 'reinitialize' -> Reinitialize
 * 'external' -> External
 * @param[in] a_str String selector. Must take one of the above strings.
 * @return The corresponding method; aborts on an unknown selector.
 */
inline ParticleMergeMethod
mergeMethodFromString(const std::string& a_str) noexcept;

/**
 * @brief Turn an input string into a kd partition rule.
 * @details Valid inputs are 'weight', 'count' and 'hybrid'.
 * @param[in] a_str Input string.
 * @return The corresponding KDPartition.
 */
inline KDPartition
kdPartitionFromString(const std::string& a_str) noexcept;

/**
 * @brief Turn an input string into a kd placement rule.
 * @details Valid inputs are 'centroid', 'sample' and 'random'.
 * @param[in] a_str Input string.
 * @return The corresponding KDPlacement.
 */
inline KDPlacement
kdPlacementFromString(const std::string& a_str) noexcept;

/**
 * @brief Turn an input string into a kd AMR boundary rule.
 * @details Valid inputs are 'carve' and 'nn'.
 * @param[in] a_str Input string.
 * @return The corresponding KDAmrBoundary.
 */
inline KDAmrBoundary
kdAmrBoundaryFromString(const std::string& a_str) noexcept;

/**
 * @brief Turn an input string into a nearest-neighbour search backend.
 * @details Valid inputs are 'tree', 'hash' and 'onecell'.
 * @param[in] a_str Input string.
 * @return The corresponding NNSearch.
 */
inline NNSearch
nnSearchFromString(const std::string& a_str) noexcept;

/**
 * @brief Concept for splitting/merging particles
 * @param[in,out] a_particles Particles to be merged/split
 * @param[in] a_cellInfo Cell info
 * @param[in] a_numTargetParticles Number of target particles
 */
template <class P, class Traits = ParticleTraits<P>>
using ParticleMerger = std::function<
  void(ParticleSoA<P, Traits>& a_particles, const CellInfo& a_cellInfo, const int a_numTargetParticles)>;

/**
 * @brief Declaration of a reconciliation function when splitting particles.
 * @details This is, for example, passed into the KD partitioning structure, whose weight median may split a
 * particle such that the weight of the two nodes differ by at most one physical particle. By default, that method will
 * call the particle copy constructor, but this function permits the user to input a reconciliation function that
 * manipulates other class members in the split particles.
 */
template <class P>
using BinaryParticleReconcile = std::function<void(P& p1, P& p2, const P& p0)>;

/**
 * @brief Internal-only implementation details of the ParticleManagement API.
 * @details These helpers have no external callers; they are the algorithms sitting behind the public
 * factory/wrapper functions (@c make*Merger, @c partitionParticleWeights) and are not part of the
 * public interface. Public callers qualify them with @c detail::.
 */
namespace detail {

/**
 * @brief Build a KD partition of a list of particles and return the leaf particle ranges.
 * @details Partitions the input particles into at most @c a_maxLeaves spatial groups ("leaves"). Each split
 * divides the node along the longest bounding-box axis, at the weight median while the node is no wider
 * than @c a_weightMedianLength and at the count median above that -- the same size-dependent rule the
 * whole-patch build uses, and the same one a fixed 'weight' or 'count' setting expresses as a crossover
 * above or below every node's size. Under the weight median the halves carry as nearly equal a weight as
 * possible (differing by at most one physical particle), which may require dividing the median particle
 * itself, in which case @c a_particleReconcile fixes up the daughters' non-weight members; the count
 * median never divides a particle. This is the partitioning primitive behind super-particle merging:
 * the caller reduces each returned leaf range to a single merged particle.
 *
 * Implementation note: the tree is built in a flat, reusable per-thread scratch arena (two ping-ponged
 * contiguous particle buffers + a flat range list), so there is no per-node heap allocation, no pointer
 * chasing, and the scratch is reused across calls on the same thread.
 *
 * Template parameters:
 *    P            -> Particle type exposing the weight/position data members below
 *    P::*weight   -> Pointer to the @c Real data member holding the particle weight
 *    P::*position -> Pointer to the @c RealVect data member holding the particle position.
 *
 * P must be an AoS-style particle type carrying weight and position as data members -- MergeParticle is
 * the one this library provides, and the one every in-tree caller uses. ParticleSoA payload types
 * (ItoParticle, ...) store those two in the container rather than the payload, so pack the columns into
 * such an intermediary AoS structure first.
 *
 * @param[in]  a_particles         Input particles (copied into internal scratch; not modified).
 * @param[in]  a_maxLeaves         Maximum number of leaf nodes (the per-cell target particle count).
 * @param[in]  a_weightMedianLength Node size, as a length, at or below which the split plane is the weight
 *                                  median rather than the count median. Negative disables the weight
 *                                  median entirely; a huge value selects it everywhere.
 * @param[in]  a_particleReconcile Reconciliation callback applied to the two daughters when the median is split.
 * @param[out] a_leaves            Filled with one {begin,end} pointer pair per leaf.
 * @warning The returned ranges alias an internal thread-local buffer that is overwritten on the next call to
 *          this function on the same thread. Consume them (reduce each leaf) before calling again.
 */
template <class P, Real P::*weight, RealVect P::*position>
inline void
buildKDCellLeaves(const std::vector<P>&                       a_particles,
                  const int                                   a_maxLeaves,
                  const Real                                  a_weightMedianLength,
                  const BinaryParticleReconcile<P>&           a_particleReconcile,
                  std::vector<std::pair<const P*, const P*>>& a_leaves) noexcept;

} // namespace detail

/**
 * @brief Create a per-cell KD-tree super-particle merger as a reusable ParticleMerger.
 * @details Returns a merger functor that partitions the particle set in a single SoA cell into at most
 * @c a_numTargetParticles KD leaves and reduces each leaf to one output particle.
 *
 * The algorithm:
 *   1. Packs each SoA slot into a lightweight intermediate (@c Packed) via @c a_gather and accumulates
 *      total weight @c W.
 *   2. Returns early if @c W < 2 or the target is non-positive.
 *   3. Calls detail::buildKDCellLeaves to partition the packed list into at most @c a_numTargetParticles
 *      leaves, using @c a_reconcile to fix up daughter particles whenever the weight median is split.
 *   4. Clears the SoA container, then reduces each leaf by calling @c a_scatterLeaf once per leaf.
 *
 * @c a_scatterLeaf receives the raw [first, last) pointer range of the leaf (pointers into the internal
 * KD scratch buffer, valid only until the next call) and the current CellInfo, so it can implement either
 * a weighted-centroid reduction or a spatially-reinitialised reduction (e.g. random point in the leaf
 * bounding box for non-cut-cells).
 *
 * Thread-local scratch vectors are reused across cells — no per-cell heap allocation after warmup.
 *
 * @tparam Packed       Lightweight AoS intermediate type (e.g. MergeParticle<Packed>). Must be copyable.
 * @tparam packWeight   Member-data pointer `Real Packed::*` naming the weight. Forwarded to
 *                      detail::buildKDCellLeaves for median splitting.
 * @tparam packPosition Member-data pointer `RealVect Packed::*` naming the position. Forwarded to
 *                      detail::buildKDCellLeaves for bounding-box splitting.
 * @tparam P            SoA payload type.
 * @tparam Traits       Column descriptor for @c P (defaults to ParticleTraits<P>).
 * @param[in] a_weightMedianCellWidths Node size, in cell widths, at or below which the split plane is the
 *                                    weight median rather than the count median. See KDPartition.
 * @param[in] a_gather      Per-slot gather: `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed`.
 * @param[in] a_reconcile   Daughter reconciliation: `(Packed& p1, Packed& p2, const Packed& p0) -> void`.
 * @param[in] a_scatterLeaf Leaf reducer: `(ParticleSoA<P,Traits>&, const Packed* first,
 *                          const Packed* last, const CellInfo&) -> void`.  Appends exactly one merged
 *                          particle per call.
 * @return A ParticleMerger<P, Traits> ready to be stored in e.g. ItoSolver::m_particleCellMerger.
 */
template <class Packed,
          Real Packed::*packWeight,
          RealVect Packed::*packPosition,
          class P,
          class Traits = ParticleTraits<P>>
inline ParticleMerger<P, Traits>
makeKDCellMerger(
  const Real                                                                                  a_weightMedianCellWidths,
  std::function<Packed(const ParticleSoA<P, Traits>&, std::size_t)>                           a_gather,
  BinaryParticleReconcile<Packed>                                                             a_reconcile,
  std::function<void(ParticleSoA<P, Traits>&, const Packed*, const Packed*, const CellInfo&)> a_scatterLeaf) noexcept;

/**
 * @brief Create a reinitialize super-particle merger as a reusable ParticleMerger.
 * @details Returns a merger functor that collapses all physical particles in a cell into at most
 * @c a_numTargetParticles new computational particles, placing them at random positions within the cell.
 *
 * The algorithm:
 *   1. Calls @c a_aggregate once on the input SoA to obtain the total physical-particle count and an
 *      arbitrary caller-defined context (e.g. the weight-averaged energy for ItoParticle).
 *   2. Returns early if the physical count is non-positive or the target is non-positive.
 *   3. Calls partitionParticleWeights to divide the physical count into at most @c a_numTargetParticles
 *      integer weights.
 *   4. Clears the SoA; for each weight draws a random position in the cell via
 *      Random::randomPosition and calls @c a_emit once.
 *
 * The cell position is computed as @c a_probLo + dx * (gridIndex + 0.5), matching the formula used
 * in ItoSolver::reinitializeParticles.
 *
 * @tparam Context  Arbitrary aggregated payload type (e.g. @c Real for average energy).  Must be
 *                  copyable and default-constructible.
 * @tparam P        SoA payload type.
 * @tparam Traits   Column descriptor for @c P (defaults to ParticleTraits<P>).
 * @param[in] a_aggregate Aggregate functor: `(const ParticleSoA<P,Traits>&) -> pair<long long, Context>`.
 *                        Returns total physical particle count and a context consumed by @c a_emit.
 * @param[in] a_emit      Per-particle emit: `(ParticleSoA<P,Traits>&, const RealVect& pos,
 *                        long long weight, const Context&) -> void`.  Appends exactly one particle.
 * @param[in] a_probLo    Lazy supplier of the domain lower-left corner (from AmrMesh::getProbLo).
 *                        Called on every merge invocation; capture the solver's AmrMesh reference
 *                        so the call is deferred until after AMR setup completes.
 * @return A ParticleMerger<P, Traits> ready to be stored in e.g. ItoSolver::m_particleCellMerger.
 */
template <class Context, class P, class Traits = ParticleTraits<P>>
inline ParticleMerger<P, Traits>
makeReinitializeMerger(std::function<std::pair<long long, Context>(const ParticleSoA<P, Traits>&)> a_aggregate,
                       std::function<void(ParticleSoA<P, Traits>&, const RealVect&, long long, const Context&)> a_emit,
                       std::function<RealVect()> a_probLo) noexcept;

/**
 * @brief Remove physical particles from an SoA container (operates on the container-owned weight column).
 * @param[in,out] a_particles           Input SoA particle container.
 * @param[in]    a_numPhysPartToRemove Number of physical particles to remove.
 * @note May fail if particle weights are not (close to) integers.
 * @tparam P      Payload type.
 * @tparam Traits Column descriptor for the payload.
 * @tparam T      Integral type for the physical-particle count.
 */
template <typename P, typename Traits, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
static inline void
removePhysicalParticles(ParticleSoA<P, Traits>& a_particles, const T a_numPhysPartToRemove) noexcept;

/**
 * @brief Remove particles from an SoA container if their weight is below the threshold (swap-and-pop).
 * @param[in,out] a_particles    Input SoA particle container.
 * @param[in]    a_weightThresh Weight threshold for removal.
 * @tparam P      Payload type.
 * @tparam Traits Column descriptor for the payload.
 */
template <typename P, typename Traits>
static inline void
deleteParticles(ParticleSoA<P, Traits>& a_particles, const Real a_weightThresh) noexcept;

/**
 * @brief Partition particle weights among a number of computational particles
 * @param[in] a_numPhysicalParticles Number of physical particles.
 * @param[in] a_maxCompParticles     Maximum number of computational particles
 * @return Return value
 */
template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
static inline std::vector<T>
partitionParticleWeights(const T a_numPhysicalParticles, const T a_maxCompParticles) noexcept;

/**
 * @brief Internal-only implementation details of the ParticleManagement API (continued).
 * @details Reopens the @c detail namespace introduced above to house @c partitionParticles alongside the
 * other internal helpers; see the first @c detail block for the full rationale. Not part of the public
 * interface -- public callers qualify with @c detail::.
 */
namespace detail {

/**
 * @brief Partition particles so that all MPI rank draw a cumulative number of particles equal to a_numParticles
 * @param[in] a_numParticles Total number of particles to be drawn.
 * @return Returns the number of particles drawn by each rank.
 */
template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
static inline T
partitionParticles(const T a_numParticles);

} // namespace detail

/**
 * @brief Draw random SoA particles into a free buffer, distributed according to a predefined distribution.
 * @details Fills a free-standing ParticleSoA leaf (a buffer) rather than an AMR container. MPI-aware --
 * it partitions a_numParticles across ranks (total = a_numParticles) and each rank draws its own share
 * (unit weight, default payload).
 * Downstream code adds the buffer to a ParticleContainer (e.g. addParticlesDestructive) to route the
 * particles to their owning patch/level/rank.
 * @param[out] a_particles    SoA buffer to fill (cleared first).
 * @param[in]  a_numParticles Total number of particles across all ranks.
 * @param[in]  a_distribution Distribution functor returning a random position.
 * @tparam P      Payload type.
 * @tparam Traits Column descriptor for the payload.
 * @tparam T      Integral count type.
 */
template <typename P, typename Traits, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
static inline void
drawRandomParticles(ParticleSoA<P, Traits>&          a_particles,
                    const T                          a_numParticles,
                    const std::function<RealVect()>& a_distribution);

/**
 * @brief Draw Gaussian-distributed particles into an SoA container (unit weight; positions only).
 * @details Each rank draws its own a_numParticles share with weight 1; routing to owners happens when
 * the buffer is added to a ParticleContainer.
 * @param[out] a_particles    Output SoA buffer (cleared on entry).
 * @param[in]  a_numParticles Number of particles.
 * @param[in]  a_center       Gaussian blob center.
 * @param[in]  a_radius       Gaussian blob radius.
 * @tparam P      Payload type.
 * @tparam Traits Column descriptor for the payload.
 * @tparam T      Integral count type.
 */
template <typename P, typename Traits, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
static inline void
drawGaussianParticles(ParticleSoA<P, Traits>& a_particles,
                      const T                 a_numParticles,
                      const RealVect&         a_center,
                      const Real              a_radius) noexcept;

/**
 * @brief Draw box-distributed particles into an SoA container (unit weight; positions only).
 * @details Each rank draws its own a_numParticles share with weight 1; routing to owners happens when
 * the buffer is added to a ParticleContainer.
 * @param[out] a_particles    Output SoA buffer (cleared on entry).
 * @param[in]  a_numParticles Number of particles.
 * @param[in]  a_loCorner     Low corner of box.
 * @param[in]  a_hiCorner     High corner of box.
 * @tparam P      Payload type.
 * @tparam Traits Column descriptor for the payload.
 * @tparam T      Integral count type.
 */
template <typename P, typename Traits, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
static inline void
drawBoxParticles(ParticleSoA<P, Traits>& a_particles,
                 const T                 a_numParticles,
                 const RealVect&         a_loCorner,
                 const RealVect&         a_hiCorner) noexcept;

/**
 * @brief Draw sphere-distributed particles into an SoA container (unit weight; positions only).
 * @details Each rank draws its own a_numParticles share with weight 1; routing to owners happens when
 * the buffer is added to a ParticleContainer.
 * @param[out] a_particles    Output SoA buffer (cleared on entry).
 * @param[in]  a_numParticles Number of particles.
 * @param[in]  a_center       Sphere center.
 * @param[in]  a_radius       Sphere radius.
 * @tparam P      Payload type.
 * @tparam Traits Column descriptor for the payload.
 * @tparam T      Integral count type.
 */
template <typename P, typename Traits, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
static inline void
drawSphereParticles(ParticleSoA<P, Traits>& a_particles,
                    const T                 a_numParticles,
                    const RealVect&         a_center,
                    const Real              a_radius) noexcept;

/**
 * @brief Minimal, payload-agnostic description of one particle as input to a distributed merge.
 * @details Shared by the kd-tree and nearest-neighbor mergers, which previously kept one of these
 * each, differing by a single field. Carries what any of them needs: position and weight, a globally
 * unique id, the owning rank, plus an opaque caller-defined payload the merge never inspects and only
 * hands to the combine callback.
 *
 * @tparam Packed Caller-defined payload type (e.g., a single Real energy for ItoParticle).
 */
template <typename Packed>
struct MergeParticle
{
  /**
   * @brief Particle position.
   */
  RealVect position = RealVect::Zero;

  /**
   * @brief Particle weight. Not required to be an integer.
   */
  Real weight = 0.0;

  /**
   * @brief Globally unique particle id, unique within one merge round across every rank.
   * @details Defaults to -1, matching ParticleSoA::s_invalidID (not spelled as that constant, which is
   * a member of a ParticleSoA instantiation this payload-agnostic struct has no business naming). The
   * default matters because the cell-based mergers have no id fields to fill: without it the KD leaf
   * builder would copy indeterminate values between daughters.
   */
  ParticleID globalID = -1;

  /**
   * @brief The rank owning this particle. Always read from the particle's own data, never inferred
   * from MPI transport metadata.
   */
  RankID ownerRank = -1;

  /**
   * @brief True iff this is a ghost copy in the patch currently holding it, i.e. not physically
   * stored there. NOT the same question as ownerRank: a particle owned by this rank via one patch is
   * still a ghost in a neighbouring patch this same rank owns. Used by the kd mergers; left false by
   * the nearest-neighbor mergers, which track the distinction separately.
   */
  bool isGhost = false;

  /**
   * @brief The AMR level this particle lives on. Used by the nearest-neighbor mergers, whose cell
   * keys are level-qualified; left at 0 by the kd mergers, which work in one patch's frame at a time.
   */
  int level = 0;

  /**
   * @brief Opaque, caller-defined payload. Never inspected by the merge logic -- only carried through
   * to the combine callback.
   */
  Packed payload;
};

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#include <CD_ParticleManagementImplem.H>

#endif
