/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleManagement.H
  @brief  Namespace containing various particle management utilities.
  @author Robert Marskar
*/

#ifndef CD_PARTICLEMANAGEMENT_H
#define CD_PARTICLEMANAGEMENT_H

// Std includes
#include <functional>
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
  @brief Namespace for various particle management tools. 
*/
namespace ParticleManagement {

  /**
    @brief Concept for splitting/merging particles
    @param[in,out] a_particles Particles to be merged/split
    @param[in] a_cellInfo Cell info
    @param[in] a_numTargetParticles Number of target particles
  */
  template <class P, class Traits = ParticleTraits<P>>
  using ParticleMerger = std::function<
    void(ParticleSoA<P, Traits>& a_particles, const CellInfo& a_cellInfo, const int a_numTargetParticles)>;

  /**
    @brief Declaration of a reconciliation function when splitting particles.
    @details This is, for example, passed into the KD-tree equal-weight partitioning structure, which may split a particle
    such that the weight of the two nodes differ by at most one physical particle. By default, that method will call the
    particle copy constructor, but this function permits the user to input a reconciliation function that manipulates other
    class members in the split particles. 
  */
  template <class P>
  using BinaryParticleReconcile = std::function<void(P& p1, P& p2, const P& p0)>;

  /**
    @brief Build an equal-weight KD partition of a list of particles and return the leaf particle ranges.
    @details Partitions the input particles into at most @c a_maxLeaves spatial groups ("leaves") whose total
    weights are as equal as possible (they differ by at most one physical particle between sibling nodes). At
    each split the particles are partitioned along the longest bounding-box axis at the weight-median; the
    median particle may itself be split so the two halves balance, in which case @c a_particleReconcile fixes
    up the daughters' non-weight members. This is the partitioning primitive behind super-particle merging:
    the caller reduces each returned leaf range to a single merged particle.

    Implementation note: the tree is built in a flat, reusable per-thread scratch arena (two ping-ponged
    contiguous particle buffers + a flat range list), so there is no per-node heap allocation, no pointer
    chasing, and the scratch is reused across calls on the same thread.

    Template parameters:
       P            -> Particle type exposing the weight()/position() member functions below
       P::*weight   -> Member function returning a mutable reference to the particle weight
       P::*position -> Member function returning a const reference to the particle position.

    P must be an AoS-style particle type with weight()/position() member functions (e.g. NonCommParticle).
    SoA payload types (ItoParticle, ...) have no such methods; pack their columns into such an intermediary
    first (as ItoSolver::makeSuperparticlesEqualWeightKD does).

    @param[in]  a_particles         Input particles (copied into internal scratch; not modified).
    @param[in]  a_maxLeaves         Maximum number of leaf nodes (the per-cell target particle count).
    @param[in]  a_particleReconcile Reconciliation callback applied to the two daughters when the median is split.
    @param[out] a_leaves            Filled with one {begin,end} pointer pair per leaf.
    @warning The returned ranges alias an internal thread-local buffer that is overwritten on the next call to
             this function on the same thread. Consume them (reduce each leaf) before calling again.
  */
  template <class P, Real& (P::*weight)(), const RealVect& (P::*position)() const>
  inline void
  buildEqualWeightKDLeaves(const std::vector<P>&                       a_particles,
                           const int                                   a_maxLeaves,
                           const BinaryParticleReconcile<P>&           a_particleReconcile,
                           std::vector<std::pair<const P*, const P*>>& a_leaves) noexcept;

  /**
    @brief Merge an already spatially-ordered particle list down to a target count by repeatedly combining
    the nearest adjacent pair.
    @details The particles are assumed to be in a spatially coherent order (e.g. sorted along a space-filling
    curve by the caller), so adjacency in the list approximates spatial proximity. The closest adjacent pair
    (by @c position()) is combined repeatedly until @c a_target particles remain. Combining the second
    particle into the first is delegated to @c a_combine (the caller decides the payload semantics, e.g. a
    weight-weighted centroid and energy). This does NOT split, so it cannot increase the count -- it is a
    no-op when the list already has @c a_target or fewer particles (the caller handles the split direction).
    Operates in place on per-thread scratch (allocation-free after warmup); on return @c a_particles holds the
    surviving merged particles.
    @param[in,out] a_particles Spatially-ordered particle list; overwritten with the merged result.
    @param[in]     a_target    Desired number of particles.
    @param[in]     a_combine   Callback that combines the second particle into the first.
  */
  template <class P, const RealVect& (P::*position)() const>
  inline void
  mergeAdjacentNearest(std::vector<P>&                          a_particles,
                       const std::size_t                        a_target,
                       const std::function<void(P&, const P&)>& a_combine) noexcept;

  /**
    @brief Create a space-filling-curve nearest-neighbor super-particle merger as a reusable ParticleMerger.
    @details Returns a merger functor that reduces (or grows) the particle set in a single SoA cell to the
    target count without forcing equal weights.

    When the cell contains more particles than the target:
      1. The intermediate list is sorted along a 21-bit Hilbert curve computed from particle positions.
      2. `mergeAdjacentNearest` merges the closest adjacent pair repeatedly until the target count is
         reached, delegating payload arithmetic to @c a_combine.

    When the cell contains fewer particles than the target: the heaviest particle (by @c packWeight) is
    repeatedly split into two co-located halves (weight divided by two, all other fields copied) until the
    target is reached or no particle has weight >= 2.

    The three user-supplied lambdas encode the particle-type-specific semantics:
      - @c a_gather packs one SoA slot into the lightweight intermediate type @c Packed.
      - @c a_combine merges the second @c Packed particle into the first (weight-weighted centroid is the
        typical choice, but the caller decides).
      - @c a_scatter appends one surviving @c Packed particle back to the (already-cleared) SoA.

    Implementation note: thread-local scratch vectors are reused across cells so there is no per-cell heap
    allocation after the first call on a given thread.

    @tparam Packed      Lightweight AoS intermediate type (e.g. NonCommParticle<M,N>). Must be copyable.
    @tparam packWeight  Member-function pointer `Real& (Packed::*)()` giving mutable access to the weight.
                        Used only for the split path.
    @tparam packPosition Member-function pointer `const RealVect& (Packed::*)() const` giving the position.
                        Used for Hilbert-key computation and forwarded to mergeAdjacentNearest.
    @tparam P           SoA payload type.
    @tparam Traits      Column descriptor for @c P (defaults to ParticleTraits<P>).
    @param[in] a_gather   Per-slot gather: `(const ParticleSoA<P,Traits>&, std::size_t) -> Packed`.
    @param[in] a_combine  In-place merge: `(Packed& first, const Packed& second) -> void`.
    @param[in] a_scatter  Append one merged particle: `(ParticleSoA<P,Traits>&, const Packed&) -> void`.
    @return A ParticleMerger<P, Traits> ready to be stored in e.g. ItoSolver::m_particleMerger.
  */
  template <class Packed,
            Real& (Packed::*packWeight)(),
            const RealVect& (Packed::*packPosition)() const,
            class P,
            class Traits = ParticleTraits<P>>
  inline ParticleMerger<P, Traits>
  makeSfcNearestNeighborMerger(std::function<Packed(const ParticleSoA<P, Traits>&, std::size_t)> a_gather,
                               std::function<void(Packed&, const Packed&)>                       a_combine,
                               std::function<void(ParticleSoA<P, Traits>&, const Packed&)>       a_scatter) noexcept;

  /**
    @brief Remove physical particles from an SoA container (operates on the container-owned weight column).
    @param[in,out] a_particles           Input SoA particle container.
    @param[in]    a_numPhysPartToRemove Number of physical particles to remove.
    @note May fail if particle weights are not (close to) integers.
    @tparam P      Payload type.
    @tparam Traits Column descriptor for the payload.
    @tparam T      Integral type for the physical-particle count.
  */
  template <typename P, typename Traits, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline void
  removePhysicalParticles(ParticleSoA<P, Traits>& a_particles, const T a_numPhysPartToRemove) noexcept;

  /**
    @brief Remove particles from an SoA container if their weight is below the threshold (swap-and-pop).
    @param[in,out] a_particles    Input SoA particle container.
    @param[in]    a_weightThresh Weight threshold for removal.
    @tparam P      Payload type.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename P, typename Traits>
  static inline void
  deleteParticles(ParticleSoA<P, Traits>& a_particles, const Real a_weightThresh) noexcept;

  /**
    @brief Partition particle weights among a number of computational particles
    @param[in] a_numPhysicalParticles Number of physical particles. 
    @param[in] a_maxCompParticles     Maximum number of computational particles
    @return Return value
  */
  template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline std::vector<T>
  partitionParticleWeights(const T a_numPhysicalParticles, const T a_maxCompParticles) noexcept;

  /**
    @brief Partition particles so that all MPI rank draw a cumulative number of particles equal to a_numParticles
    @param[in] a_numParticles Total number of particles to be drawn. 
    @return Returns the number of particles drawn by each rank. 
  */
  template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline T
  partitionParticles(const T a_numParticles);

  /**
    @brief Draw random SoA particles into a free buffer, distributed according to a predefined distribution.
    @details Fills a free-standing ParticleSoA leaf (a buffer) rather than an AMR container. MPI-aware --
    it partitions a_numParticles across ranks (total = a_numParticles) and each rank draws its own share
    (unit weight, default payload).
    Downstream code adds the buffer to a ParticleContainer (e.g. addParticlesDestructive) to route the
    particles to their owning patch/level/rank.
    @param[out] a_particles    SoA buffer to fill (cleared first).
    @param[in]  a_numParticles Total number of particles across all ranks.
    @param[in]  a_distribution Distribution functor returning a random position.
    @tparam P      Payload type.
    @tparam Traits Column descriptor for the payload.
    @tparam T      Integral count type.
  */
  template <typename P, typename Traits, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline void
  drawRandomParticles(ParticleSoA<P, Traits>&          a_particles,
                      const T                          a_numParticles,
                      const std::function<RealVect()>& a_distribution);

  /**
    @brief Draw Gaussian-distributed particles into an SoA container (unit weight; positions only).
    @details Each rank draws its own a_numParticles share with weight 1; routing to owners happens when
    the buffer is added to a ParticleContainer.
    @param[out] a_particles    Output SoA buffer (cleared on entry).
    @param[in]  a_numParticles Number of particles.
    @param[in]  a_center       Gaussian blob center.
    @param[in]  a_radius       Gaussian blob radius.
    @tparam P      Payload type.
    @tparam Traits Column descriptor for the payload.
    @tparam T      Integral count type.
  */
  template <typename P, typename Traits, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline void
  drawGaussianParticles(ParticleSoA<P, Traits>& a_particles,
                        const T                 a_numParticles,
                        const RealVect&         a_center,
                        const Real              a_radius) noexcept;

  /**
    @brief Draw box-distributed particles into an SoA container (unit weight; positions only).
    @details Each rank draws its own a_numParticles share with weight 1; routing to owners happens when
    the buffer is added to a ParticleContainer.
    @param[out] a_particles    Output SoA buffer (cleared on entry).
    @param[in]  a_numParticles Number of particles.
    @param[in]  a_loCorner     Low corner of box.
    @param[in]  a_hiCorner     High corner of box.
    @tparam P      Payload type.
    @tparam Traits Column descriptor for the payload.
    @tparam T      Integral count type.
  */
  template <typename P, typename Traits, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline void
  drawBoxParticles(ParticleSoA<P, Traits>& a_particles,
                   const T                 a_numParticles,
                   const RealVect&         a_loCorner,
                   const RealVect&         a_hiCorner) noexcept;

  /**
    @brief Draw sphere-distributed particles into an SoA container (unit weight; positions only).
    @details Each rank draws its own a_numParticles share with weight 1; routing to owners happens when
    the buffer is added to a ParticleContainer.
    @param[out] a_particles    Output SoA buffer (cleared on entry).
    @param[in]  a_numParticles Number of particles.
    @param[in]  a_center       Sphere center.
    @param[in]  a_radius       Sphere radius.
    @tparam P      Payload type.
    @tparam Traits Column descriptor for the payload.
    @tparam T      Integral count type.
  */
  template <typename P, typename Traits, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline void
  drawSphereParticles(ParticleSoA<P, Traits>& a_particles,
                      const T                 a_numParticles,
                      const RealVect&         a_center,
                      const Real              a_radius) noexcept;

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#include <CD_ParticleManagementImplem.H>

#endif
