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
#include <memory>
#include <vector>

// Chombo includes
#include <RealVect.H>
#include <List.H>

// Our includes
#include <CD_KDNode.H>
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
    @brief Partition and split a KD-node with particles so that the weights are approximately the same. 

    @param[in,out] a_node Node to split. This MUST be a leaf node and it MUST be possible to split it. 
    @param[in] a_particleReconcile Optional reconciliation function when splitting a particle into two new particles.
    @details This will partition particles on two sides of a plane parallel to one of the coordinate axis. The middle
    particle closest to this plane may be split, and the weights of the resulting particles are assigned such that the
    total weight of the two halves differ by at most one. The template parameters are 

       P            -> Particle type
       P::*weight   -> Function pointer to particle weight
       P::*position -> Function pointer to particle position.

    A possible call signature is e.g. recursivePartitionAndSplitEqualWeightKD<P, &P::weight, &P::position>.

    The user can input a particle reconciliation function that manipulates the particle properties of the split particles. By default,
    the split particles will use the copy constructor and thus inherit class members from p0, with the exception of the particle weights.
    The reconcile function lets the user manipulate other particle properties, e.g. ones that are not properly captured by the particle copy
    constructor, or that need some other form of reconciliation. 
  */
  template <class P, Real& (P::*weight)(), const RealVect& (P::*position)() const>
  static inline void
  partitionAndSplitEqualWeightKD(
    KDNode<P>&                       a_node,
    const BinaryParticleReconcile<P> a_particleReconcile = [](P& p1, P& p2, const P& p0) -> void {
    }) noexcept;

  /**
    @brief Recursively build a KD-tree following the "equal weight" principle when partitioning nodes.
    @param[in,out] a_inputParticles Input particles. These are destroyed on output.
    @param[in] a_maxLeaves Maximum number of leaves in the tree.
    @return Returns leaf nodes containing the particles
    @details If the number of leaves is a factor of two, leaves exist on the same level and the weight in each node will differ by
    @param[in] a_particleReconcile Description
    at most one physical particle. The template parameters are 

       P            -> Particle type
       P::*weight   -> Function pointer to particle weight
       P::*position -> Function pointer to particle position.

    A possible call signature is e.g. recursivePartitionAndSplitEqualWeightKD<P, &P::weight, &P::position>.

    The user can input a particle reconciliation function that manipulates the particle properties of the split particles. By default,
    the split particles will use the copy constructor and thus inherit class members from p0, with the exception of the particle weights.
    The reconcile function lets the user manipulate other particle properties, e.g. ones that are not properly captured by the particle copy
    constructor, or that need some other form of reconciliation. 

  */
  template <class P, Real& (P::*weight)(), const RealVect& (P::*position)() const>
  static inline std::vector<std::shared_ptr<KDNode<P>>>
  recursivePartitionAndSplitEqualWeightKD(
    typename KDNode<P>::ParticleList& a_inputParticles,
    const int                         a_maxLeaves,
    const BinaryParticleReconcile<P>  a_particleReconcile = [](P& p1, P& p2, const P& p0) -> void {
    }) noexcept;

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
