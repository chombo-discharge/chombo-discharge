/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParticleManagement.H
  @brief  Namespace containing various particle management utilities.
  @author Robert Marskar
*/

#ifndef CD_ParticleManagement_H
#define CD_ParticleManagement_H

// Std includes
#include <memory>
#include <vector>

// Chombo includes
#include <RealVect.H>
#include <List.H>

// Our includes
#include <CD_KDNode.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Namespace for various particle management tools. 
*/
namespace ParticleManagement {

  /*!
    @brief Concept for splitting a particle into more particles
  */
  template<class P, size_t N>
  using ParticleSplitter = std::function<std::array<P, N>(const P& p)>;

  /*!
    @brief Partition and split a KD-node with particles so that the weights are approximately the same. 
    @details This will partition particles on two sides of a plane parallel to one of the coordinate axis. The middle
    particle closest to this plane may be split, and the weights of the resulting particles are assigned such that the
    total weight of the two halves differ by at most one. The template parameters are 

       P            -> Particle type
       P::*weight   -> Function pointer to particle weight
       P::*position -> Function pointer to particle position.

    A possible call signature is e.g. recursivePartitionAndSplitEqualWeightKD<P, &P::weight, &P::position>.
  */
  template <class P,  Real& (P::*weight)(), const RealVect& (P::*position)() const>
  static inline void
  partitionAndSplitEqualWeightKD(KDNode<P>& a_node) noexcept;

  /*!
    @brief Recursively build a KD-tree following the "equal weight" principle when partitioning nodes.
    @details If the number of leaves is a factor of two, leaves exist on the same level and the weight in each node will differ by
    at most one physical particle. The template parameters are 

       P            -> Particle type
       P::*weight   -> Function pointer to particle weight
       P::*position -> Function pointer to particle position.

    A possible call signature is e.g. recursivePartitionAndSplitEqualWeightKD<P, &P::weight, &P::position>.
    @param[inout] a_inputParticles Input particles. These are destroyed on output.
    @param[in] a_maxLeaves Maximum number of leaves in the tree.
    @return Returns leaf nodes with particles. 
  */
  template <class P,  Real& (P::*weight)(), const RealVect& (P::*position)() const>  
  static inline std::vector<std::shared_ptr<KDNode<P>>>
  recursivePartitionAndSplitEqualWeightKD(typename KDNode<P>::ParticleList& a_inputParticles, const int a_maxLeaves) noexcept;

  /*!
    @brief Remove physical particles from the input particles.
    @param[inout] a_particles           Input list of particles. Must have a weight function. 
    @param[in]    a_numPhysPartToRemove Number of physical particles to remove
    @note May fail if particle weights are not (close to) integers. 
  */
  template <typename P, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline void
  removePhysicalParticles(List<P>& a_particles, const T a_numPhysPartToRemove) noexcept;

  /*!
    @brief Remove particles if their weight is below the input weight
    @param[inout] a_particles    Input list of particles. Must have a weight function. 
    @param[in]    a_weightThresh Weight threshold for removal
  */
  template <typename P>
  static inline void
  deleteParticles(List<P>& a_particles, const Real a_weightThresh) noexcept;

  /*!
    @brief Partition particle weights among a number of computational particles
    @param[in] a_numPhysicalParticles Number of physical particles. 
    @param[in] a_maxCompParticles     Maximum number of computational particles
  */
  template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline std::vector<T>
  partitionParticleWeights(const T a_numPhysicalParticles, const T a_maxCompParticles) noexcept;

  /*!
    @brief Partition particles so that all MPI rank draw a cumulative number of particles equal to a_numParticles
    @param[in] a_numParticles Total number of particles to be drawn. 
    @return Returns the number of particles drawn by each rank. 
  */
  template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline T
  partitionParticles(const T a_numParticles);

  /*!
    @brief Draw a specified number of random particles distributed according to a predefined distribution. 
    @details The implementation will partition a_numParticles among the ranks so the total number of particles drawn is a_numParticles.
    P is the particle type (must derived from BinItem) and T is an integer type. 
    @param[out] a_particles    Output particles (different for each rank)
    @param[in]  a_numParticles Number of particles
    @param[in]  a_distribution Particle distribution. Must return a random position when calling the function. 
  */
  template <typename P, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline void
  drawRandomParticles(List<P>& a_particles, const T a_numParticles, const std::function<RealVect()>& a_distribution);

  /*!
    @brief Draw particles in a sphere
    @details The implementation will partition a_numParticles among the ranks so the total number of particles drawn is a_numParticles.
    P is the particle type (must have a RealVect& position() function) and T is an integer type. 
    @param[out] a_particles    Output particles (different for each rank)
    @param[in]  a_numParticles Number of particles
    @param[in]  a_center       Sphere center
    @param[in]  a_radius       Sphere radius
  */
  template <typename P, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline void
  drawSphereParticles(List<P>&       a_particles,
                      const T        a_numParticles,
                      const RealVect a_center,
                      const Real     a_radius) noexcept;

  /*!
    @brief Draw particles in a box
    @details The implementation will partition a_numParticles among the ranks so the total number of particles drawn is a_numParticles.
    P is the particle type (must have a RealVect& position() function) and T is an integer type. 
    @param[out] a_particles    Output particles (different for each rank)
    @param[in]  a_numParticles Number of particles
    @param[in]  a_loCorner     Low corner of box
    @param[in]  a_hiCorner     High corner of box
  */
  template <typename P, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline void
  drawBoxParticles(List<P>&       a_particles,
                   const T        a_numParticles,
                   const RealVect a_loCorner,
                   const RealVect a_hiCorner) noexcept;

  /*!
    @brief Draw particles from a Gaussian distribution
    @details The implementation will partition a_numParticles among the ranks so the total number of particles drawn is a_numParticles.
    P is the particle type (must have a RealVect& position() function) and T is an integer type. 
    @param[out] a_particles    Output particles (different for each rank)
    @param[in]  a_numParticles Number of particles
    @param[in]  a_center       Gaussian blob Center
    @param[in]  a_radius       Gaussian blob radius
  */
  template <typename P, typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
  static inline void
  drawGaussianParticles(List<P>&       a_particles,
                        const T        a_numParticles,
                        const RealVect a_center,
                        const Real     a_radius) noexcept;

} // namespace ParticleManagement

#include <CD_NamespaceFooter.H>

#include <CD_ParticleManagementImplem.H>

#endif
