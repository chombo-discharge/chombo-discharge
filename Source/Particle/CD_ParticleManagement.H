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
#include <CD_NamespaceHeader.H>

/*!
  @brief Namespace for various particle management tools. 
*/
namespace ParticleManagement {

  /*!
    @brief Node in a particle-merging KD-tree. 
    @details This node type is used for partitioning particles into spatial subsets and merging/splitting them. 
    @details The template argument P is the particle type - the users decides how to merge/split/partition. 
  */
  template <class P>
  class KDNode
  {
  public:
    /*!
      @brief List of particles. This is aliased because the KD-tree construction may require both random access
      and something sortable. 
    */
    using ParticleList = std::vector<P>;

    /*!
      @brief Node partitioner for KD-tree construction.
    */
    using Partitioner = std::function<void(KDNode&)>;

    /*!
      @brief Stop function for KD-tree construction. 
    */
    using StopFunction = std::function<bool(const ParticleList&)>;

    /*!
      @brief Disallowed constructor
    */
    KDNode(const KDNode&) = delete;

    /*!
      @brief Disallowed constructor
    */
    KDNode(const KDNode&&) = delete;

    /*!
      @brief Disallowed assignment
    */
    KDNode&
    operator=(const KDNode&) = delete;

    /*!
      @brief Disallowed assignment
    */
    KDNode&
    operator=(const KDNode&&) = delete;

    /*!
      @brief Default constructor
    */
    KDNode();

    /*!
      @brief Valid constructor. Takes list of particles.
      @param[inout] a_particles List of particles to be partitioned later. 
      @note The input particle list is transferred to m_particles. 
    */
    KDNode(ParticleList& a_particles);

    /*!
      @brief Destructor. Does nothing.
    */
    virtual ~KDNode();

    /*!
      @brief Get the node weight
    */
    inline const Real&
    weight() const noexcept;

    /*!
      @brief Get the node weight
    */
    inline Real&
    weight() noexcept;

    /*!
      @brief Get particles in this node. 
    */
    inline const ParticleList&
    getParticles() const noexcept;

    /*!
      @brief Get particles in this node. 
    */
    inline ParticleList&
    getParticles() noexcept;

    /*!
      @brief Is leaf node or not.
    */
    inline bool
    isLeafNode() const noexcept;

    /*!
      @brief Is leaf node or not.
    */
    inline bool
    isInteriorNode() const noexcept;

    /*!
      @brief Gather particles further down in the subtree and return all particles (in the leaf nodes)
    */
    inline ParticleList
    gatherParticles() const noexcept;

    /*!
      @brief Move the particles list further down in the subtree into this vector.
      @note This clears the particles from the leaf nodes. Use gatherParticles if you do not
      want to change the tree contents. 
    */
    inline ParticleList
    moveParticles() noexcept;

    /*!
      @brief Get the left node.
    */
    inline std::shared_ptr<KDNode<P>>&
    getLeft() noexcept;

    /*!
      @brief Get the right node.
    */
    inline std::shared_ptr<KDNode<P>>&
    getRight() noexcept;

  protected:
    /*!
      @brief Left KD-node
    */
    std::shared_ptr<KDNode> m_left;

    /*!
      @brief Right KD-node
    */
    std::shared_ptr<KDNode> m_right;

    /*!
      @brief Node weight
    */
    Real m_weight;

    /*!
      @brief List of particles
    */
    ParticleList m_particles;

    /*!
      @brief Gather particles further down in the subtree and return all particles (in the leaf nodes)
      @param[inout] a_particles List of particles
    */
    inline void
    gatherParticles(ParticleList& a_particles) const noexcept;

    /*!
      @brief Move particles further down in the subtree and return all particles (in the leaf nodes)
      @param[inout] a_particles List of particles
      @note This clears the particles from the leaf nodes. Use gatherParticles if you do not
      want to change the tree contents. 
    */
    inline void
    moveParticles(ParticleList& a_particles) noexcept;
  };

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
