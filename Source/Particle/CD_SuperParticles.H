/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SuperParticles.H
  @brief  Namespace containing various particle merging utilities.
  @author Robert Marskar
*/

#ifndef CD_SuperParticles_H
#define CD_SuperParticles_H

// Std includes
#include <memory>
#include <vector>

// Chombo includes
#include <RealVect.H>
#include <List.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Namespace for various particle merging tools. 
*/
namespace SuperParticles {

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
    @brief Partitioning method for KDNode<P> which equilibriates the weight across the two halves when splitting a node
    @note This only adjusts the weight of the median particle. The particle class P should have a copy constructor for the other
    fields; note that these fields are simply copied over when the particle is split and are not adjusted. 
  */
  template <class P, Real& (P::*weight)(), const RealVect& (P::*position)() const>
  typename KDNode<P>::Partitioner PartitionEqualWeight;

  /*!
    @brief Stop function which returns true if weight in node < 2
    @note The particle type is P and P::*weight is a pointer to a member function providing const-reference access to the particle weight.
  */
  template <class P, const Real& (P::*weight)() const>
  typename KDNode<P>::StopFunction InsufficientWeightForSplitting;
} // namespace SuperParticles

#include <CD_NamespaceFooter.H>

#include <CD_SuperParticlesImplem.H>

#endif
