/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_KDNodeImplem.H
  @brief  Implementation of CD_KDNode.H
  @author Robert Marskar
*/

#ifndef CD_KDNodeImplem_H
#define CD_KDNodeImplem_H

// Std includes
#include <utility>
#include <type_traits>

// Our includes
#include <CD_KDNode.H>
#include <CD_NamespaceHeader.H>

template <class P>
inline KDNode<P>::KDNode() : m_left(nullptr), m_right(nullptr), m_weight(0.0)
{}

template <class P>
inline KDNode<P>::KDNode(ParticleList& a_particles)
  : m_left(nullptr), m_right(nullptr), m_weight(0.0), m_particles(std::move(a_particles))
{}

template <class P>
inline KDNode<P>::~KDNode()
{}

template <class P>
inline const typename KDNode<P>::ParticleList&
KDNode<P>::getParticles() const noexcept
{
  return m_particles;
}

template <class P>
inline typename KDNode<P>::ParticleList&
KDNode<P>::getParticles() noexcept
{
  return m_particles;
}

template <class P>
inline const Real&
KDNode<P>::weight() const noexcept
{
  return m_weight;
}

template <class P>
inline Real&
KDNode<P>::weight() noexcept
{
  return m_weight;
}

template <class P>
inline bool
KDNode<P>::isLeafNode() const noexcept
{
  return (m_left == nullptr) && (m_right == nullptr);
}

template <class P>
inline bool
KDNode<P>::isInteriorNode() const noexcept
{
  return !(this->isLeafNode());
}

template <class P>
inline typename KDNode<P>::ParticleList
KDNode<P>::gatherParticles() const noexcept
{
  ParticleList primitives;

  this->gatherParticles(primitives);

  return primitives;
}

template <class P>
inline void
KDNode<P>::gatherParticles(ParticleList& a_particles) const noexcept
{
  if (this->isLeafNode()) {
    a_particles.reserve(a_particles.size() + m_particles.size());
    a_particles.insert(a_particles.end(), m_particles.begin(), m_particles.end());
  }
  else {
    m_left->gatherParticles(a_particles);
    m_right->gatherParticles(a_particles);
  }
}

template <class P>
inline typename KDNode<P>::ParticleList
KDNode<P>::moveParticles() noexcept
{
  ParticleList primitives;

  this->moveParticles(primitives);

  return primitives;
}

template <class P>
inline void
KDNode<P>::moveParticles(ParticleList& a_particles) noexcept
{
  if (this->isLeafNode()) {
    a_particles.reserve(a_particles.size() + m_particles.size());

    std::move(m_particles.begin(), m_particles.end(), std::back_inserter(a_particles));
  }
  else {
    m_left->moveParticles(a_particles);
    m_right->moveParticles(a_particles);
  }
}

template <class P>
inline std::shared_ptr<KDNode<P>>&
KDNode<P>::getLeft() noexcept
{
  return m_left;
}

template <class P>
inline std::shared_ptr<KDNode<P>>&
KDNode<P>::getRight() noexcept
{
  return m_right;
}

#include <CD_NamespaceFooter.H>

#endif
