/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaPhotoReactionImplem.H
  @brief  Implementation of CD_ItoPlasmaPhotoReaction.H
  @author Robert Marskar
*/

// Our includes
#ifndef CD_ItoPlasmaPhotoReactionImplem_H
#define CD_ItoPlasmaPhotoReactionImplem_H

#include <CD_ItoPlasmaPhotoReaction.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

inline ItoPlasmaPhotoReaction::ItoPlasmaPhotoReaction(const size_t a_source, const std::list<size_t> a_targets) noexcept
{
  this->define(a_source, a_targets);
}

inline ItoPlasmaPhotoReaction::~ItoPlasmaPhotoReaction() noexcept {}

inline void
ItoPlasmaPhotoReaction::define(const size_t a_source, const std::list<size_t> a_targets) noexcept
{

  m_source  = a_source;
  m_targets = a_targets;
}

inline const size_t&
ItoPlasmaPhotoReaction::getSourcePhoton() const noexcept
{
  return m_source;
}

inline const std::list<size_t>&
ItoPlasmaPhotoReaction::getTargetSpecies() const noexcept
{
  return m_targets;
}

#include <CD_NamespaceFooter.H>

#endif
