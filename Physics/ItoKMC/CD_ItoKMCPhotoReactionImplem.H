/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCPhotoReactionImplem.H
  @brief  Implementation of CD_ItoKMCPhotoReaction.H
  @author Robert Marskar
*/

// Our includes
#ifndef CD_ItoKMCPhotoReactionImplem_H
#define CD_ItoKMCPhotoReactionImplem_H

#include <CD_ItoKMCPhotoReaction.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

inline ItoKMCPhotoReaction::ItoKMCPhotoReaction(const size_t            a_source,
                                                const std::list<size_t> a_targets,
                                                const Real              a_efficiency) noexcept
{
  this->define(a_source, a_targets, a_efficiency);
}

inline ItoKMCPhotoReaction::~ItoKMCPhotoReaction() noexcept
{}

inline void
ItoKMCPhotoReaction::define(const size_t a_source, const std::list<size_t> a_targets, const Real a_efficiency) noexcept
{
  m_source     = a_source;
  m_targets    = a_targets;
  m_efficiency = a_efficiency;

  const std::string baseError = "CD_ItoKMCPhotoReactionImplem.H in function ";
  if (m_efficiency < 0.0) {
    const std::string derivError = "'ItoKMCPhotoReaction::define' - efficiency must be >= 0.0";

    MayDay::Error((baseError + derivError).c_str());
  }
}

inline const size_t&
ItoKMCPhotoReaction::getSourcePhoton() const noexcept
{
  return m_source;
}

inline const std::list<size_t>&
ItoKMCPhotoReaction::getTargetSpecies() const noexcept
{
  return m_targets;
}

inline const Real&
ItoKMCPhotoReaction::getEfficiency() const noexcept
{
  return m_efficiency;
}

#include <CD_NamespaceFooter.H>

#endif
