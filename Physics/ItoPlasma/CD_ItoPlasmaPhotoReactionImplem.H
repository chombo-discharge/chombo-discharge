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

inline const int&
ItoPlasmaPhotoReaction::source() const
{
  return m_source;
}

inline const std::deque<int>&
ItoPlasmaPhotoReaction::targets() const
{
  return m_targets;
}

inline const std::deque<std::pair<int, Real>>&
ItoPlasmaPhotoReaction::getTargetEnergies() const
{
  return m_energies;
}

#include <CD_NamespaceFooter.H>

#endif
