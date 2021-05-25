/*!
  @file   ItoPlasmaPhotoReactionI.H
  @brief  Implementation of ItoPlasmaPhotoReaction.H
  @author Robert Marskar
  @date   Aug. 2020
*/

#ifndef _PHOTO_REACTIONI_H_
#define _PHOTO_REACTIONI_H_

#include <CD_ItoPlasmaPhotoReaction.H>

#include "CD_NamespaceHeader.H"
using namespace Physics::ItoPlasma;

inline
const int& ItoPlasmaPhotoReaction::source() const {
  return m_source;
}

inline
const std::deque<int>& ItoPlasmaPhotoReaction::targets() const {
  return m_targets;
}

inline
const std::deque<std::pair<int, Real> >& ItoPlasmaPhotoReaction::get_target_energies() const{
  return m_energies;
}
#include "CD_NamespaceFooter.H"

#endif
