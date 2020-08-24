/*!
  @file   photo_reactionI.H
  @brief  Implementation of photo_reaction.H
  @author Robert Marskar
  @date   Aug. 2020
*/

#ifndef _PHOTO_REACTIONI_H_
#define _PHOTO_REACTIONI_H_

#include "photo_reaction.H"

using namespace physics::ito_plasma;

inline const int& photo_reaction::source() const {
  return m_source;
}

inline const std::vector<int>& photo_reaction::targets() const {
  return m_targets;
}

#endif
