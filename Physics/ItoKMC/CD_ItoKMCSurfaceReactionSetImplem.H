/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCSurfaceReactionSetImplem.H
  @brief  Declaration of a class that encapsulates secondary emission at EBs and domain edges/faces
  @author Robert Marskar
*/

#ifndef CD_ItoKMCSurfaceReactionSetImplem_H
#define CD_ItoKMCSurfaceReactionSetImplem_H

// Our includes
#include <CD_ItoKMCSurfaceReactionSet.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

inline ItoKMCSurfaceReactionSet::~ItoKMCSurfaceReactionSet() noexcept {}

inline void
ItoKMCSurfaceReactionSet::addSurfaceReaction(const ItoKMCSurfaceReaction&            a_reaction,
                                             const ItoKMCSurfaceReactionSet::Surface a_surface,
                                             const ItoKMCSurfaceReactionSet::Species a_species) noexcept
{
  if (a_surface == Surface::Dielectric) {
    if (a_species == Species::Plasma) {
      m_dielectricPlasmaReactions.push_back(a_reaction);
    }
    else if (a_species == Species::Photon) {
      m_dielectricPhotonReactions.push_back(a_reaction);
    }
    else {
      MayDay::Error("ItoKMCSurfaceReaction::addSurfaceReaction - logic bust 1");
    }
  }
  else if (a_surface == Surface::Electrode) {
    if (a_species == Species::Plasma) {
      m_electrodePlasmaReactions.push_back(a_reaction);
    }
    else if (a_species == Species::Photon) {
      m_electrodePhotonReactions.push_back(a_reaction);
    }
    else {
      MayDay::Error("ItoKMCSurfaceReaction::addSurfaceReaction - logic bust 2");
    }
  }
  else {
    MayDay::Error("ItoKMCSurfaceReaction::addSurfaceReaction - logic bust 3");
  }
}

inline const std::vector<ItoKMCSurfaceReaction>&
ItoKMCSurfaceReactionSet::getDielectricPlasmaReactions() const noexcept
{
  return m_dielectricPlasmaReactions;
}

inline const std::vector<ItoKMCSurfaceReaction>&
ItoKMCSurfaceReactionSet::getDielectricPhotonReactions() const noexcept
{
  return m_dielectricPhotonReactions;
}

inline const std::vector<ItoKMCSurfaceReaction>&
ItoKMCSurfaceReactionSet::getElectrodePlasmaReactions() const noexcept
{
  return m_electrodePlasmaReactions;
}

inline const std::vector<ItoKMCSurfaceReaction>&
ItoKMCSurfaceReactionSet::getElectrodePhotonReactions() const noexcept
{
  return m_electrodePhotonReactions;
}

#include <CD_NamespaceFooter.H>

#endif
