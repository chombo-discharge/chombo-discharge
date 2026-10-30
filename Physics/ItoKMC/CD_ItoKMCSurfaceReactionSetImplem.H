/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ItoKMCSurfaceReactionSetImplem.H
 * @brief  Implementation of CD_ItoKMCSurfaceReactionSet.H
 * @author Robert Marskar
 */

#ifndef CD_ITOKMCSURFACEREACTIONSETIMPLEM_H
#define CD_ITOKMCSURFACEREACTIONSETIMPLEM_H

// Our includes
#include <CD_ItoKMCSurfaceReactionSet.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

inline ItoKMCSurfaceReactionSet::~ItoKMCSurfaceReactionSet() noexcept
{}

inline void
ItoKMCSurfaceReactionSet::add(const size_t                            a_speciesIndex,
                              const ItoKMCSurfaceReaction&            a_reaction,
                              const ItoKMCSurfaceReactionSet::Surface a_surface,
                              const ItoKMCSurfaceReactionSet::Species a_speciesType) noexcept
{
  if (a_surface == Surface::Dielectric) {
    if (a_speciesType == Species::Plasma) {
      m_dielectricPlasmaReactions[a_speciesIndex].add(a_reaction);
    }
    else if (a_speciesType == Species::Photon) {
      m_dielectricPhotonReactions[a_speciesIndex].add(a_reaction);
    }
    else {
      MayDay::Error("ItoKMCSurfaceReaction::add - logic bust 1");
    }
  }
  else if (a_surface == Surface::Electrode) {
    if (a_speciesType == Species::Plasma) {
      m_electrodePlasmaReactions[a_speciesIndex].add(a_reaction);
    }
    else if (a_speciesType == Species::Photon) {
      m_electrodePhotonReactions[a_speciesIndex].add(a_reaction);
    }
    else {
      MayDay::Error("ItoKMCSurfaceReaction::add - logic bust 2");
    }
  }
  else {
    MayDay::Error("ItoKMCSurfaceReaction::add - logic bust 3");
  }
}

inline void
ItoKMCSurfaceReactionSet::addDomain(const size_t                            a_speciesIndex,
                                    const ItoKMCSurfaceReaction&            a_reaction,
                                    const int                               a_dir,
                                    const Side::LoHiSide                    a_side,
                                    const ItoKMCSurfaceReactionSet::Species a_speciesType) noexcept
{
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  const int face = domainFaceIndex(a_dir, a_side);

  if (a_speciesType == Species::Plasma) {
    m_domainPlasmaReactions[face][a_speciesIndex].add(a_reaction);
  }
  else if (a_speciesType == Species::Photon) {
    m_domainPhotonReactions[face][a_speciesIndex].add(a_reaction);
  }
  else {
    MayDay::Error("ItoKMCSurfaceReaction::addDomain - logic bust");
  }
}

inline const std::map<size_t, ItoKMCSurfaceReactions>&
ItoKMCSurfaceReactionSet::getDielectricPlasmaReactions() const noexcept
{
  return m_dielectricPlasmaReactions;
}

inline const std::map<size_t, ItoKMCSurfaceReactions>&
ItoKMCSurfaceReactionSet::getDielectricPhotonReactions() const noexcept
{
  return m_dielectricPhotonReactions;
}

inline const std::map<size_t, ItoKMCSurfaceReactions>&
ItoKMCSurfaceReactionSet::getElectrodePlasmaReactions() const noexcept
{
  return m_electrodePlasmaReactions;
}

inline const std::map<size_t, ItoKMCSurfaceReactions>&
ItoKMCSurfaceReactionSet::getElectrodePhotonReactions() const noexcept
{
  return m_electrodePhotonReactions;
}

inline const std::map<size_t, ItoKMCSurfaceReactions>&
ItoKMCSurfaceReactionSet::getDomainPlasmaReactions(const int a_dir, const Side::LoHiSide a_side) const noexcept
{
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  return m_domainPlasmaReactions[domainFaceIndex(a_dir, a_side)];
}

inline const std::map<size_t, ItoKMCSurfaceReactions>&
ItoKMCSurfaceReactionSet::getDomainPhotonReactions(const int a_dir, const Side::LoHiSide a_side) const noexcept
{
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  return m_domainPhotonReactions[domainFaceIndex(a_dir, a_side)];
}

inline int
ItoKMCSurfaceReactionSet::domainFaceIndex(const int a_dir, const Side::LoHiSide a_side) noexcept
{
  return 2 * a_dir + ((a_side == Side::Lo) ? 0 : 1);
}

#include <CD_NamespaceFooter.H>

#endif
