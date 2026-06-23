/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_PhotonImplem.H
  @brief  Implementation of the PhotonSeed bridge class in CD_Photon.H
  @author Robert Marskar
*/

#ifndef CD_PHOTONIMPLEM_H
#define CD_PHOTONIMPLEM_H

// Our includes
#include <CD_Photon.H>
#include <CD_NamespaceHeader.H>

inline PhotonSeed::PhotonSeed()
{}

inline PhotonSeed::~PhotonSeed()
{}

inline PhotonSeed::PhotonSeed(const RealVect& a_position,
                              const RealVect& a_velocity,
                              const Real&     a_kappa,
                              const Real      a_weight)
{
  this->define(a_position, a_velocity, a_kappa, a_weight);
}

inline void
PhotonSeed::define(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_weight)
{
  this->position() = a_position;
  this->velocity() = a_velocity;
  this->kappa()    = a_kappa;
  this->weight()   = a_weight;
}

inline Real&
PhotonSeed::weight()
{
  return this->real<0>();
}

inline const Real&
PhotonSeed::weight() const
{
  return this->real<0>();
}

inline Real&
PhotonSeed::kappa()
{
  return this->real<1>();
}

inline const Real&
PhotonSeed::kappa() const
{
  return this->real<1>();
}

inline RealVect&
PhotonSeed::velocity()
{
  return this->vect<0>();
}

inline const RealVect&
PhotonSeed::velocity() const
{
  return this->vect<0>();
}

#include <CD_NamespaceFooter.H>

#endif
