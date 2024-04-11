/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PhotonImplem.H
  @brief  Implementation of CD_Photon.H
  @author Robert Marskar
*/

#ifndef CD_PhotonImplem_H
#define CD_PhotonImplem_H

// Our includes
#include <CD_Photon.H>
#include <CD_NamespaceHeader.H>

inline Photon::Photon()
{}

inline Photon::~Photon()
{}

inline Photon::Photon(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_weight)
{
  this->define(a_position, a_velocity, a_kappa, a_weight);
}

inline void
Photon::define(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_weight)
{
  this->position() = a_position;
  this->velocity() = a_velocity;
  this->kappa()    = a_kappa;
  this->weight()   = a_weight;
}

inline Real&
Photon::weight()
{
  return this->real<0>();
}

inline const Real&
Photon::weight() const
{
  return this->real<0>();
}

inline Real&
Photon::kappa()
{
  return this->real<1>();
}

inline const Real&
Photon::kappa() const
{
  return this->real<1>();
}

inline RealVect&
Photon::velocity()
{
  return this->vect<0>();
}

inline const RealVect&
Photon::velocity() const
{
  return this->vect<0>();
}

#include <CD_NamespaceFooter.H>

#endif
