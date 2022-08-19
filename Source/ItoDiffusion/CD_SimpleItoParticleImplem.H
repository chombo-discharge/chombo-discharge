/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SimpleItoParticleImplem.H
  @brief  Implementation of CD_SimpleItoParticle.H
  @author Robert Marskar
*/

#ifndef CD_SimpleItoParticleParticleImplem_H
#define CD_SimpleItoParticleParticleImplem_H

// Our includes
#include <CD_SimpleItoParticle.H>
#include <CD_NamespaceHeader.H>

inline SimpleItoParticle::SimpleItoParticle() {}

inline SimpleItoParticle::SimpleItoParticle(const Real a_weight, const RealVect a_position, const Real a_energy)
{
  this->define(a_weight, a_position, a_energy);
}

inline SimpleItoParticle::~SimpleItoParticle() {}

inline void
SimpleItoParticle::define(const Real a_weight, const RealVect a_position, const Real a_energy)
{
  this->weight()   = a_weight;
  this->position() = a_position;
  this->energy()   = a_energy;
}

inline Real&
SimpleItoParticle::weight()
{
  return this->real<0>();
}

inline const Real&
SimpleItoParticle::weight() const
{
  return this->real<0>();
}

inline Real&
SimpleItoParticle::energy()
{
  return this->real<1>();
}

inline const Real&
SimpleItoParticle::energy() const
{
  return this->real<1>();
}

#include <CD_NamespaceFooter.H>

#endif
