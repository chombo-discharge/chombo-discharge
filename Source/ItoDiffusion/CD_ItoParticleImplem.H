/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoParticleImplem.H
  @brief  Implementation of CD_ItoParticle.H
  @author Robert Marskar
*/

#ifndef CD_ItoParticleImplem_H
#define CD_ItoParticleImplem_H

// Std includes
#include <string>
#include <typeinfo>
#include <cstdio>
#include <cmath>

// Our includes
#include <CD_ItoParticle.H>
#include <CD_NamespaceHeader.H>

inline ItoParticle::ItoParticle() {}

inline ItoParticle::ItoParticle(const Real      a_weight,
                                const RealVect& a_position,
                                const RealVect& a_velocity,
                                const Real      a_diffusion,
                                const Real      a_mobility,
                                const Real      a_energy)
{
  this->define(a_weight, a_position, a_velocity, a_diffusion, a_mobility, a_energy);
}

inline ItoParticle::ItoParticle(const ItoParticle& a_other)
{
  this->weight()      = a_other.weight();
  this->position()    = a_other.position();
  this->oldPosition() = a_other.oldPosition();
  this->velocity()    = a_other.velocity();
  this->diffusion()   = a_other.diffusion();
  this->mobility()    = a_other.mobility();
  this->energy()      = a_other.energy();
}

inline ItoParticle::~ItoParticle() {}

inline void
ItoParticle::define(const Real      a_weight,
                    const RealVect& a_position,
                    const RealVect& a_velocity,
                    const Real      a_diffusion,
                    const Real      a_mobility,
                    const Real      a_energy)
{
  this->weight()    = a_weight;
  this->position()  = a_position;
  this->velocity()  = a_velocity;
  this->diffusion() = a_diffusion;
  this->mobility()  = a_mobility;
  this->energy()    = a_energy;
}

inline Real&
ItoParticle::weight()
{
  return this->real<0>();
}

inline const Real&
ItoParticle::weight() const
{
  return this->real<0>();
}

inline Real&
ItoParticle::mobility()
{
  return this->real<1>();
}

inline const Real&
ItoParticle::mobility() const
{
  return this->real<1>();
}

inline Real&
ItoParticle::diffusion()
{
  return this->real<2>();
}

inline const Real&
ItoParticle::diffusion() const
{
  return this->real<2>();
}

inline Real&
ItoParticle::energy()
{
  return this->real<3>();
}

inline const Real&
ItoParticle::energy() const
{
  return this->real<3>();
}

inline RealVect&
ItoParticle::oldPosition()
{
  return this->vect<0>();
}

inline const RealVect&
ItoParticle::oldPosition() const
{
  return this->vect<0>();
}

inline RealVect&
ItoParticle::velocity()
{
  return this->vect<1>();
}

inline const RealVect&
ItoParticle::velocity() const
{
  return this->vect<1>();
}

inline Real
ItoParticle::conductivity() const
{
  return this->weight() * this->mobility();
}

inline Real
ItoParticle::diffusivity() const
{
  return this->weight() * this->diffusion();
}

inline Real
ItoParticle::totalEnergy() const
{
  return this->weight() * this->energy();
}

inline Real&
ItoParticle::tmpReal()
{
  return this->real<4>();
}

inline const Real&
ItoParticle::tmpReal() const
{
  return this->real<4>();
}

inline RealVect&
ItoParticle::tmpVect()
{
  return this->vect<2>();
}

inline const RealVect&
ItoParticle::tmpVect() const
{
  return this->vect<2>();
}

#include <CD_NamespaceFooter.H>

#include <CD_ItoParticleImplem.H>

#endif
