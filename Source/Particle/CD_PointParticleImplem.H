/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_PointParticleImplem.H
  @brief  Implementation of CD_PointParticle.H
  @author Robert Marskar
*/

#ifndef CD_POINTPARTICLEIMPLEM_H
#define CD_POINTPARTICLEIMPLEM_H

// Our includes
#include <CD_PointParticle.H>
#include <CD_NamespaceHeader.H>

PointParticle::PointParticle() = default;

inline PointParticle::PointParticle(const RealVect& a_position, const Real a_weight)
{
  this->define(a_position, a_weight);
}

PointParticle::~PointParticle() = default;

inline void
PointParticle::define(const RealVect& a_position, const Real a_weight)
{
  this->position() = a_position;
  this->weight()   = a_weight;
}

inline Real&
PointParticle::weight()
{
  return this->real<0>();
}

inline const Real&
PointParticle::weight() const
{
  return this->real<0>();
}

#include <CD_NamespaceFooter.H>

#endif
