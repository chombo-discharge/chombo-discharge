/* chombo-discharge
  * Copyright © 2021 SINTEF Energy Research.
  * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
  */

/*!
  @file   CD_ItoPlasmaGodunovParticle.cpp
  @brief  Implementation of CD_ItoPlasmaGodunovParticle.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ItoPlasmaGodunovParticle.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaGodunovParticle::ItoPlasmaGodunovParticle() : BinItem() {}

ItoPlasmaGodunovParticle::ItoPlasmaGodunovParticle(const RealVect a_position, const Real a_mass)
{
  m_mass     = a_mass;
  m_position = a_position;
}

ItoPlasmaGodunovParticle::~ItoPlasmaGodunovParticle() {}

void
ItoPlasmaGodunovParticle::define(const RealVect a_position, const Real a_mass)
{
  setMass(a_mass);
  setPosition(a_position);
}

void
ItoPlasmaGodunovParticle::setMass(const Real a_mass)
{
  m_mass = a_mass;
}

Real&
ItoPlasmaGodunovParticle::mass()
{
  return m_mass;
}

const Real&
ItoPlasmaGodunovParticle::mass() const
{
  return m_mass;
}

int
ItoPlasmaGodunovParticle::size() const
{
  return (BinItem::size() + sizeof(m_mass));
}

void
ItoPlasmaGodunovParticle::linearOut(void* buf) const
{
  Real* buffer = (Real*)buf;
  D_TERM6(*buffer++ = m_position[0];, *buffer++ = m_position[1];, *buffer++ = m_position[2];, *buffer++ = m_position[3];
          , *buffer++ = m_position[4];
          , *buffer++ = m_position[5];);

  *buffer = m_mass;
}

void
ItoPlasmaGodunovParticle::linearIn(void* buf)
{
  Real* buffer = (Real*)buf;
  D_TERM6(m_position[0] = *buffer++;, m_position[1] = *buffer++;, m_position[2] = *buffer++;, m_position[3] = *buffer++;
          , m_position[4] = *buffer++;
          , m_position[5] = *buffer++;);

  m_mass = *buffer;
}

#include <CD_NamespaceFooter.H>
