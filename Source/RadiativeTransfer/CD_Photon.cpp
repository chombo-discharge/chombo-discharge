/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Photon.cpp
  @brief  Implementation of CD_Photon.H
  @author Robert Marskar
*/

// Our includes
#include <CD_Photon.H>
#include <CD_NamespaceHeader.H>

Photon::Photon() : BinItem() {}

Photon::~Photon() {}

Photon::Photon(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_mass)
{
  this->define(a_position, a_velocity, a_kappa, a_mass);
}

void
Photon::define(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_mass)
{
  this->setPosition(a_position);
  this->setVelocity(a_velocity);
  this->setKappa(a_kappa);
  this->setMass(a_mass);
}

void
Photon::setKappa(const Real a_kappa)
{
  m_kappa = a_kappa;
}

void
Photon::setMass(const Real a_mass)
{
  m_mass = a_mass;
}

void
Photon::setVelocity(const RealVect& a_velocity)
{
  m_velocity = a_velocity;
}

Real&
Photon::mass()
{
  return m_mass;
}
Real&
Photon::kappa()
{
  return m_kappa;
}

RealVect&
Photon::velocity()
{
  return m_velocity;
}

const Real&
Photon::kappa() const
{
  return m_kappa;
}

const Real&
Photon::mass() const
{
  return m_mass;
}

const RealVect&
Photon::velocity() const
{
  return m_velocity;
}

bool
Photon::operator==(const Photon* a_p) const
{
  return (*this == *a_p);
}

bool
Photon::operator!=(const Photon& a_p) const
{
  return !(*this == a_p);
}

bool
Photon::operator==(const Photon& a_p) const
{
  return (m_position == a_p.m_position && m_velocity == a_p.m_velocity && m_kappa == a_p.m_kappa &&
          m_mass == a_p.m_mass);
}

int
Photon::size() const
{
  return (BinItem::size() + sizeof(m_velocity) + sizeof(m_kappa) + sizeof(m_mass));
}

void
Photon::linearOut(void* buf) const
{
  Real* buffer = (Real*)buf;
  D_TERM6(*buffer++ = m_position[0];, *buffer++ = m_position[1];, *buffer++ = m_position[2];, *buffer++ = m_position[3];
          , *buffer++ = m_position[4];
          , *buffer++ = m_position[5];);

  D_TERM6(*buffer++ = m_velocity[0];, *buffer++ = m_velocity[1];, *buffer++ = m_velocity[2];, *buffer++ = m_velocity[3];
          , *buffer++ = m_velocity[4];
          , *buffer++ = m_velocity[5];);

  *buffer = m_kappa;
  buffer++;
  *buffer = m_mass;
}

void
Photon::linearIn(void* buf)
{
  Real* buffer = (Real*)buf;
  D_TERM6(m_position[0] = *buffer++;, m_position[1] = *buffer++;, m_position[2] = *buffer++;, m_position[3] = *buffer++;
          , m_position[4] = *buffer++;
          , m_position[5] = *buffer++;);

  D_TERM6(m_velocity[0] = *buffer++;, m_velocity[1] = *buffer++;, m_velocity[2] = *buffer++;, m_velocity[3] = *buffer++;
          , m_velocity[4] = *buffer++;
          , m_velocity[5] = *buffer++;);

  m_kappa = *buffer;
  buffer++;
  m_mass = *buffer;
}

#include <CD_NamespaceFooter.H>
