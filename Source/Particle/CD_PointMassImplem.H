/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PointMassImplem.H
  @brief  Implementation of CD_PointMass.H
  @author Robert Marskar
*/

#ifndef CD_PointMassImplem_H
#define CD_PointMassImplem_H

// Our includes
#include <CD_PointMass.H>
#include <CD_NamespaceHeader.H>

inline PointMass::PointMass()
{
  m_pos    = RealVect::Zero;
  m_mass   = 0.0;
  m_energy = 0.0;
}

inline PointMass::PointMass(const RealVect a_pos, const Real a_mass, const Real a_energy)
{
  m_pos    = a_pos;
  m_mass   = a_mass;
  m_energy = a_energy;
}

inline PointMass::PointMass(const std::vector<PointMass>& a_pointMasses)
{
  m_pos    = RealVect::Zero;
  m_mass   = 0.0;
  m_energy = 0.0;

  for (const auto& pm : a_pointMasses) {
    const RealVect& p = pm.pos();
    const Real&     m = pm.mass();
    const Real&     e = pm.energy();

    m_mass += m;
    m_pos += m * p;
    m_energy += m * e;
  }

  m_pos    = m_pos / m_mass;
  m_energy = m_energy / m_mass;
}

inline PointMass::~PointMass() {}

inline void
PointMass::define(const RealVect& a_pos, const Real& a_mass, const Real& a_energy)
{
  m_pos    = a_pos;
  m_mass   = a_mass;
  m_energy = a_energy;
}

inline Real&
PointMass::operator[](const int a_dir)
{
  return m_pos[a_dir];
}

inline const Real&
PointMass::operator[](const int a_dir) const
{
  return m_pos[a_dir];
}

inline RealVect&
PointMass::pos()
{
  return m_pos;
}

inline const RealVect&
PointMass::pos() const
{
  return m_pos;
}

inline Real&
PointMass::mass()
{
  return m_mass;
}

inline const Real&
PointMass::mass() const
{
  return m_mass;
}

inline Real&
PointMass::energy()
{
  return m_energy;
}

inline const Real&
PointMass::energy() const
{
  return m_energy;
}

inline bool
PointMass::canSplit() const
{
  return m_mass * (1.0 + 1.E-6) >= 2.0;
}

#include <CD_NamespaceFooter.H>

#endif
