/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Dielectric.cpp
  @brief  Implementation of CD_Dielectric.H
  @author Robert marskar
  @date   Nov. 2017
*/

// Our includes
#include <CD_Dielectric.H>
#include <CD_NamespaceHeader.H>

Dielectric::Dielectric()
{
  CH_TIME("Dielectric::Dielectric()");

  m_isDefined = false;
}

Dielectric::Dielectric(const RefCountedPtr<BaseIF>& a_baseIF, const Real a_permittivity) : Dielectric()
{
  CH_TIME("Dielectric::Dielectric(RefCountedPtr<BaseIF>, Real)");

  CH_assert(!a_baseIF.isNull());

  this->define(a_baseIF, a_permittivity);
}

Dielectric::Dielectric(const RefCountedPtr<BaseIF>&                     a_baseIF,
                       const std::function<Real(const RealVect a_pos)>& a_permittivity)
  : Dielectric()
{
  CH_TIME("Dielectric::Dielectric(RefCountedPtr<BaseIF>, std::function<Real(const RealVect a_pos)>)");

  CH_assert(!a_baseIF.isNull());

  this->define(a_baseIF, a_permittivity);
}

Dielectric::~Dielectric()
{}

void
Dielectric::define(const RefCountedPtr<BaseIF>& a_baseIF, const Real a_permittivity)
{
  CH_TIME("Dielectric::define(RefCountedPtr<BaseIF>, Real)");

  CH_assert(!a_baseIF.isNull());

  m_baseIF               = a_baseIF;
  m_constantPermittivity = a_permittivity;
  m_useConstant          = true;
  m_isDefined            = true;
}

void
Dielectric::define(const RefCountedPtr<BaseIF>&                     a_baseIF,
                   const std::function<Real(const RealVect a_pos)>& a_permittivity)
{
  CH_TIME("Dielectric::define(RefCountedPtr<BaseIF>, std::function<Real(const RealVect a_pos)>)");

  CH_assert(!a_baseIF.isNull());

  m_baseIF               = a_baseIF;
  m_variablePermittivity = a_permittivity;
  m_useConstant          = false;
  m_isDefined            = true;
}

const RefCountedPtr<BaseIF>&
Dielectric::getImplicitFunction() const
{
  CH_TIME("Dielectric::getImplicitFunction()");

  CH_assert(m_isDefined);

  return (m_baseIF);
}

Real
Dielectric::getPermittivity(const RealVect a_pos) const
{
  CH_TIME("Dielectric::getPermittivity(RealVect)");

  CH_assert(m_isDefined);

  Real ret = 1.0;

  if (m_useConstant) {
    ret = m_constantPermittivity;
  }
  else {
    ret = m_variablePermittivity(a_pos);
  }

  return ret;
}

#include <CD_NamespaceFooter.H>
