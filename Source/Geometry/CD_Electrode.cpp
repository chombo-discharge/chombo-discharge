/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Electrode.cpp
  @brief  Implementation of CD_Electrode.H
  @author Robert marskar
*/

// Our includes
#include <CD_Electrode.H>
#include <CD_NamespaceHeader.H>

Electrode::Electrode()
{
  CH_TIME("Electrode::Electrode()");

  m_isDefined = false;
}

Electrode::Electrode(const RefCountedPtr<BaseIF>& a_baseIF, const bool a_live, const Real a_voltageFraction)
  : Electrode()
{
  CH_TIME("Electrode::Electrode(RefCountedPtr<BaseIF>, bool, Real");

  CH_assert(!a_baseIF.isNull());

  this->define(a_baseIF, a_live, a_voltageFraction);
}

Electrode::~Electrode()
{}

void
Electrode::define(const RefCountedPtr<BaseIF>& a_baseIF, const bool a_live, const Real a_voltageFraction)
{
  CH_TIME("Electrode::define(RefCountedPtr<BaseIF>, bool, Real");

  CH_assert(!a_baseIF.isNull());

  m_baseIF          = a_baseIF;
  m_isLive          = a_live;
  m_voltageFraction = a_voltageFraction;
  m_isDefined       = true;
}

const RefCountedPtr<BaseIF>&
Electrode::getImplicitFunction() const
{
  CH_TIME("Electrode::getImplicitFunction()");

  CH_assert(m_isDefined);

  return (m_baseIF);
}

const bool&
Electrode::isLive() const
{
  CH_TIME("Electrode::isLive()");

  CH_assert(m_isDefined);

  return (m_isLive);
}

const Real&
Electrode::getFraction() const
{
  CH_TIME("Electrode::getFraction()");

  CH_assert(m_isDefined);

  return (m_voltageFraction);
}

#include <CD_NamespaceFooter.H>
