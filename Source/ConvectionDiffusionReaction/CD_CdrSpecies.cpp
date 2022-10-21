/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrSpecies.cpp
  @brief  Implementation of CD_CdrSpecies.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_CdrSpecies.H>
#include <CD_NamespaceHeader.H>

CdrSpecies::CdrSpecies()
{
  CH_TIME("CdrSpecies::CdrSpecies()");

  // Default settings
  m_name         = "CdrSpecies";
  m_chargeNumber = 0;
  m_isDiffusive  = true;
  m_isMobile     = true;

  m_initialParticles.clear();
}

CdrSpecies::CdrSpecies(const std::string a_name,
                       const int         a_chargeNumber,
                       const bool        a_isMobile,
                       const bool        a_isDiffusive)
{
  CH_TIME("CdrSpecies::CdrSpecies(string, int, bool, bool)");

  m_name         = a_name;
  m_chargeNumber = a_chargeNumber;
  m_isMobile     = a_isMobile;
  m_isDiffusive  = a_isDiffusive;

  m_initialParticles.clear();
}

CdrSpecies::~CdrSpecies() {}

Real
CdrSpecies::initialData(const RealVect a_pos, const Real a_time) const
{
  CH_TIME("CdrSpecies::initialData(RealVect, Real)");

  return 0.;
}

std::string
CdrSpecies::getName() const
{
  CH_TIME("CdrSpecies::getName()");

  return m_name;
}

int
CdrSpecies::getChargeNumber() const
{
  CH_TIME("CdrSpecies::getChargeNumber()");

  return m_chargeNumber;
}

bool
CdrSpecies::isDiffusive() const
{
  CH_TIME("CdrSpecies::isDiffusive()");

  return m_isDiffusive;
}

bool
CdrSpecies::isMobile() const
{
  CH_TIME("CdrSpecies::isMobile()");

  return m_isMobile;
}

const List<PointParticle>&
CdrSpecies::getInitialParticles() const
{
  CH_TIME("CdrSpecies::getInitialParticles()");

  return m_initialParticles;
}

List<PointParticle>&
CdrSpecies::getInitialParticles()
{
  CH_TIME("CdrSpecies::getInitialParticles()");

  return m_initialParticles;
}

#include <CD_NamespaceFooter.H>
