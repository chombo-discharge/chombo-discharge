/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoSpecies.cpp
  @brief  Implementation of CD_ItoSpecies.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ItoSpecies.H>
#include <CD_NamespaceHeader.H>

ItoSpecies::ItoSpecies()
{
  m_name         = "ItoSpecies";
  m_isMobile     = false;
  m_isDiffusive  = false;
  m_chargeNumber = 0;
}

ItoSpecies::ItoSpecies(const std::string a_name, const int a_chargeNumber, const bool a_mobile, const bool a_diffusive)
{
  m_name         = a_name;
  m_chargeNumber = a_chargeNumber;
  m_isMobile     = a_mobile;
  m_isDiffusive  = a_diffusive;
}

ItoSpecies::~ItoSpecies()
{}

std::string
ItoSpecies::getName() const
{
  return m_name;
}

int
ItoSpecies::getChargeNumber() const
{
  return m_chargeNumber;
}

bool
ItoSpecies::isDiffusive() const
{
  return m_isDiffusive;
}

bool
ItoSpecies::isMobile() const
{
  return m_isMobile;
}

List<ItoParticle>&
ItoSpecies::getInitialParticles()
{
  return m_initialParticles;
}

const List<ItoParticle>&
ItoSpecies::getInitialParticles() const
{
  return m_initialParticles;
}

Real
ItoSpecies::mobility(const Real a_energy) const
{
  return 0.0;
}

Real
ItoSpecies::diffusion(const Real a_energy) const
{
  return 0.0;
}

#include <CD_NamespaceFooter.H>
