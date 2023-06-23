/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCBackgroundSpecies.H
  @brief  Implementation of ItoKMCBackgroundSpecies.H
  @author Robert Marskar
*/

#ifndef CD_ItoKMCBackgroundSpeciesImplem_H
#define CD_ItoKMCBackgroundSpeciesImplem_H

// Our includes
#include <CD_ItoKMCBackgroundSpecies.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

ItoKMCBackgroundSpecies::ItoKMCBackgroundSpecies() noexcept { m_isDefined = false; }

ItoKMCBackgroundSpecies::ItoKMCBackgroundSpecies(const std::string&   a_name,
                                                 const NumberDensity& a_numberDensity) noexcept
{
  this->define(a_name, a_numberDensity);
}

ItoKMCBackgroundSpecies::~ItoKMCBackgroundSpecies() noexcept {

}

void
ItoKMCBackgroundSpecies::define(const std::string& a_name, const NumberDensity& a_numberDensity) noexcept
{
  m_name          = a_name;
  m_numberDensity = a_numberDensity;
  m_isDefined     = true;
}

Real
ItoKMCBackgroundSpecies::operator()(const RealVect a_pos) const noexcept
{
  CH_assert(m_isDefined);

  return m_numberDensity(a_pos);
}

std::string
ItoKMCBackgroundSpecies::getName() const noexcept
{
  CH_assert(m_isDefined);

  return m_name;
}

#include <CD_NamespaceFooter.H>

#endif
