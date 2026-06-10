/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoKMCBackgroundSpecies.cpp
  @brief  Implementation of CD_ItoKMCBackgroundSpecies.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ItoKMCBackgroundSpecies.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

ItoKMCBackgroundSpecies::ItoKMCBackgroundSpecies() noexcept : m_isDefined(false)
{}

ItoKMCBackgroundSpecies::ItoKMCBackgroundSpecies(const std::string&   a_name,
                                                 const MolarFraction& a_molarFraction) noexcept
{
  this->define(a_name, a_molarFraction);
}

ItoKMCBackgroundSpecies::~ItoKMCBackgroundSpecies() noexcept = default;

void
ItoKMCBackgroundSpecies::define(const std::string& a_name, const MolarFraction& a_molarFraction) noexcept
{
  m_name          = a_name;
  m_molarFraction = a_molarFraction;
  m_isDefined     = true;
}

Real
ItoKMCBackgroundSpecies::molarFraction(const RealVect& a_pos) const noexcept
{
  CH_assert(m_isDefined);

  return m_molarFraction(a_pos);
}

std::string
ItoKMCBackgroundSpecies::getName() const noexcept
{
  CH_assert(m_isDefined);

  return m_name;
}

#include <CD_NamespaceFooter.H>
