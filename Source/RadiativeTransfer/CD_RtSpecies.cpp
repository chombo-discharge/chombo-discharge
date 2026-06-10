/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RtSpecies.cpp
  @brief  Implementation of CD_RtSpecies.H
  @author Robert Marskar
*/

// Our includes
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

RtSpecies::RtSpecies() : m_name("DefaultRtSpecies")
{
  // Default settings
}

RtSpecies::~RtSpecies() = default;

std::string
RtSpecies::getName() const
{
  return m_name;
}

Real
RtSpecies::getScatteringCoefficient(const RealVect /*a_pos*/) const
{
  return 0.0;
}

#include <CD_NamespaceFooter.H>
