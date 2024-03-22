/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtSpecies.cpp
  @brief  Implementation of CD_RtSpecies.H
  @author Robert Marskar
*/

// Our includes
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

RtSpecies::RtSpecies()
{
  // Default settings
  m_name = "DefaultRtSpecies";
}

RtSpecies::~RtSpecies()
{}

std::string
RtSpecies::getName() const
{
  return m_name;
}

Real
RtSpecies::getScatteringCoefficient(const RealVect a_pos) const
{
  return 0.0;
}

#include <CD_NamespaceFooter.H>
