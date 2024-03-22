/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NeutralSpeciesJSON.cpp
  @brief  Implementation of CD_NeutralSpeciesJSON.H
  @author Robert Marskar
*/

// Our includes
#include <CD_NeutralSpeciesJSON.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

NeutralSpeciesJSON::NeutralSpeciesJSON()
{
  m_isDefined = false;
}

NeutralSpeciesJSON::NeutralSpeciesJSON(const std::string                               a_name,
                                       const Real                                      a_molarFraction,
                                       const NeutralSpeciesJSON::NumberDensityFunction a_function)
{
  this->define(a_name, a_molarFraction, a_function);
}

NeutralSpeciesJSON::~NeutralSpeciesJSON()
{}

void
NeutralSpeciesJSON::define(const std::string                               a_name,
                           const Real                                      a_molarFraction,
                           const NeutralSpeciesJSON::NumberDensityFunction a_function)
{
  m_name          = a_name;
  m_molarFraction = a_molarFraction, m_function = a_function;
  m_isDefined = true;
}

std::string
NeutralSpeciesJSON::getName() const
{
  return m_name;
}

Real
NeutralSpeciesJSON::getMolarFraction() const
{
  return m_molarFraction;
}

Real
NeutralSpeciesJSON::operator()(const RealVect a_pos) const
{
  return m_function(a_pos);
}

#include <CD_NamespaceFooter.H>
