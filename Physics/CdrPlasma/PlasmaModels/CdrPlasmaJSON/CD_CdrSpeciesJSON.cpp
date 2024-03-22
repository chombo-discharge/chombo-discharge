/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrSpeciesJSON.cpp
  @brief  Implementation of CD_CdrSpeciesJSON.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_assert.H>

// Our includes
#include <CD_CdrSpeciesJSON.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrSpeciesJSON::CdrSpeciesJSON(const std::string                  a_name,
                               const int                          a_Z,
                               const bool                         a_diffusive,
                               const bool                         a_mobile,
                               const CdrSpeciesJSON::InitFunction a_initialData)
{
  this->define(a_name, a_Z, a_diffusive, a_mobile, a_initialData);
}

CdrSpeciesJSON::~CdrSpeciesJSON()
{}

void
CdrSpeciesJSON::define(const std::string                  a_name,
                       const int                          a_Z,
                       const bool                         a_diffusive,
                       const bool                         a_mobile,
                       const CdrSpeciesJSON::InitFunction a_initialData)
{
  m_name         = a_name;
  m_chargeNumber = a_Z;
  m_isDiffusive  = a_diffusive;
  m_isMobile     = a_mobile;
  m_initFunction = a_initialData;
  m_isDefined    = true;
}

Real
CdrSpeciesJSON::initialData(const RealVect a_pos, const Real a_time) const
{
  CH_assert(m_isDefined);

  return m_initFunction(a_pos, a_time);
}

#include <CD_NamespaceFooter.H>
