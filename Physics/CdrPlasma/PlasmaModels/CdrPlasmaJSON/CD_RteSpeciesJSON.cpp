/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RteSpeciesJSON.cpp
  @brief  Implementation of CD_RteSpeciesJSON.H
  @author Robert Marskar
*/
// Our includes
#include <CD_RteSpeciesJSON.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

RteSpeciesJSON::RteSpeciesJSON(const std::string& a_name, const RteSpeciesJSON::KappaFunction& a_kappaFunction) {
  this->define(a_name, a_kappaFunction);
}

RteSpeciesJSON::~RteSpeciesJSON() {
  // Nothing to do here. 
}

void RteSpeciesJSON::define(const std::string& a_name, const RteSpeciesJSON::KappaFunction& a_kappaFunction) {
  m_name               = a_name;
  m_absorptionFunction = a_kappaFunction;
}

Real RteSpeciesJSON::getAbsorptionCoefficient(const RealVect a_pos) const {
  return m_absorptionFunction(a_pos);
}
  
#include <CD_NamespaceFooter.H>
