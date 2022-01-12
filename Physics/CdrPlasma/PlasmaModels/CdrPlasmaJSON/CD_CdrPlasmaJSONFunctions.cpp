/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaJSONFunctions.cpp
  @brief  Implementation of CD_CdrPlasmaJSONFunctions.H
  @author Robert Marskar
*/

// Our includes
#include <CD_CdrPlasmaJSONFunctions.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaJSONFunctionENA::CdrPlasmaJSONFunctionENA(const Real a_A, const Real a_B, const Real a_C) {
  m_A = a_A;
  m_B = a_B;
  m_C = a_C;  
}

CdrPlasmaJSONFunctionENA::~CdrPlasmaJSONFunctionENA(){

}

Real CdrPlasmaJSONFunctionENA::operator()(const RealVect a_E, const Real a_N) const {
  const Real E = a_E.vectorLength();

  return m_A * std::pow(E, m_B) / std::pow(a_N, m_C);
}

#include <CD_NamespaceFooter.H>
