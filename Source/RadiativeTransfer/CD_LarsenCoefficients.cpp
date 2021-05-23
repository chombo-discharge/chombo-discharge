/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LarsenCoefficients.cpp
  @brief  Implementation of CD_LarsenCoefficients.H
  @author Robert Marskar
*/

// Our includes
#include <CD_LarsenCoefficients.H>
#include <CD_NamespaceHeader.H>

LarsenCoefficients::LarsenCoefficients(const RefCountedPtr<RtSpecies>& a_RtSpecies,
				       const Real                         a_r1,
				       const Real                         a_r2){
  m_RtSpecies = a_RtSpecies;
  m_reflectionCoefficientOne = a_r1;
  m_reflectionCoefficientTwo = a_r2;
}

LarsenCoefficients::~LarsenCoefficients(){

}

Real LarsenCoefficients::aco(const RealVect a_pos) const {

  Real val = m_RtSpecies->getKappa(a_pos)*m_RtSpecies->getKappa(a_pos);
  val *= 3.0/2.0;
  val *= (1 + 3*m_reflectionCoefficientTwo)/(1 - 2*m_reflectionCoefficientOne);

  return val;
}

Real LarsenCoefficients::bco(const RealVect a_pos) const {
  return -m_RtSpecies->getKappa(a_pos);
}

Real LarsenCoefficients::rhs(const RealVect a_pos) const {
  return 0.;
}

#include <CD_NamespaceFooter.H>
