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
				       const Real                      a_r1,
				       const Real                      a_r2){
  m_RtSpecies                = a_RtSpecies;
  m_reflectionCoefficientOne = a_r1;
  m_reflectionCoefficientTwo = a_r2;

  m_bcFunction = [](const RealVect a_pos, const Real a_time){
    return 0.0;
  };
}

LarsenCoefficients::LarsenCoefficients(const RefCountedPtr<RtSpecies>& a_RtSpecies,
				       const Real                      a_r1,
				       const Real                      a_r2,
				       const std::function<Real(const RealVect a_pos, const Real a_time)>& a_bcFunction){
  m_RtSpecies                = a_RtSpecies;
  m_reflectionCoefficientOne = a_r1;
  m_reflectionCoefficientTwo = a_r2;

  m_bcFunction = a_bcFunction;
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
  const Real dummyDt = 0.0;

  return m_bcFunction(a_pos, dummyDt);
}

void LarsenCoefficients::setRhsFunction(const std::function<Real(const RealVect a_pos, const Real a_time)>& a_bcFunction){
  m_bcFunction = a_bcFunction;
}

#include <CD_NamespaceFooter.H>
