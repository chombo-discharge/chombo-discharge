/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CdrSpecies.cpp
  @brief  Implementation of CD_CdrSpecies.H
  @author Robert Marskar
*/


// Our includes
#include <CD_CdrSpecies.H>
#include <CD_NamespaceHeader.H>
  
CdrSpecies::CdrSpecies(){
  m_name         = "default_CdrSpecies";
  m_unit         = "default_unit";
  m_chargeNumber       = 0;
  m_isDiffusive    = true;
  m_isMobile       = true;
  m_forceOutput = false;

  m_initializeWithFunction  = true;
  m_initializeWithParticles = true;
  
  m_deposition = DepositionType::NGP;

  m_initialParticles.clear();
}

CdrSpecies::CdrSpecies(const std::string a_name, const int a_chargeNumber, const bool a_mobile, const bool a_diffusive){
  m_name         = a_name;
  m_chargeNumber = a_chargeNumber;
  m_isMobile     = a_mobile;
  m_isDiffusive  = a_diffusive;

  m_initializeWithFunction  = true;
  m_initializeWithParticles = true;
  
  m_deposition = DepositionType::NGP;
  m_initialParticles.clear();
}

CdrSpecies::~CdrSpecies(){

}

Real CdrSpecies::initialData(const RealVect a_pos, const Real a_time) const{
  return 0.;
}

std::string CdrSpecies::getName() const {
  return m_name;
}

std::string CdrSpecies::get_unit() const {
  return m_unit;
}

int CdrSpecies::getChargeNumber() const {
  return m_chargeNumber;
}

bool CdrSpecies::isDiffusive() const {
  return m_isDiffusive;
}

bool CdrSpecies::isMobile() const {
  return m_isMobile;
}

bool CdrSpecies::forceOutput() const {
  return m_forceOutput;
}

bool CdrSpecies::initializeWithParticles() const {
  return m_initializeWithParticles;
}

bool CdrSpecies::initializeWithFunction() const{
  return m_initializeWithFunction;
}

DepositionType::Which CdrSpecies::getDeposition() {
  return m_deposition;
}

List<Particle>& CdrSpecies::getInitialParticles() {
  return m_initialParticles;
}

#include <CD_NamespaceFooter.H>
