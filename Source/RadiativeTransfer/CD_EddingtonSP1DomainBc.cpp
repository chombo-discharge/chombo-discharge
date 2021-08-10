/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EddingtonSP1DomainBc.cpp
  @brief  Implementation of CD_EddingtonSP1DomainBc.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EddingtonSP1DomainBc.H>
#include <CD_NamespaceHeader.H>

EddingtonSP1DomainBc::EddingtonSP1DomainBc() {
  m_bcFunctions.clear();
}

EddingtonSP1DomainBc::~EddingtonSP1DomainBc() {
  m_bcFunctions.clear();
}

void EddingtonSP1DomainBc::setBc(const DomainSide a_domainSide, const Bc a_func){
  m_bcFunctions.emplace(a_domainSide, a_func);
}

EddingtonSP1DomainBc::Bc& EddingtonSP1DomainBc::getBc(const DomainSide a_domainSide) {
  if(m_bcFunctions.find(a_domainSide) == m_bcFunctions.end()){
    MayDay::Abort("EddingtonSP1DomainBc::getBc -- BC not found. Perhaps you've forgotten to set it...?");
  }

  return m_bcFunctions.at(a_domainSide);
}

const EddingtonSP1DomainBc::Bc& EddingtonSP1DomainBc::getBc(const DomainSide a_domainSide) const{
  if(m_bcFunctions.find(a_domainSide) == m_bcFunctions.end()){
    MayDay::Abort("EddingtonSP1DomainBc::getBc -- BC not found. Perhaps you've forgotten to set it...?");
  }

  return m_bcFunctions.at(a_domainSide);
}

#include <CD_NamespaceFooter.H>
