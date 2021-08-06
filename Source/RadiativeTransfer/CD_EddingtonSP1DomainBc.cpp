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

void EddingtonSP1DomainBc::setBc(const Wall a_wall, const Bc a_func){
  m_bcFunctions.emplace(a_wall, a_func);
}

EddingtonSP1DomainBc::Bc& EddingtonSP1DomainBc::getBc(const Wall a_wall) {
  if(m_bcFunctions.find(a_wall) == m_bcFunctions.end()){
    MayDay::Abort("EddingtonSP1DomainBc::getBc -- BC not found. Perhaps you've forgotten to set it...?");
  }

  return m_bcFunctions.at(a_wall);
}

const EddingtonSP1DomainBc::Bc& EddingtonSP1DomainBc::getBc(const Wall a_wall) const{
  if(m_bcFunctions.find(a_wall) == m_bcFunctions.end()){
    MayDay::Abort("EddingtonSP1DomainBc::getBc -- BC not found. Perhaps you've forgotten to set it...?");
  }

  return m_bcFunctions.at(a_wall);
}

#include <CD_NamespaceFooter.H>
