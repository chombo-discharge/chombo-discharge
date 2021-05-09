/*!
  @file   CD_ElectrostaticDomainBc.cpp
  @brief  Implementation of CD_ElectrostaticDomainBc.H
  @author Robert Marskar
  @date   May 2021
*/

#include "CD_ElectrostaticDomainBc.H"

#include "CD_NamespaceHeader.H"

ElectrostaticDomainBc::ElectrostaticDomainBc() {
  m_bcFunctions.clear();
}

ElectrostaticDomainBc::~ElectrostaticDomainBc() {
  m_bcFunctions.clear();
}


void ElectrostaticDomainBc::setBc(const Wall a_wall, const Bc a_func){
  m_bcFunctions.emplace(a_wall, a_func);
}

ElectrostaticDomainBc::Bc ElectrostaticDomainBc::getBc(const Wall a_wall) const{
  if(m_bcFunctions.find(a_wall) != m_bcFunctions.end()){
    MayDay::Abort("ElectrostaticDomainBc::getBc -- BC not found. Perhaps you've forgotten to set it...?");
  }

  return m_bcFunctions.at(a_wall);
}

#include "CD_NamespaceFooter.H"
