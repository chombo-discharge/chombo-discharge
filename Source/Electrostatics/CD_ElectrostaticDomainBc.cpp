/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ElectrostaticDomainBc.cpp
  @brief  Implementation of CD_ElectrostaticDomainBc.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ElectrostaticDomainBc.H>
#include <CD_NamespaceHeader.H>

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
  if(m_bcFunctions.find(a_wall) == m_bcFunctions.end()){
    MayDay::Abort("ElectrostaticDomainBc::getBc -- BC not found. Perhaps you've forgotten to set it...?");
  }

  return m_bcFunctions.at(a_wall);
}

#include <CD_NamespaceFooter.H>
