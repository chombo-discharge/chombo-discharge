/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ElectrostaticEbBc.cpp
  @brief  Implementation of CD_ElectrostaticEbBc.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ElectrostaticEbBc.H>
#include <CD_NamespaceHeader.H>

ElectrostaticEbBc::ElectrostaticEbBc(){

}

ElectrostaticEbBc::~ElectrostaticEbBc(){

}

void ElectrostaticEbBc::clear(){
  m_bcFunctions.resize(0);
}

void ElectrostaticEbBc::addEbBc(const Electrode& a_electrode, const BcFunction& a_bcFunction){
  m_bcFunctions.emplace_back(std::make_pair(a_electrode, a_bcFunction));
}

void ElectrostaticEbBc::setEbBc(const int a_electrode, const BcFunction& a_bcFunction){
  if(a_electrode < 0 || a_electrode > m_bcFunctions.size()){
    MayDay::Abort("ElectrostaticEbBc::setEbBc -- index is out of range!");
  }
  
  m_bcFunctions[a_electrode].second = a_bcFunction;
}


ElectrostaticEbBc::BcFunction& ElectrostaticEbBc::getBc(const int a_electrode){
  if(a_electrode < 0 || a_electrode > m_bcFunctions.size()){
    MayDay::Abort("ElectrostaticEbBc::getBc -- index is out of range!");
  }

  return m_bcFunctions[a_electrode].second;
}

const ElectrostaticEbBc::BcFunction& ElectrostaticEbBc::getBc(const int a_electrode) const {
  if(a_electrode < 0 || a_electrode > m_bcFunctions.size()){
    MayDay::Abort("ElectrostaticEbBc::getBc -- index is out of range!");
  }

  return m_bcFunctions[a_electrode].second;
}

std::vector<std::pair<Electrode, ElectrostaticEbBc::BcFunction> >& ElectrostaticEbBc::getBcs(){
  return m_bcFunctions;
}

const std::vector<std::pair<Electrode, ElectrostaticEbBc::BcFunction> >&  ElectrostaticEbBc::getBcs() const {
  return m_bcFunctions;
}

#include <CD_NamespaceFooter.H>
