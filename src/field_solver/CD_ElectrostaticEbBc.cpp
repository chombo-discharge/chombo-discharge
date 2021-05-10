/*!
  @file   CD_ElectrostaticEbBc.cpp
  @brief  Implementation of CD_ElectrostaticEbBc.H
  @author Robert Marskar
  @date   May 2021
*/

#include "CD_ElectrostaticEbBc.H"
#include "CD_NamespaceHeader.H"

ElectrostaticEbBc::ElectrostaticEbBc(){

}

ElectrostaticEbBc::~ElectrostaticEbBc(){

}

void ElectrostaticEbBc::addEbBc(const std::shared_ptr<electrode>& a_electrode, const BcFunction& a_bcFunction){
  if(m_bcFunctions.find(a_electrode) == m_bcFunctions.end()){ // 
    this->setEbBc(a_electrode, a_bcFunction);
  }
  else{
    m_bcFunctions.emplace(a_electrode, a_bcFunction);
  }
}

void ElectrostaticEbBc::setEbBc(const std::shared_ptr<electrode>& a_electrode, const BcFunction& a_bcFunction){
  if(m_bcFunctions.find(a_electrode) != m_bcFunctions.end()){
    m_bcFunctions.at(a_electrode) = a_bcFunction;
  }
  else{
    MayDay::Abort("ElectrostaticEbBc::setEbBc -- electrode was not found! You need to add it first!");
  }
}


ElectrostaticEbBc::BcFunction& ElectrostaticEbBc::getBc(const std::shared_ptr<electrode>& a_electrode){
  if(m_bcFunctions.find(a_electrode) == m_bcFunctions.end()){
    MayDay::Abort("ElectrostaticEbBc::getBc -- electrode was not found! You need to add it first!");
  }

  return m_bcFunctions.at(a_electrode);
}

const ElectrostaticEbBc::BcFunction& ElectrostaticEbBc::getBc(const std::shared_ptr<electrode>& a_electrode) const {
  if(m_bcFunctions.find(a_electrode) == m_bcFunctions.end()){
    MayDay::Abort("ElectrostaticEbBc::getBc -- electrode was not found! You need to add it first!");
  }

  return m_bcFunctions.at(a_electrode);
}

std::map<std::shared_ptr<electrode>, ElectrostaticEbBc::BcFunction>& ElectrostaticEbBc::getBcs(){
  return m_bcFunctions;
}

const std::map<std::shared_ptr<electrode>, ElectrostaticEbBc::BcFunction>& ElectrostaticEbBc::getBcs() const {
  return m_bcFunctions;
}

#include "CD_NamespaceFooter.H"
