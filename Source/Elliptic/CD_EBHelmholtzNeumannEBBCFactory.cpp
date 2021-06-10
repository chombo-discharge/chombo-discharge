/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannEBBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzNeumannEBBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzNeumannEBBC.H>
#include <CD_EBHelmholtzNeumannEBBCFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzNeumannEBBCFactory::EBHelmholtzNeumannEBBCFactory(){
  m_multByBco   = true;
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzNeumannEBBCFactory::EBHelmholtzNeumannEBBCFactory(const Real a_DphiDn){
  this->setDphiDn(a_DphiDn);
}

EBHelmholtzNeumannEBBCFactory::EBHelmholtzNeumannEBBCFactory(const std::function<Real(const RealVect& a_pos)>& a_DphiDn){
  this->setDphiDn(a_DphiDn);
}

EBHelmholtzNeumannEBBCFactory::~EBHelmholtzNeumannEBBCFactory(){

}

void EBHelmholtzNeumannEBBCFactory::setDphiDn(const int a_DphiDn){
  m_multByBco   = false;
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantDphiDn = a_DphiDn;
}

void EBHelmholtzNeumannEBBCFactory::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn){
  m_multByBco   = false;
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionDphiDn = a_DphiDn;
}

void EBHelmholtzNeumannEBBCFactory::setBxDphiDn(const int a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);
  m_multByBco = false;
}

void EBHelmholtzNeumannEBBCFactory::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);
  m_multByBco = false;
}

RefCountedPtr<EBHelmholtzEBBC> EBHelmholtzNeumannEBBCFactory::create() {
  
  auto bc = new EBHelmholtzNeumannEBBC();

  if(m_multByBco){
    if(m_useConstant) {
      bc->setDphiDn(m_constantDphiDn);
    }
    else if(m_useFunction){
      bc->setDphiDn(m_functionDphiDn);
    }
    else{
      MayDay::Error("EBHelmholtzEBBC::create - logic bust");
    }
  }
  else{
    if(m_useConstant) {
      bc->setBxDphiDn(m_constantDphiDn);
    }
    else if(m_useFunction){
      bc->setBxDphiDn(m_functionDphiDn);
    }
    else{
      MayDay::Error("EBHelmholtzEBBC::create - logic bust");
    }
  }

  return RefCountedPtr<EBHelmholtzEBBC>(bc);
}

#include <CD_NamespaceFooter.H>
