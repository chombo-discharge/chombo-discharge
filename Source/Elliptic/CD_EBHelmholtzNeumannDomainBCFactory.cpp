/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannDomainBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzNeumannDomainBCFactory.H
  @author Robert Marskar
*/

// Our includes
//#include <CD_EBHelmholtzNeumannDomainBC.H>
#include <CD_EBHelmholtzNeumannDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzNeumannDomainBCFactory::EBHelmholtzNeumannDomainBCFactory(){
  m_multByBco   = true;
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzNeumannDomainBCFactory::EBHelmholtzNeumannDomainBCFactory(const Real a_DphiDn){
  this->setDphiDn(a_DphiDn);
}

EBHelmholtzNeumannDomainBCFactory::EBHelmholtzNeumannDomainBCFactory(const std::function<Real(const RealVect& a_pos)>& a_DphiDn){
  this->setDphiDn(a_DphiDn);
}


EBHelmholtzNeumannDomainBCFactory::~EBHelmholtzNeumannDomainBCFactory(){

}

void EBHelmholtzNeumannDomainBCFactory::setDphiDn(const int a_DphiDn){
  m_multByBco   = false;
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantDphiDn = a_DphiDn;
}

void EBHelmholtzNeumannDomainBCFactory::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn){
  m_multByBco   = false;
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionDphiDn = a_DphiDn;
}

void EBHelmholtzNeumannDomainBCFactory::setBxDphiDn(const int a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);
  m_multByBco = false;
}

void EBHelmholtzNeumannDomainBCFactory::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);
  m_multByBco = false;
}

RefCountedPtr<EBHelmholtzDomainBC> EBHelmholtzNeumannDomainBCFactory::create() {

}

#include <CD_NamespaceFooter.H>
