/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzNeumannEBBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzNeumannEBBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzNeumannEBBCFactory.H>
#include <CD_MFHelmholtzNeumannEBBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzNeumannEBBCFactory::MFHelmholtzNeumannEBBCFactory(){
  m_order       = -1;
  m_weight      = -1;
  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzNeumannEBBCFactory::MFHelmholtzNeumannEBBCFactory(const int a_order, const int a_weight, const Real a_DphiDn){
  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setDphiDn(a_DphiDn);
}

MFHelmholtzNeumannEBBCFactory::MFHelmholtzNeumannEBBCFactory(const int a_order, const int a_weight, const std::function<Real(const RealVect& a_pos)>& a_DphiDn){
  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setDphiDn(a_DphiDn);
}

MFHelmholtzNeumannEBBCFactory::~MFHelmholtzNeumannEBBCFactory(){
}

void MFHelmholtzNeumannEBBCFactory::setDphiDn(const Real a_DphiDn){
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantDphiDn = a_DphiDn;
}

void MFHelmholtzNeumannEBBCFactory::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn){
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionDphiDn = a_DphiDn;
}

void MFHelmholtzNeumannEBBCFactory::setBxDphiDn(const Real a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);
  m_multByBco = false;
}

void MFHelmholtzNeumannEBBCFactory::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn){
  this->setDphiDn(a_BxDphiDn);
  m_multByBco = false;
}

RefCountedPtr<EBHelmholtzEBBC> MFHelmholtzNeumannEBBCFactory::create(const int a_iphase, const RefCountedPtr<JumpBC>& a_jumpBC) const {
  auto bc = new MFHelmholtzNeumannEBBC(a_iphase, a_jumpBC);

  bc->setOrder(m_order);
  bc->setWeight(m_weight);
  
  if(m_multByBco){
    if(m_useConstant){
      bc->setDphiDn(m_constantDphiDn);
    }
    else if(m_useFunction){
      bc->setDphiDn(m_functionDphiDn);
    }
    else{
      MayDay::Error("MFHelmholtzNeumannEBBCFactory::create() - logic bust. Not using constant or function");
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
      MayDay::Error("MFHelmholtzNeumannEBBCFactory::create() - logic bust");
    }
  }

  return RefCountedPtr<EBHelmholtzEBBC> (bc);
}

#include <CD_NamespaceFooter.H>
