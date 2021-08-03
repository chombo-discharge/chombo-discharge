/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzDirichletEBBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzDirichletEBBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzDirichletEBBCFactory.H>
#include <CD_MFHelmholtzDirichletEBBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzDirichletEBBCFactory::MFHelmholtzDirichletEBBCFactory(){
  m_order       = -1;
  m_weight      = -1;
  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzDirichletEBBCFactory::MFHelmholtzDirichletEBBCFactory(const int a_order, const int a_weight, const Real a_value){
  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setValue(a_value);
}

MFHelmholtzDirichletEBBCFactory::MFHelmholtzDirichletEBBCFactory(const int a_order, const int a_weight, const std::function<Real(const RealVect& a_pos)>& a_value){
  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setValue(a_value);
}

MFHelmholtzDirichletEBBCFactory::~MFHelmholtzDirichletEBBCFactory(){
}

void MFHelmholtzDirichletEBBCFactory::setValue(const Real a_value){
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantValue = a_value;
}

void MFHelmholtzDirichletEBBCFactory::setValue(const std::function<Real(const RealVect& a_pos)>& a_value){
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionValue = a_value;
}

RefCountedPtr<EBHelmholtzEBBC> MFHelmholtzDirichletEBBCFactory::create(const int a_iphase, const RefCountedPtr<JumpBC>& a_jumpBC) const {
  auto bc = new MFHelmholtzDirichletEBBC(a_iphase, a_jumpBC);

  bc->setOrder(m_order);
  bc->setWeight(m_weight);
  if(m_useConstant){
    bc->setValue(m_constantValue);
  }
  else if(m_useFunction){
    bc->setValue(m_functionValue);
  }
  else{
    MayDay::Error("MFHelmholtzDirichletEBBCFactory::create() - logic bust. Not using constant or function");
  }

  return RefCountedPtr<EBHelmholtzEBBC> (bc);
}

#include <CD_NamespaceFooter.H>
