/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletEBBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzDirichletEBBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzDirichletEBBC.H>
#include <CD_EBHelmholtzDirichletEBBCFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzDirichletEBBCFactory::EBHelmholtzDirichletEBBCFactory(){
  m_order       = -1;
  m_weight      = -1;
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzDirichletEBBCFactory::EBHelmholtzDirichletEBBCFactory(const int a_order, const int a_weight, const Real a_value){
  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setValue(a_value);
}


EBHelmholtzDirichletEBBCFactory::EBHelmholtzDirichletEBBCFactory(const int a_order, const int a_weight, const std::function<Real(const RealVect& a_pos)>& a_value){
  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setValue(a_value);
}

EBHelmholtzDirichletEBBCFactory::~EBHelmholtzDirichletEBBCFactory(){

}

void EBHelmholtzDirichletEBBCFactory::setOrder(const int a_order){
  m_order = a_order;
}

void EBHelmholtzDirichletEBBCFactory::setWeight(const int a_weight){
  m_weight = a_weight;
}

void EBHelmholtzDirichletEBBCFactory::setValue(const Real a_value){
  m_useConstant = true;
  m_useFunction = false;
  
  m_constantValue = a_value;
}

void EBHelmholtzDirichletEBBCFactory::setValue(const std::function<Real(const RealVect& a_pos)>& a_value){
  m_useConstant = false;
  m_useFunction = true;
  
  m_functionValue = a_value;
}

RefCountedPtr<EBHelmholtzEBBC> EBHelmholtzDirichletEBBCFactory::create() {

  auto bc = new EBHelmholtzDirichletEBBC();

  bc->setOrder(m_order);
  bc->setWeight(m_weight);
  if(m_useConstant){
    bc->setValue(m_constantValue);
  }
  else if(m_useFunction){
    bc->setValue(m_functionValue);
  }

  return RefCountedPtr<EBHelmholtzEBBC>(bc);
}

#include <CD_NamespaceFooter.H>
