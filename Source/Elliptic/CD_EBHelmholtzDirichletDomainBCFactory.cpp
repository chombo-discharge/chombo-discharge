/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletDomainBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzDirichletDomainBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzDirichletDomainBC.H>
#include <CD_EBHelmholtzDirichletDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzDirichletDomainBCFactory::EBHelmholtzDirichletDomainBCFactory(){
  CH_TIME("EBHelmholtzDirichletDomainBCFactory::EBHelmholtzDirichletDomainBCFactory()");
  
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzDirichletDomainBCFactory::EBHelmholtzDirichletDomainBCFactory(const Real a_value){
  CH_TIME("EBHelmholtzDirichletDomainBCFactory::EBHelmholtzDirichletDomainBCFactory(Real)");
  
  this->setValue(a_value);
}

EBHelmholtzDirichletDomainBCFactory::EBHelmholtzDirichletDomainBCFactory(const std::function<Real(const RealVect& a_pos)>& a_value){
  CH_TIME("EBHelmholtzDirichletDomainBCFactory::EBHelmholtzDirichletDomainBCFactory(std::function<Real(const RealVect)>)");
  
  this->setValue(a_value);
}

EBHelmholtzDirichletDomainBCFactory::~EBHelmholtzDirichletDomainBCFactory(){
  CH_TIME("EBHelmholtzDirichletDomainBCFactory::~EBHelmholtzDirichletDomainBCFactory()");
}

void EBHelmholtzDirichletDomainBCFactory::setValue(const Real a_value){
  CH_TIME("EBHelmholtzDirichletDomainBCFactory::setValue(Real)");
  
  m_useConstant = true;
  m_useFunction = false;

  m_constantValue = a_value;
}

void EBHelmholtzDirichletDomainBCFactory::setValue(const std::function<Real(const RealVect& a_pos)>& a_value){
  CH_TIME("EBHelmholtzDirichletDomainBCFactory::setValue(std::function<Real(const RealVect)>)");
  
  m_useConstant = false;
  m_useFunction = true;

  m_functionValue = a_value;
}

RefCountedPtr<EBHelmholtzDomainBC> EBHelmholtzDirichletDomainBCFactory::create() const {
  CH_TIME("EBHelmholtzDirichletDomainBCFactory::create()");

  CH_assert(m_useConstant || m_useFunction);

  // Also issue a run-time error if user forgot to set the BC. 
  if(!(m_useConstant || m_useFunction)) MayDay::Error("EBHelmholtzDirichletDomainBCFactory::create - logic bust, not using function or constant!");
  
  auto bc = new EBHelmholtzDirichletDomainBC();

  if(m_useConstant) {
    bc->setValue(m_constantValue);
  }
  else if(m_useFunction){
    bc->setValue(m_functionValue);
  }

  return RefCountedPtr<EBHelmholtzDomainBC>(bc);
}

#include <CD_NamespaceFooter.H>
