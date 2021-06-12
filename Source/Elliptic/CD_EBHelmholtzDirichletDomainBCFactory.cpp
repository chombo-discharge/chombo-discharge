/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletDomainBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzDirichletDomainBCFactory.H
  @author Robert Marskar
*/

// Our includes
//#include <CD_EBHelmholtzDirichletDomainBC.H>
#include <CD_EBHelmholtzDirichletDomainBCFactory.H>
#include <CD_NamespaceHeader.H>


EBHelmholtzDirichletDomainBCFactory::EBHelmholtzDirichletDomainBCFactory(){
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzDirichletDomainBCFactory::~EBHelmholtzDirichletDomainBCFactory(){

}

void EBHelmholtzDirichletDomainBCFactory::setValue(const Real a_value){
  m_useConstant = true;
  m_useFunction = false;

  m_constantValue = a_value;
}

void EBHelmholtzDirichletDomainBCFactory::setValue(const std::function<Real(const RealVect& a_pos)>& a_value){
  m_useConstant = false;
  m_useFunction = true;

  m_functionValue = a_value;
}

RefCountedPtr<EBHelmholtzDomainBC> EBHelmholtzDirichletDomainBCFactory::create() {

  //  auto bc = new EBHelmholtzDirichletDomainBC();
}

#include <CD_NamespaceFooter.H>
