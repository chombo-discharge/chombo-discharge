/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzRobinEBBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzRobinEBBCFactory.cpp
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzRobinEBBC.H>
#include <CD_EBHelmholtzRobinEBBCFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzRobinEBBCFactory::EBHelmholtzRobinEBBCFactory(){
  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzRobinEBBCFactory::EBHelmholtzRobinEBBCFactory(const Real a_A, const Real a_B, const Real a_C){
  this->setCoefficients(a_A, a_B, a_C);
}

EBHelmholtzRobinEBBCFactory::EBHelmholtzRobinEBBCFactory(const std::function<Real(const RealVect& a_pos) >& a_A,
							 const std::function<Real(const RealVect& a_pos) >& a_B,
							 const std::function<Real(const RealVect& a_pos) >& a_C){
  this->setCoefficients(a_A, a_B, a_C);
}

EBHelmholtzRobinEBBCFactory::~EBHelmholtzRobinEBBCFactory(){

}

void EBHelmholtzRobinEBBCFactory::setCoefficients(const Real a_A, const Real a_B, const Real a_C){
  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void EBHelmholtzRobinEBBCFactory::setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
						  const std::function<Real(const RealVect& a_pos) >& a_B,
						  const std::function<Real(const RealVect& a_pos) >& a_C){
  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

RefCountedPtr<EBHelmholtzEBBC> EBHelmholtzRobinEBBCFactory::create() {
  EBHelmholtzRobinEBBC* bc = new EBHelmholtzRobinEBBC();

  if(m_useConstant){
    bc->setCoefficients(m_constantA, m_constantB, m_constantC);
  }
  else if(m_useFunction){
    bc->setCoefficients(m_functionA, m_functionB, m_functionC);
  }
  else {
    MayDay::Error("EBHelmholtzRobinEBBCFactory::create - logic bust, you must set the BC values!");
  }

  return RefCountedPtr<EBHelmholtzEBBC> (bc);
}


#include <CD_NamespaceFooter.H>
