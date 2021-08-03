/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzRobinDomainBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzRobinDomainBCFactory.cpp
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzRobinDomainBC.H>
#include <CD_MFHelmholtzRobinDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzRobinDomainBCFactory::MFHelmholtzRobinDomainBCFactory(){
  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzRobinDomainBCFactory::MFHelmholtzRobinDomainBCFactory(const Real a_A, const Real a_B, const Real a_C){
  this->setCoefficients(a_A, a_B, a_C);
}

MFHelmholtzRobinDomainBCFactory::MFHelmholtzRobinDomainBCFactory(const std::function<Real(const RealVect& a_pos) >& a_A,
								 const std::function<Real(const RealVect& a_pos) >& a_B,
								 const std::function<Real(const RealVect& a_pos) >& a_C){
  this->setCoefficients(a_A, a_B, a_C);
}

MFHelmholtzRobinDomainBCFactory::~MFHelmholtzRobinDomainBCFactory(){

}

void MFHelmholtzRobinDomainBCFactory::setCoefficients(const Real a_A, const Real a_B, const Real a_C){
  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void MFHelmholtzRobinDomainBCFactory::setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
						      const std::function<Real(const RealVect& a_pos) >& a_B,
						      const std::function<Real(const RealVect& a_pos) >& a_C){
  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

RefCountedPtr<EBHelmholtzDomainBC> MFHelmholtzRobinDomainBCFactory::create(const int a_iphase) const {
  if(!(m_useConstant || m_useFunction)) MayDay::Abort("MFHelmholtzRobinDomaniBCFactory::create -- not using constant or function. Did you forget to set coefficients?");
  
  EBHelmholtzRobinDomainBC* bc = new EBHelmholtzRobinDomainBC();

  if(m_useConstant){
    bc->setCoefficients(m_constantA, m_constantB, m_constantC);
  }
  else if(m_useFunction){
    bc->setCoefficients(m_functionA, m_functionB, m_functionC);
  }
  else {
    MayDay::Error("MFHelmholtzRobinDomainBCFactory::create - logic bust, you must set the Robin coefficients!!");
  }

  return RefCountedPtr<EBHelmholtzDomainBC> (bc);
}


#include <CD_NamespaceFooter.H>
