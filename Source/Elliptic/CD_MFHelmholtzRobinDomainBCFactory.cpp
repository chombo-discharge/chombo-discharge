/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzRobinDomainBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzRobinDomainBCFactory.cpp
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzRobinDomainBC.H>
#include <CD_MFHelmholtzRobinDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzRobinDomainBCFactory::MFHelmholtzRobinDomainBCFactory(){
  CH_TIME("MFHelmholtzRobinDomainBCFactory::MFHelmholtzRobinDomainBCFactory()");
  
  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzRobinDomainBCFactory::MFHelmholtzRobinDomainBCFactory(const Real a_A, const Real a_B, const Real a_C){
  CH_TIME("MFHelmholtzRobinDomainBCFactory::MFHelmholtzRobinDomainBCFactory(Real, Real, Real)");
  
  this->setCoefficients(a_A, a_B, a_C);
}

MFHelmholtzRobinDomainBCFactory::MFHelmholtzRobinDomainBCFactory(const std::function<Real(const RealVect& a_pos) >& a_A,
								 const std::function<Real(const RealVect& a_pos) >& a_B,
								 const std::function<Real(const RealVect& a_pos) >& a_C){
  CH_TIME("MFHelmholtzRobinDomainBCFactory::MFHelmholtzRobinDomainBCFactory(3x std::function<Real(RealVect)>)");
  
  this->setCoefficients(a_A, a_B, a_C);
}

MFHelmholtzRobinDomainBCFactory::~MFHelmholtzRobinDomainBCFactory(){
  CH_TIME("MFHelmholtzRobinDomainBCFactory::~MFHelmholtzRobinDomainBCFactory()");
}

void MFHelmholtzRobinDomainBCFactory::setCoefficients(const Real a_A, const Real a_B, const Real a_C){
  CH_TIME("MFHelmholtzRobinDomainBCFactory::setCoefficients(Real, Real, Real)");
  
  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void MFHelmholtzRobinDomainBCFactory::setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
						      const std::function<Real(const RealVect& a_pos) >& a_B,
						      const std::function<Real(const RealVect& a_pos) >& a_C){
  CH_TIME("MFHelmholtzRobinDomainBCFactory::setCoefficients(3x std::function<Real(RealVect)>)");
  
  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

RefCountedPtr<EBHelmholtzDomainBC> MFHelmholtzRobinDomainBCFactory::create(const int a_iphase) const {
  CH_TIME("MFHelmholtzRobinDomainBCFactory::create(int)");

  CH_assert(m_useFunction || m_useConstant);
  
  // Also issue run-time rror
  if(!(m_useConstant || m_useFunction)) {
    MayDay::Error("MFHelmholtzRobinDomaniBCFactory::create -- not using constant or function. Did you forget to set coefficients?");
  }
  
  EBHelmholtzRobinDomainBC* bc = new EBHelmholtzRobinDomainBC();

  if(m_useConstant){
    bc->setCoefficients(m_constantA, m_constantB, m_constantC);
  }
  else if(m_useFunction){
    bc->setCoefficients(m_functionA, m_functionB, m_functionC);
  }

  return RefCountedPtr<EBHelmholtzDomainBC> (bc);
}

#include <CD_NamespaceFooter.H>
