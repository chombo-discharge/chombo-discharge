/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzRobinEBBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzRobinEBBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzRobinEBBCFactory.H>
#include <CD_MFHelmholtzRobinEBBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzRobinEBBCFactory::MFHelmholtzRobinEBBCFactory(const Real a_A,
							 const Real a_B,
							 const Real a_C){
  this->setCoefficients(a_A, a_B, a_C);
}

MFHelmholtzRobinEBBCFactory::MFHelmholtzRobinEBBCFactory(const std::function<Real(const RealVect& a_pos)>& a_A,
							 const std::function<Real(const RealVect& a_pos)>& a_B,
							 const std::function<Real(const RealVect& a_pos)>& a_C){
  this->setCoefficients(a_A, a_B, a_C);
}

MFHelmholtzRobinEBBCFactory::~MFHelmholtzRobinEBBCFactory(){
}

void MFHelmholtzRobinEBBCFactory::setCoefficients(const Real a_A, const Real a_B, const Real a_C){
  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void MFHelmholtzRobinEBBCFactory::setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
						  const std::function<Real(const RealVect& a_pos) >& a_B,
						  const std::function<Real(const RealVect& a_pos) >& a_C){
  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

RefCountedPtr<EBHelmholtzEBBC> MFHelmholtzRobinEBBCFactory::create(const int a_iphase, const RefCountedPtr<JumpBC>& a_jumpBC) const {
  auto bc = new MFHelmholtzRobinEBBC(a_iphase, a_jumpBC);

  if(m_useConstant){
    bc->setCoefficients(m_constantA, m_constantB, m_constantC);
  }
  else if(m_useFunction){
    bc->setCoefficients(m_functionA, m_functionB, m_functionC);
  }
  else {
    MayDay::Error("MFHelmholtzRobinEBBCFactory::create - logic bust, you must set the Robin coefficients!");
  }

  return RefCountedPtr<EBHelmholtzEBBC> (bc);
}

#include <CD_NamespaceFooter.H>
