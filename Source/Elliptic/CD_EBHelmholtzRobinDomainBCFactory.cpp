/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzRobinDomainBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzRobinDomainBCFactory.cpp
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzRobinDomainBC.H>
#include <CD_EBHelmholtzRobinDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzRobinDomainBCFactory::EBHelmholtzRobinDomainBCFactory()
{
  CH_TIME("EBHelmholtzRobinDomainBCFactory::EBHelmholtzRobinDomainBCFactory()");

  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzRobinDomainBCFactory::EBHelmholtzRobinDomainBCFactory(const Real a_A, const Real a_B, const Real a_C)
{
  CH_TIME("EBHelmholtzRobinDomainBCFactory::EBHelmholtzRobinDomainBCFactory(Real, Real, Real)");

  this->setCoefficients(a_A, a_B, a_C);
}

EBHelmholtzRobinDomainBCFactory::EBHelmholtzRobinDomainBCFactory(const std::function<Real(const RealVect& a_pos)>& a_A,
                                                                 const std::function<Real(const RealVect& a_pos)>& a_B,
                                                                 const std::function<Real(const RealVect& a_pos)>& a_C)
{
  CH_TIME("EBHelmholtzRobinDomainBCFactory::EBHelmholtzRobinDomainBCFactory(3x std::function<Real(RealVect)>)");

  this->setCoefficients(a_A, a_B, a_C);
}

EBHelmholtzRobinDomainBCFactory::~EBHelmholtzRobinDomainBCFactory()
{}

void
EBHelmholtzRobinDomainBCFactory::setCoefficients(const Real a_A, const Real a_B, const Real a_C)
{
  CH_TIME("EBHelmholtzRobinDomainBCFactory::setCoefficients(Real, Real, Real)");

  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void
EBHelmholtzRobinDomainBCFactory::setCoefficients(const std::function<Real(const RealVect& a_pos)>& a_A,
                                                 const std::function<Real(const RealVect& a_pos)>& a_B,
                                                 const std::function<Real(const RealVect& a_pos)>& a_C)
{
  CH_TIME("EBHelmholtzRobinDomainBCFactory::setCoefficients(3x std::function<Real(RealVect)>)");

  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

RefCountedPtr<EBHelmholtzDomainBC>
EBHelmholtzRobinDomainBCFactory::create() const
{
  CH_TIME("EBHelmholtzRobinDomainBCFactory::create()");

  CH_assert(m_useConstant || m_useFunction);

  if (!(m_useConstant || m_useFunction)) {
    MayDay::Error(
      "EBHelmholtzRobinDomainBCFactory::create -- not using constant or function. Did you forget to set coefficients?");
  }

  EBHelmholtzRobinDomainBC* bc = new EBHelmholtzRobinDomainBC();

  if (m_useConstant) {
    bc->setCoefficients(m_constantA, m_constantB, m_constantC);
  }
  else if (m_useFunction) {
    bc->setCoefficients(m_functionA, m_functionB, m_functionC);
  }
  else {
    MayDay::Error("EBHelmholtzRobinDomainBCFactory::create - logic bust, you must set the Robin coefficients!!");
  }

  return RefCountedPtr<EBHelmholtzDomainBC>(bc);
}

#include <CD_NamespaceFooter.H>
