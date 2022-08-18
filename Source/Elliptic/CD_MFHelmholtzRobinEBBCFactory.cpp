/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzRobinEBBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzRobinEBBCFactory.H
  @author Robert Marskar
*/

// Chombo inculdes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzRobinEBBCFactory.H>
#include <CD_MFHelmholtzRobinEBBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzRobinEBBCFactory::MFHelmholtzRobinEBBCFactory(const int  a_order,
                                                         const int  a_weight,
                                                         const Real a_A,
                                                         const Real a_B,
                                                         const Real a_C)
{
  CH_TIME("MFHelmholtzRobinEBBCFactory::MFHelmholtzRobinEBBCFactory(int, int, Real, Real, Real)");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);

  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setCoefficients(a_A, a_B, a_C);
  this->setDomainDropOrder(-1);
}

MFHelmholtzRobinEBBCFactory::MFHelmholtzRobinEBBCFactory(const int                                         a_order,
                                                         const int                                         a_weight,
                                                         const std::function<Real(const RealVect& a_pos)>& a_A,
                                                         const std::function<Real(const RealVect& a_pos)>& a_B,
                                                         const std::function<Real(const RealVect& a_pos)>& a_C)
{
  CH_TIME("MFHelmholtzRobinEBBCFactory::MFHelmholtzRobinEBBCFactory(int, int, Real, Real, Real)");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);

  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setCoefficients(a_A, a_B, a_C);
  this->setDomainDropOrder(-1);  
}

MFHelmholtzRobinEBBCFactory::~MFHelmholtzRobinEBBCFactory()
{
  CH_TIME("MFHelmholtzRobinEBBCFactory::~MFHelmholtzRobinEBBCFactory()");
}

void
MFHelmholtzRobinEBBCFactory::setOrder(const int a_order)
{
  CH_TIME("MFHelmholtzRobinEBBCFactory::setOrder(int)");

  CH_assert(a_order > 0);

  m_order = a_order;
}

void
MFHelmholtzRobinEBBCFactory::setWeight(const int a_weight)
{
  CH_TIME("MFHelmholtzRobinEBBCFactory::setWeight(int)");

  CH_assert(a_weight > 0);

  m_weight = a_weight;
}

void
MFHelmholtzRobinEBBC::setDomainDropOrder(const int a_domainSize) {
  CH_TIME("MFHelmholtzRobinEBBC::setDomainDropOrder()");

  m_domainDropOrder = a_domainSize;
}

void
MFHelmholtzRobinEBBCFactory::setCoefficients(const Real a_A, const Real a_B, const Real a_C)
{
  CH_TIME("MFHelmholtzRobinEBBCFactory::setCoefficients(Real, Real, Real)");

  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void
MFHelmholtzRobinEBBCFactory::setCoefficients(const std::function<Real(const RealVect& a_pos)>& a_A,
                                             const std::function<Real(const RealVect& a_pos)>& a_B,
                                             const std::function<Real(const RealVect& a_pos)>& a_C)
{
  CH_TIME("MFHelmholtzRobinEBBCFactory::setCoefficients(3xstd::function<Real(RealVect)>)");

  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

RefCountedPtr<EBHelmholtzEBBC>
MFHelmholtzRobinEBBCFactory::create(const int a_iphase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC) const
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::create()");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);

  auto bc = new MFHelmholtzRobinEBBC(a_iphase, a_jumpBC);

  bc->setOrder(m_order);
  bc->setWeight(m_weight);
  bc->setDomainDropOrder(m_domainDropOrder);
  
  if (m_useConstant) {
    bc->setCoefficients(m_constantA, m_constantB, m_constantC);
  }
  else if (m_useFunction) {
    bc->setCoefficients(m_functionA, m_functionB, m_functionC);
  }
  else {
    MayDay::Error("MFHelmholtzRobinEBBCFactory::create - logic bust, you must set the Robin coefficients!");
  }

  return RefCountedPtr<EBHelmholtzEBBC>(bc);
}

#include <CD_NamespaceFooter.H>
