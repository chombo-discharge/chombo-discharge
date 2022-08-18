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

EBHelmholtzRobinEBBCFactory::EBHelmholtzRobinEBBCFactory()
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::EBHelmholtzRobinEBBCFactory()");

  m_order           = -1;
  m_weight          = -1;
  m_domainDropOrder = -1;
  m_useConstant     = false;
  m_useFunction     = false;
}

EBHelmholtzRobinEBBCFactory::EBHelmholtzRobinEBBCFactory(const int  a_order,
                                                         const int  a_weight,
                                                         const Real a_A,
                                                         const Real a_B,
                                                         const Real a_C)
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::EBHelmholtzRobinEBBCFactory(int, int, Real, Real, Real)");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);

  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setCoefficients(a_A, a_B, a_C);
}

EBHelmholtzRobinEBBCFactory::EBHelmholtzRobinEBBCFactory(const int                                         a_order,
                                                         const int                                         a_weight,
                                                         const std::function<Real(const RealVect& a_pos)>& a_A,
                                                         const std::function<Real(const RealVect& a_pos)>& a_B,
                                                         const std::function<Real(const RealVect& a_pos)>& a_C)
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::EBHelmholtzRobinEBBCFactory(int, int, 3xstd::function<Real(RealVect)>)");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);

  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setCoefficients(a_A, a_B, a_C);
}

EBHelmholtzRobinEBBCFactory::~EBHelmholtzRobinEBBCFactory()
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::~EBHelmholtzRobinEBBCFactory()");
}

void
EBHelmholtzRobinEBBCFactory::setOrder(const int a_order)
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::setOrder(int)");

  CH_assert(a_order > 0);

  m_order = a_order;
}

void
EBHelmholtzRobinEBBCFactory::setWeight(const int a_weight)
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::setWeight(int)");

  CH_assert(a_weight > 0);

  m_weight = a_weight;
}

void
EBHelmholtzRobinEBBCFactory::setDomainDropOrder(const int a_domainSize)
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::setDomainDropOrder()");

  m_domainDropOrder = a_domainSize;
}

void
EBHelmholtzRobinEBBCFactory::setCoefficients(const Real a_A, const Real a_B, const Real a_C)
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::setCoefficients(Real, Real, Real)");

  m_constantA = a_A;
  m_constantB = a_B;
  m_constantC = a_C;

  m_useConstant = true;
  m_useFunction = false;
}

void
EBHelmholtzRobinEBBCFactory::setCoefficients(const std::function<Real(const RealVect& a_pos)>& a_A,
                                             const std::function<Real(const RealVect& a_pos)>& a_B,
                                             const std::function<Real(const RealVect& a_pos)>& a_C)
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::setCoefficients(3xstd::function<Real(RealVect)>)");

  m_functionA = a_A;
  m_functionB = a_B;
  m_functionC = a_C;

  m_useConstant = false;
  m_useFunction = true;
}

RefCountedPtr<EBHelmholtzEBBC>
EBHelmholtzRobinEBBCFactory::create()
{
  CH_TIME("EBHelmholtzRobinEBBCFactory::create()");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);

  EBHelmholtzRobinEBBC* bc = new EBHelmholtzRobinEBBC();

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
    MayDay::Error("EBHelmholtzRobinEBBCFactory::create - logic bust, you must set the Robin coefficients!");
  }

  return RefCountedPtr<EBHelmholtzEBBC>(bc);
}

#include <CD_NamespaceFooter.H>
