/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzNeumannEBBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzNeumannEBBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzNeumannEBBCFactory.H>
#include <CD_MFHelmholtzNeumannEBBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzNeumannEBBCFactory::MFHelmholtzNeumannEBBCFactory()
{
  CH_TIME("MFHelmholtzNeumannEBBCFactory::MFHelmholtzNeumannEBBCFactory()");

  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzNeumannEBBCFactory::MFHelmholtzNeumannEBBCFactory(const Real a_DphiDn)
{
  CH_TIME("MFHelmholtzNeumannEBBCFactory::MFHelmholtzNeumannEBBCFactory(Real)");

  this->setDphiDn(a_DphiDn);
}

MFHelmholtzNeumannEBBCFactory::MFHelmholtzNeumannEBBCFactory(const std::function<Real(const RealVect& a_pos)>& a_DphiDn)
{
  CH_TIME("MFHelmholtzNeumannEBBCFactory::MFHelmholtzNeumannEBBCFactory(std::function<Real(RealVect)>)");

  this->setDphiDn(a_DphiDn);
}

MFHelmholtzNeumannEBBCFactory::~MFHelmholtzNeumannEBBCFactory()
{
  CH_TIME("MFHelmholtzNeumannEBBCFactory::~MFHelmholtzNeumannEBBCFactory()");
}

void
MFHelmholtzNeumannEBBCFactory::setDphiDn(const Real a_DphiDn)
{
  CH_TIME("MFHelmholtzNeumannEBBCFactory::setDphiDn(Real)");

  m_multByBco = true;

  m_useConstant = true;
  m_useFunction = false;

  m_constantDphiDn = a_DphiDn;
}

void
MFHelmholtzNeumannEBBCFactory::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn)
{
  CH_TIME("MFHelmholtzNeumannEBBCFactory::setDphiDn(std::function<Real(RealVect)>)");

  m_multByBco = true;

  m_useConstant = false;
  m_useFunction = true;

  m_functionDphiDn = a_DphiDn;
}

void
MFHelmholtzNeumannEBBCFactory::setBxDphiDn(const Real a_BxDphiDn)
{
  CH_TIME("MFHelmholtzNeumannEBBCFactory::setBxDphiDn(Real)");

  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

void
MFHelmholtzNeumannEBBCFactory::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn)
{
  CH_TIME("MFHelmholtzNeumannEBBCFactory::setBxDphiDn(std::function<Real(RealVect)>)");

  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

RefCountedPtr<EBHelmholtzEBBC>
MFHelmholtzNeumannEBBCFactory::create(const int a_iphase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC) const
{
  CH_TIME("MFHelmholtzNeumannEBBCFactory::create(int, RefCountedPtr<MFHelmholtzJumpBC>)");

  CH_assert(m_useFunction || m_useConstant);

  auto bc = new MFHelmholtzNeumannEBBC(a_iphase, a_jumpBC);

  if (m_multByBco) {
    if (m_useConstant) {
      bc->setDphiDn(m_constantDphiDn);
    }
    else if (m_useFunction) {
      bc->setDphiDn(m_functionDphiDn);
    }
    else {
      MayDay::Error("MFHelmholtzNeumannEBBCFactory::create() - logic bust. Not using constant or function");
    }
  }
  else {
    if (m_useConstant) {
      bc->setBxDphiDn(m_constantDphiDn);
    }
    else if (m_useFunction) {
      bc->setBxDphiDn(m_functionDphiDn);
    }
    else {
      MayDay::Error("MFHelmholtzNeumannEBBCFactory::create() - logic bust");
    }
  }

  return RefCountedPtr<EBHelmholtzEBBC>(bc);
}

#include <CD_NamespaceFooter.H>
