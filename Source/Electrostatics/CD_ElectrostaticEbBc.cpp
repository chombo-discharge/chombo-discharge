/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ElectrostaticEbBc.cpp
  @brief  Implementation of CD_ElectrostaticEbBc.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ElectrostaticEbBc.H>
#include <CD_NamespaceHeader.H>

ElectrostaticEbBc::ElectrostaticEbBc()
{
  CH_TIME("ElectrostaticEbBc::ElectrostaticEbBc()");
}

ElectrostaticEbBc::~ElectrostaticEbBc()
{
  CH_TIME("ElectrostaticEbBc::~ElectrostaticEbBc()");
}

void
ElectrostaticEbBc::clear()
{
  CH_TIME("ElectrostaticEbBc::clear()");

  m_bcFunctions.resize(0);
}

void
ElectrostaticEbBc::addEbBc(const Electrode& a_electrode, const BcFunction& a_bcFunction)
{
  CH_TIME("ElectrostaticEbBc::addEbBc(Electrode, BcFunction)");

  m_bcFunctions.emplace_back(a_electrode, a_bcFunction);
}

void
ElectrostaticEbBc::setEbBc(const int a_electrode, const BcFunction& a_bcFunction)
{
  CH_TIME("ElectrostaticEbBc::setEbBc(int, BcFunction)");

  if (a_electrode < 0 || a_electrode > m_bcFunctions.size()) {
    MayDay::Error("ElectrostaticEbBc::setEbBc -- index is out of range!");
  }

  m_bcFunctions[a_electrode].second = a_bcFunction;
}

ElectrostaticEbBc::BcFunction&
ElectrostaticEbBc::getBc(const int a_electrode)
{
  CH_TIME("ElectrostaticEbBc::getBc(int)");

  if (a_electrode < 0 || a_electrode > m_bcFunctions.size()) {
    MayDay::Error("ElectrostaticEbBc::getBc -- index is out of range!");
  }

  return m_bcFunctions[a_electrode].second;
}

const ElectrostaticEbBc::BcFunction&
ElectrostaticEbBc::getBc(const int a_electrode) const
{
  CH_TIME("ElectrostaticEbBc::getBc(int)");

  if (a_electrode < 0 || a_electrode > m_bcFunctions.size()) {
    MayDay::Error("ElectrostaticEbBc::getBc -- index is out of range!");
  }

  return m_bcFunctions[a_electrode].second;
}

std::vector<std::pair<Electrode, ElectrostaticEbBc::BcFunction>>&
ElectrostaticEbBc::getBcs()
{
  CH_TIME("ElectrostaticEbBc::getBcs()");

  return m_bcFunctions;
}

const std::vector<std::pair<Electrode, ElectrostaticEbBc::BcFunction>>&
ElectrostaticEbBc::getBcs() const
{
  CH_TIME("ElectrostaticEbBc::getBcs()");

  return m_bcFunctions;
}

#include <CD_NamespaceFooter.H>
