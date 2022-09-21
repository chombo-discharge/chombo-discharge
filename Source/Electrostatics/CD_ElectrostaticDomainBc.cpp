/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ElectrostaticDomainBc.cpp
  @brief  Implementation of CD_ElectrostaticDomainBc.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_ElectrostaticDomainBc.H>
#include <CD_NamespaceHeader.H>

ElectrostaticDomainBc::ElectrostaticDomainBc()
{
  CH_TIME("ElectrostaticDomainBc::ElectrostaticDomainBc()");

  m_bcFunctions.clear();

  // Make a lambda which returns zero everywhere.
  auto zero = [](const RealVect a_pos, const Real a_time) {
    return 0.0;
  };

  // Set the boundary condition functions on each side to be homogeneous Neumann boundary conditions.
  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const DomainSide domainSide = std::make_pair(dir, sit());

      m_bcFunctions.emplace(domainSide, std::make_pair(BcType::Neumann, zero));
    }
  }
}

ElectrostaticDomainBc::~ElectrostaticDomainBc()
{
  CH_TIME("ElectrostaticDomainBc::~ElectrostaticDomainBc()");

  m_bcFunctions.clear();
}

void
ElectrostaticDomainBc::setBc(const DomainSide a_domainSide, const Bc a_bc)
{
  CH_TIME("ElectrostaticDomainBc::setBc(DomainSide, Bc)");

  m_bcFunctions.at(a_domainSide) = a_bc;
}

ElectrostaticDomainBc::Bc
ElectrostaticDomainBc::getBc(const DomainSide a_domainSide) const
{
  CH_TIME("ElectrostaticDomainBc::getBc(DomainSide)");

  return m_bcFunctions.at(a_domainSide);
}

#include <CD_NamespaceFooter.H>
