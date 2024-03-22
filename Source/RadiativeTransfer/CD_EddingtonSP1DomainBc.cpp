/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EddingtonSP1DomainBc.cpp
  @brief  Implementation of CD_EddingtonSP1DomainBc.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EddingtonSP1DomainBc.H>
#include <CD_NamespaceHeader.H>

EddingtonSP1DomainBc::EddingtonSP1DomainBc()
{
  m_bcFunctions.clear();

  auto zero = [](const RealVect a_pos, const Real a_time) {
    return 0.0;
  };

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const DomainSide domainSide = std::make_pair(dir, sit());

      m_bcFunctions.emplace(domainSide, std::make_pair(BcType::Neumann, zero));
    }
  }
}

EddingtonSP1DomainBc::~EddingtonSP1DomainBc()
{
  m_bcFunctions.clear();
}

void
EddingtonSP1DomainBc::setBc(const DomainSide a_domainSide, const Bc a_bc)
{
  m_bcFunctions.at(a_domainSide) = a_bc;
}

EddingtonSP1DomainBc::Bc&
EddingtonSP1DomainBc::getBc(const DomainSide a_domainSide)
{
  return m_bcFunctions.at(a_domainSide);
}

const EddingtonSP1DomainBc::Bc&
EddingtonSP1DomainBc::getBc(const DomainSide a_domainSide) const
{
  return m_bcFunctions.at(a_domainSide);
}

#include <CD_NamespaceFooter.H>
