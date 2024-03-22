/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrDomainBCImplem.H
  @brief  Implementation of CD_CdrDomainBC.H
  @author Robert Marskar
*/

#ifndef CD_CdrDomainBCImplem_H
#define CD_CdrDomainBCImplem_H

// Our includes
#include <CD_CdrDomainBC.H>
#include <CD_NamespaceHeader.H>

CdrDomainBC::CdrDomainBC()
{
  // TLDR: The default constructor populates the m_bcTypes and m_bcFunctions with wall functions and types. Users
  //       can later override these.

  // Flux function which returns zero everywhere.
  auto wallBc = [](const RealVect a_position, const Real a_time) -> Real {
    return 0.0;
  };

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const DomainSide domainSide = this->domainSide(dir, sit());

      m_bcTypes.emplace(domainSide, BcType::Wall);
      m_bcFunctions.emplace(domainSide, wallBc);
    }
  }
}

CdrDomainBC::~CdrDomainBC()
{}

CdrDomainBC::DomainSide
CdrDomainBC::domainSide(const int a_dir, const Side::LoHiSide a_side) const
{
  return std::make_pair(a_dir, a_side);
}

const CdrDomainBC::BcType&
CdrDomainBC::getBcType(const DomainSide& a_domainSide) const
{
  return m_bcTypes.at(a_domainSide);
}

const CdrDomainBC::FluxFunction&
CdrDomainBC::getBcFunction(const DomainSide& a_domainSide) const
{
  return m_bcFunctions.at(a_domainSide);
}

void
CdrDomainBC::setBcType(const DomainSide& a_domainSide, const BcType& a_bcType)
{
  m_bcTypes.at(a_domainSide) = a_bcType;
}

void
CdrDomainBC::setBcFunction(const DomainSide& a_domainSide, const FluxFunction& a_function)
{
  m_bcFunctions.at(a_domainSide) = a_function;
}

#include <CD_NamespaceFooter.H>

#endif
