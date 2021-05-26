/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ConductivityEddingtonSP1DomainBcFactory.cpp
  @brief  Implementation of ConductivityEddingtonSP1DomainBcFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ConductivityEddingtonSP1DomainBcFactory.H>
#include <CD_NamespaceHeader.H>

ConductivityEddingtonSP1DomainBcFactory::ConductivityEddingtonSP1DomainBcFactory(const EddingtonSP1DomainBc& a_domainBc,
										 const LarsenMap&            a_larsenCoefficients,
										 const RealVect a_probLo){
  m_domainBc = a_domainBc;
  m_probLo   = a_probLo;
  m_larsen   = a_larsenCoefficients;
}

ConductivityEddingtonSP1DomainBcFactory::~ConductivityEddingtonSP1DomainBcFactory(){

}

ConductivityEddingtonSP1DomainBc* ConductivityEddingtonSP1DomainBcFactory::create(const ProblemDomain& a_domain, const EBISLayout& a_ebisl, const RealVect& a_dx) {
  return new ConductivityEddingtonSP1DomainBc(m_domainBc, m_larsen, m_probLo);
}

#include <CD_NamespaceFooter.H>
