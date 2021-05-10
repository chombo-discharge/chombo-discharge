/*!
  @file   CD_ConductivityElectrostaticDomainBcFactory.cpp
  @brief  Implementation of ConductivityElectrostaticDomainBcFactory.H
  @author Robert Marskar
  @date   May 2021
*/

#include "CD_ConductivityElectrostaticDomainBcFactory.H"
#include "CD_NamespaceHeader.H"

ConductivityElectrostaticDomainBcFactory::ConductivityElectrostaticDomainBcFactory(const WallBcFuncs& a_bcFunctions, const WallBcTypes& a_bcTypes){
  m_bcFunctions = a_bcFunctions;
  m_bcTypes     = a_bcTypes;
}

ConductivityElectrostaticDomainBcFactory::~ConductivityElectrostaticDomainBcFactory(){

}

ConductivityElectrostaticDomainBc* ConductivityElectrostaticDomainBcFactory::create(const ProblemDomain& a_domain, const EBISLayout& a_ebisl, const RealVect& a_dx) {

  return new ConductivityElectrostaticDomainBc(m_bcFunctions, m_bcTypes);
}

#include "CD_NamespaceFooter.H"
