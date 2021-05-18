/*!
  @file   ConductivityDomainBcWrapperFactory.cpp
  @brief  Implementation of ConductivityDomainBcWrapperFactory
  @author Robert Marskar
  @date   June 2018
*/

#include <CD_ConductivityDomainBcWrapperFactory.H>

#include "CD_NamespaceHeader.H"

ConductivityDomainBcWrapperFactory::ConductivityDomainBcWrapperFactory(){
  m_hasbc = false;
}

ConductivityDomainBcWrapperFactory::~ConductivityDomainBcWrapperFactory(){

}

void ConductivityDomainBcWrapperFactory::setWallBc(const Vector<RefCountedPtr<WallBc> >& a_wallbc){
  m_wallbc = a_wallbc;
  m_hasbc = true;
}

void ConductivityDomainBcWrapperFactory::setPotentials(const Vector<RefCountedPtr<BaseBCFuncEval> >& a_potentials){
  m_potentials = a_potentials;
}

void ConductivityDomainBcWrapperFactory::setRobinCoefficients(const Vector<RefCountedPtr<RobinCoefficients> >& a_robinco){
  m_robinco = a_robinco;
}

ConductivityDomainBcWrapper* ConductivityDomainBcWrapperFactory::create(const ProblemDomain& a_domain,
									   const EBISLayout&    a_ebisl,
									   const RealVect&      a_dx){
  CH_assert(m_hasbc);
    
  ConductivityDomainBcWrapper* fresh = new ConductivityDomainBcWrapper();

  fresh->setPotentials(m_potentials);
  fresh->setRobinCoefficients(m_robinco);
  fresh->setWallBc(m_wallbc);

    
  return fresh;
}
#include "CD_NamespaceFooter.H"
