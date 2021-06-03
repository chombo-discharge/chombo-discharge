/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */
/*!
  @file   CD_EBHelmholtzOpFactory.cpp
  @brief  Implementation of CD_EBHelmholtzOpFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzOpFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzOpFactory::EBHelmholtzOpFactory(const Real&              a_alpha,
					   const Real&              a_beta,
					   const AmrLevelGrids&     a_amrLevelGrids,
					   const AmrInterpolators&  a_amrInterpolators,
					   const AmrFluxRegisters&  a_amrFluxRegisters,
					   const AmrCoarseners&     a_amrCoarseners,
					   const AmrRefRatios&      a_amrRefRatios,
					   const AmrResolutions&    a_amrResolutions,					   
					   const AmrCellData&       a_amrAcoef,
					   const AmrFluxData&       a_amrBcoef,
					   const AmrIrreData&       a_amrBcoefIrreg,
					   const DomainBCFactory&   a_domainBCFactory,
					   const EBBCFactory&       a_ebbcFactory,
					   const IntVect&           a_ghostPhi,
					   const IntVect&           a_ghostRHS,
					   const RelaxationMethod&  a_relaxationMethod,
					   const ProblemDomain&     a_bottomDomain,
					   const AmrLevelGrids&     a_deeperLevelGrids){
  m_alpha = a_alpha;
  m_beta  = a_beta;  

  m_amrLevelGrids       = a_amrLevelGrids;
  m_amrInterpolators    = a_amrInterpolators;
  m_amrFluxRegisters    = a_amrFluxRegisters;
  m_amrCoarseners       = a_amrCoarseners;
  m_amrResolutions      = a_amrResolutions;
  m_amrRefRatios        = a_amrRefRatios;

  m_amrAcoef      = a_amrAcoef;
  m_amrBcoef      = a_amrBcoef;
  m_amrBcoefIrreg = a_amrBcoefIrreg;




}

  
EBHelmholtzOpFactory::~EBHelmholtzOpFactory(){

}

EBHelmholtzOp* EBHelmholtzOpFactory::MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly) {
  return NULL; // Not implemented (yet);
}

EBHelmholtzOp* EBHelmholtzOpFactory::AMRnewOp(const ProblemDomain& a_domain) {
  return NULL; // Not implemented (yet);
}

int EBHelmholtzOpFactory::refToFiner(const ProblemDomain& a_domain) const {
  int ref    = -1;
  bool found = false;
  
  for (int ilev = 0; ilev < m_amrLevelGrids.size(); ilev++){
    if(m_amrLevelGrids[ilev]->getDomain() == a_domain){
      ref   = m_amrRefRatios[ilev];
      found = true;
    }
  }

  // I will call this an error
  if(!found) MayDay::Abort("EBHelmholtzOpFactory::refToFiner - Domain not found in the AMR hierarchy");

  return ref;
}


#include <CD_NamespaceFooter.H>
