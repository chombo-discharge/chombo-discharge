/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBHelmholtzOpFactory.cpp
  @brief  Implementation of CD_EBHelmholtzOpFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <BRMeshRefine.H>

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
					   const DomainBcFactory&   a_domainBcFactory,
					   const EbBcFactory&       a_ebBcFactory,
					   const IntVect&           a_ghostPhi,
					   const IntVect&           a_ghostRhs,
					   const RelaxationMethod&  a_relaxationMethod,
					   const ProblemDomain&     a_bottomDomain,
					   const int&               a_mgBlockingFactor,
					   const AmrLevelGrids&     a_deeperLevelGrids){

  // Define constructor arguments. 
  m_alpha               = a_alpha;
  m_beta                = a_beta;
  
  m_amrLevelGrids       = a_amrLevelGrids;
  m_amrInterpolators    = a_amrInterpolators;
  m_amrFluxRegisters    = a_amrFluxRegisters;
  m_amrCoarseners       = a_amrCoarseners;
  m_amrResolutions      = a_amrResolutions;
  m_amrRefRatios        = a_amrRefRatios;
  
  m_amrAcoef            = a_amrAcoef;
  m_amrBcoef            = a_amrBcoef;
  m_amrBcoefIrreg       = a_amrBcoefIrreg;
  
  m_domainBcFactory     = a_domainBcFactory;
  m_ebBcFactory         = a_ebBcFactory;
  
  m_ghostPhi            = a_ghostPhi;
  m_ghostRhs            = a_ghostRhs;
  
  m_relaxMethod         = a_relaxationMethod;
  m_bottomDomain        = a_bottomDomain;
  m_mgBlockingFactor    = a_mgBlockingFactor;
  m_deeperLevelGrids    = a_deeperLevelGrids;

  m_numAmrLevels = m_amrLevelGrids.size();

  // Asking multigrid to do the bottom solve at a refined AMR level is classified as bad input.
  if(this->isFiner(m_bottomDomain, m_amrLevelGrids[0]->getDomain())){
    MayDay::Abort("EBHelmholtzOpFactory -- bottomsolver domain can't be larger than the base AMR domain!");
  }
  
  this->defineMultigridLevels();
}
  
EBHelmholtzOpFactory::~EBHelmholtzOpFactory(){

}

void EBHelmholtzOpFactory::defineMultigridLevels(){
  // TLDR: This routine defines what is needed for making the multigrid levels. This includes the intermediate
  // levels (if you run with refinement factor 4) as well as the deeper multigrid levels that are coarsenings of the base AMR level. 
  
  m_mgLevelGrids.resize(m_numAmrLevels);
  m_mgAcoef.resize(m_numAmrLevels);
  m_mgBcoef.resize(m_numAmrLevels);
  m_mgBcoefIrreg.resize(m_numAmrLevels);
  m_hasMgLevels.resize(m_numAmrLevels);

  for (int amrLevel = 0; amrLevel < m_numAmrLevels; amrLevel++){
    if(amrLevel == 0 || m_amrRefRatios[amrLevel] > 2){
      m_hasMgLevels[amrLevel] = true;
      
      const int mgRefRatio = 2;

      m_mgLevelGrids[amrLevel].resize(0);
      m_mgAcoef[amrLevel].resize(0);
      m_mgBcoef[amrLevel].resize(0);
      m_mgBcoefIrreg[amrLevel].resize(0);

      m_mgLevelGrids[amrLevel].push_back(m_amrLevelGrids[amrLevel]);
      m_mgAcoef[amrLevel].push_back(m_amrAcoef[amrLevel]);
      m_mgBcoef[amrLevel].push_back(m_amrBcoef[amrLevel]);
      m_mgBcoefIrreg[amrLevel].push_back(m_amrBcoefIrreg[amrLevel]);


      bool hasCoarser = true;
      while(hasCoarser){
	
	// Current number of multigrid levels and the multigrid level which we will coarsen. Note the inverse order here (first entry is finest level)
	const int curMgLevels        =  m_mgLevelGrids[amrLevel].size();          
	const EBLevelGrid mgEblgFine = *m_mgLevelGrids[amrLevel][curMgLevels-1]; 

	// BoxLayout and domains for coarsening.
	EBLevelGrid mgEblgCoar;

	// This is an overriding option where we use the pre-defined coarsenings in m_deeperMultigridLevels. This is only valid for coarsenings of
	// the base AMR level, hence amrLevel == 0. Once those levels are exhausted we begin with direct coarsening. 
	if(amrLevel == 0 && curMgLevels < m_deeperLevelGrids.size()){    
	  hasCoarser = true;                                // Note that m_deeperLevelGrids[0] should be a factor 2 coarsening of the 
	  mgEblgCoar = *m_deeperLevelGrids[curMgLevels-1];  // coarsest AMR level. So curMgLevels-1 is correct. 
	}
	else{
	  // Let the operator factory do the coarsening this time. 
	  EBLevelGrid coarEblg;
	  hasCoarser = this->getCoarserLayout(mgEblgCoar, mgEblgFine, mgRefRatio, 16);
	}

	// Do not coarsen further if we end up with a domain smaller than m_bottomDomain. In this case
	// we will terminate the coarsening and let AMRMultiGrid do the bottom solve. 
	if(this->isCoarser(mgEblgCoar.getDomain(), m_bottomDomain)) hasCoarser = false;

	// Ok, we have a valid coarser domain, and it is given by mgCoarGrids and mgCoarDomain. 
	if(hasCoarser){

	}
      }
    }
    else{
      m_hasMgLevels[amrLevel] = false;
    }
  }
}

bool EBHelmholtzOpFactory::isCoarser(const ProblemDomain& A, const ProblemDomain& B) const{
  return A.domainBox().numPts() < B.domainBox().numPts();
}


bool EBHelmholtzOpFactory::isFiner(const ProblemDomain& A, const ProblemDomain& B) const{
  return A.domainBox().numPts() > B.domainBox().numPts();
}

bool EBHelmholtzOpFactory::getCoarserLayout(EBLevelGrid& a_coarEblg, const EBLevelGrid& a_fineEblg, const int a_refRat, const int a_blockingFactor) const {
  // TLDR: This creates a coarsening of a_fineGrid with refinement factor 2. We first try to split the grid using a_blockingFactor. If that does not work we
  //       coarsen directly. If that does not work we are out of ideas.

  bool hasCoarser;
  

  // This returns true if the fine grid fully covers the domain. The nature of this makes it
  // always true for the "deeper" multigridlevels,  but not so for the intermediate levels. 
  auto isFullyCovered = [&] (const EBLevelGrid& a_eblg) -> bool {
    unsigned long long numPtsLeft = a_eblg.getDomain().domainBox().numPts(); // Number of grid points in
    const DisjointBoxLayout& dbl  = a_eblg.getDBL();

    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit){
      numPtsLeft -= dbl[lit()].numPts();
    }

    return (numPtsLeft == 0ULL);
  };


  // Check if we can get a coarsenable domain. 
  if(a_fineEblg.coarsenable(a_refRat)){
    DisjointBoxLayout coarDbl;
    ProblemDomain     coarDomain = coarsen(a_fineEblg.getDomain(), a_refRat);

    bool doAggregation;
    bool doCoarsen;

    // Check if we are coarsening to an intermediate or "deep" level. If we are coarsening to a deep level we know that
    // the grids will cover the domain, so the aggregation check is just a matter of testing for valid domain decompositions. 
    if(this->isCoarser(coarDomain, m_amrLevelGrids[0]->getDomain())) {
      doAggregation = (refine(coarsen(coarDomain, a_blockingFactor), a_blockingFactor) == coarDomain);
    }
    else{ // The "simplest" way we can run with aggregation is if the fine grid covers the entire domain. There are probably generalizations of this. 
      doAggregation = isFullyCovered(a_fineEblg);
    }

    // Prefer aggregation over direct coarsening
    if(doAggregation){ 
      Vector<Box> boxes;
      Vector<int> procs;

      // We could have use our load balancing here, but I don't see why. 
      domainSplit(coarDomain, boxes, a_blockingFactor);
      mortonOrdering(boxes);
      LoadBalance(procs, boxes);

      coarDbl.define(boxes, procs, coarDomain);
      a_coarEblg.define(coarDbl, coarDomain, m_ghostPhi.max(), a_fineEblg.getEBIS());

      hasCoarser = true;
    }
    else if(doCoarsen){ // But use coarsening if we must.
      coarsen(coarDbl, a_fineEblg.getDBL(), a_refRat);
      a_coarEblg.define(coarDbl, coarDomain, m_ghostPhi.max(), a_fineEblg.getEBIS());

      hasCoarser = true;
    }
    else{ // Out of ideas. 
      hasCoarser = false;
    }
  }
  else{ // Nothing we can do. 
    hasCoarser = false;
  }

  return hasCoarser;
}

EBHelmholtzOp* EBHelmholtzOpFactory::MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly) {

  // Look through m_amrLevelGrids to see if we find a domain corresponding to a_fineDomain.
  // If we do, great. If we don't, then something has gone wrong and this function should never have been called.
  int amrLevel = -1;
  for (int ilev = 0; ilev < m_numAmrLevels; ilev++){
    if(a_fineDomain == m_amrLevelGrids[ilev]->getDomain()){
      amrLevel= ilev;
    }
  }

  if(amrLevel < 0) MayDay::Abort("EBHelmholtzOpFactory::MGnewOp - logic bust in MGnewOp, no corresponding amr level found");
  
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
