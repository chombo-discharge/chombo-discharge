/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOpFactory.cpp
  @brief  Implementation of CD_MFHelmholtzOpFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <BRMeshRefine.H>

// Our includes
#include <CD_MFHelmholtzOpFactory.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

constexpr int MFHelmholtzOpFactory::m_comp;
constexpr int MFHelmholtzOpFactory::m_nComp;

MFHelmholtzOpFactory::MFHelmholtzOpFactory(const Real&             a_alpha,
					   const Real&             a_beta,
					   const RealVect&         a_probLo,
					   const AmrLevelGrids&    a_amrLevelGrids,
					   const AmrInterpolators& a_amrInterpolators,
					   const AmrFluxRegisters& a_amrFluxRegisters,
					   const AmrCoarseners&    a_amrCoarseners,
					   const AmrRefRatios&     a_amrRefRatios,
					   const AmrResolutions&   a_amrResolutions,
					   const AmrCellData&      a_amrAcoef,
					   const AmrFluxData&      a_amrBcoef,
					   const AmrIrreData&      a_amrBcoefIrreg,
					   const DomainBCFactory&  a_domainBcFactory,
					   const IntVect&          a_ghostPhi,
					   const IntVect&          a_ghostRhs,
					   const RelaxType&        a_relaxationMethod,
					   const ProblemDomain&    a_bottomDomain,
					   const int&              a_ebbcOrder,
					   const int&              a_jumpOrder,
					   const int&              a_blockingFactor,
					   const AmrLevelGrids&    a_deeperLevelGrids){

  m_alpha  = a_alpha;
  m_beta   = a_beta;
  m_probLo = a_probLo;

  m_amrLevelGrids    = a_amrLevelGrids;
  m_amrInterpolators = a_amrInterpolators;
  m_amrFluxRegisters = a_amrFluxRegisters;
  m_amrCoarseners    = a_amrCoarseners;
  m_amrRefRatios     = a_amrRefRatios;  
  m_amrResolutions   = a_amrResolutions;

  m_amrAcoef      = a_amrAcoef;
  m_amrBcoef      = a_amrBcoef;
  m_amrBcoefIrreg = a_amrBcoefIrreg;

  m_domainBcFactory = a_domainBcFactory;

  m_ghostPhi = a_ghostPhi;
  m_ghostRhs = a_ghostRhs;

  m_relaxMethod  = a_relaxationMethod;
  m_bottomDomain = a_bottomDomain;

  m_ebbcOrder = a_ebbcOrder;
  m_jumpOrder = a_jumpOrder;

  m_mgBlockingFactor = a_blockingFactor;

  m_deeperLevelGrids = a_deeperLevelGrids;

  m_numAmrLevels = m_amrLevelGrids.size();

  // Asking multigrid to do the bottom solve at a refined AMR level is classified as bad input.
  if(this->isFiner(m_bottomDomain, m_amrLevelGrids[0].getDomain())){
    MayDay::Abort("EBHelmholtzOpFactory -- bottomsolver domain can't be larger than the base AMR domain!");
  }

  // Define the jump data and the multigrid levels. 
  this->defineJump();
  this->defineMultigridLevels();
}

MFHelmholtzOpFactory::~MFHelmholtzOpFactory(){

}

void MFHelmholtzOpFactory::setJump(const EBAMRIVData& a_sigma, const Real& a_scale){
  MayDay::Warning("MFHelmholtzOp::setJump -- can't do data-based jump just yet");
}

void MFHelmholtzOpFactory::setJump(const Real& a_sigma, const Real& a_scale){
  DataOps::setValue(m_amrJump, a_sigma*a_scale);
  for (int i = 0; i < m_mgJump.size(); i++){
    DataOps::setValue(m_mgJump[i], a_sigma*a_scale);
  }
}

void MFHelmholtzOpFactory::defineJump(){
  // TLDR: This defines m_amrJump on the first phase (gas phase). This is irregular data intended to be interfaced into the 
  //       boundary condition class. When we match the BC we get the data from here (the gas phase). Note that we define
  //       m_amrJump on all irregular cells, but the operators will do matching on a subset of them. 
  
  m_amrJump.resize(1 + m_numAmrLevels);

  constexpr int mainPhase = 0;

  for (int lvl = 0; lvl < m_numAmrLevels; lvl++){
    const DisjointBoxLayout& dbl = m_amrLevelGrids[lvl].getGrids();
    const EBLevelGrid& eblg      = m_amrLevelGrids[lvl].getEBLevelGrid(mainPhase);
    const EBISLayout& ebisl      = eblg.getEBISL();

    LayoutData<IntVectSet> irregCells(dbl);
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box box          = dbl[dit()];
      const EBISBox& ebisbox = ebisl[dit()];

      irregCells[dit()] = ebisbox.getIrregIVS(box);
    }

    BaseIVFactory<Real> fact(ebisl, irregCells);
    m_amrJump[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real> > > (new LevelData<BaseIVFAB<Real> >(dbl, m_comp, IntVect::Zero, fact));
  }
}

void MFHelmholtzOpFactory::defineMultigridLevels(){
  MayDay::Warning("MFHelmholtzOpFactory::defineMultigridLevels -- not implemented");

  m_mgLevelGrids.resize(m_numAmrLevels);
  m_mgAcoef.resize(m_numAmrLevels);
  m_mgBcoef.resize(m_numAmrLevels);
  m_mgBcoefIrreg.resize(m_numAmrLevels);
  m_mgJump.resize(m_numAmrLevels);
  m_hasMgLevels.resize(m_numAmrLevels);

  for (int amrLevel = 0; amrLevel < m_numAmrLevels; amrLevel++){
    m_hasMgLevels[amrLevel] = false;

    if(amrLevel == 0 && this->isCoarser(m_bottomDomain, m_amrLevelGrids[amrLevel].getDomain())){
      m_hasMgLevels[amrLevel] = true;
    }

    if(amrLevel > 0){
      if(m_amrRefRatios[amrLevel-1] > 2){
	m_hasMgLevels[amrLevel] = true;
      }
    }

    if(m_hasMgLevels[amrLevel]){
      
      constexpr int mgRefRatio = 2;

      m_mgLevelGrids[amrLevel].resize(0);
      m_mgAcoef[amrLevel].resize(0);
      m_mgBcoef[amrLevel].resize(0);
      m_mgBcoefIrreg[amrLevel].resize(0);
      m_mgJump[amrLevel].resize(0);

      m_mgLevelGrids[amrLevel].push_back(m_amrLevelGrids[amrLevel]);
      m_mgAcoef[amrLevel].push_back(m_amrAcoef[amrLevel]);
      m_mgBcoef[amrLevel].push_back(m_amrBcoef[amrLevel]);
      m_mgBcoefIrreg[amrLevel].push_back(m_amrBcoefIrreg[amrLevel]);
      m_mgJump[amrLevel].push_back(m_amrJump[amrLevel]);

      bool hasCoarser = true;

      while(hasCoarser){

	const int curMgLevels         = m_mgLevelGrids[amrLevel].size();
	const MFLevelGrid& mgMflgFine = m_mgLevelGrids[amrLevel].back();

	// This is the one we will define
	MFLevelGrid mgMflgCoar;

	// This is an overriding option where we use the pre-defined coarsenings in m_deeperMultigridLevels. This is only valid for coarsenings of
	// the base AMR level, hence amrLevel == 0. Once those levels are exhausted we begin with direct coarsening. 
	if(amrLevel == 0 && curMgLevels < m_deeperLevelGrids.size()){    
	  hasCoarser = true;                                // Note that m_deeperLevelGrids[0] should be a factor 2 coarsening of the 
	  mgMflgCoar = m_deeperLevelGrids[curMgLevels-1];  // coarsest AMR level. So curMgLevels-1 is correct.
	}
	else{
	  // Let the operator factory do the coarsening this time. 
	  hasCoarser = this->getCoarserLayout(mgMflgCoar, mgMflgFine, mgRefRatio, m_mgBlockingFactor);
	}

	// Ok, we can coarsen the domain define by mgMflgCoar. Set up that domain as a multigrid level. 
	if(hasCoarser){
	  
	}

	

	hasCoarser = false;
      }
    }
  }
}

bool MFHelmholtzOpFactory::getCoarserLayout(MFLevelGrid& a_coarMflg, const MFLevelGrid& a_fineMflg, const int a_refRat, const int a_blockingFactor) const {
  MayDay::Abort("MFHelmholtzOpFactory::getCoarserLayout -- not implemented");

  // This returns true if the fine grid fully covers the domain. The nature of this makes it
  // always true for the "deeper" multigridlevels,  but not so for the intermediate levels. 
  auto isFullyCovered = [&] (const MFLevelGrid& a_mflg) -> bool {
    unsigned long long numPtsLeft = a_mflg.getDomain().domainBox().numPts(); // Number of grid points in
    const DisjointBoxLayout& dbl  = a_mflg.getGrids();

    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit){
      numPtsLeft -= dbl[lit()].numPts();
    }

    return (numPtsLeft == 0);
  };

  const ProblemDomain fineDomain   = a_fineMflg.getDomain();
  const ProblemDomain coarDomain   = coarsen(fineDomain, a_refRat);
  const DisjointBoxLayout& fineDbl = a_fineMflg.getGrids();
  DisjointBoxLayout coarDbl;

  // Check if we can get a coarsenable domain. Don't want to coarsen to 1x1 so hence the factor of 2 in the test here. 
  ProblemDomain test = fineDomain;
  if(refine(coarsen(test, 2*a_refRat), 2*a_refRat) == fineDomain){
    const bool doCoarsen = a_fineMflg.getGrids().coarsenable(2*a_refRat);

    // Use coarsening if we can
    if(doCoarsen){

    }
    else{ // Check if we can use box aggregation 
      if(isFullyCovered(a_fineMflg)){
	Vector<Box> boxes;
	Vector<int> procs;

	// We could have use our load balancing here, but I don't see why. 
	domainSplit(coarDomain, boxes, a_blockingFactor);
	mortonOrdering(boxes);
	LoadBalance(procs, boxes);

	coarDbl.define(boxes, procs, coarDomain);
      }
    }
  }
}

MFHelmholtzOp* MFHelmholtzOpFactory::MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly) {
  return nullptr; // Multigrid turned off, for now. 
}

MFHelmholtzOp* MFHelmholtzOpFactory::AMRnewOp(const ProblemDomain& a_domain) {
  return new MFHelmholtzOp();
}

bool MFHelmholtzOpFactory::isCoarser(const ProblemDomain& A, const ProblemDomain& B) const{
  return A.domainBox().numPts() < B.domainBox().numPts();
}

bool MFHelmholtzOpFactory::isFiner(const ProblemDomain& A, const ProblemDomain& B) const{
  return A.domainBox().numPts() > B.domainBox().numPts();
}

int MFHelmholtzOpFactory::refToFiner(const ProblemDomain& a_domain) const {
  int ref    = -1;
  bool found = false;

  for (int ilev = 0; ilev < m_amrLevelGrids.size(); ilev++){
    if(m_amrLevelGrids[ilev].getDomain() == a_domain){
      found = true;
      ref   = m_amrRefRatios[ilev];
    }
  }

  if(!found) MayDay::Abort("EBHelmholtzOpFactory::refToFiner - domain not found in AMR hiearchy.");
}

int MFHelmholtzOpFactory::findAmrLevel(const ProblemDomain& a_domain) const{
  int amrLevel = -1;
  for (int lvl = 0; lvl < m_amrLevelGrids.size(); lvl++){
    if(m_amrLevelGrids[lvl].getDomain() == a_domain){
      amrLevel = lvl;
      break;
    }
  }

  if(amrLevel < 0) MayDay::Abort("EBHelmholtzOpFactory::findAmrLevel - no corresponding amr level found!");

  return amrLevel;
}

#include <CD_NamespaceFooter.H>

