/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBHelmholtzOpFactory.cpp
  @brief  Implementation of CD_EBHelmholtzOpFactory.H
  @author Robert Marskar
*/

#define EBHelmholtzOpFactoryDebug 0
#if EBHelmholtzOpFactoryDebug
#include <EBAMRPoissonOp.H>
#endif

// Chombo includes
#include <ParmParse.H>
#include <BRMeshRefine.H>

// Our includes
#include <CD_EBHelmholtzOpFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzOpFactory::EBHelmholtzOpFactory(const Location::Cell    a_dataLocation,
					   const Real&             a_alpha,
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
					   const EBBCFactory&      a_ebbcFactory,
					   const IntVect&          a_ghostPhi,
					   const IntVect&          a_ghostRhs,
					   const Smoother&         a_smoother,
					   const ProblemDomain&    a_bottomDomain,
					   const int&              a_mgBlockingFactor,
					   const AmrLevelGrids&    a_deeperLevelGrids){
  // Define constructor arguments.
  m_dataLocation        = a_dataLocation;
  m_alpha               = a_alpha;
  m_beta                = a_beta;

  m_probLo              = a_probLo;
  
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
  m_ebBcFactory         = a_ebbcFactory;
  
  m_ghostPhi            = a_ghostPhi;
  m_ghostRhs            = a_ghostRhs;
  
  m_smoother            = a_smoother;
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
    m_hasMgLevels[amrLevel] = false;

    if(amrLevel == 0 && this->isCoarser(m_bottomDomain, m_amrLevelGrids[amrLevel]->getDomain())){
      m_hasMgLevels[amrLevel] = true;
    }

    if(amrLevel > 0){
      if(m_amrRefRatios[amrLevel-1] > 2){
	m_hasMgLevels[amrLevel] = true;
      }
    }

    if(m_hasMgLevels[amrLevel]){
      
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
	const int curMgLevels         =  m_mgLevelGrids[amrLevel].size();          
	const EBLevelGrid& mgEblgFine = *m_mgLevelGrids[amrLevel].back();

	// BoxLayout and domains for coarsening.
	RefCountedPtr<EBLevelGrid> mgEblgCoar = RefCountedPtr<EBLevelGrid>(new EBLevelGrid());
	// This is an overriding option where we use the pre-defined coarsenings in m_deeperMultigridLevels. This is only valid for coarsenings of
	// the base AMR level, hence amrLevel == 0. Once those levels are exhausted we begin with direct coarsening. 
	if(amrLevel == 0 && curMgLevels < m_deeperLevelGrids.size()){    
	  hasCoarser = true;                                // Note that m_deeperLevelGrids[0] should be a factor 2 coarsening of the 
	  mgEblgCoar = m_deeperLevelGrids[curMgLevels-1];  // coarsest AMR level. So curMgLevels-1 is correct.
	}
	else{
	  // Let the operator factory do the coarsening this time. 
	  hasCoarser = this->getCoarserLayout(*mgEblgCoar, mgEblgFine, mgRefRatio, m_mgBlockingFactor);
	}

	// Do not coarsen further if we end up with a domain smaller than m_bottomDomain. In this case
	// we will terminate the coarsening and let AMRMultiGrid do the bottom solve. 
	if(hasCoarser){
	  if(this->isCoarser(mgEblgCoar->getDomain(), m_bottomDomain)){
	    hasCoarser = false;
	  }
	  else{
	    // Not so sure about this one, will we ever be asked to make an coarsened MG level which is also an AMR level? If not, this code
	    // will reduce the coarsening efforts.
	    for (int iamr = 0; iamr < m_numAmrLevels; iamr++){
	      if(mgEblgCoar->getDomain() == m_amrLevelGrids[iamr]->getDomain()){
	    	hasCoarser = false;
	      }
	    }
	  }
	}

	// Ok, we have a valid coarser domain which is given by mgEblgCoar. Use that domain to make the coefficients. 
	if(hasCoarser){
	  const EBLevelGrid& eblgCoar        = *mgEblgCoar;
	  const EBLevelGrid& eblgFine        =  mgEblgFine;

	  const EBISLayout& ebislCoar        = eblgCoar.getEBISL();
	  const EBISLayout& ebislFine        = eblgFine.getEBISL();

	  const DisjointBoxLayout& gridsCoar = eblgCoar.getDBL();
	  const DisjointBoxLayout& gridsFine = eblgFine.getDBL();

	  const ProblemDomain& domainCoar    = eblgCoar.getDomain();
	  const ProblemDomain& domainFine    = eblgFine.getDomain();

	  // Make the irregular sets
	  const int nghost = 1;
	  LayoutData<IntVectSet> irregSets(gridsCoar);
	  for (DataIterator dit(gridsCoar); dit.ok(); ++dit){
	    Box bx = gridsCoar[dit()];
	    bx.grow(nghost);
	    bx &= domainCoar;

	    const EBISBox& ebisbox = ebislCoar[dit()];

	    irregSets[dit()] = ebisbox.getIrregIVS(bx);
	  }
	  
	  // Factories
	  EBCellFactory       cellFactory(mgEblgCoar->getEBISL());
	  EBFluxFactory       fluxFactory(mgEblgCoar->getEBISL());
	  BaseIVFactory<Real> irreFactory(mgEblgCoar->getEBISL(), irregSets);

	  // Define coarsened coefficients
	  RefCountedPtr<LevelData<EBCellFAB> >        coarAcoef     (new LevelData<EBCellFAB>       (gridsCoar, 1, nghost*IntVect::Unit, cellFactory));
	  RefCountedPtr<LevelData<EBFluxFAB> >        coarBcoef     (new LevelData<EBFluxFAB>       (gridsCoar, 1, nghost*IntVect::Unit, fluxFactory));
	  RefCountedPtr<LevelData<BaseIVFAB<Real> > > coarBcoefIrreg(new LevelData<BaseIVFAB<Real> >(gridsCoar, 1, nghost*IntVect::Unit, irreFactory));

	  // These are the fine coefficients which will be coarsened.
	  const LevelData<EBCellFAB>&        fineAcoef      = *m_mgAcoef[amrLevel].back();
	  const LevelData<EBFluxFAB>&        fineBcoef      = *m_mgBcoef[amrLevel].back();
	  const LevelData<BaseIVFAB<Real> >& fineBcoefIrreg = *m_mgBcoefIrreg[amrLevel].back();

	  // Coarsen the coefficients
	  this->coarsenCoefficients(*coarAcoef,
	  			    *coarBcoef,
	  			    *coarBcoefIrreg,
	  			    fineAcoef,
	  			    fineBcoef,
	  			    fineBcoefIrreg,
	  			    eblgCoar,
	  			    eblgFine,
	  			    mgRefRatio);


	  m_mgLevelGrids[amrLevel].push_back(mgEblgCoar);
	  m_mgAcoef[amrLevel].push_back(coarAcoef);
	  m_mgBcoef[amrLevel].push_back(coarBcoef);
	  m_mgBcoefIrreg[amrLevel].push_back(coarBcoefIrreg);
	}
      }
    }

    if(m_mgLevelGrids[amrLevel].size() <= 1) {
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
  bool hasCoarser = false;

#if EBHelmholtzOpFactoryDebug
  DisjointBoxLayout dblCoFi;
  ProblemDomain domainCoFi;
  bool dummy;
  hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoFi,
						 domainCoFi,
						 a_fineEblg.getDBL(),
						 a_fineEblg.getEBISL(),
						 a_fineEblg.getDomain(),
						 a_refRat,
						 a_fineEblg.getEBIS(),
						 a_blockingFactor,
						 dummy,
						 2);

  if(hasCoarser){
    a_coarEblg.define(dblCoFi, domainCoFi, 4, a_fineEblg.getEBIS());
  }
#else
  
  // This returns true if the fine grid fully covers the domain. The nature of this makes it
  // always true for the "deeper" multigridlevels,  but not so for the intermediate levels. 
  auto isFullyCovered = [&] (const EBLevelGrid& a_eblg) -> bool {
    unsigned long long numPtsLeft = a_eblg.getDomain().domainBox().numPts(); // Number of grid points in
    const DisjointBoxLayout& dbl  = a_eblg.getDBL();

    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit){
      numPtsLeft -= dbl[lit()].numPts();
    }

    return (numPtsLeft == 0);
  };

  const ProblemDomain fineDomain   = a_fineEblg.getDomain();
  const ProblemDomain coarDomain   = coarsen(fineDomain, a_refRat);
  const DisjointBoxLayout& fineDbl = a_fineEblg.getDBL();
  DisjointBoxLayout coarDbl;

  // Check if we can get a coarsenable domain. Don't want to coarsen to 1x1 so hence the factor of 2 in the test here. 
  ProblemDomain test = fineDomain;
  if(refine(coarsen(test, 2*a_refRat), 2*a_refRat) == fineDomain){
    const bool doCoarsen = a_fineEblg.getDBL().coarsenable(2*a_refRat);

    // Use coarsening if we can. 
    if(a_fineEblg.getDBL().coarsenable(2*a_refRat)){
      coarsen(coarDbl, a_fineEblg.getDBL(), a_refRat);
      a_coarEblg.define(coarDbl, coarDomain, m_ghostPhi.max(), a_fineEblg.getEBIS());

      hasCoarser = true;
    }
    else { // Check if we can use box aggregation
      if(isFullyCovered(a_fineEblg)){
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
      else{ // out of ideas
	hasCoarser = false;
      }
    }
  }
  else{ // Nothing we can do. 
    hasCoarser = false;
  }
#endif
  return hasCoarser;
}

void EBHelmholtzOpFactory::coarsenCoefficients(LevelData<EBCellFAB>&              a_coarAcoef,
					       LevelData<EBFluxFAB>&              a_coarBcoef,
					       LevelData<BaseIVFAB<Real> >&       a_coarBcoefIrreg,
					       const LevelData<EBCellFAB>&        a_fineAcoef,
					       const LevelData<EBFluxFAB>&        a_fineBcoef,
					       const LevelData<BaseIVFAB<Real> >& a_fineBcoefIrreg,
					       const EBLevelGrid&                 a_eblgCoar,
					       const EBLevelGrid&                 a_eblgFine,
					       const int                          a_refRat){

  const Interval interv(0,0);
  
  if(a_refRat == 1){
    a_fineAcoef.copyTo(a_coarAcoef);
    a_fineBcoef.copyTo(a_coarBcoef);
    a_fineBcoefIrreg.copyTo(a_coarBcoefIrreg);
  }
  else{
    EbCoarAve averageOp(a_eblgFine.getDBL(),
			a_eblgCoar.getDBL(),
			a_eblgFine.getEBISL(),
			a_eblgCoar.getEBISL(),
			a_eblgCoar.getDomain(),
			a_refRat,
			1,
			a_eblgCoar.getEBIS());

    averageOp.average(a_coarAcoef,      a_fineAcoef,      interv);
    averageOp.average(a_coarBcoef,      a_fineBcoef,      interv);
    averageOp.average(a_coarBcoefIrreg, a_fineBcoefIrreg, interv);

    a_coarAcoef.exchange();
    a_coarBcoef.exchange();
    a_coarBcoefIrreg.exchange();
  }
}

EBHelmholtzOp* EBHelmholtzOpFactory::MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly) {
  EBHelmholtzOp* mgOp = nullptr;

  ParmParse pp;
  bool turnOffMG = false;
  pp.query("turn_off_multigrid", turnOffMG);
  if(turnOffMG){
    pout() << "turning off multigrid for EBHelmholtzOpFactory because 'turn_off_multigrid' = true" << endl;
    return mgOp;
  }
  
  const int amrLevel = this->findAmrLevel(a_fineDomain); // Run-time abort if a_fineDomain is not found in anhy amr level. 
  
  const AmrLevelGrids& mgLevelGrids = m_mgLevelGrids[amrLevel];

  const int mgRefRat = 2;

  // Things that are needed for defining the operator. This might seem weird but
  // the multigrid operators only do relaxation and there's no fine-to-coar stuff. So hasFine = hasCoar = false.
  // But we might need to restrict and interpolate another, even coarser multigrid level so in that case we need
  // to define that. 
  EBLevelGrid eblg;
  EBLevelGrid eblgMgCoar;
    
  RefCountedPtr<EBMultigridInterpolator>             interpolator;  // Only if defined on an AMR level
  RefCountedPtr<EBFluxRegister>                      fluxReg;       // Only if defined on an AMR level
  RefCountedPtr<EbCoarAve>                           coarsener;     // Only if defined on an AMR level

  bool hasMGObjects;
    
  RefCountedPtr<LevelData<EBCellFAB > >       Acoef;
  RefCountedPtr<LevelData<EBFluxFAB > >       Bcoef;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > > BcoefIrreg;

  bool foundMgLevel = false;

  if(a_depth == 0){ // Asking for the AMR level.
    eblg       = *m_amrLevelGrids[amrLevel];
    Acoef      =  m_amrAcoef[amrLevel];
    Bcoef      =  m_amrBcoef[amrLevel];
    BcoefIrreg =  m_amrBcoefIrreg[amrLevel];

    interpolator = m_amrInterpolators[amrLevel];
    fluxReg      = m_amrFluxRegisters[amrLevel];
    coarsener    = m_amrCoarseners[amrLevel];

    hasMGObjects = m_hasMgLevels[amrLevel];

    if(hasMGObjects){
      eblgMgCoar = *m_mgLevelGrids[amrLevel][1];
    }
      
    foundMgLevel = true;
  }
  else{ // Asking for a coarsening. No interp or flux reg object here.
    // TLDR: Go through the coarsened levels for the specified amr level and see if we find a coarsening at the
    //       specified depth.
    const ProblemDomain coarDomain = coarsen(a_fineDomain, std::pow(mgRefRat, a_depth));

    // These are the things that live below the AMR level corresponding to a_fineDomain. 
    const AmrLevelGrids& mgLevelGrids = m_mgLevelGrids[amrLevel];
    const AmrCellData&   mgAcoef      = m_mgAcoef[amrLevel];
    const AmrFluxData&   mgBcoef      = m_mgBcoef[amrLevel];
    const AmrIrreData&   mgBcoefIrreg = m_mgBcoefIrreg[amrLevel];

    // See if we have a corresponding multigrid level.
    int mgLevel;
    for (int img = 0; img < mgLevelGrids.size(); img++){
      if(mgLevelGrids[img]->getDomain() == coarDomain){
	mgLevel      = img;
	foundMgLevel = true;
	break;
      }
    }

    // Found the multigrid level. We can define the operator. 
    if(foundMgLevel){
      Acoef      =  mgAcoef[mgLevel];
      Bcoef      =  mgBcoef[mgLevel];
      BcoefIrreg =  mgBcoefIrreg[mgLevel];
      eblg       = *mgLevelGrids[mgLevel];

      hasMGObjects = (mgLevel < mgLevelGrids.size() - 1); // This just means that mgLevel was not the last entry in mgLevelGrids so there's even coarser stuff below.
      if(hasMGObjects){
	eblgMgCoar = *mgLevelGrids[mgLevel+1];
      }
    }
  }

  // Make the operator
  if(foundMgLevel){
    const Real dx     = m_amrResolutions[amrLevel]*std::pow(mgRefRat, a_depth); // 

    auto domBC = m_domainBcFactory->create();
    auto ebBC  = m_ebBcFactory->create();
    
    mgOp = new EBHelmholtzOp(m_dataLocation,
			     EBLevelGrid(), // Multigrid operator, so no fine. 
			     eblg,
			     EBLevelGrid(), // Multigrid operator, so no cofi. 
			     EBLevelGrid(), // Multigrid operator, so no coarse. 
			     eblgMgCoar,
			     interpolator,  // Defined if an amr level
			     fluxReg,       // Defined if an amr level
			     coarsener,     // Defined if an amr level
			     domBC,
			     ebBC,
			     m_probLo,
			     dx,            // Set from depth
			     1,             // Multigrid operator. Set to 1 in operator anyways.
			     1,             // Multigrid operator. Set to 1 in operator anyways. 
			     false,         // Multigrid operator, so false.
			     false,         // Multigrid operator, so false.
			     hasMGObjects,
			     m_alpha,   
			     m_beta,
			     Acoef,
			     Bcoef,
			     BcoefIrreg,
			     m_ghostPhi,
			     m_ghostRhs,
			     m_smoother);


  }
  
  return mgOp;
}

EBHelmholtzOp* EBHelmholtzOpFactory::AMRnewOp(const ProblemDomain& a_domain) {
  EBHelmholtzOp* op;

  const int amrLevel = this->findAmrLevel(a_domain);

  const bool hasFine = amrLevel < m_numAmrLevels - 1;
  const bool hasCoar = amrLevel > 0;

  EBLevelGrid eblgFine;
  EBLevelGrid eblg;
  EBLevelGrid eblgCoFi;
  EBLevelGrid eblgCoar;
  EBLevelGrid eblgCoarMG;

  Real dx;

  eblg = *m_amrLevelGrids[amrLevel];
  dx   = m_amrResolutions[amrLevel];

  int refToCoar = 1;
  int refToFine = 1;
  
  if(hasCoar){
    eblgCoar  = *m_amrLevelGrids[amrLevel-1];
    refToCoar = m_amrRefRatios[amrLevel-1];
  }

  if(hasFine){
    refToFine  = m_amrRefRatios[amrLevel];
    eblgFine   = *m_amrLevelGrids[amrLevel+1];
  }

  const bool hasMGObjects = m_hasMgLevels[amrLevel];
  if(hasMGObjects){
    eblgCoarMG = *m_mgLevelGrids[amrLevel][1];
    CH_assert(eblgCoarMG.isDefined());
  }

  if(hasCoar){
    this->getCoarserLayout(eblgCoFi, eblg, refToCoar, m_mgBlockingFactor);
    CH_assert(eblgCoFi.isDefined());
  }

  op = new EBHelmholtzOp(m_dataLocation,
			 eblgFine,
			 eblg,
			 eblgCoFi,
			 eblgCoar,
			 eblgCoarMG,
			 m_amrInterpolators[amrLevel],
			 m_amrFluxRegisters[amrLevel],
			 m_amrCoarseners[amrLevel],
			 m_domainBcFactory->create(),
			 m_ebBcFactory->create(),
			 m_probLo,			 
			 dx,
			 refToFine,
			 refToCoar,
			 hasFine,
			 hasCoar,
			 hasMGObjects,
			 m_alpha,
			 m_beta,
			 m_amrAcoef[amrLevel],
			 m_amrBcoef[amrLevel],
			 m_amrBcoefIrreg[amrLevel],
			 m_ghostPhi,
			 m_ghostRhs,
			 m_smoother);

  return op;
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

int EBHelmholtzOpFactory::findAmrLevel(const ProblemDomain& a_domain) const{
  int amrLevel = -1;
  for (int lvl = 0; lvl < m_amrLevelGrids.size(); lvl++){
    if(m_amrLevelGrids[lvl]->getDomain() == a_domain){
      amrLevel = lvl;
      break;
    }
  }

  if(amrLevel < 0) MayDay::Abort("EBHelmholtzOpFactory::findAmrLevel - no corresponding amr level found!");

  return amrLevel;
}


#include <CD_NamespaceFooter.H>
