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

EBHelmholtzOpFactory::EBHelmholtzOpFactory(const Real&                       a_alpha,
					   const Real&                       a_beta,
					   const AmrLevelGrids&              a_amrLevelGrids,
					   const AmrInterpolators&           a_amrInterpolators,
					   const AmrFluxRegisters&           a_amrFluxRegisters,
					   const AmrCoarseners&              a_amrCoarseners,
					   const AmrRefRatios&               a_amrRefRatios,
					   const AmrResolutions&             a_amrResolutions,					   
					   const AmrCellData&                a_amrAcoef,
					   const AmrFluxData&                a_amrBcoef,
					   const AmrIrreData&                a_amrBcoefIrreg,
					   const EBHelmholtzDomainBcFactory& a_domainBcFactory,
					   const EBHelmholtzEbBcFactory&     a_ebBcFactory,
					   const IntVect&                    a_ghostPhi,
					   const IntVect&                    a_ghostRhs,
					   const RelaxType&                  a_relaxationMethod,
					   const ProblemDomain&              a_bottomDomain,
					   const int&                        a_mgBlockingFactor,
					   const AmrLevelGrids&              a_deeperLevelGrids){

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
    m_hasMgLevels[amrLevel] = false;
    
    if(amrLevel == 0 && this->isCoarser(m_bottomDomain, m_amrLevelGrids[amrLevel]->getDomain())){
      m_hasMgLevels[amrLevel] = true;
    }
    else if(m_amrRefRatios[amrLevel-1] > 2){ // There must be one intermediate level 
      m_hasMgLevels[amrLevel] = true;
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
	const EBLevelGrid& mgEblgFine = *m_mgLevelGrids[amrLevel][curMgLevels-1]; 

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
	  hasCoarser = this->getCoarserLayout(mgEblgCoar, mgEblgFine, mgRefRatio, m_mgBlockingFactor);
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
    else{
      m_hasMgLevels[amrLevel] = false;
    }

#if 1 // Debug
    if(procID() == 0){
      std::cout << "amrLevel = " << amrLevel << "\t domain = " << m_amrLevelGrids[amrLevel]->getDomain() << ":\n";
      for (int mglevel = 0; mglevel < m_mgLevelGrids[amrLevel].size(); mglevel++){
	std::cout << "\t mg level = " << mglevel << "\t mg domain = " << m_mgLevelGrids[amrLevel][mglevel]->getDomain() << "\n";
      }
    }
#endif
  }
}

bool EBHelmholtzOpFactory::isCoarser(const ProblemDomain& A, const ProblemDomain& B) const{
  return A.domainBox().numPts() < B.domainBox().numPts();
}


bool EBHelmholtzOpFactory::isFiner(const ProblemDomain& A, const ProblemDomain& B) const{
  return A.domainBox().numPts() > B.domainBox().numPts();
}

bool EBHelmholtzOpFactory::getCoarserLayout(RefCountedPtr<EBLevelGrid>& a_coarEblg, const EBLevelGrid& a_fineEblg, const int a_refRat, const int a_blockingFactor) const {
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

    doCoarsen = a_fineEblg.getDBL().coarsenable(a_refRat);

    // Prefer aggregation over direct coarsening
    if(doAggregation){ 
      Vector<Box> boxes;
      Vector<int> procs;

      // We could have use our load balancing here, but I don't see why. 
      domainSplit(coarDomain, boxes, a_blockingFactor);
      mortonOrdering(boxes);
      LoadBalance(procs, boxes);

      coarDbl.define(boxes, procs, coarDomain);
      a_coarEblg->define(coarDbl, coarDomain, m_ghostPhi.max(), a_fineEblg.getEBIS());

      hasCoarser = true;
    }
    else if(doCoarsen){ // But use coarsening if we must.
      coarsen(coarDbl, a_fineEblg.getDBL(), a_refRat);
      a_coarEblg->define(coarDbl, coarDomain, m_ghostPhi.max(), a_fineEblg.getEBIS());

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
  }
}

EBHelmholtzOp* EBHelmholtzOpFactory::MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly) {
  EBHelmholtzOp* mgOp = nullptr;

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

  bool hasCoarMg;
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

    hasCoarMg = m_hasMgLevels[amrLevel];

    if(hasCoarMg){
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
      if(mgLevelGrids[mgLevel]->getDomain() == coarDomain){
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

      hasCoarMg = (mgLevel < mgLevelGrids.size() - 1); // This just means that mgLevel was not the last entry in mgLevelGrids so there's even coarser stuff below.
      if(hasCoarMg){
	eblgMgCoar = *mgLevelGrids[mgLevel+1];
      }
    }
  }

  // Make the operator
  if(foundMgLevel){

    const Real dx     = m_amrResolutions[amrLevel]*std::pow(mgRefRat, a_depth); // 
    const Real dxCoar = (amrLevel > 0) ? m_amrResolutions[amrLevel-1] : -1.0;

    auto dobc = this->makeDomainBcObject(eblg, dx);
    auto ebbc = this->makeEbBcObject    (eblg, dx);
									
    mgOp = new EBHelmholtzOp(EBLevelGrid(), // Multigrid operator, so no fine. 
			     eblg,
			     EBLevelGrid(), // Multigrid operator, so no coarse. 
			     eblgMgCoar,
			     interpolator,  // Defined if an amr level
			     fluxReg,       // Defined if an amr level
			     coarsener,     // Defined if an amr level
			     dobc,
			     ebbc,          
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
			     m_relaxMethod);
  }
  
  return mgOp;
 
}

EBHelmholtzOp* EBHelmholtzOpFactory::AMRnewOp(const ProblemDomain& a_domain) {
  EBHelmholtzOp* op;

  const int amrLevel = this->findAmrLevel(a_domain);

  const bool hasFine = amrLevel < m_numAmrLevels - 1;
  const bool hasCoar = amrLevel > 0;
  
  EBLevelGrid eblgFine;
  EBLevelGrid eblgCoar;
  EBLevelGrid eblgCoarMG;

  int refToFine = 1;
  int refToCoar = 1;
  
  if(hasCoar){
    eblgCoar  = *m_amrLevelGrids[amrLevel-1];
    refToCoar = m_amrRefRatios[amrLevel-1];
  }

  if(hasFine){
    eblgFine   = *m_amrLevelGrids[amrLevel];
    refToFine  = m_amrRefRatios[amrLevel];
  }

  const bool hasMGObjects = m_hasMgLevels[amrLevel];
  if(hasMGObjects){
    eblgCoarMG = *m_mgLevelGrids[amrLevel][1];
  }

  auto dobc = this->makeDomainBcObject(*m_amrLevelGrids[amrLevel], m_amrResolutions[amrLevel]);
  auto ebbc = this->makeEbBcObject    (*m_amrLevelGrids[amrLevel], m_amrResolutions[amrLevel]);

  op = new EBHelmholtzOp(eblgFine,
			 *m_amrLevelGrids[amrLevel],
			 eblgCoar,
			 eblgCoarMG,
			 m_amrInterpolators[amrLevel],
			 m_amrFluxRegisters[amrLevel],
			 m_amrCoarseners[amrLevel],
			 dobc,
			 ebbc,
			 m_amrResolutions[amrLevel],
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
			 m_relaxMethod);

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

RefCountedPtr<EBHelmholtzOp::EBHelmholtzEbBc> EBHelmholtzOpFactory::makeEbBcObject(const EBLevelGrid& a_eblg, const Real& a_dx) const {
  EBHelmholtzOp::EBHelmholtzEbBc* bc = (EBHelmholtzOp::EBHelmholtzEbBc*) m_ebBcFactory->create(a_eblg.getDomain(),
											       a_eblg.getEBISL(),
											       a_dx*RealVect::Unit,
											       &m_ghostPhi,
											       &m_ghostRhs);

  return RefCountedPtr<EBHelmholtzOp::EBHelmholtzEbBc>(bc);
}

RefCountedPtr<EBHelmholtzOp::EBHelmholtzDomainBc> EBHelmholtzOpFactory::makeDomainBcObject(const EBLevelGrid& a_eblg, const Real& a_dx) const {
  EBHelmholtzOp::EBHelmholtzDomainBc* bc = (EBHelmholtzOp::EBHelmholtzDomainBc*) m_domainBcFactory->create(a_eblg.getDomain(),
													   a_eblg.getEBISL(),
													   a_dx*RealVect::Unit);
  
  return RefCountedPtr<EBHelmholtzOp::EBHelmholtzDomainBc>(bc);
}

int EBHelmholtzOpFactory::findAmrLevel(const ProblemDomain& a_domain) const{
  int amrLevel = -1;
  for (int lvl = 0; lvl < m_amrLevelGrids.size(); lvl++){
    if(m_amrLevelGrids[lvl]->getDomain() == a_domain){
      amrLevel = lvl;
      break;
    }
  }

  if(amrLevel < 0) MayDay::Abort("EBHelmholtzOpFactory::AMRnewOp - no corresponding amr level found!");

  return amrLevel;
}

#include <CD_NamespaceFooter.H>
