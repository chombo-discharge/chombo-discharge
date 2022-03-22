/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrGodunov.cpp
  @brief  Implementation of CD_CdrGodunov.H
  @author Robert Marskar
*/

// Chombo includes
#include <ExtrapAdvectBC.H>
#include <EBArith.H>
#include <ParmParse.H>

// Our includes
#include <CD_CdrGodunov.H>
#include <CD_DataOps.H>
#include <CD_ParallelOps.H>
#include "CD_NamespaceHeader.H"

CdrGodunov::CdrGodunov() : CdrMultigrid() {
  CH_TIME("CdrGodunov::CdrGodunov()");

  // Class name and instantiation name. 
  m_className = "CdrGodunov";
  m_name      = "CdrGodunov";
}

CdrGodunov::~CdrGodunov(){
  CH_TIME("CdrGodunov::~CdrGodunov()");
}

void CdrGodunov::parseOptions(){
  CH_TIME("CdrGodunov::parseOptions()");
  if(m_verbosity > 5){
    pout() << m_name + "::parseOptions()" << endl;
  }
  
  this->parseDomainBc();               // Parses domain BC options
  this->parseSlopeLimiter();           // Parses slope limiter settings
  this->parsePlotVariables();          // Parses plot variables
  this->parsePlotMode();               // Parses plot mode
  this->parseMultigridSettings();      // Parses multigrid settings
  this->parseExtrapolateSourceTerm();  // Parses source term extrapolation for Godunov time extrapolation. 
  this->parseDivergenceComputation();  // Parses non-conservative divergence blending
}

void CdrGodunov::parseRuntimeOptions(){
  CH_TIME("CdrGodunov::parseRuntimeOptions()");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRuntimeOptions()" << endl;
  }

  this->parseSlopeLimiter();          // Parses slope limiter (on/off)
  this->parsePlotVariables();         // Parses plot variables
  this->parsePlotMode();              // Parses plot mode
  this->parseMultigridSettings();     // Parses multigrid settings
  this->parseDomainBc();              // Parses domain BCs
  this->parseExtrapolateSourceTerm(); // Parses source term extrapolation
  this->parseDivergenceComputation(); // Parses non-conservative divergence blending. 
}

Real CdrGodunov::computeAdvectionDt(){
  CH_TIME("CdrGodunov::computeAdvectionDt()");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAdvectionDt()" << endl;
  }

  // TLDR: For advection, Bell, Collela, and Glaz says we must have dt <= dx/max(|vx|, |vy|, |vz|). See these two papers for details:
  //
  //       Bell, Colella, Glaz, J. Comp. Phys 85 (257), 1989
  //       Minion, J. Comp. Phys 123 (435), 1996

  Real minDt = std::numeric_limits<Real>::max();

  if(m_isMobile){
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)              [lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real               dx    = m_amr->getDx()                        [lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit){
	const Box        cellBox = dbl  [dit()];	
	const EBCellFAB& velo    = (*m_cellVelocity[lvl])[dit()];
	const EBISBox&   ebisBox = ebisl[dit()];

	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

	// Regular grid data.
	const BaseFab<Real>& veloReg = velo.getSingleValuedFAB();	

	// Compute dt = dx/(|vx|+|vy|+|vz|) and check if it's smaller than the smallest so far. 
	auto regularKernel = [&](const IntVect& iv) -> void {
	  Real velMax = 0.0;	  
	  if(!ebisBox.isCovered(iv)){
	    for (int dir = 0; dir < SpaceDim; dir++){
	      velMax = std::max(velMax, std::abs(veloReg(iv,dir)));
	    }
	  }

	  if(velMax > 0.0){
	    minDt = std::min(dx/velMax, minDt);
	  }	  
	};

	// Same kernel, but for cut-cells.
	auto irregularKernel = [&](const VolIndex& vof) -> void {
	  Real velMax = 0.0;
	  for (int dir = 0; dir < SpaceDim; dir++){
	    velMax = std::max(velMax, std::abs(velo(vof, dir)));
	  }

	  if(velMax > 0.0){
	    minDt = std::min(dx/velMax, minDt);
	  }
	};


	// Execute the kernels.
	BoxLoops::loop(cellBox, regularKernel  );
	BoxLoops::loop(vofit,   irregularKernel);
      }
    }

    // If we are using MPI then ranks need to know of each other's time steps.
    minDt = ParallelOps::min(minDt);
  }

  return minDt;
}

void CdrGodunov::parseExtrapolateSourceTerm(){
  CH_TIME("CdrGodunov::parseExtrapolateSourceTerm()");
  if(m_verbosity > 5){
    pout() << m_name + "::parseExtrapolateSourceTerm()" << endl;
  }
  
  ParmParse pp(m_className.c_str());

  pp.get("extrap_source", m_extrapolateSourceTerm);
}

void CdrGodunov::parseSlopeLimiter(){
  CH_TIME("CdrGodunov::parseSlopeLimiter()");
  if(m_verbosity > 5){
    pout() << m_name + "::parseSlopeLimiter()" << endl;
  }
  
  ParmParse pp(m_className.c_str());

  pp.get("limit_slopes", m_limitSlopes);
}

void CdrGodunov::allocateInternals(){
  CH_TIME("CdrSolver::allocateInternals()");
  if(m_verbosity > 5){
    pout() << m_name + "::allocateInternals()" << endl;
  }

  // CdrMultigrid allocates everything except storage needed for the advection object. 
  CdrMultigrid::allocateInternals();


  // Allocate levelAdvect only if the solver is mobile. See Chombo design docs for how the EBAdvectLevelIntegrator operates. 
  if(m_isMobile){
    const Vector<RefCountedPtr<EBLevelGrid> >& eblgs = m_amr->getEBLevelGrid(m_realm, m_phase);
    const Vector<int>& refRatios                     = m_amr->getRefinementRatios();
    const Vector<Real>& dx                           = m_amr->getDx();
    const int finestLevel                            = m_amr->getFinestLevel();
    
    m_levelAdvect.resize(1 + finestLevel);
  
    for (int lvl = 0; lvl <= finestLevel; lvl++){

      const bool hasCoar = lvl > 0;
      const bool hasFine = lvl < finestLevel;

      int refRat = 1;
      EBLevelGrid coarEblg;
      if(hasCoar){
	coarEblg = *eblgs[lvl-1];
	refRat   = refRatios[lvl-1];
      }

      // Note: There is a "bug" in the function signature in Chombo. The second-to-last argument is the slope limiter and not the EBCF things. 
      m_levelAdvect[lvl] = RefCountedPtr<EBAdvectLevelIntegrator> (new EBAdvectLevelIntegrator(*eblgs[lvl],
												coarEblg,
												refRat,
												dx[lvl]*RealVect::Unit,
												hasCoar,
												hasFine,
												false,
												m_limitSlopes,
												m_ebis));
    }
  }
}
  
void CdrGodunov::advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi, const Real a_extrapDt){
  CH_TIME("CdrGodunov::advectToFaces(EBAMRFluxDat, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::advectToFaces(EBAMRFluxDat, EBAMRCellData, Real)" << endl;
  }

  CH_assert(a_facePhi[0]->nComp() == 1);
  CH_assert(a_cellPhi[0]->nComp() == 1);

    // If we are extrapolating in time the source term will yield different states on the face centers. The source term is the source
    // S^k + Div(D*Grad(Phi)), which we add to the solver below. It is stored on m_scratch (which won't be used elsewhere, I think).
  DataOps::setValue(m_scratch, 0.0);

  // Compute viscous source term for the advection.
#if 0
  if(m_isDiffusive && a_extrapDt > 0.0) {
    this->setupDiffusionSolver();

    this->computeKappaLphi(m_scratch, a_cellPhi);
  }
#endif

  // If user asks for it, add the source term to the extrapolation. 
  if(m_extrapolateSourceTerm && a_extrapDt > 0.0){
    DataOps::incr(m_scratch, m_source, 1.0);
  }

  m_amr->averageDown(m_scratch, m_realm, m_phase);
  m_amr->interpGhost(m_scratch, m_realm, m_phase);  

  // This code extrapolates the cell-centered state to face centers on every grid level, in both space and time. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      EBFluxFAB&       facePhi = (*a_facePhi     [lvl])[dit()];
      const EBCellFAB& cellPhi = (*a_cellPhi     [lvl])[dit()];
      const EBCellFAB& cellVel = (*m_cellVelocity[lvl])[dit()];
      const EBFluxFAB& faceVel = (*m_faceVelocity[lvl])[dit()];
      const EBCellFAB& source  = (*m_scratch     [lvl])[dit()]; // Note: the source consists of the source term and the diffusive contribution. 
      const Real time          = 0.0;

      EBAdvectPatchIntegrator& ebAdvectPatch = m_levelAdvect[lvl]->getPatchAdvect(dit());

      // These are settings for EBAdvectPatchIntegrator -- it's not a very pretty design but the object has settings
      // that permits it to run advection code (through setDoingVel(0)). 
      ebAdvectPatch.setVelocities(cellVel, faceVel);       // Set cell/face velocities
      ebAdvectPatch.setDoingVel(0);                        // If setDoingVel(0) EBAdvectLevelIntegrator advects a scalar. 
      ebAdvectPatch.setEBPhysIBC(ExtrapAdvectBCFactory()); // Set the BC object. It won't matter what we use here because CdrSolver runs its own BC routines. 
      ebAdvectPatch.setCurComp(m_comp);                    // Solving for m_comp = 0

      // Extrapolate to face-centers. The face-centered states are Godunov-style extrapolated in time to a_extrapDt. 
      ebAdvectPatch.extrapolateBCG(facePhi, cellPhi, source, dit(), time, a_extrapDt);
    }
  }
}

void CdrGodunov::extrapFluxToEB(EBAMRIVData& a_fluxEB, const EBAMRCellData& a_cellPhi, const Real a_extrapDt) {
  CH_TIME("CdrGodunov::extrapFluxToEB(EBAMRIVData, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::extrapFluxToEB(EBAMRIVData, EBAMRCellData, Real)" << endl;
  }

  CH_assert(a_fluxEB [0]->nComp() == 1);
  CH_assert(a_cellPhi[0]->nComp() == 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain&     domain = m_amr->getDomains()[lvl];    
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx     = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box      cellBox = dbl  [dit()];
      const EBISBox& ebisBox = ebisl[dit()];

      BaseIVFAB<Real>& flux    = (*a_fluxEB [lvl])[dit()];
      const EBCellFAB& phi     = (*a_cellPhi[lvl])[dit()];

      const bool isAllRegular = ebisBox.isAllRegular();
      const bool isAllCovered = ebisBox.isAllCovered();
      const bool isIrregular  = !isAllRegular && !isAllCovered;

      if(isIrregular) {

	EBCellFAB slopes[SpaceDim];
	if(m_isMobile){ // Compute slopes
	  const EBCellFAB& cellVel = (*m_cellVelocity[lvl])[dit()];
	  const EBFluxFAB& faceVel = (*m_faceVelocity[lvl])[dit()];	
	
	  EBAdvectPatchIntegrator& ebAdvectPatch = m_levelAdvect[lvl]->getPatchAdvect(dit());

	  // These are settings for EBAdvectPatchIntegrator -- it's not a very pretty design but the object has settings
	  // that permits it to run advection code (through setDoingVel(0)). 
	  ebAdvectPatch.setVelocities(cellVel, faceVel);       // Set cell/face velocities
	  ebAdvectPatch.setDoingVel(0);                        // If setDoingVel(0) EBAdvectLevelIntegrator advects a scalar. 
	  ebAdvectPatch.setEBPhysIBC(ExtrapAdvectBCFactory()); // Set the BC object. It won't matter what we use here because CdrSolver runs its own BC routines. 
	  ebAdvectPatch.setCurComp(m_comp);                    // Solving for m_comp = 0

	  // Define and compute slopes. 
	  for (int dir = 0; dir < SpaceDim; dir++){
	    slopes[dir].clone(phi);

	    ebAdvectPatch.slope(slopes[dir],
				phi,
				dir, cellBox);
	  }
	}	

	// Extrapolation stencils for velocity and diffusion coefficient.
	const BaseIVFAB<VoFStencil>& ebInterpStencils  = m_amr->getEbCentroidInterpolationStencils(m_realm, m_phase)[lvl][dit()];
	const BaseIVFAB<VoFStencil>& normDerivStencils = (*(m_stencilsDphiDn[lvl]))[dit()];

	auto extrapKernel = [&] (const VolIndex& vof) -> void {
	  flux(vof, m_comp) = 0;

	  // Interpolation stencil and stencil for computing d(phi)/dn
	  const VoFStencil& interpSten    = ebInterpStencils (vof, m_comp);
	  const VoFStencil& normDerivSten = normDerivStencils(vof, m_comp);

	  const RealVect bndryCentroid = ebisBox.bndryCentroid(vof);

	  // If mobile, add extrapolated advective flux. 
	  if(m_isMobile){
	    const EBCellFAB& cellVel = (*m_cellVelocity[lvl])[dit()];
	    
	    Real phiEB  = 0.0; // Phi on EB
	    Real velEB  = 0.0; // Velocity on EB

	    // Compute normal velocity on EB face.
	    RealVect v = RealVect::Zero;
	    for (int i = 0; i < interpSten.size(); i++){
	      for (int dir = 0; dir < SpaceDim; dir++){
		v[dir] += interpSten.weight(i) * cellVel(interpSten.vof(i), dir);
	      }
	    }
	    velEB = v.dotProduct(ebisBox.normal(vof)); // Negative sign because normal points into the computational region. 

	    // Extrapolate phi to EB centroid.
	    phiEB = phi(vof, m_comp);
	    for (int dir = 0; dir < SpaceDim; dir++){
	      phiEB += slopes[dir](vof, m_comp) * bndryCentroid[dir];
	    }

	    flux(vof, m_comp) += phiEB * velEB;
	  }

	  // If diffusive, add diffusion flux D * d(phi)/dn. 
	  if(m_isDiffusive) {
	    const BaseIVFAB<Real>& dcoEB = (*m_ebCenteredDiffusionCoefficient[lvl])[dit()];      	    

	    // Compute d(phi)/dn
	    Real gradEB = 0.0; // grad(phi) on EB.	    	    
	    for (int i = 0; i < normDerivSten.size(); i++){
	      gradEB += normDerivSten.weight(i) * phi(normDerivSten.vof(i), m_comp);
	    }

	    flux(vof, m_comp) -= dcoEB(vof, m_comp) * gradEB;
	  }

	  // Flip sign because everything was computed with the EB normal point into the cell but we want
	  // to have it outwards.
	  flux(vof, m_comp) = -flux(vof, m_comp);
	};

	// Kernel region.
	VoFIterator& vofit = (*(m_amr->getVofIterator(m_realm, m_phase)[lvl]))[dit()];

	// Execute kernel.
	BoxLoops::loop(vofit, extrapKernel);
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
