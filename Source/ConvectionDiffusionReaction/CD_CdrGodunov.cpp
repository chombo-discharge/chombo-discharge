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
#include "CD_NamespaceHeader.H"

CdrGodunov::CdrGodunov() : CdrTGA() {
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
  this->parseRngSeed();                // Creates a seed (if one runs with FHD). 
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

  // CdrTGA allocates everything except storage needed for the advection object. 
  CdrTGA::allocateInternals();


  // Allocate levelAdvect only if the solver is mobile. See Chombo design docs for how the EBAdvectLevelIntegrator operates. 
  if(m_isMobile){
    const Vector<RefCountedPtr<EBLevelGrid> >& eblgs = m_amr->getEBLevelGrid(m_realm, m_phase);
    const Vector<DisjointBoxLayout>& grids           = m_amr->getGrids(m_realm);
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
  if(m_extrapolateSourceTerm && a_extrapDt > 0.0){
    this->computeDivD(m_scratch, (EBAMRCellData&) a_cellPhi, false, false); // This will touch ghost cells in a_cellPhi, but those should be linearly interpolated anyways. 

    DataOps::incr(m_scratch, m_source, 1.0);
    
    m_amr->averageDown(m_scratch, m_realm, m_phase);
    m_amr->interpGhost(m_scratch, m_realm, m_phase);
  }
  else{
    DataOps::setValue(m_scratch, 0.0);
  }

  // This code extrapolates the cell-centered state to face centers on every grid level. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const Real dx                = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator();dit.ok(); ++dit){
      
      EBFluxFAB& facePhi       = (*a_facePhi     [lvl])[dit()];
      const EBCellFAB& cellPhi = (*a_cellPhi     [lvl])[dit()];
      const EBCellFAB& cellVel = (*m_cellVelocity[lvl])[dit()];
      const EBFluxFAB& faceVel = (*m_faceVelocity[lvl])[dit()];
      const EBCellFAB& source  = (*m_scratch     [lvl])[dit()];
      const EBISBox& ebisbox   = ebisl[dit()];
      const Real time          = 0.0;

      EBAdvectPatchIntegrator& ebpatchad = m_levelAdvect[lvl]->getPatchAdvect(dit());

      // These are settings for EBAdvectPatchIntegrator -- it's not a very pretty design but the object has settings
      // that permits it to run advection code (through setDoingVel(0)). 
      ebpatchad.setVelocities(cellVel, faceVel);       // Set cell/face velocities
      ebpatchad.setDoingVel(0);                        // If setDoingVel(0) EBAdvectLevelIntegrator advects a scalar. 
      ebpatchad.setEBPhysIBC(ExtrapAdvectBCFactory()); // Set the BC object. It won't matter what we use here because CdrSolver runs its own BC routines. 
      ebpatchad.setCurComp(m_comp);                    // Solving for m_comp = 0

      // Extrapolate to face-centers. The face-centered states are Godunov-style extrapolated in time to a_extrapDt. 
      ebpatchad.extrapolateBCG(facePhi, cellPhi, source, dit(), time, a_extrapDt);
    }
  }
}

#include <CD_NamespaceFooter.H>
