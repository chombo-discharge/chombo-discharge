/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrSolver.cpp
  @brief  Implementation of CD_CdrSolver.H
  @author Robert Marskar
  @todo   Need to implement parseDomainBc
  @todo   Should computeAdvectionFlux and computeDiffusionFlux be ghosted... filling fluxes on ghost faces...?
*/

// Chombo includes
#include <ParmParse.H>
#include <EBAMRIO.H>
#include <EBArith.H>

// Our includes
#include <CD_CdrSolver.H>
#include <CD_CdrSolverF_F.H>
#include <CD_CdrFhdF_F.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

constexpr int CdrSolver::m_comp;
constexpr int CdrSolver::m_nComp;

CdrSolver::CdrSolver(){

  // Default options.
  m_verbosity  = -1;
  m_name       = "CdrSolver";
  m_className  = "CdrSolver";

  this->setRealm(Realm::Primal);
  this->setDefaultDomainBC();    // Set default domain BCs (wall)
}

CdrSolver::~CdrSolver(){

}

void CdrSolver::setDefaultDomainBC(){
  CH_TIME("CdrSolver::setDefaultDomainBC()");
  if(m_verbosity > 5){
    pout() << m_name + "::setDefaultDomainBC()" << endl;
  }

  // TLDR: This sets the domain boundary condition to be a wall BC (no incoming/outgoing mass). 

  // Lambda function for wall bc -- mostly left in place so I can remind myself how to do this.
  auto zero = [](const RealVect a_position, const Real a_time){
    return 0.0;
  };

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const CdrDomainBC::DomainSide domainSide = std::make_pair(dir, sit());

      this->setDomainBcType    (domainSide, CdrDomainBC::BcType::Wall);
      this->setDomainBcFunction(domainSide, zero);
    }
  }
}

void CdrSolver::setDomainBcType(const CdrDomainBC::DomainSide a_domainSide, const CdrDomainBC::BcType a_bcType){
  CH_TIME("CdrSolver::setDomainBcType(CdrDomainBC::DomainSide, CdrDomainBC::BcType)");
  if(m_verbosity > 5){
    pout() << m_name + "::setDomainBcType(CdrDomainBC::DomainSide, CdrDomainBC::BcType)" << endl;
  }

  m_domainBC.setBcType(a_domainSide, a_bcType);
}


void CdrSolver::setDomainBcFunction(const CdrDomainBC::DomainSide a_domainSide, const CdrDomainBC::FluxFunction a_fluxFunction){
  CH_TIME("CdrSolver::setDomainBcFunction(CdrDomainBC::DomainSide, CdrDomainBC::FluxFunction)");
  if(m_verbosity > 5){
    pout() << m_name + "::setDomainBcFunction(CdrDomainBC::DomainSide, CdrDomainBC::FluxFunction)" << endl;
  }

  m_domainBC.setBcFunction(a_domainSide, a_fluxFunction);
}

std::string CdrSolver::getName() const {
  CH_TIME("CdrSolver::getName()");
  if(m_verbosity > 5){
    pout() << m_name + "::getName()" << endl;
  }
  
  return m_name;
}

std::string CdrSolver::getRealm() const {
  CH_TIME("CdrSolver::getRealm()");
  if(m_verbosity > 5){
    pout() << m_name + "::getRealm()" << endl;
  }
  
  return m_realm;
}

void CdrSolver::setRealm(const std::string a_realm) {
  CH_TIME("CdrSolver::setRealm(std::string)");
  if(m_verbosity > 5){
    pout() << m_name + "::setRealm(std::string)" << endl;
  }
  
  m_realm.assign(a_realm);
}

Vector<std::string> CdrSolver::getPlotVariableNames() const {
  CH_TIME("CdrSolver::getPlotVariableNames()");
  if(m_verbosity > 5){
    pout() << m_name + "::getPlotVariableNames()" << endl;
  }

  // TLDR: Possible plot variables is the density (m_phi), diffusion coefficient, source term, velocity, and eb flux. This
  //       function returns the associated plot variable names, and will be used in the plot files. 
  
  Vector<std::string> plotVarNames(0);
  
  if(m_plotPhi)                                     plotVarNames.push_back(m_name + " phi");
  if(m_plotDiffusionCoefficient && m_isDiffusive)   plotVarNames.push_back(m_name + " diffusion_coefficient");
  if(m_plotSource)                                  plotVarNames.push_back(m_name + " source");
  if(m_plotVelocity && m_isMobile)                  plotVarNames.push_back("x-Velocity " + m_name);
  if(m_plotVelocity && m_isMobile)                  plotVarNames.push_back("y-Velocity " + m_name);
  if(m_plotVelocity && m_isMobile && SpaceDim == 3) plotVarNames.push_back("z-Velocity " + m_name);
  if(m_plotEbFlux && (m_isMobile || m_isDiffusive)) plotVarNames.push_back(m_name + " eb_flux");
  
  return plotVarNames;
}

int CdrSolver::getNumberOfPlotVariables() const {
  CH_TIME("CdrSolver::getNumberOfPlotVariables()");
  if(m_verbosity > 5){
    pout() << m_name + "::getNumberOfPlotVariables()" << endl;
  }

  // TLDR: Possible plot variables is the density (m_phi), diffusion coefficient, source term, velocity, and eb flux. This
  //       function returns the number of variables that this sum up to (a vector field is 2/3 variables in 2D/3D)

  int numPlotVars = 0;

  if(m_plotPhi)                                   numPlotVars = numPlotVars + 1;
  if(m_plotDiffusionCoefficient && m_isDiffusive) numPlotVars = numPlotVars + 1;
  if(m_plotSource)                                numPlotVars = numPlotVars + 1;
  if(m_plotVelocity && m_isMobile)                numPlotVars = numPlotVars + SpaceDim;
  if(m_plotEbFlux && m_isMobile)                  numPlotVars = numPlotVars + 1;

  return numPlotVars;
}

void CdrSolver::advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt){
  CH_TIME("CdrSolver::advanceEuler(EBAMRCellData, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceEuler(EBAMRCellData, EBAMRCellData, Real)" << endl;
  }

  CH_assert(a_newPhi[0]->nComp() == 1);
  CH_assert(a_oldPhi[0]->nComp() == 1);

  // TLDR: We are solving phi^(k+1) - phi^k - dt*Div(D*Grad(phi^(k+1))) = 0.0. We create a source term = 0 and call the implementation
  //       version.
  if(m_isDiffusive){
    EBAMRCellData src;
    m_amr->allocate(src, m_realm, m_phase, m_nComp);
    DataOps::setValue(src, 0.0);

    this->advanceEuler(a_newPhi, a_oldPhi, src, a_dt);
  }
  else{
    DataOps::copy(a_newPhi, a_oldPhi);
  }
}

void CdrSolver::advanceTGA(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt){
  CH_TIME("CdrSolver::advanceTGA(EBAMRCellData, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceTGA(EBAMRCellData, EBAMRCellData, Real)" << endl;
  }

  CH_assert(a_newPhi[0]->nComp() == 1);
  CH_assert(a_oldPhi[0]->nComp() == 1);

  // TLDR: We are solving the TGA equation (a second order implicit Runge-Kutta scheme). We create a source term = 0 and call the
  //       implementation version. 

  if(m_isDiffusive){
    EBAMRCellData src;
    m_amr->allocate(src, m_realm, m_phase, m_nComp);
    DataOps::setValue(src, 0.0);

    this->advanceTGA(a_newPhi, a_oldPhi, src, a_dt);
  }
  else{
    DataOps::copy(a_newPhi, a_oldPhi);
  }
}

void CdrSolver::allocateInternals(){
  CH_TIME("CdrSolver::allocateInternals()");
  if(m_verbosity > 5){
    pout() << m_name + "::allocateInternals()" << endl;
  }

  // TLDR: This allocates a number of fields for this solver. We wish to trim memory if we can, so we don't allocate
  //       storage for velocities/diffusion coefficients if those fields are not going to be used. Somewhat confusingly,
  //       we allocate pointers for these fields but no memory blocks (for interface reasons). 

  // These three are allocated no matter what -- we will always need the state (and probably also the source term)
  m_amr->allocate(m_phi,     m_realm, m_phase, m_nComp);
  m_amr->allocate(m_source,  m_realm, m_phase, m_nComp);
  m_amr->allocate(m_scratch, m_realm, m_phase, m_nComp);
  
  DataOps::setValue(m_phi,     0.0);
  DataOps::setValue(m_source,  0.0);
  DataOps::setValue(m_scratch, 0.0);

  // Only allocate memory for cell-centered and face-centered velocities if the solver is mobile. Otherwise, allocate
  // a NULL pointer that we can pass around in TimeStepper in order to handle special cases
  if(m_isMobile){
    m_amr->allocate(m_cellVelocity, m_realm, m_phase, SpaceDim);
    m_amr->allocate(m_faceVelocity, m_realm, m_phase, m_nComp );      
    m_amr->allocate(m_faceStates,   m_realm, m_phase, m_nComp );

    DataOps::setValue(m_cellVelocity, 0.0);    
    DataOps::setValue(m_faceVelocity, 0.0);
    DataOps::setValue(m_faceStates,   0.0);
  }
  else{
    m_amr->allocatePointer(m_faceVelocity);
    m_amr->allocatePointer(m_cellVelocity);
    m_amr->allocatePointer(m_faceStates  );
  }

  // Only allocate memory for diffusion coefficients if we need it. Otherwise, allocate a NULL pointer that we can
  // pass around pointers in case we need it. This might seem confusing but its easier to resize those vectors
  // and set nullpointers when we iterate through grid levels and patches. 
  if(m_isDiffusive){
    m_amr->allocate(m_faceCenteredDiffusionCoefficient, m_realm, m_phase, m_nComp);
    m_amr->allocate(m_ebCenteredDiffusionCoefficient,   m_realm, m_phase, m_nComp);
    
    DataOps::setValue(m_faceCenteredDiffusionCoefficient, 0.0);
    DataOps::setValue(m_ebCenteredDiffusionCoefficient,   0.0);
  }
  else{
    m_amr->allocatePointer(m_faceCenteredDiffusionCoefficient);
    m_amr->allocatePointer(m_ebCenteredDiffusionCoefficient  );
  }

  // Allocate stuff for holding fluxes -- this data is used when computing advection and diffusion fluxes. 
  if(m_isDiffusive || m_isMobile){
    m_amr->allocate(m_scratchFluxOne, m_realm, m_phase, m_nComp);
    m_amr->allocate(m_scratchFluxTwo, m_realm, m_phase, m_nComp);
  }

  // These don't consume (much) memory so we always allocate them. 
  m_amr->allocate(m_ebFlux,              m_realm, m_phase, m_nComp);
  m_amr->allocate(m_ebZero,              m_realm, m_phase, m_nComp);
  m_amr->allocate(m_domainFlux,          m_realm, m_phase, m_nComp);
  m_amr->allocate(m_massDifference,      m_realm, m_phase, m_nComp);
  m_amr->allocate(m_nonConservativeDivG, m_realm, m_phase, m_nComp);
  
  DataOps::setValue(m_ebFlux,              0.0);
  DataOps::setValue(m_ebZero,              0.0);
  DataOps::setValue(m_domainFlux,          0.0);
  DataOps::setValue(m_massDifference,      0.0);
  DataOps::setValue(m_nonConservativeDivG, 0.0);

  // Define interpolation stencils
  this->defineInterpolationStencils();
}

void CdrSolver::deallocateInternals(){
  CH_TIME("CdrSolver::deallocateInternals()");
  if(m_verbosity > 5){
    pout() << m_name + "::deallocateInternals()" << endl;
  }

  // TLDR: This deallocates a bunch of storage. This can be used during regrids to trim memory (because the Berger-Rigoutsous algorithm eats memory). 
  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_source);
  m_amr->deallocate(m_faceVelocity);
  m_amr->deallocate(m_cellVelocity);
  m_amr->deallocate(m_ebFlux);
  m_amr->deallocate(m_faceCenteredDiffusionCoefficient);
  m_amr->deallocate(m_ebCenteredDiffusionCoefficient);
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_scratchFluxOne);
  m_amr->deallocate(m_scratchFluxTwo);
  m_amr->deallocate(m_faceStates);
}

void CdrSolver::averageVelocityToFaces(){
  CH_TIME("CdrSolver::averageVelocityToFaces()");
  if(m_verbosity > 5){
    pout() << m_name + "::averageVelocityToFaces()" << endl;
  }

  this->averageVelocityToFaces(m_faceVelocity, m_cellVelocity); // Average velocities to face centers for all levels
}

void CdrSolver::averageVelocityToFaces(EBAMRFluxData& a_faceVelocity, const EBAMRCellData& a_cellVelocity){
  CH_TIME("CdrSolver::averageVelocityToFaces(EBAMRFluxData, EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << m_name + "::averageVelocityToFaces(EBAMRFluxData, EBAMRCellData)" << endl;
  }

  CH_assert(a_faceVelocity[0]->nComp() == 1       );
  CH_assert(a_cellVelocity[0]->nComp() == SpaceDim);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    DataOps::averageCellVectorToFaceScalar(*a_faceVelocity[lvl], *a_cellVelocity[lvl], m_amr->getDomains()[lvl]);
    
    a_faceVelocity[lvl]->exchange();
  }
}

void CdrSolver::preRegrid(const int a_lmin, const int a_oldFinestLevel){
  CH_TIME("CdrSolver::preRegrid(int, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::preRegrid(int, int)" << endl;
  }
  
  m_amr->allocate(m_cachePhi,    m_realm, m_phase, m_nComp);
  m_amr->allocate(m_cacheSource, m_realm, m_phase, m_nComp);
  
  for (int lvl = 0; lvl <= a_oldFinestLevel; lvl++){
    m_phi   [lvl]->localCopyTo(*m_cachePhi   [lvl]);
    m_source[lvl]->localCopyTo(*m_cacheSource[lvl]);
  }
}

void CdrSolver::coarseFineIncrement(const EBAMRIVData& a_massDifference){
  CH_TIME("CdrSolver::coarseFineIncrement(EBAMRIVData)");
  if(m_verbosity > 5){
    pout() << m_name + "::coarseFineIncrement(EBAMRIVData)" << endl;
  }

  CH_assert(a_massDifference[0]->nComp() == 1);

  // TLDR: This routine increments masses for the coarse-fine redistribution. In the algorithm we can have an EBCF situation
  //       where we 1) redistribute mass to ghost cells, 2) redistribute mass into the part under the fine grid. Fortunately,
  //       there are registers in Chombo that let us do this kind of redistribution magic. 

  const Interval interv(m_comp, m_comp);

  const int finestLevel = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finestLevel; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    RefCountedPtr<EBFineToCoarRedist>& fine2coarRedist = m_amr->getFineToCoarRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fineRedist = m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coarRedist = m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];

    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < finestLevel;

    if(hasCoar){
      fine2coarRedist->setToZero();

    }
    if(hasFine){
      coar2fineRedist->setToZero();
      coar2coarRedist->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      if(hasCoar){
	fine2coarRedist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
      }

      if(hasFine){
	coar2fineRedist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
	coar2coarRedist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
      }
    }
  }
}

void CdrSolver::coarseFineRedistribution(EBAMRCellData& a_phi){
  CH_TIME("CdrSolver::coarseFineRedistribution(EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << m_name + "::coarseFineRedistribution(EBAMRCellData)" << endl;
  }
  
  CH_assert(a_phi[0]->nComp() == 1);

  // TLDR: This routine does the coarse-fine redistribution. In the algorithm we can have an EBCF situation
  //       where we 1) redistribute mass to ghost cells, 2) redistribute mass into the part under the fine grid. Fortunately,
  //       there are registers in Chombo that let us do precisely this kind of redistribution magic. 

  const Interval interv(m_comp, m_comp);
  
  const int finestLevel = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finestLevel; lvl++){
    
    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < finestLevel;

    RefCountedPtr<EBCoarToFineRedist>& coar2fineRedist = m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coarRedist = m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBFineToCoarRedist>& fine2coarRedist = m_amr->getFineToCoarRedist(m_realm, m_phase)[lvl];
    
    if(hasCoar){
      fine2coarRedist->redistribute(*a_phi[lvl-1], interv);
      fine2coarRedist->setToZero();
    }

    if(hasFine){
      coar2fineRedist->redistribute(*a_phi[lvl+1], interv);
      coar2coarRedist->redistribute(*a_phi[lvl  ], interv);

      coar2fineRedist->setToZero();
      coar2coarRedist->setToZero();
    }
  }
}

void CdrSolver::computeDivG(EBAMRCellData& a_divG, EBAMRFluxData& a_G, const EBAMRIVData& a_ebFlux){
  CH_TIME("CdrSolver::computeDivG(EBAMRCellData, EBAMRFluxData, EBAMRIVData)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDivG(EBAMRCellData, EBAMRFluxData, EBAMRIVData)" << endl;
  }

  CH_assert(a_divG  [0]->nComp() == 1);
  CH_assert(a_G     [0]->nComp() == 1);
  CH_assert(a_ebFlux[0]->nComp() == 1);

  // TLDR: This routine computes a finite volume approximation to Div(G) where G is a flux, stored on face centers and eb faces. The routine uses
  //       flux matching on refinement boundaries and the so-called hybrid divergence in the cut-cells. The mass which is missed by the hybrid
  //       divergence is smooshed back in through Chombo's redistribution registers. 

  DataOps::setValue(a_divG, 0.0);

  // Interpolate fluxes to centroids and then coarsen them so we have consistent fluxes. 
  this->interpolateFluxToFaceCentroids(a_G);
  m_amr->averageDown(a_G, m_realm, m_phase);
  
  this->conservativeDivergenceNoKappaDivision(a_divG, a_G, a_ebFlux);      // Make the conservative divergence without kappa-division.
  this->nonConservativeDivergence(m_nonConservativeDivG, a_divG);          // Compute the non-conservative divergence
  this->hybridDivergence(a_divG, m_massDifference, m_nonConservativeDivG); // a_divG becomes hybrid divergence. Mass diff computed. 

  if(m_whichRedistribution != Redistribution::None){
    // Increment redistribution registers with the mass difference.
    this->incrementRedist    (m_massDifference);
    this->coarseFineIncrement(m_massDifference);

    // Redistribute mass on the level and across the coarse-fine boundaries.
    this->hyperbolicRedistribution(a_divG);      // Level redistribution. 
    this->coarseFineRedistribution(a_divG);      // Coarse-fine redistribution
  }
}

void CdrSolver::computeAdvectionFlux(EBAMRFluxData&       a_flux,
				     const EBAMRFluxData& a_facePhi,
				     const EBAMRFluxData& a_faceVelocity,
				     const bool           a_addDomainFlux){
  CH_TIME("CdrSolver::computeAdvectionFlux(EBAMRFluxData, EBAMRFluxData, EBAMRFluxData, bool)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAdvectionFlux(EBAMRFluxData, EBAMRFluxData, EBAMRFluxData, bool)" << endl;
  }

  CH_assert(a_flux        [0]->nComp() == 1);
  CH_assert(a_facePhi     [0]->nComp() == 1);
  CH_assert(a_faceVelocity[0]->nComp() == 1);

  // Computes the advection flux F = phi*v. This includes domain boundary faces. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    this->computeAdvectionFlux(*a_flux[lvl], *a_facePhi[lvl], *a_faceVelocity[lvl], lvl);
  }

  // Set domain faces to zero or enforce domain BCs (depending on flag)
  if(a_addDomainFlux){
    this->fillDomainFlux(a_flux);
  }
  else{
    this->resetDomainFlux(a_flux);
  }
}

void CdrSolver::computeAdvectionFlux(LevelData<EBFluxFAB>&       a_flux,
				     const LevelData<EBFluxFAB>& a_facePhi,
				     const LevelData<EBFluxFAB>& a_faceVelocity,
				     const int                   a_lvl){
  CH_TIME("CdrSolver::computeAdvectionFlux(LD<EBFluxFAB>, LD<EBFluxFAB>, LD<EBFluxFAB>, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAdvectionFlux(LD<EBFluxFAB>, LD<EBFluxFAB>, LD<EBFluxFAB>, int)" << endl;
  }
  
  CH_assert(a_flux.        nComp() == 1);
  CH_assert(a_facePhi.     nComp() == 1);
  CH_assert(a_faceVelocity.nComp() == 1);

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();

    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux      = a_flux        [dit()][dir];
      const EBFaceFAB& phi = a_facePhi     [dit()][dir];
      const EBFaceFAB& vel = a_faceVelocity[dit()][dir];

      flux.setVal(0.0, m_comp);
      flux += phi;
      flux *= vel;

      // Irregular faces
      const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;
      for (FaceIterator faceit(ebisbox.getIrregIVS(box), ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();
	
	flux(face, m_comp) = vel(face, m_comp)*phi(face, m_comp);
      }
    }
  }
}

void CdrSolver::computeDiffusionFlux(EBAMRFluxData& a_flux, const EBAMRCellData& a_phi, const bool a_addDomainFlux){
  CH_TIME("CdrSolver::computeDiffusionFlux(EBAMRFluxData, EBAMRCellData, bool) ");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDiffusionFlux(EBAMRFluxData, EBAMRCellData, bool)" << endl;
  }

  CH_assert(a_flux[0]->nComp() == 1);
  CH_assert(a_phi [0]->nComp() == 1);

  // Computes the diffusion flux F = D*Grad(phi). This includes domain boundary faces. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    this->computeDiffusionFlux(*a_flux[lvl], *a_phi[lvl], lvl);
  }

  // Set domain faces to zero or enforce domain BCs (depending on flag)  
  if(a_addDomainFlux){
    this->fillDomainFlux(a_flux);
  }
  else{
    this->resetDomainFlux(a_flux);
  }
}

void CdrSolver::computeDiffusionFlux(LevelData<EBFluxFAB>& a_flux, const LevelData<EBCellFAB>& a_phi, const int a_lvl){
  CH_TIME("CdrSolver::computeDiffusionFlux(LD<EBFluxFAB>, LD<EBCellFAB>, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDiffusionFlux(LD<EBFluxFAB>, LD<EBCellFAB>, int)" << endl;
  }

  CH_assert(a_flux.nComp() == 1);
  CH_assert(a_phi. nComp() == 1);

  // TLDR: This routine computes the diffusion flux F = D*Grad(phi) on face centers. Since this uses centered differencing
  //       and we don't have valid data outside the computational domain we only do the differencing for interior faces,
  //       setting the flux to zero on domain faces. 

  const Real dx                = m_amr->getDx()[a_lvl];
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];
  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const EBCellFAB& phi     = a_phi[dit()];
    const Box&       cellBox = dbl  [dit()];    
    const EBISBox&   ebisbox = ebisl[dit()];
    const EBGraph&   ebgraph = ebisbox.getEBGraph();

    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux      = a_flux[dit()][dir];
      const EBFaceFAB& dco = (*m_faceCenteredDiffusionCoefficient[a_lvl])[dit()][dir];

      flux.setVal(0.0);
      
      // Only want interior faces -- the domain flux will be set to zero (what else would it be...?)
      Box compBox = cellBox;
      compBox.grow(dir, 1);
      compBox &= domain;
      compBox.grow(dir, -1);
      compBox.surroundingNodes(dir);

      // Fortran kernel -- this computes f = D(face)*(phi(iv_high) - phi(iv-low))/dx
      BaseFab<Real>& regFlux      = flux.getSingleValuedFAB();
      const BaseFab<Real>& regPhi = phi. getSingleValuedFAB();
      const BaseFab<Real>& regDco = dco. getSingleValuedFAB();
      
      FORT_DFLUX_REG(CHF_FRA1(regFlux, m_comp),
		     CHF_CONST_FRA1(regPhi, m_comp),
		     CHF_CONST_FRA1(regDco, m_comp),
		     CHF_CONST_INT(dir),
		     CHF_CONST_REAL(dx),
		     CHF_BOX(compBox));


      // Irregular faces need to be redone. 
      for (FaceIterator faceit(ebisbox.getIrregIVS(cellBox), ebgraph, dir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();

	if(!face.isBoundary()){ 
	  const VolIndex hiVoF = face.getVoF(Side::Hi);
	  const VolIndex loVoF = face.getVoF(Side::Lo);

	  flux(face, m_comp) = dco(face,m_comp) * (phi(hiVoF,m_comp) - phi(loVoF, m_comp))/dx;
	}
      }
    }
  }
}

void CdrSolver::resetDomainFlux(EBAMRFluxData& a_flux){
  CH_TIME("CdrSolver::resetDomainFlux(EBAMRFluxData)");
  if(m_verbosity > 5){
    pout() << m_name + "::resetDomainFlux(EBAMRFluxData)" << endl;
  }

  CH_assert(a_flux[0]->nComp() == 1);

  // TLDR: This routine iterates through all faces in a_flux which are boundary faces and sets the flux there to zero. No
  //       Fortran kernel here because the box we iterate over has a smaller dimension (it's a slice). 

  constexpr Real zero = 0.0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box      cellBox = dbl  [dit()];
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      
      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& flux        = (*a_flux[lvl])[dit()][dir];
	BaseFab<Real>& regFlux = flux.getSingleValuedFAB();
	
	for (SideIterator sit; sit.ok(); ++sit){
	  const Side::LoHiSide side = sit();
	  
	  // Create a cell box which lies next to the domain 
	  Box boundaryCellBox;
	  boundaryCellBox  = adjCellBox(domain.domainBox(), dir, side, -1); 
	  boundaryCellBox &= cellBox;

	  // Do the regular cells -- note the shift to make sure that cell indices map to face indices. 
	  const int     sign  = (side == Side::Lo) ? 0 : 1;
	  const IntVect shift = sign*BASISV(dir);

	  for (BoxIterator bit(boundaryCellBox); bit.ok(); ++bit){
	    regFlux(bit() + shift, m_comp) = zero;
	  }

	  // Do the irregular cells
	  const IntVectSet irregIVS = ebisbox.getIrregIVS(boundaryCellBox);
	  for (FaceIterator faceIt(irregIVS, ebgraph, dir, FaceStop::AllBoundaryOnly); faceIt.ok(); ++faceIt){
	    flux(faceIt(), m_comp) = zero;
	  }
	}
      }
    }
  }
}

void CdrSolver::fillDomainFlux(EBAMRFluxData& a_flux){
  CH_TIME("CdrSolver::fillDomainFlux(EBAMRFluxData)");
  if(m_verbosity > 5){
    pout() << m_name + "::fillDomainFlux(EBAMRFluxData)" << endl;
  }

  CH_assert(a_flux[0]->nComp() == 1);
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    this->fillDomainFlux(*a_flux[lvl], lvl);
  }
}

void CdrSolver::fillDomainFlux(LevelData<EBFluxFAB>& a_flux, const int a_level) {
  CH_TIME("CdrSolver::fillDomainFlux(LD<EBFluxFAB>, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::fillDomainFlux(LD<EBFluxFAB>, int)" << endl;
  }

  CH_assert(a_flux.nComp() == 1);

  // TLDR: This iterates through all domain faces and sets the BC flux to either data-based, function-based, wall bc, or extrapolated outflow. This
  //       routine uses a face-iterator (because of BaseIFFAB<T>), but performance should be acceptable since we go through a very small number of faces. 

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  const ProblemDomain& domain  = m_amr->getDomains()[a_level];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const Box      cellBox = dbl[dit()];
    
    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux = a_flux[dit()][dir];

      // Iterate over domain cells. I'm REALLY not sure about if this is performant...
      for (SideIterator sit; sit.ok(); ++sit){
	const Side::LoHiSide side = sit();
	
	const CdrDomainBC::DomainSide    domainSide = std::make_pair(dir, side);
	const CdrDomainBC::BcType&       bcType     = m_domainBC.getBcType(domainSide);
	const CdrDomainBC::FluxFunction& bcFunction = m_domainBC.getBcFunction(domainSide);

	// Data-based BC holders -- needed if bcType is CdrDomainBC::BcType::DataBased. 
	const BaseIFFAB<Real>& dataBasedFlux = (*m_domainFlux[a_level])[dit()](dir, side);

	// Create a box which abuts the current domain side
	Box boundaryCellBox;
	boundaryCellBox  = adjCellBox(domain.domainBox(), dir, side, -1); 
	boundaryCellBox &= cellBox;

	// Use a face-iterator here to go through domain faces -- sort of have to do this because of (the odd) design decision to use a BaseIFFAB for holding
	// the domain fluxes.
	//
	// A FaceIterator may be inefficient, but we are REALLY going through a very small number of faces so there shouldn't be a performance hit here. 
	const IntVectSet ivs(boundaryCellBox);
	for (FaceIterator faceit(ivs, ebgraph, dir, FaceStop::AllBoundaryOnly); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  // Set the flux on the BC. This can be data-based, wall, function-based, outflow (extrapolated). 
	  switch(bcType){
	  case CdrDomainBC::BcType::DataBased:{
	    flux(face, m_comp) = dataBasedFlux(face, m_comp);
	    break;
	  }
	  case CdrDomainBC::BcType::Wall:{	    
	    flux(face, m_comp) = 0.0;
	    break;
	  }
	  case CdrDomainBC::BcType::Function:{	    	    
	    const RealVect pos = EBArith::getFaceLocation(face, m_amr->getDx()[a_level],m_amr->getProbLo());
	    flux(face, m_comp)   = -sign(sit()) * bcFunction(pos, m_time);
	    break;
	  }
	  case CdrDomainBC::BcType::Outflow:{
	    flux(face, m_comp) = 0.0;
	    
	    // Get the next interior face(s) and extrapolate from these. 
	    const VolIndex interiorVof = face.getVoF(flip(sit()));
	    
	    const std::vector<FaceIndex> neighborFaces = ebisbox.getFaces(interiorVof, dir, flip(sit())).stdVector();

	    if(neighborFaces.size() > 0){
	      Real sumArea = 0.0;
	      
	      for (const auto& f : neighborFaces){
		const Real areaFrac = ebisbox.areaFrac(f);
		
		sumArea            += areaFrac;
		flux(face, m_comp) += areaFrac * flux(f, m_comp);
		
	      }
	      
	      flux(face, m_comp) = sign(sit()) * std::max(0.0, sign(sit())*flux(face, m_comp))/sumArea;
	    }
	    
	    break;
	  }
	  case CdrDomainBC::BcType::Solver:{ // Don't do anything.
	    break; 
	  }
	  default:
	    MayDay::Error("CdrSolver::fillDomainFlux(LD<EBFluxFAB>, int) - trying to fill unsupported domain bc flux type");
	    break;
	  }
	}
      }
    }
  }
}

void CdrSolver::conservativeDivergenceNoKappaDivision(EBAMRCellData& a_conservativeDivergence, EBAMRFluxData& a_flux, const EBAMRIVData& a_ebFlux){
  CH_TIME("CdrSolver::conservativeDivergenceNoKappaDivision(EBAMRCellData, EBAMRFluxData, EBAMRIVData)");
  if(m_verbosity > 5){
    pout() << m_name + "::conservativeDivergenceNoKappaDivision(EBAMRCellData, EBAMRFluxData, EBAMRIVData)" << endl;
  }

  CH_assert(a_conservativeDivergence[0]->nComp() == 1);
  CH_assert(a_flux                  [0]->nComp() == 1);
  CH_assert(a_ebFlux                [0]->nComp() == 1);

  // This routine computes the "regular" finite volume conservative divergence Sum(fluxes) but does not divide by the volume fraction (which a naive implementation would).
  // Rather, this no-kappa-divided divergence is used for computing another divergence (a hybrid) divergence which represents a stable update to
  // the hyperbolic equation.
  //
  // This routine computes just that: a_conservativeDivergence = kappa*div(F) = Sum(fluxes).
  //
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    a_flux[lvl]->exchange();
    
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    this->conservativeDivergenceRegular(*a_conservativeDivergence[lvl], *a_flux[lvl], lvl);              // Compute kappa*div(F) in regular cells
    this->computeDivergenceIrregular(*a_conservativeDivergence[lvl], *a_flux[lvl], *a_ebFlux[lvl], lvl); // Recompute divergence on irregular cells

    a_conservativeDivergence[lvl]->exchange();
  }

  m_amr->averageDown(a_conservativeDivergence, m_realm, m_phase);
  m_amr->interpGhost(a_conservativeDivergence, m_realm, m_phase);
}

void CdrSolver::computeDivergenceIrregular(LevelData<EBCellFAB>&              a_divG,
					   const LevelData<EBFluxFAB>&        a_centroidFlux,
					   const LevelData<BaseIVFAB<Real> >& a_ebFlux,
					   const int                          a_lvl){
  CH_TIME("CdrSolver::computeDivergenceIrregular(LD<EBCellFAB>, LD<EBFluxFAB>, LD<BaseIVFAB<Real> >, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDivergenceIrregular(LD<EBCellFAB>, LD<EBFluxFAB>, LD<BaseIVFAB<Real> >, int)" << endl;
  }

  CH_assert(a_divG.        nComp() == 1);
  CH_assert(a_centroidFlux.nComp() == 1);
  CH_assert(a_ebFlux.      nComp() == 1);

  // This computes the conservative divergence kappa*div(F) = sum(fluxes) in cut-cells. Fluxes
  // must be face centroid centered!

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];
  const Real dx                = m_amr->getDx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const EBISBox& ebisbox = ebisl[dit()];

    EBCellFAB& divG               = a_divG[dit()];
    const BaseIVFAB<Real>& ebflux = a_ebFlux[dit()];

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof    = vofit();
      const Real      ebArea = ebisbox.bndryArea(vof);

      // Add the EB flux contribution to our sum(fluxes)
      divG(vof, m_comp) = ebflux(vof,m_comp)*ebArea;

      // Add face flux contributions to our sum(fluxes). 
      for (int dir = 0; dir < SpaceDim; dir++){
	const EBFaceFAB& flux = a_centroidFlux[dit()][dir];

	for (SideIterator sit; sit.ok(); ++sit){
	  const Side::LoHiSide side = sit();
	  
	  const int isign               = sign(side);
	  const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, side);

	  for (int iface = 0; iface < faces.size(); iface++){
	    const FaceIndex face     = faces[iface];
	    const Real      faceArea = ebisbox.areaFrac(face);
	    
	    divG(vof, m_comp) += isign*faceArea*flux(face, m_comp);
	  }
	}
      }

      // Scale divG by dx but not by kappa.
      divG(vof, m_comp) *= 1./dx;
    }
  }
}

void CdrSolver::conservativeDivergenceRegular(LevelData<EBCellFAB>& a_divJ, const LevelData<EBFluxFAB>& a_flux, const int a_lvl){
  CH_TIME("CdrSolver::conservativeDivergenceRegular(LD<EBCellFAB>, LD<EBFluxFAB>, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::conservativeDivergenceRegular(LD<EBCellFAB>, LD<EBFluxFAB>, int)" << endl;
  }

  CH_assert(a_divJ.nComp() == 1);
  CH_assert(a_flux.nComp() == 1);

  // We compute the conservative divergence in regular grid cells. The divergence is in the form
  // kappa*div(J) = sum(fluxes)/dx


  const DisjointBoxLayout dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain domain  = m_amr->getDomains()[a_lvl];
  const Real dx               = m_amr->getDx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box cellBox       = dbl.get(dit());
    
    EBCellFAB&     divJ     = a_divJ[dit()];
    BaseFab<Real>& divJReg  = divJ.getSingleValuedFAB();

    divJ.setVal(0.0);
    
    for (int dir = 0; dir < SpaceDim; dir++){
      const EBFaceFAB&     flux    = a_flux[dit()][dir];
      const BaseFab<Real>& fluxReg = flux.getSingleValuedFAB();

      // Fortran kernel -- this does divJ(iv) += (fluxReg(iv+1)-fluxReg(iv))/dx. 
      FORT_CONSDIV_REG(CHF_FRA1(divJReg, m_comp),
		       CHF_CONST_FRA1(fluxReg, m_comp),
		       CHF_CONST_INT(dir),
		       CHF_CONST_REAL(dx),
		       CHF_BOX(cellBox));
    }

    // Reset irregular grid cells -- these are set by interpolating fluxes to centroids and calling computeDivergenceIrregular(....)
    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      divJ(vof, m_comp) = 0.0;
    }
  }

  a_divJ.exchange();
}

void CdrSolver::defineInterpolationStencils(){
  CH_TIME("CdrSolver::defineInterpolationStencils()");
  if(m_verbosity > 5){
    pout() << m_name + "::defineInterpolationStencils()" << endl;
  }

  // This routine defines interpolation stencils for going from face-centered fluxes to face-centroid fluxes. 

  const int finestLevel    = m_amr->getFinestLevel();
  
  for (int dir = 0; dir < SpaceDim; dir++){
    (m_interpStencils[dir]).resize(1 + finestLevel);

    for (int lvl = 0; lvl <= finestLevel; lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const ProblemDomain& domain  = m_amr->getDomains()[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

      m_interpStencils[dir][lvl] = RefCountedPtr<LayoutData<BaseIFFAB<FaceStencil> > > (new LayoutData<BaseIFFAB<FaceStencil> >(dbl));

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	BaseIFFAB<FaceStencil>& sten = (*m_interpStencils[dir][lvl])[dit()];
	const Box        cellBox  = dbl  [dit()];
	const EBISBox&   ebisbox  = ebisl[dit()];
	const EBGraph&   ebgraph  = ebisbox.getEBGraph();
	const IntVectSet irregIVS = ebisbox.getIrregIVS(cellBox);

	sten.define(irregIVS, ebisbox.getEBGraph(), dir, m_nComp);
	
	for (FaceIterator faceIt(irregIVS, ebgraph, dir, FaceStop::SurroundingWithBoundary); faceIt.ok(); ++faceIt){
	  const FaceIndex& face = faceIt();

	  sten(face, m_comp) = EBArith::getInterpStencil(face, IntVectSet(), ebisbox, domain.domainBox());
	}
      }
    }
  }
}

void CdrSolver::incrementRedistFlux(){
  CH_TIME("CdrSolver::incrementRedistFlux()");
  if(m_verbosity > 5){
    pout() << m_name + "::incrementRedistFlux()" << endl;
  }

  // This routine lets the flux register know what is going on with the cut-cell mass redistribution
  // across the coarse-fine boundary. 

  const int finestLevel = m_amr->getFinestLevel();
  
  const Interval interv(m_comp, m_comp);

  for (int lvl = 0; lvl <= finestLevel; lvl++){
    const Real dx       = m_amr->getDx()[lvl];
    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < finestLevel;
    
    if(hasFine){
      RefCountedPtr<EBFluxRegister>&     fluxReg         = m_amr->getFluxRegister    (m_realm, m_phase)[lvl];
      RefCountedPtr<EBCoarToFineRedist>& coar2fineRedist = m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];
      RefCountedPtr<EBCoarToCoarRedist>& coar2coarRedist = m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];
      
      const Real scale = -dx;

      fluxReg->incrementRedistRegister(*coar2fineRedist, interv, scale);
      fluxReg->incrementRedistRegister(*coar2coarRedist, interv, scale);
    }
  }
}

void CdrSolver::initialData(){
  CH_TIME("CdrSolver::initialData()");
  if(m_verbosity > 5){
    pout() << m_name + "::initialData()" << endl;
  }

  // CdrSolver can be initialized in several way -- we can fill the initial data with analytic data from a mesh, or we can
  // deposit particles (with an NGP scheme) on the mesh. This function does both -- first particles if we have them and
  // then we increment with the function. 

  const bool depositParticles = m_species->getInitialParticles().length() > 0;

  DataOps::setValue(m_phi, 0.0);

  // Deposit particles if we have them. 
  if(depositParticles){
    this->initialDataParticles();
  }

  // Increment with function values if this is also called for. Note that this increments,
  // why is why initialDataParticles is called first!
  this->initialDataDistribution();

  m_amr->averageDown(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);
}

void CdrSolver::initialDataDistribution(){
  CH_TIME("CdrSolver::initialDataDistribution()");
  if(m_verbosity > 5){
    pout() << m_name + "::initialDataDistribution()" << endl;
  }

  // TLDR: We just run through every cell in the grid and increment by m_species->initialData
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const Real               dx     = m_amr->getDx()[lvl];
    const RealVect           probLo = m_amr->getProbLo();
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box        cellBox = dbl[dit()];
      const EBISBox&   ebisbox = ebisl[dit()];
      
      EBCellFAB&       phi     = (*m_phi[lvl])[dit()];
      BaseFab<Real>&   regPhi  = phi.getSingleValuedFAB();

      // Regular cells
      for (BoxIterator bit(cellBox); bit.ok(); ++bit){
	const IntVect iv   = bit();
	if(ebisbox.isRegular(iv)){	  
	  const RealVect pos = probLo + RealVect(iv)*dx + 0.5*dx*RealVect::Unit;

	  regPhi(iv, m_comp) = regPhi(iv, m_comp) + m_species->initialData(pos, m_time);
	}
      }

      // Irreg and multicells
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = EBArith::getVofLocation(vof, m_amr->getDx()[lvl]*RealVect::Unit, probLo);
	
	phi(vof, m_comp) = phi(vof, m_comp) + m_species->initialData(pos, m_time);
      }
    }
  }

  m_amr->averageDown(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);

  DataOps::setCoveredValue(m_phi, 0, 0.0);
}

void CdrSolver::initialDataParticles(){
  CH_TIME("CdrSolver::initialDataParticles");
  if(m_verbosity > 5){
    pout() << m_name + "::initialDataParticles" << endl;
  }

  // TLDR: This function deposits a list of initial particles on the mesh using
  //       an NGP scheme, ignoring conservation in cut-cells. Please refers to
  //       the particle API to see the details of the ParticleContainer<T> type and
  //       deposition methods.

  const List<Particle>& initialParticles = m_species->getInitialParticles();

  if(initialParticles.length() > 0){

    constexpr int pvrBuffer = 0; // No PVR buffer means that particles live in their natural grid cells. 

    // Make a ParticleContainer<T> and redistribute particles over the AMR hierarchy. 
    ParticleContainer<Particle> particles;
    m_amr->allocate(particles, pvrBuffer, m_realm);
    particles.addParticles(m_species->getInitialParticles());

    // This function will be called BEFORE initialDataFunction, we it is safe to set m_phi to zero
    DataOps::setValue(m_phi, 0.0);
  
    // Deposit onto mesh. 
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const RealVect           dx     = m_amr->getDx()[lvl]*RealVect::Unit;
      const RealVect           probLo = m_amr->getProbLo();
      const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    
      // 2. Deposit this levels particles. We use an NGP scheme to do this. 
      for (DataIterator dit(dbl); dit.ok(); ++dit){
	const Box      cellBox = dbl  [dit()];
	const EBISBox& ebisbox = ebisl[dit()];

	// Make the deposition object and put the particles on the grid. 
	const bool forceIrregNGP = true;
	EbParticleInterp interp(cellBox, ebisbox, dx, probLo, forceIrregNGP);
	
	interp.deposit(particles[lvl][dit()].listItems(), (*m_phi[lvl])[dit()].getFArrayBox(), DepositionType::NGP);
      }

#if CH_SPACEDIM==2 // Scale for 2D Cartesian. We do this because the 2D deposition object will normalize by 1/(dx*dx), but we want 1/(dx*dx*dx) in both 2D and 3D
      DataOps::scale(*m_phi[lvl], 1./dx[0]);
#endif      
    }
  }
}

void CdrSolver::hybridDivergence(EBAMRCellData&     a_hybridDivergence,
				 EBAMRIVData&       a_massDifference,
				 const EBAMRIVData& a_nonConservativeDivergence){
  CH_TIME("CdrSolver::hybridDivergence(EBAMRCellData, EBAMRIVData, EBAMRIVData)");
  if(m_verbosity > 5){
    pout() << m_name + "::hybridDivergence(EBAMRCellData, EBAMRIVData, EBAMRIVData)" << endl;
  }

  CH_assert(a_hybridDivergence         [0]->nComp() == 1);
  CH_assert(a_massDifference           [0]->nComp() == 1);
  CH_assert(a_nonConservativeDivergence[0]->nComp() == 1);

  // On input, a_hybridDivergence must contain kappa*div(F). We also have the non-conservative divergence computed, so we just compute
  // the regular hybrid divergence and the mass loss/gain. 

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    this->hybridDivergence(*a_hybridDivergence[lvl], *a_massDifference[lvl], *a_nonConservativeDivergence[lvl], lvl);
  }
}

void CdrSolver::hybridDivergence(LevelData<EBCellFAB>&              a_hybridDivergence,
				 LevelData<BaseIVFAB<Real> >&       a_massDifference,
				 const LevelData<BaseIVFAB<Real> >& a_nonConservativeDivergence,
				 const int                          a_lvl){
  CH_TIME("CdrSolver::hybridDivergence(LD<EBCellFAB>, LD<BaseIVFAB<Real> >, LD<BaseIVFAB<Real> >, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::hybridDivergence(LD<EBCellFAB>, LD<BaseIVFAB<Real> >, LD<BaseIVFAB<Real> >, int)" << endl;
  }

  CH_assert(a_hybridDivergence.         nComp() == 1);
  CH_assert(a_massDifference.           nComp() == 1);
  CH_assert(a_nonConservativeDivergence.nComp() == 1);  


  // TLDR: We want to compute a stable approximation to the divergence using Dh = kappa*Dc + (1-kappa)*Dnc where
  //       Dh is the hybrid divergence and Dc/Dnc are the conservative and non-conservative divergences. On the way in
  //       we had a_hybridDivergence = kappa*Dc. 
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];
    
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB&             divH   = a_hybridDivergence         [dit()];  // On input, this contains kappa*div(F)
    BaseIVFAB<Real>&       deltaM = a_massDifference           [dit()];
    const BaseIVFAB<Real>& divNC  = a_nonConservativeDivergence[dit()];  // On input, this contains the non-conservative divergence.
    
    const EBISBox& ebisbox = ebisl[dit()];

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real kappa    = ebisbox.volFrac(vof);
      const Real dc       = divH (vof, m_comp);
      const Real dnc      = divNC(vof, m_comp);

      // Note to self: deltaM = (1-kappa)*(dc - kappa*dnc) because dc was not divided by kappa,
      // which it would be otherwise. 
      divH  (vof, m_comp) = dc + (1-kappa)*dnc;          // On output, contains hybrid divergence
      deltaM(vof, m_comp) = (1-kappa)*(dc - kappa*dnc);  // On output, contains mass loss/gain. 
    }
  }
}

void CdrSolver::setRedistWeights(const EBAMRCellData& a_weights){
  CH_TIME("CdrSolver::setRedistWeights(EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << m_name + "::setRedistWeights(EBAMRCellData)" << endl;
  }

  CH_assert(a_weights[0]->nComp() == 1);

  // This function sets the weights for redistribution (the default is to use volume-weighted). This applies
  // to level-redistribution, fine-to-coar redistribution, coarse-to-fine redistribution, and the heinous
  // re-redistribution. 

  const int finestLevel = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finestLevel; lvl++){

    // Level redistribution
    EBLevelRedist& redist = *m_amr->getLevelRedist(m_realm, m_phase)[lvl];
    redist.resetWeights(*a_weights[lvl], m_comp);

    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < finestLevel;

    if(hasCoar){
      EBFineToCoarRedist& fine2coarRedist = *m_amr->getFineToCoarRedist(m_realm, m_phase)[lvl];
      fine2coarRedist.resetWeights(*a_weights[lvl-1], m_comp);
    }
    if(hasFine){
      EBCoarToCoarRedist& coar2coarRedist = *m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];
      EBCoarToFineRedist& coar2fineRedist = *m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];

      coar2coarRedist.resetWeights(*a_weights[lvl], m_comp);
      coar2fineRedist.resetWeights(*a_weights[lvl], m_comp);
    }
  }
}

void CdrSolver::hyperbolicRedistribution(EBAMRCellData& a_divF){
  CH_TIME("CdrSolver::hyberbolicRedistribution(EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << m_name + "::hyperbolicRedistribution(EBAMRCellData)" << endl;
  }

  CH_assert(a_divF[0]->nComp() == 1);

  // This function ONLY redistributes on a grid level. If you don't have EBxCF crossings, thats the end of the story but in
  // general we also need to account for mass that redistributes over refinement boundaries. 

  const Interval interv(m_comp, m_comp);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    EBLevelRedist& levelRedist = *(m_amr->getLevelRedist(m_realm, m_phase)[lvl]);
    
    levelRedist.redistribute(*a_divF[lvl], interv);
    levelRedist.setToZero();
  }
}

void CdrSolver::interpolateFluxToFaceCentroids(EBAMRFluxData& a_flux){
  CH_TIME("CdrSolver::interpolateFluxToFaceCentroids(EBAMRFluxData)");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateFluxToFaceCentroids(EBAMRFluxData)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    a_flux[lvl]->exchange();

    this->interpolateFluxToFaceCentroids(*a_flux[lvl], lvl);;
  }
}

void CdrSolver::interpolateFluxToFaceCentroids(LevelData<EBFluxFAB>& a_flux, const int a_lvl){
  CH_TIME("CdrSolver::interpolateFluxToFaceCentroids(LD<EBFluxFAB>, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateFluxToFaceCentroids(LD<EBFluxFAB>, int)" << endl;
  }

  CH_assert(a_flux.nComp() == 1);

  // TLDR: We are given face-centered fluxes which we want to put on face centroids. To do this we store face-centered
  //       fluxes on temporary storage and just interpolate them to face centroids.
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box        cellBox  = dbl.get(dit());
    const EBISBox&   ebisbox  = ebisl[dit()];
    const EBGraph&   ebgraph  = ebisbox.getEBGraph();
    const IntVectSet irregIVS = ebisbox.getIrregIVS(cellBox);

    const bool isRegular   = ebisbox.isRegular(cellBox);
    const bool isCovered   = ebisbox.isCovered(cellBox);
    const bool isIrregular = !isRegular && !isCovered;

    if(isIrregular){
      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& faceFlux = a_flux[dit()][dir];

	// Make a copy of faceFlux which we use as interpolant. EBFaceFAB
	// needs a cell box, so here you go. 
	EBFaceFAB clone(ebisbox, faceFlux.getCellRegion(), dir, m_nComp);
	clone.setVal(0.0);
	clone += faceFlux;

	// Compute face centroid flux on cut-cell face centroids. Since a_flux enforces boundary conditions
	// we include domain boundary cut-cell faces in the interpolation. 
	for (FaceIterator faceit(irregIVS, ebgraph, dir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit){
	  const FaceIndex& face   = faceit();
	  const FaceStencil& sten = (*m_interpStencils[dir][a_lvl])[dit()](face, m_comp);

	  faceFlux(face, m_comp) = 0.;

	  for (int i = 0; i < sten.size(); i++){
	    const FaceIndex& iface = sten.face(i);
	    const Real iweight     = sten.weight(i);
	  
	    faceFlux(face, m_comp) += iweight * clone(iface, m_comp);
	  }
	}
      }
    }
  }
}

void CdrSolver::resetFluxRegister(){
  CH_TIME("CdrSolver::resetFluxRegister()");
  if(m_verbosity > 5){
    pout() << m_name + "::resetFluxRegister()" << endl;
  }

  // Function just sets the flux register to zero (so we can increment it later). 
  const int finestLevel = m_amr->getFinestLevel();
  
  for (int lvl = 0; lvl <= finestLevel; lvl++){
    if(lvl < finestLevel){ // Because flux registers between level l and level l+1 live on level l. 
      RefCountedPtr<EBFluxRegister>& fluxReg = m_amr->getFluxRegister(m_realm, m_phase)[lvl];
      
      fluxReg->setToZero();
    }
  }
}

void CdrSolver::incrementFluxRegister(const EBAMRFluxData& a_flux){
  CH_TIME("CdrSolver::incrementFluxRegister(EBAMRFluxData)");
  if(m_verbosity > 5){
    pout() << m_name + "::incrementFluxRegister(EBAMRFluxData)" << endl;
  }

  CH_assert(a_flux[0]->nComp() == 1);

  // TLDR: To enforce flux matching on refinement boundaries we must have coarse flux = sum(fine fluxes). Fortunately
  //       Chombo has operators to handle that choreography, and this function provides the necessary fluxes to that operator.
  
  const int finestLevel = m_amr->getFinestLevel();
  
  const Interval interv(m_comp, m_comp);

  Vector<RefCountedPtr<EBFluxRegister> >& fluxReg = m_amr->getFluxRegister(m_realm, m_phase);
  
  for (int lvl = 0; lvl <= finestLevel; lvl++){
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    
    const bool hasCoar = lvl > 0;
    const bool hasFine = lvl < finestLevel;

    // Reset flux register first. 
    if(hasFine){
      fluxReg[lvl]->setToZero();
    }

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      for (int dir = 0; dir < SpaceDim; dir++){
	const Real       scale = 1.0;
	const EBFaceFAB& flux  = (*a_flux[lvl])[dit()][dir];

	// Increment flux register for irregular/regular. Add both from coarse to fine and from fine to coarse. 
	if(hasFine){ 
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxReg[lvl]->incrementCoarseBoth(flux, scale, dit(), interv, dir, sit()); // Register between lvl and lvl+1
	  }
	}
	
	if(hasCoar){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxReg[lvl-1]->incrementFineBoth(flux, scale, dit(), interv, dir, sit()); // Register between lvl-1 and lvl. 
	  }
	}
      }
    }
  }
}

void CdrSolver::incrementRedist(const EBAMRIVData& a_massDifference){
  CH_TIME("CdrSolver::incrementRedist(EBAMRIVData)");
  if(m_verbosity > 5){
    pout() << m_name + "::incrementRedist(EBAMRIVData)" << endl;
  }

  CH_assert(a_massDifference[0]->nComp() == 1);

  // a_massDifference is the mass to be smooshed back onto the grid through redistribution. Here, we increment the level-only
  // redistribution object with that mass. 

  const Interval interv(m_comp, m_comp);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    
    EBLevelRedist& levelRedist = *(m_amr->getLevelRedist(m_realm, m_phase)[lvl]);
    
    levelRedist.setToZero();

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      levelRedist.increment((*a_massDifference[lvl])[dit()], dit(), interv);
    }
  }
}

void CdrSolver::nonConservativeDivergence(EBAMRIVData& a_nonConservativeDivergence, const EBAMRCellData& a_divG){
  CH_TIME("CdrSolver::nonConservativeDivergence(EBAMRIVData, EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << m_name + "::nonConservativeDivergence(EBAMRIVData, EBAMRCellData)" << endl;
  }

  CH_assert(a_nonConservativeDivergence[0]->nComp() == 1);
  CH_assert(a_divG                     [0]->nComp() == 1);

  // TLDR: This routine computes the non-conservative divergence divNC(G) = sum(kappa*div(G))/sum(kappa). This is done through
  //       AmrMesh's stencil in cut cells. The neighborhood consists of cells that can be reached with a monotone path. 

  if(m_blendConservation){ 
    const IrregAmrStencil<NonConservativeDivergenceStencil>& stencils = m_amr->getNonConservativeDivergenceStencils(m_realm, m_phase);
    
    stencils.apply(a_nonConservativeDivergence, a_divG);
  }
  else{// Mostly for debugging. 
    DataOps::setValue(a_nonConservativeDivergence, 0.0);
  }
}

void CdrSolver::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("CdrSolver::regrid(int, int, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid(int, int, int)" << endl;
  }

  const Interval interv(m_comp, m_comp);

  // Allocate internal storage. This will fill m_phi with invalid data, but preRegrid(...) should have been called
  // prior to this routine so the old m_phi (before the regrid) is stored in m_cachedPhi. 
  this->allocateInternals();

  // These levels have not changed, but load balancing could have changed the processor ownership of data. 
  for (int lvl = 0; lvl <= Max(0,a_lmin-1); lvl++){
    m_cachePhi   [lvl]->copyTo(*m_phi   [lvl]); // Base level should never change, but ownership can.
    m_cacheSource[lvl]->copyTo(*m_source[lvl]); // Base level should never change, but ownership can.
  }

  // These levels have changed
  for (int lvl = Max(1,a_lmin); lvl <= a_newFinestLevel; lvl++){
    RefCountedPtr<EBPWLFineInterp>& interpolator = m_amr->getPwlInterpolator(m_realm, m_phase)[lvl]; // The interpolator for interpolator from lvl-1 to lvl lives on lvl

    // Linearly interpolate (with limiters) from level lvl-1 to lvl
    interpolator->interpolate(*m_phi   [lvl], *m_phi   [lvl-1], interv);
    interpolator->interpolate(*m_source[lvl], *m_source[lvl-1], interv);

    // There could be parts of the new grid that overlapped with the old grid (on level lvl) -- we don't want
    // to pollute the solution with interpolation there since we already have valid data. 
    if(lvl <= Min(a_oldFinestLevel, a_newFinestLevel)){
      m_cachePhi   [lvl]->copyTo(*m_phi[lvl]);
      m_cacheSource[lvl]->copyTo(*m_source[lvl]);
    }
  }

  m_amr->averageDown(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);

  m_amr->averageDown(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void CdrSolver::reflux(EBAMRCellData& a_phi){
  CH_TIME("CdrSolver::reflux(EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << m_name + "::reflux(EBAMRCellData)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);

  const int finestLevel = m_amr->getFinestLevel();
  
  const Interval interv(m_comp, m_comp);

  Vector<RefCountedPtr<EBFluxRegister > >& fluxreg = m_amr->getFluxRegister(m_realm, m_phase);

  for (int lvl = 0; lvl <= finestLevel; lvl++){
    if(lvl < finestLevel){
      RefCountedPtr<EBFluxRegister>& fluxReg = m_amr->getFluxRegister(m_realm, m_phase)[lvl]; // Flux matching between lvl and lvl-1 lives on lvl-1
      
      const Real dx    = m_amr->getDx()[lvl];
      const Real scale = 1.0/dx;

      fluxReg->reflux(*a_phi[lvl], interv, scale);
      fluxReg->setToZero();
    }
  }
}

void CdrSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("CdrSolver::setAmr(AmrMesh)");
  if(m_verbosity > 5){
    pout() << m_name + "::setAmr(AmrMesh)" << endl;
  }

  CH_assert(!a_amr.isNull());

  m_amr = a_amr;
}

void CdrSolver::registerOperators(){
  CH_TIME("CdrSolver::registerOperators()");
  if(m_verbosity > 5){
    pout() << m_name + "::registerOperators()" << endl;
  }

  CH_assert(!m_amr.isNull());

  m_amr->registerOperator(s_eb_coar_ave,     m_realm, m_phase);
  m_amr->registerOperator(s_eb_fill_patch,   m_realm, m_phase);
  m_amr->registerOperator(s_eb_pwl_interp,   m_realm, m_phase);
  m_amr->registerOperator(s_eb_flux_reg,     m_realm, m_phase);
  m_amr->registerOperator(s_eb_redist,       m_realm, m_phase);
  m_amr->registerOperator(s_eb_irreg_interp, m_realm, m_phase);
  m_amr->registerOperator(s_noncons_div,     m_realm, m_phase);
}

void CdrSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry){
  CH_TIME("CdrSolver::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)");
  if(m_verbosity > 5){
    pout() << m_name + "::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)" << endl;
  }

  CH_assert(!a_computationalGeometry.isNull());
  
  m_computationalGeometry = a_computationalGeometry;

  const RefCountedPtr<MultiFluidIndexSpace> MultiFluidIndexSpace = m_computationalGeometry->getMfIndexSpace();
  
  this->setEbIndexSpace(MultiFluidIndexSpace->getEBIndexSpace(m_phase));
}

void CdrSolver::setDiffusionCoefficient(const EBAMRFluxData& a_diffusionCoefficient, const EBAMRIVData& a_ebDiffusionCoefficient){
  CH_TIME("CdrSolver::setDiffusionCoefficient(EBAMRFluxData, EBAMRIVData)");
  if(m_verbosity > 5){
    pout() << m_name + "::setDiffusionCoefficient(EBAMRFluxData, EBAMRIVData)" << endl;
  }

  CH_assert(a_diffusionCoefficient  [0]->nComp() == 1);
  CH_assert(a_ebDiffusionCoefficient[0]->nComp() == 1);

  // Do a copy -- realms do not have to be the same. 
  m_faceCenteredDiffusionCoefficient.copy(a_diffusionCoefficient);
  m_ebCenteredDiffusionCoefficient.  copy(a_ebDiffusionCoefficient);

  m_amr->averageDown(m_faceCenteredDiffusionCoefficient, m_realm, m_phase);
  m_amr->averageDown(m_ebCenteredDiffusionCoefficient,   m_realm, m_phase);
}

void CdrSolver::setDiffusionCoefficient(const Real a_diffusionCoefficient){
  CH_TIME("CdrSolver::setDiffusionCoefficient(Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::setDiffusionCoefficient(Real)" << endl;
  }

  DataOps::setValue(m_faceCenteredDiffusionCoefficient, a_diffusionCoefficient);
  DataOps::setValue(m_ebCenteredDiffusionCoefficient,   a_diffusionCoefficient);

  m_amr->averageDown(m_faceCenteredDiffusionCoefficient, m_realm, m_phase);
  m_amr->averageDown(m_ebCenteredDiffusionCoefficient,   m_realm, m_phase);
}

void CdrSolver::setDiffusionCoefficient(const std::function<Real(const RealVect a_position)>& a_diffCo){
  CH_TIME("CdrSolver::setDiffusionCoefficient(std::function<Real(const RealVect a_position)>))");
  if(m_verbosity > 5){
    pout() << m_name + "::setDiffusionCoefficient(std::function<Real(const RealVect a_position)>)" << endl;
  }
  
  DataOps::setValue(m_faceCenteredDiffusionCoefficient, a_diffCo, m_amr->getProbLo(), m_amr->getDx(), m_comp);
  DataOps::setValue(m_ebCenteredDiffusionCoefficient,   a_diffCo, m_amr->getProbLo(), m_amr->getDx(), m_comp);
}

void CdrSolver::setEbFlux(const EBAMRIVData& a_ebFlux){
  CH_TIME("CdrSolver::setEbFlux(a_ebFlux)");
  if(m_verbosity > 5){
    pout() << m_name + "::setEbFlux(a_ebFlux)" << endl;
  }

  CH_assert(a_ebFlux[0]->nComp() == 1);

  m_ebFlux.copy(a_ebFlux);
}

void CdrSolver::setEbFlux(const Real a_ebFlux){
  CH_TIME("CdrSolver::setEbFlux(Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::setEbFlux(Real)" << endl;
  }

  DataOps::setValue(m_ebFlux, a_ebFlux);
}

void CdrSolver::setEbIndexSpace(const RefCountedPtr<EBIndexSpace>& a_ebis){
  CH_TIME("CdrSolver::setEbIndexSpace(RefCountedPtr<EBIndexSpace>)");
  if(m_verbosity > 5){
    pout() << m_name + "::setEbIndexSpace(RefCountedPtr<EBIndexSpace>)" << endl;
  }

  CH_assert(!a_ebis.isNull());

  m_ebis = a_ebis;
}

void CdrSolver::setSpecies(const RefCountedPtr<CdrSpecies>& a_species){
  CH_TIME("CdrSolver::setSpecies(RefCountedPtr<CdrSpecies>)");
  if(m_verbosity > 5){
    pout() << m_name + "::setSpecies(RefCountedPtr<CdrSpecies>)" << endl;
  }

  CH_assert(!a_species.isNull());

  m_species     = a_species;
  m_name        = m_species->getName();
  m_isDiffusive = m_species->isDiffusive();
  m_isMobile    = m_species->isMobile();
}

void CdrSolver::setSource(const EBAMRCellData& a_source){
  CH_TIME("CdrSolver::setSource(EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << m_name + "::setSource(EBAMRCellData)" << endl;
  }

  CH_assert(a_source[0]->nComp() == 1);

  m_source.copy(a_source);

  m_amr->averageDown(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void CdrSolver::setSource(const Real a_source){
  CH_TIME("CdrSolver::setSource(Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::setSource(Real)" << endl;
  }

  DataOps::setValue(m_source, a_source);

  m_amr->averageDown(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void CdrSolver::setSource(const std::function<Real(const RealVect a_position)> a_source){
  CH_TIME("CdrSolver::setSource(std::function<Real(const RealVect a_position)>)");
  if(m_verbosity > 5){
    pout() << m_name + "::setSource(std::function<Real(const RealVect a_position)>)" << endl;
  }
  
  DataOps::setValue(m_source, a_source, m_amr->getProbLo(), m_amr->getDx(), m_comp);

  m_amr->averageDown(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void CdrSolver::setTime(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("CdrSolver::setTime(int, Real, Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::setTime(int, Real, Real)" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;
}

void CdrSolver::setVelocity(const EBAMRCellData& a_velo){
  CH_TIME("CdrSolver::setVelocity");
  if(m_verbosity > 5){
    pout() << m_name + "::setVelocity" << endl;
  }

  CH_assert(a_velo[0]->nComp() == SpaceDim);

  m_cellVelocity.copy(a_velo);

  m_amr->averageDown(m_cellVelocity, m_realm, m_phase);
  m_amr->interpGhost(m_cellVelocity, m_realm, m_phase);
}

void CdrSolver::setVelocity(const RealVect a_velo){
  CH_TIME("CdrSolver::setVelocity(RealVect)");
  if(m_verbosity > 5){
    pout() << m_name + "::setVelocity(RealVect)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    for (int dir = 0; dir < SpaceDim; dir++){
      DataOps::setValue(*m_cellVelocity[lvl], a_velo[dir], dir);
    }
  }

  m_amr->averageDown(m_cellVelocity, m_realm, m_phase);
  m_amr->interpGhost(m_cellVelocity, m_realm, m_phase);
}

void CdrSolver::setVelocity(const std::function<RealVect(const RealVect a_pos)>& a_velo){
  CH_TIME("CdrSolver::setVelocity(std::function<RealVect(const RealVect a_pos)>)");
  if(m_verbosity > 5){
    pout() << m_name + "::setVelocity(std::function<RealVect(const RealVect a_pos)>)" << endl;
  }
  
  DataOps::setValue(m_cellVelocity, a_velo, m_amr->getProbLo(), m_amr->getDx());
}

void CdrSolver::setPhase(const phase::which_phase a_phase){
  CH_TIME("CdrSolver::setPhase(phase::which_phase)");
  if(m_verbosity > 5){
    pout() << m_name + "::setPhase(phase::which_phase)" << endl;
  }

  m_phase = a_phase;
}

void CdrSolver::setVerbosity(const int a_verbosity){
  CH_TIME("CdrSolver::setVerbosity(int)");
  m_verbosity = a_verbosity;
  
  if(m_verbosity > 5){
    pout() << m_name + "::setVerbosity(int)" << endl;
  }
}

void CdrSolver::writePlotFile(){
  CH_TIME("CdrSolver::writePlotFile()");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotFile()" << endl;
  }

  // We recycle the writePlotData(...) routine to write a plot file. That routine is used for external
  // calls to this solver when we want to populate an external plot file with data from this class. Nothing
  // wrong with reusing it here. 

  // Number of output components and their names
  const int numPlotVars                  = this->getNumberOfPlotVariables();
  const Vector<std::string> plotVarNames = this->getPlotVariableNames();

  // Allocate storage
  EBAMRCellData output;
  m_amr->allocate(output, m_realm, m_phase, numPlotVars);
  DataOps::setValue(output, 0.0);

  // Copy internal data to be plotted over to 'output'
  int icomp = 0;
  this->writePlotData(output, icomp);

  // Filename
  char filename[100];
  sprintf(filename, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_timeStep, SpaceDim);
  std::string fname(filename);
  
  // Alias, because Chombo's EBAMRIO wants raw pointers (but we stick to smart pointers, like God intended). 
  Vector<LevelData<EBCellFAB>* > outputPtr;
  m_amr->alias(outputPtr, output);

#ifdef CH_USE_HDF5
  writeEBHDF5(fname,
	      m_amr->getGrids(m_realm),
	      outputPtr,
	      plotVarNames,
	      m_amr->getDomains()[0].domainBox(),
	      m_amr->getDx()[0],
	      m_dt,
	      m_time,
	      m_amr->getRefinementRatios(),
	      m_amr->getFinestLevel() + 1,
	      false,
	      Vector<Real>(numPlotVars, 0.0), //coveredValues,
	      IntVect::Unit);
#endif
}

void CdrSolver::writePlotData(EBAMRCellData& a_output, int& a_icomp){
  CH_TIME("CdrSolver::writePlotData(EBAMRCellData, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotData(EBAMRCellData, int)" << endl;
  }

  // TLDR: This routine writes plot data to an "external" data holder a_output, starting on a_icomp. The intention
  //       behind that is that a_output is a plot file which accumulates data over many solvers. Note that because
  //       this solver is cell-centered, we interpolate m_phi to centroids when outputting it. 

  // Plot state
  if(m_plotPhi) {
    this->writeData(a_output, a_icomp, m_phi, true);
  }

  // Plot diffusion coefficients. These are stored on face centers but we need cell-centered data for output. 
  if(m_plotDiffusionCoefficient && m_isDiffusive) { 
    DataOps::setValue(m_scratch, 0.0);
    DataOps::averageFaceToCell(m_scratch, m_faceCenteredDiffusionCoefficient, m_amr->getDomains());
    this->writeData(a_output, a_icomp, m_scratch, false);
  }

  // Plot source terms
  if(m_plotSource) {
    this->writeData(a_output, a_icomp, m_source, false);
  }

  // Plot velocities
  if(m_plotVelocity && m_isMobile) {
    this->writeData(a_output, a_icomp, m_cellVelocity, false);
  }

  // Plot EB fluxes. These are stored on sparse data structures but we need them on cell centers. So copy them to scratch and write data. 
  if(m_plotEbFlux && m_isMobile){
    DataOps::setValue(m_scratch, 0.0);
    DataOps::incr(m_scratch, m_ebFlux, 1.0);
    this->writeData(a_output, a_icomp, m_scratch, false);
  }
}

void CdrSolver::writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp){
  CH_TIME("CdrSolver::writeData(EBAMRCellData, int, EBAMRCellData, bool)");
  if(m_verbosity > 5){
    pout() << m_name + "::writeData(EBAMRCellData, int, EBAMRCellData, bool)" << endl;
  }

  // TLDR: This routine copies from a data holder to another data, but in a general way where components don't align (which prevents
  //       us from using one-line methods). A special flag (a_interp) tells us to interpolate to cell centroids or not.

  // Number of components we are working with. 
  const int comp  = 0;
  const int ncomp = a_data[0]->nComp();

  const Interval srcInterv(0, ncomp-1);                 // This is the component range in scratch which we want to copy from
  const Interval dstInterv(a_comp, a_comp + ncomp - 1); // This is the component range in a_output which we want to copy into

  // Copy data to scratch and interpolate scratch to cell centroids if we are asked to. 
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, ncomp);
  DataOps::copy(scratch, a_data);

  if(a_interp){
    m_amr->interpToCentroids(scratch, m_realm, phase::gas);
  }

  m_amr->averageDown(scratch, m_realm, m_phase);
  m_amr->interpGhost(scratch, m_realm, m_phase);

  // Need a more general copy method because we can't call DataOps::copy (because realms might not be the same) and
  // we can't call EBAMRData<T>::copy either (because components don't align). So -- general type of copy here.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(m_realm == a_output.getRealm()){
      scratch[lvl]->localCopyTo(srcInterv, *a_output[lvl], dstInterv);
    }
    else{
      scratch[lvl]->copyTo(srcInterv, *a_output[lvl], dstInterv);
    }
  }

  DataOps::setCoveredValue(a_output, a_comp, 0.0);

  a_comp += ncomp;
}

#ifdef CH_USE_HDF5
void CdrSolver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("CdrSolver::writeCheckpointLevel(HDF5Handle, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::writeCheckpointLevel(HDF5Handle, int)" << endl;
  }

  // Write state vector
  write(a_handle, *m_phi   [a_level], m_name       );
  write(a_handle, *m_source[a_level], m_name+"_src");
}
#endif

#ifdef CH_USE_HDF5
void CdrSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level){
  CH_TIME("CdrSolver::readCheckpointLevel(HDF5Handle, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::readCheckpointLevel(HDF5Handle, int)" << endl;
  }
  
  const Interval interv(m_comp, m_comp);

  read<EBCellFAB>(a_handle, *m_phi   [a_level], m_name,        m_amr->getGrids(m_realm)[a_level], interv, false);
  read<EBCellFAB>(a_handle, *m_source[a_level], m_name+"_src", m_amr->getGrids(m_realm)[a_level], interv, false);
}
#endif

Real CdrSolver::computeAdvectionDt(){
  CH_TIME("CdrSolver::computeAdvectionDt()");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAdvectionDt()" << endl;
  }

  // TLDR: For advection we must have dt <= dx/(|vx|+|vy|+|vz|). E.g., with first order upwind phi^(k+1)_i = phi^k_i - (v*dt) * (phi^k_i - phi^k_(i-1))/dx so
  //       if phi^k_(i-1) == 0 then (1 - v*dt/dx) > 0.0 yields a positive definite solution (more general analysis certainly possible..)

  Real minDt = std::numeric_limits<Real>::max();

  if(m_isMobile){

    // We use m_scratch for holding dx/(|vx| + |vy| + |vz|)
    DataOps::setValue(m_scratch, minDt);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real dx                = m_amr->getDx()[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit){
	EBCellFAB& dt          = (*m_scratch     [lvl])[dit()];
	const EBCellFAB& velo  = (*m_cellVelocity[lvl])[dit()];
	const Box cellBox      = dbl.get(dit());

	// Regular cells -- the Fortran kernels computes dtReg = dx/(|vx|+|vy|+|vz|)
	BaseFab<Real>&       dtReg   = dt.getSingleValuedFAB();
	const BaseFab<Real>& veloReg = velo.getSingleValuedFAB();
	FORT_ADVECTION_DT(CHF_FRA1(dtReg, m_comp),
			  CHF_CONST_FRA(veloReg),
			  CHF_CONST_REAL(dx),
			  CHF_BOX(cellBox));


	// Same kernel as above, but for irregular cells. 
	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	
	  Real vel = 0.0;
	  for (int dir = 0; dir < SpaceDim; dir++){
	    vel += std::abs(velo(vof, dir));
	  }

	  dt(vof, m_comp) = dx/vel;
	}
      }
    }

    Real maxVal = std::numeric_limits<Real>::max();
    
    DataOps::setCoveredValue(m_scratch, m_comp, maxVal);  // Covered cells don't matter -- make sure they don't trigger anything (because MaxMin might reach into them). 
    DataOps::getMaxMin(maxVal, minDt, m_scratch, m_comp); // Get maximum and minimum. 
  }


  return minDt;
}

Real CdrSolver::computeDiffusionDt(){
  CH_TIME("CdrSolver::computeDiffusionDt()");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDiffusionDt()" << endl;
  }

  // TLDR: For advection we must have dt <= (dx*dx)/(2*d*D) where D is diffusion coefficient and d is spatial dimensions.
  //       E.g. in 1D, centered differencing yields
  //
  //          phi^(k+1)_i = phi^k_i - dt*D*(phi^k_(i+1) - 2*phi^k_i + phi^k_(i-1))/(dx*dx).
  //
  //       So for phi^k_(i-1) = phi^k_(i+1) = 0 we have phi^(k+1)_i = phi^k * (1 - 2*dt*D/(dx*dx)) which is positive definite only for dt < (dx*dx)/(2*D). Again,
  //       a more general analysis could be made.

  Real minDt = std::numeric_limits<Real>::max();

  if(m_isDiffusive){

    // We use scratch for holding dx*dx/(2*d*D)
    DataOps::setValue(m_scratch, minDt);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real dx                = m_amr->getDx()[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit){
	EBCellFAB&             dt         = (*m_scratch[lvl])[dit()];
	const Box              cellBox    = dbl  [dit()];
	const EBISBox&         ebisbox    = ebisl[dit()];
	const EBFluxFAB&       diffCoFace = (*m_faceCenteredDiffusionCoefficient[lvl])[dit()];
	const BaseIVFAB<Real>& diffcoEB   = (*m_ebCenteredDiffusionCoefficient[lvl])[dit()];

	// Fortran kernel. Strictly speaking, we should have a kernel which increments with the diffusion coefficients on each face since the finite volume
	// approximation to the Laplacian becomes (in 1D) Div*(D*Grad(phi)) = -D_(i-1/2)*(phi_i - phi_(i-1)) + D_(i+1/2)*(phi_(i+1)-phi_i). A good kernel
	// would do just that, but a lazy programmer just find the largest diffusion coefficient and uses that as an approximation. 
	for (int dir = 0; dir < SpaceDim; dir++){
	  BaseFab<Real>&       dtReg     = dt.getSingleValuedFAB();
	  const EBFaceFAB&     diffCo    = diffCoFace[dir];
	  const BaseFab<Real>& diffCoReg = diffCo.getSingleValuedFAB();
	  
	  FORT_DIFFUSION_DT(CHF_FRA1(dtReg, m_comp),
			    CHF_CONST_FRA1(diffCoReg, m_comp),
			    CHF_CONST_REAL(dx),
			    CHF_CONST_INT(dir),
			    CHF_BOX(cellBox));
	}

	// Same kernel as above, but including the EB face. 
	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex vof = vofit();

	  // Go through faces in the cut cell and fine the face which has the largest diffusion coefficient. We just
	  // use that as an approximation for computing dt. 
	  Real maxD  = diffcoEB(vof, m_comp);
	  
	  for (int dir = 0; dir < SpaceDim; dir++){
	    const EBFaceFAB& diffcoFace = diffCoFace[dir];
	    for (SideIterator sit; sit.ok(); ++sit){
	      const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());

	      for (int iface = 0; iface < faces.size(); iface++){
		const FaceIndex& face = faces[iface];
		
		maxD = std::max(maxD, diffcoFace(face, m_comp));
	      }
	    }
	  }

	  dt(vof, m_comp) = (dx*dx) / (2.0*SpaceDim*maxD);
	}
      }

      Real maxVal = std::numeric_limits<Real>::max();
    
      DataOps::setCoveredValue(m_scratch, m_comp, maxVal);  // Set covered to something large (in case MaxMin triggers covered cells)
      DataOps::getMaxMin(maxVal, minDt, m_scratch, m_comp);
    }
  }

  return minDt;
}

Real CdrSolver::computeAdvectionDiffusionDt(){
  CH_TIME("CdrSolver::computeAdvectionDiffusionDt()");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAdvectionDiffusionDt()" << endl;
  }

  // In 1D we have, e.g. d(phi)/dt = -d/dx(v*phi) + D*d^2(phi)/dx^2. Discretizing it with e.g. first order upwind and centered differencing yields
  //
  //    phi^(k+1)_i = phi^k_i - dt*v*[phi^k_i - phi^k_(i-1)]/dx + dt*D*[phi^k_(i+1) - 2*phi^k_i + phi^k_(i-1)]/(dx*dx).
  //
  // Setting phi_(i-1) = phi_(i+1) = 0 yields phi^(k+1)_i = phi^k_i - (dt*v/dx)*phi^k_i - 2*D*dt/(dx*dx)*phi_i^k. This is positive definite only if
  // 
  //    (1 - dt*v/dx - 2*D*dt/(dx*dx)) > 0.0
  //
  // which yields
  //
  //    dt <= 1/(v/dx + 2*D/(dx*dx).
  //
  // Generalization to 3D yields dt <= 1/[(|vx|+|vy|+|vz)/dx + 2*d*D/(dx*dx)]. The below kernels compute the factor in all grid cells, and we return the
  // smallest value. 
  

  Real minDt = std::numeric_limits<Real>::max();

  // Default to advection or diffusion time steps if the solver is only advective/diffusive. 
  if(m_isMobile && !m_isDiffusive){
    minDt = this->computeAdvectionDt();
  }
  else if(!m_isMobile && m_isDiffusive){
    minDt = this->computeDiffusionDt();
  }
  else if(m_isMobile && m_isDiffusive){

    // scratch is used for holding dt in the form above. 
    DataOps::setValue(m_scratch, 0.0);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real dx                = m_amr->getDx()[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit){
	EBCellFAB&             dt         = (*m_scratch                       [lvl])[dit()];
	const EBCellFAB&       velo       = (*m_cellVelocity                  [lvl])[dit()];
	const BaseIVFAB<Real>& diffCoEB   = (*m_ebCenteredDiffusionCoefficient[lvl])[dit()];
	const EBFluxFAB&       diffCoFace = (*m_faceCenteredDiffusionCoefficient[lvl])[dit()];
	
	const EBISBox& ebisbox          = ebisl[dit()];	
	const Box      cellBox          = dbl  [dit()];

	dt.setVal(0.0);

	// Regular faces. The first kernel finds the diffusion contribution only, using the maximum diffusion coefficient
	// on one of the faces. After the dir-loop dtReg will contain 2*d*D/(dx*dx) where D is the largest coefficient on one
	// of the faces of the cell. 
	BaseFab<Real>& dtReg = dt.getSingleValuedFAB();
	
	for (int dir = 0; dir < SpaceDim; dir++){
	  const BaseFab<Real>& diffCoReg = diffCoFace[dir].getSingleValuedFAB();
	  FORT_ADVECTION_DIFFUSION_DT_ONE(CHF_FRA1(dtReg, m_comp),
					  CHF_CONST_FRA1(diffCoReg, m_comp),
					  CHF_CONST_REAL(dx),
					  CHF_CONST_INT(dir),
					  CHF_BOX(cellBox));
	}

	// The second kernel increments with the advective contribution so that dtReg = (|Vx|+|Vy|+|Vz|)/dx + 2*d*D/(dx*dx)
	const BaseFab<Real>& veloReg = velo.getSingleValuedFAB();
	FORT_ADVECTION_DIFFUSION_DT_TWO(CHF_FRA1(dtReg, m_comp),  
					CHF_CONST_FRA(veloReg),
					CHF_CONST_REAL(dx),
					CHF_BOX(cellBox));

	// Invert the result so dtReg = 1/(|Vx|+|Vy|+|Vz|)/dx + 2*d*D/(dx*dx). 
	FORT_ADVECTION_DIFFUSION_DT_INVERT(CHF_FRA1(dtReg, m_comp),
					   CHF_BOX(cellBox));


	// Same kernel as above, but for cut-cells only. 
	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex vof = vofit();

	  // Compute |vx|+|vy|+|vz| in the cut-cell. 
	  Real vel = 0.0;
	  for (int dir = 0; dir < SpaceDim; dir++){
	    vel += std::abs(velo(vof, dir));
	  }

	  // Find largest diffusion coefficient in one of the faces of the cut-cell. 
	  Real maxD  = diffCoEB(vof, m_comp);
	  
	  for (int dir = 0; dir < SpaceDim; dir++){
	    const EBFaceFAB& diffCo = diffCoFace[dir];
	    
	    for (SideIterator sit; sit.ok(); ++sit){
	      const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());

	      for (int iface = 0; iface < faces.size(); iface++){
		const FaceIndex& face = faces[iface];
		
		maxD = std::max(maxD, diffCo(face, m_comp));
	      }
	    }
	  }

	  const Real inverseDtA  = vel/dx;
	  const Real inverseDtD  = (2*SpaceDim*maxD)/(dx*dx);

	  dt(vof, m_comp) = 1.0/(inverseDtA + inverseDtD);
	}
      }
    }

    Real maxVal = std::numeric_limits<Real>::max();
    
    DataOps::setCoveredValue(m_scratch, m_comp, maxVal);  // Covered cells are bogus -- set them to a large value for safety. 
    DataOps::getMaxMin(maxVal, minDt, m_scratch, m_comp);
  }

  return minDt;
}

Real CdrSolver::computeSourceDt(const Real a_max, const Real a_tolerance){
  CH_TIME("CdrSolver::computeSourceDt(Real, Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeSourceDt(Real, Real)" << endl;
  }

  // This routine computes the time step dt = phi/source but only for cells where phi lies within a_tolerance*a_max. It's an old
  // routine -- and one that is hardly ever used. 
  
  Real minDt = std::numeric_limits<Real>::max();

  if(a_max > 0.0){
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real dx                = m_amr->getDx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const EBCellFAB& phi     = (*m_phi   [lvl])[dit()];
	const EBCellFAB& source  = (*m_source[lvl])[dit()];
	const Box        cellBox = dbl             [dit()];

	// Fortran kernel -- computes dt = phi/source but only if phi > a_max*a_tolerance and source > 0.0
	const BaseFab<Real>& phiReg = phi.getSingleValuedFAB();
	const BaseFab<Real>& srcReg = source.getSingleValuedFAB();

	FORT_SOURCE_DT(CHF_REAL(minDt),
		       CHF_CONST_FRA1(phiReg, m_comp),
		       CHF_CONST_FRA1(srcReg, m_comp),
		       CHF_CONST_REAL(a_tolerance),
		       CHF_CONST_REAL(a_max),
		       CHF_BOX(cellBox));

	// Same as the Fortran kernel, but for cut-cells only. 
	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex vof = vofit();

	  const Real curPhi = phi   (vof, m_comp);
	  const Real curSrc = source(vof, m_comp);

	  Real curDt = std::numeric_limits<Real>::max();
	  
	  if(std::abs(curPhi) > a_tolerance*a_max && curSrc > 0.0){
	    curDt = std::abs(curPhi/curSrc);
	  }
	  
	  minDt = std::min(minDt, curDt);
	}
      }
    }
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&minDt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("CdrSolver::computeSourceDt() - communication error on norm");
  }
  minDt = tmp;
#endif
  
  return minDt;
}

void CdrSolver::weightedUpwind(EBAMRCellData& a_weightedUpwindPhi, const int a_pow) {
  CH_TIME("CdrSolver::weightedUpwind()");
  if(m_verbosity > 5){
    pout() << m_name + "::weightedUpwind()" << endl;
  }

  if(m_isMobile){
    m_amr->averageDown(m_cellVelocity, m_realm, m_phase);
    m_amr->interpGhost(m_cellVelocity, m_realm, m_phase);

    m_amr->averageDown(m_phi, m_realm, m_phase);
    m_amr->interpGhost(m_phi, m_realm, m_phase);        

    // Compute velocity on faces and EBs. The data is currently face-centered.
    this->averageVelocityToFaces(m_faceVelocity, m_cellVelocity);
    this->advectToFaces(m_faceStates, m_phi, 0.0); 

    DataOps::setValue(m_scratch, 0.0); // Used to store sum(alpha*v)

    // Grid loop. 
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){

      const DisjointBoxLayout& dbl    = m_amr->getGrids     (m_realm         )[lvl];
      const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const ProblemDomain&     domain = m_amr->getDomains()[lvl];
      const Real&              dx     = m_amr->getDx()[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit){
	const Box&     cellBox = dbl  [dit()];
	const EBISBox& ebisBox = ebisl[dit()];
	
	EBCellFAB&       sumPhi    = (*a_weightedUpwindPhi[lvl])[dit()];
	EBCellFAB&       sumWeight = (*m_scratch          [lvl])[dit()];
	const EBFluxFAB& facePhi   = (*m_faceStates       [lvl])[dit()];
	const EBFluxFAB& faceVel   = (*m_faceVelocity     [lvl])[dit()];

	sumPhi.   setVal(0.0);
	sumWeight.setVal(0.0);


	// Regular cells. Note that in the WEIGHTED_UPWIND Fortran kernel we compute the weighted sum of the upwinded value of phi. This means that the kernel
	// will reach out of the domain boundary. We fill the first ghost layer outside the domain boundary with the value in the valid cell immediately inside
	// the domain so that the kernel does, in fact, use the correct math.
	BaseFab<Real>& regSumPhi    = sumPhi.   getSingleValuedFAB();
	BaseFab<Real>& regSumWeight = sumWeight.getSingleValuedFAB();
	
	for (int dir = 0; dir < SpaceDim; dir++){
	  const BaseFab<Real>& regFacePhi = facePhi[dir].getSingleValuedFAB();
	  const BaseFab<Real>& regFaceVel = faceVel[dir].getSingleValuedFAB();

	  FORT_WEIGHTED_UPWIND(CHF_FRA1      (regSumPhi,    0),
			       CHF_FRA1      (regSumWeight, 0),
			       CHF_CONST_FRA1(regFacePhi,   0),
			       CHF_CONST_FRA1(regFaceVel,   0),
			       CHF_CONST_INT (dir),
			       CHF_CONST_INT (a_pow),
			       CHF_BOX       (cellBox));
	}

	// Irregular cells. This is a bit more involved. Note that we do want to compute everything at face centroids since
	// that is the flux that makes its way into the cells anyways when we advect. 
	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();

	  sumPhi   (vof, m_comp) = 0.0;
	  sumWeight(vof, m_comp) = 0.0;

	  for (int dir = 0; dir < SpaceDim; dir++){
	    const Vector<FaceIndex> facesLo = ebisBox.getFaces(vof, dir, Side::Lo);
	    const Vector<FaceIndex> facesHi = ebisBox.getFaces(vof, dir, Side::Hi);	    

	    // Add contribution from cut-cell faces in the low side. Observe that if the upwind value of phi is outside the domain
	    // then the "upwinded" value is just the cell-centered value. 
	    for (int iface = 0; iface < facesLo.size(); iface++){
	      const FaceIndex& faceLo   = facesLo[iface];
	      const VolIndex&  vofLo    = faceLo.getVoF(Side::Lo);
	      const Real       areaFrac = ebisBox.areaFrac(faceLo);

	      const FaceStencil& interpSten = (*m_interpStencils[dir][lvl])[dit()](faceLo, m_comp);

	      Real phiLo = 0.0;
	      Real velLo = 0.0;

	      for (int i = 0; i < interpSten.size(); i++){
		const FaceIndex& interpFace   = interpSten.face(i);
		const Real&      interpWeight = interpSten.weight(i);

		phiLo += interpWeight * facePhi[dir](interpFace, m_comp);
		velLo += interpWeight * faceVel[dir](interpFace, m_comp);		
	      }

	      if(velLo > 0.0){
		sumWeight(vof, m_comp) += areaFrac * std::abs(std::pow(velLo, a_pow))        ;
		sumPhi   (vof, m_comp) += areaFrac * std::abs(std::pow(velLo, a_pow)) * phiLo;
	      }
	    }

	    // Add contribution from cut-cell faces on the high side. 
	    for (int iface = 0; iface < facesHi.size(); iface++){
	      const FaceIndex& faceHi   = facesHi[iface];
	      const VolIndex&  vofHi    = faceHi.getVoF(Side::Hi);
	      const Real       areaFrac = ebisBox.areaFrac(faceHi);
	      
	      const FaceStencil& interpSten = (*m_interpStencils[dir][lvl])[dit()](faceHi, m_comp);

	      Real phiHi = 0.0;
	      Real velHi = 0.0;

	      for (int i = 0; i < interpSten.size(); i++){
		const FaceIndex& interpFace   = interpSten.face(i);
		const Real&      interpWeight = interpSten.weight(i);

		phiHi += interpWeight * facePhi[dir](interpFace, m_comp);
		velHi += interpWeight * faceVel[dir](interpFace, m_comp);		
	      }

	      if(velHi < 0.0){
		sumWeight(vof, m_comp) += areaFrac * std::abs(std::pow(velHi, a_pow))        ;
		sumPhi   (vof, m_comp) += areaFrac * std::abs(std::pow(velHi, a_pow)) * phiHi;
	      }
	    }
	  }

	  // If we don't have an inflow face set the upwind stuff to zero. 
	  if(!(sumWeight(vof, m_comp) > 0.0)) {
	    sumPhi(vof, m_comp)    = 0.0;
	    sumWeight(vof, m_comp) = 1.0;
	  }
	}
      }
    }

    // Divide. Set to zero m_phi if there are no inflow faces.
    EBAMRCellData zero;
    m_amr->allocate(zero, m_realm, m_phase, m_nComp);
    DataOps::setValue(zero, 0.0);
    DataOps::divideFallback(a_weightedUpwindPhi, m_scratch, zero);

    m_amr->averageDown(a_weightedUpwindPhi, m_realm, m_phase);
    m_amr->interpGhost(a_weightedUpwindPhi, m_realm, m_phase);
  }
  else{
    DataOps::copy(a_weightedUpwindPhi, m_phi);
  }  
}

Real CdrSolver::computeMass(){
  CH_TIME("CdrSolver::computeMass()");
  if(m_verbosity > 5){
    pout() << m_name + "::computeMass()" << endl;
  }

  // TLDR: We have conservative coarsening so we just coarsen the solution and compute the mass on the coarsest level only.
  
  m_amr->averageDown(m_phi, m_realm, m_phase);

  const Real dx  = m_amr->getDx()[0];

  Real mass = 0.;
  DataOps::kappaSum(mass, *m_phi[0]);

  mass *= pow(dx, SpaceDim);
  
  return mass;
}

Real CdrSolver::computeCharge(){
  CH_TIME("CdrSolver::computeCharge()");
  if(m_verbosity > 5){
    pout() << m_name + "::computeCharge()" << endl;
  }

  // TLDR: Compute the mass on the coarsest level and multiply by charge number. 

  const Real Q = this->computeMass()*m_species->getChargeNumber();

  return Q;
}

bool CdrSolver::isDiffusive(){
  CH_TIME("CdrSolver::isDiffusive()");
  if(m_verbosity > 5){
    pout() << m_name + "::isDiffusive()" << endl;
  }
  
  return m_isDiffusive;
}

bool CdrSolver::isMobile(){
  CH_TIME("CdrSolver::isMobile()");
  if(m_verbosity > 5){
    pout() << m_name + "::isMobile()" << endl;
  }
  
  return m_isMobile;
}

EBAMRCellData& CdrSolver::getPhi(){
  CH_TIME("CdrSolver::getPhi()");
  if(m_verbosity > 5){
    pout() << m_name + "::getPhi()" << endl;
  }
  
  return m_phi;
}

EBAMRCellData& CdrSolver::getSource(){
  CH_TIME("CdrSolver::getSource()");
  if(m_verbosity > 5){
    pout() << m_name + "::getSource()" << endl;
  }
  
  return m_source;
}

EBAMRCellData& CdrSolver::getCellCenteredVelocity(){
  CH_TIME("CdrSolver::getCellCenteredVelocity()");
  if(m_verbosity > 5){
    pout() << m_name + "::getCellCenteredVelocity()" << endl;
  }
  
  return m_cellVelocity;
}

EBAMRFluxData& CdrSolver::getFaceCenteredVelocity(){
  CH_TIME("CdrSolver::getFaceCenteredVelocity()");
  if(m_verbosity > 5){
    pout() << m_name + "::getFaceCenteredVelocity()" << endl;
  }
  
  return m_faceVelocity;
}

EBAMRIVData& CdrSolver::getEbCenteredVelocity(){
  CH_TIME("CdrSolver::getEbCenteredVelocity()");
  if(m_verbosity > 5){
    pout() << m_name + "::getEbCenteredVelocity()" << endl;
  }
  
  return m_ebVelocity;
}

EBAMRFluxData& CdrSolver::getFaceCenteredDiffusionCoefficient(){
  CH_TIME("CdrSolver::getFaceCenteredDiffusionCoefficient()");
  if(m_verbosity > 5){
    pout() << m_name + "::getFaceCenteredDiffusionCoefficient()" << endl;
  }
  
  return m_faceCenteredDiffusionCoefficient;
}

EBAMRIVData& CdrSolver::getEbCenteredDiffusionCoefficient(){
  CH_TIME("CdrSolver::getEbCenteredDiffusionCoefficient()");
  if(m_verbosity > 5){
    pout() << m_name + "::getEbCenteredDiffusionCoefficient()" << endl;
  }
  
  return m_ebCenteredDiffusionCoefficient;
}

EBAMRIVData& CdrSolver::getEbFlux(){
  CH_TIME("CdrSolver::getEbFlux()");
  if(m_verbosity > 5){
    pout() << m_name + "::getEbFlux()" << endl;
  }
  
  return m_ebFlux;
}

EBAMRIFData& CdrSolver::getDomainFlux(){
  CH_TIME("CdrSolver::getDomainFlux()");
  if(m_verbosity > 5){
    pout() << m_name + "::getDomainFlux()" << endl;
  }
  
  return m_domainFlux;
}

std::string CdrSolver::makeBcString(const int a_dir, const Side::LoHiSide a_side) const {
  CH_TIME("FieldSolver::makeBcString(int, Side::LoHiSide)");
  if(m_verbosity > 5){
    pout() << "FieldSolver::makeBcString(int, Side::LoHiSide)" << endl;
  }

  std::string strDir;
  std::string strSide;
  
  if(a_dir == 0){
    strDir = "x";
  }
  else if(a_dir == 1){
    strDir = "y";
  }
  else if(a_dir == 2){
    strDir = "z";
  }

  if(a_side == Side::Lo){
    strSide = "lo";
  }
  else if(a_side == Side::Hi){
    strSide = "hi";
  }

  const std::string ret = std::string("bc.") + strDir + std::string(".") + strSide;

  return ret;
}

void CdrSolver::parseDomainBc(){
  CH_TIME("CdrSolver::parseDomainBc()");
  if(m_verbosity > 5){
    pout() << m_name + "::parseDomainBc()" << endl;
  }
  
  ParmParse pp(m_className.c_str());

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Side::LoHiSide side = sit();
      
      const CdrDomainBC::DomainSide domainSide = std::make_pair(dir, side);
      const std::string             bcString   = this->makeBcString(dir, side);

      std::string str;
      pp.get(bcString.c_str(), str);

      // Set domain bc from input script evaluations.
      if(str == "data"){
	this->setDomainBcType(domainSide, CdrDomainBC::BcType::DataBased);
      }
      else if (str == "function"){
	this->setDomainBcType(domainSide, CdrDomainBC::BcType::Function);
      }
      else if (str == "wall"){
	this->setDomainBcType(domainSide, CdrDomainBC::BcType::Wall);
      }
      else if (str == "outflow"){
	this->setDomainBcType(domainSide, CdrDomainBC::BcType::Outflow);
      }
      else if (str == "solver"){
	this->setDomainBcType(domainSide, CdrDomainBC::BcType::Solver);
      }
      else{
	const std::string errorString = "CdrSolver::parseDomain BC - bad argument '" + bcString + "'";
	MayDay::Error(errorString.c_str());
      }
    }
  }
}


void CdrSolver::parseDivergenceComputation(){
  CH_TIME("CdrSolver::parseDivergenceComputation()");
  if(m_verbosity > 5){
    pout() << m_name + "::parseDivergenceComputation()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;

  // Get redistribution type. 
  pp.get("which_redistribution", str);
  if(str == "volume"){
    m_whichRedistribution = Redistribution::VolumeWeighted;
  }
  else if(str == "mass"){
    m_whichRedistribution = Redistribution::MassWeighted;
  }
  else if(str == "none"){
    m_whichRedistribution = Redistribution::None;
  }
  else{
    MayDay::Error("CdrSolver::parseDivergenceComputation -- logic bust. Must specify 'volume', 'mass', or 'none'");
  }

  // Blend conservation or not. 
  pp.get("blend_conservation",   m_blendConservation);
}

void CdrSolver::parsePlotVariables(){
  CH_TIME("CdrSolver::parsePlotVariables()");
  if(m_verbosity > 5){
    pout() << m_name + "::parsePlotVariables()" << endl;
  }

  m_plotPhi                  = false;
  m_plotVelocity             = false;
  m_plotDiffusionCoefficient = false;
  m_plotSource               = false;
  m_plotEbFlux               = false;

  // Parse plot variables from command line
  ParmParse pp(m_className.c_str());
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  // Set plot variables
  for (int i = 0; i < num; i++){
    if(     str[i] == "phi")    m_plotPhi                  = true;
    else if(str[i] == "vel")    m_plotVelocity             = true;
    else if(str[i] == "dco")    m_plotDiffusionCoefficient = true; 
    else if(str[i] == "src")    m_plotSource               = true;
    else if(str[i] == "ebflux") m_plotEbFlux               = true;
  }
}

void CdrSolver::gwnDiffusionSource(EBAMRCellData& a_noiseSource, const EBAMRCellData& a_cellPhi){
  CH_TIME("CdrSolver::gwnDiffusionSource(EBAMRCellData, EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << m_name + "::gwnDiffusionSource(EBAMRCellData, EBAMRCellData)" << endl;
  }

  CH_assert(a_noiseSource[0]->nComp() == 1);
  CH_assert(a_cellPhi    [0]->nComp() == 1);

  // This routine computes a Gaussian white noise source term for the fluctuating convection-diffusion equation. The source
  // term enters the convection-diffusion equation as
  //
  //    d(phi)/dt = Div(D*Grad(phi) - Z*sqrt(2*D*phi)) where Z is a Gaussian random field.
  //
  // In this routine we compute the finite volume approximation to Div(Z*sqrt(2*D*phi)) using standard finite volume routines. The only difference
  // is that the flux is defined as Z*sqrt(2*D*phi) but we don't have phi on face centers. Worse, setting phi to be the arithmetic average can easily
  // lead to negative density, so we use the Heaviside smoothing from Kim et. al. The paper is
  // "Stochastic Simulation of Reaction-Diffusion Systems: A Fluctuating-Hydrodynamics Approach", 10.1063/1.4978775.
  //
  // This does not guarantee against negative densities, but its surely better than the naive approach. Note that once we have the fluctuating term, we simply
  // compute Div(G). 

  if(m_isDiffusive){
    this->smoothHeavisideFaces(m_scratchFluxOne, a_cellPhi);                  // m_scratchFluxOne = phis on faces (smoothing as to avoid negative densities)
    DataOps::multiply(m_scratchFluxOne, m_faceCenteredDiffusionCoefficient);  // m_scratchFluxOne = D*phis
    DataOps::scale(m_scratchFluxOne, 2.0);                                    // m_scratchFluxOne = 2*D*phis
    DataOps::squareRoot(m_scratchFluxOne);                                    // m_scratchFluxOne = sqrt(2*D*phis)

#if 0 // Debug, check if we have negative face values
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      for (int dir = 0; dir <SpaceDim; dir++){
	Real max, min;
	EBLevelDataOps::getMaxMin(max, min, *m_scratchFluxOne[lvl], 0, dir);
	if(min < 0.0 || max < 0.0){
	  MayDay::Abort("CdrSolver::gwnDiffusionSource - negative face value");
	}
      }
    }
#endif

    this->fillGwn(m_scratchFluxTwo, 1.0);                         // Gaussian White Noise = W/sqrt(dV)
    DataOps::multiply(m_scratchFluxOne, m_scratchFluxTwo);        // Now m_scratchFluxOne holds the fluctuating cell-centered flux Z*sqrt(2*D*phi). 
    this->computeDivG(a_noiseSource, m_scratchFluxOne, m_ebZero); // Compute the finite volume approximation and make it into a source term. 
  
#if 0 // Debug, check if we have NaNs/Infs
    m_amr->averageDown(a_noiseSource, m_realm, m_phase);
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      if(EBLevelDataOps::checkNANINF(*a_noiseSource[lvl])){
	MayDay::Abort("CdrSolver::gwnDiffusionSource - something is wrong");
      }
    }
#endif
  }
  else{
    DataOps::setValue(a_noiseSource, 0.0);
  }
}

void CdrSolver::smoothHeavisideFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi){
  CH_TIME("CdrSolver::smoothHeavisideFaces(EBAMRFluxData, EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << m_name + "::smoothHeavisideFaces(EBAMRFluxData, EBAMRCellData)" << endl;
  }

  CH_assert(a_facePhi[0]->nComp() == 1);
  CH_assert(a_cellPhi[0]->nComp() == 1);

  // This routine is taken from  Kim et. al. "Stochastic Simulation of Reaction-Diffusion Systems: A Fluctuating-Hydrodynamics Approach" (10.1063/1.4978775).
  // It is inspired by a desire to gradually turn off fluctuations as the number of particles in a grid become small, as to avoid negative densities. So, we
  // compute the value of phi on faces with the following rules:
  //
  //    1. If theres more than one particle in the cells, we take the arithmetic average as usual. 
  //    2. If there's between zero and one particle in the cells, we use an averaging function phis = 0.5*(phiLo + phiHi) * loFactor * hiFactor
  //       where loFactor and hiFactor are the number of particles in the grid cell. 

  // Loop over levels. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    const Real vol               = pow(dx, SpaceDim);
    
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box&        cellBox  = dbl  [dit()];
      const EBISBox&    ebisbox  = ebisl[dit()];
      const EBGraph&    ebgraph  = ebisbox.getEBGraph();
      const IntVectSet& irregIVS = ebisbox.getIrregIVS(cellBox);      

      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB&       facePhi = (*a_facePhi[lvl])[dit()][dir];
	const EBCellFAB& cellPhi = (*a_cellPhi[lvl])[dit()];

	BaseFab<Real>&       faceReg = facePhi.getSingleValuedFAB();
	const BaseFab<Real>& cellReg = cellPhi.getSingleValuedFAB();

	// Regular faces (this will also do irregular faces but we fix those below)
	Box faceBox = cellBox;
	faceBox.surroundingNodes(dir);	
	FORT_HEAVISIDE_MEAN(CHF_FRA1(faceReg, m_comp),  
			    CHF_CONST_FRA1(cellReg, m_comp),
			    CHF_CONST_INT(dir),
			    CHF_CONST_REAL(dx),
			    CHF_BOX(faceBox));

	// Fix up irregular cell faces. This is the same kernel as above but irregular faces only...
	const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingNoBoundary;
	for (FaceIterator faceit(irregIVS, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  const VolIndex loVof = face.getVoF(Side::Lo);
	  const VolIndex hiVof = face.getVoF(Side::Hi);

	  const Real loVal = std::max(0.0, cellPhi(loVof, m_comp));
	  const Real hiVal = std::max(0.0, cellPhi(hiVof, m_comp));

	  Real Hlo;
	  if(loVal*vol <= 0.0){
	    Hlo = 0.0;
	  }
	  else if(loVal*vol >= 1.0){
	    Hlo = 1.0;
	  }
	  else{
	    Hlo = loVal*vol;
	  }

	  Real Hhi;
	  if(hiVal*vol <= 0.0){
	    Hhi = 0.0;
	  }
	  else if(hiVal*vol >= 1.0){
	    Hhi = 1.0;
	  }
	  else{
	    Hhi = hiVal*vol;
	  }

	  facePhi(face, m_comp) = 0.5*(hiVal + loVal)*Hlo*Hhi;
	}
      }
    }

    // Covered faces are bogus. 
    EBLevelDataOps::setCoveredVal(*a_facePhi[lvl], 0.0);
  }

  // No random flux on domain faces. 
  this->resetDomainFlux(a_facePhi);
}

void CdrSolver::fillGwn(EBAMRFluxData& a_noise, const Real a_sigma){
  CH_TIME("CdrSolver::fillGwn(EBAMRFluxData, Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::fillGwn(EBAMRFluxData, Real)" << endl;
  }

  CH_assert(a_noise[0]->nComp() == 1);

  // TLDR: This routine draws a random number from a Gaussian distribution at each cell face. This is used
  //       in FHD routines where we want to compute a finite volume approximation to the FHD diffusion noise
  //       term. 

  // Gaussian white noise distribution -- this should be centered at 0 and we will usually have a_sigma=1 
  std::normal_distribution<double> GWN(0.0, a_sigma);

  DataOps::setValue(a_noise, 0.0);
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    const Real vol               = pow(dx, SpaceDim);
    const Real ivol              = sqrt(1./vol);
    
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box&     cellBox  = dbl  [dit()];
      const EBISBox& ebisbox  = ebisl[dit()];
      const EBGraph& ebgraph  = ebisbox.getEBGraph();
      const IntVectSet& irreg = ebisbox.getIrregIVS(cellBox);      

      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& noise = (*a_noise[lvl])[dit()][dir];

	// Regular faces
	Box facebox = cellBox;
	facebox.surroundingNodes(dir);
	
	BaseFab<Real>& noiseReg = noise.getSingleValuedFAB();
	for (BoxIterator bit(facebox); bit.ok(); ++bit){
	  const IntVect iv = bit();
	  noiseReg(iv, m_comp) = GWN(m_rng)*ivol;
	}

	// Irregular faces
	const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingNoBoundary;
	for (FaceIterator faceit(irreg, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  noise(face, m_comp) = GWN(m_rng)*ivol;
	}
      }
    }
  }
}

void CdrSolver::parseRngSeed(){
  CH_TIME("CdrSolver::parseRngSeed()");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRngSeed()" << endl;
  }

  // Parse a random seed. If a seed < 0 is used we generate one from the system clock. 
  
  ParmParse pp(m_className.c_str());
  pp.get("seed", m_seed);
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();

  m_rng = std::mt19937_64(m_seed);
}

void CdrSolver::parsePlotMode(){
  CH_TIME("CdrSolver::parsePlotMode()");
  if(m_verbosity > 5){
    pout() << m_name + "::parsePlotMode()" << endl;
  }

  // Parses plot mode. 
  
  ParmParse pp(m_className.c_str());

  m_plotNumbers = false;
  
  std::string str;
  pp.get("plot_mode", str);
  if(str == "density"){
    m_plotNumbers = false;
  }
  else if(str == "numbers"){
    m_plotNumbers = true;
  }
}

#include <CD_NamespaceFooter.H>
