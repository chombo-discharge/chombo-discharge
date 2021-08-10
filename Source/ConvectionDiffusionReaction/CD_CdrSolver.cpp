/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrSolver.cpp
  @brief  Implementation of CD_CdrSolver.H
  @author Robert Marskar
  @todo   Need to parse domain BCs!
  @todo   Code walk starts on line 859 ::conservativeDivergenceRegular(...)
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
  m_name       = "CdrSolver";
  m_className  = "CdrSolver";

  this->setRealm(Realm::Primal); // Always primal as default.
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
  
  m_realm = a_realm;
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
    DataOps::averageCellToFace(*a_faceVelocity[lvl], *a_cellVelocity[lvl], m_amr->getDomains()[lvl]);
    
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
  
  this->conservativeDivergenceNoKappaDivision(a_divG, a_G, a_ebFlux);      // Make the conservative divergence without kappa-division.
  this->nonConservativeDivergence(m_nonConservativeDivG, a_divG);          // Compute the non-conservative divergence
  this->hybridDivergence(a_divG, m_massDifference, m_nonConservativeDivG); // a_divG becomes hybrid divergence. Mass diff computed. 
  this->incrementFluxRegister(a_G);                                        // Increment flux register
  this->incrementRedist(m_massDifference);                                 // Increment level redistribution register

  // Redistribution and reflux magic. 
  this->coarseFineIncrement(m_massDifference); // Compute C2F, F2C, and C2C mass transfers
  this->incrementRedistFlux();                 // Tell flux register about whats going on on refinement boundaries. 
  this->hyperbolicRedistribution(a_divG);      // Level redistribution. 
  this->coarseFineRedistribution(a_divG);      // Do the coarse-fine redistribution
  this->reflux(a_divG);                        // Reflux
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
	    VolIndex interiorVof = face.getVoF(flip(sit()));
	    Vector<FaceIndex> neighborFaces = ebisbox.getFaces(interiorVof, dir, flip(sit()));

	    if(neighborFaces.size() > 0){
	      for (const auto& f : neighborFaces.stdVector()){
		flux(face, m_comp) += flux(f, m_comp);
	      }
	      flux(face, m_comp) = std::max(0.0, sign(sit())*flux(face, m_comp))/neighborFaces.size();
	    }
	    
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
    this->interpolateFluxToFaceCentroids(*a_flux[lvl], lvl);                                             // Interpolate fluxes to centroids. 
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

  const bool depositFunction  = m_species->initializeWithFunction();
  const bool depositParticles = m_species->initializeWithParticles();

  DataOps::setValue(m_phi, 0.0);

  // Deposit particles if we have them. 
  if(depositParticles){
    this->initialDataParticles();
  }

  // Increment with function values if this is also called for. 
  if(depositFunction){
    this->initialDataDistribution();
  }

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
	const RealVect pos = probLo + RealVect(iv)*dx + 0.5*dx*RealVect::Unit;

	regPhi(iv, m_comp) = regPhi(iv, m_comp) + m_species->initialData(pos, m_time);
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

void CdrSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry> a_computationalGeometry){
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

void CdrSolver::setSpecies(const RefCountedPtr<CdrSpecies> a_species){
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
  CH_TIME("CdrSolver::writeData");
  if(m_verbosity > 5){
    pout() << m_name + "::writeData" << endl;
  }

  const int comp = 0;
  const int ncomp = a_data[0]->nComp();

  const Interval src_interv(0, ncomp-1);
  const Interval dst_interv(a_comp, a_comp + ncomp - 1);

  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, ncomp);
  DataOps::copy(scratch, a_data);

  if(a_interp){
    m_amr->interpToCentroids(scratch, m_realm, phase::gas);
  }

  m_amr->averageDown(scratch, m_realm, m_phase);
  m_amr->interpGhost(scratch, m_realm, m_phase);


  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(m_realm == a_output.getRealm()){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else{
      scratch[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }

  DataOps::setCoveredValue(a_output, a_comp, 0.0);

  a_comp += ncomp;
}

void CdrSolver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("CdrSolver::writeCheckpointLevel");
  if(m_verbosity > 5){
    pout() << m_name + "::writeCheckpointLevel" << endl;
  }

  // Write state vector
  write(a_handle, *m_phi[a_level], m_name);
  write(a_handle, *m_source[a_level],   m_name+"_src");
}

void CdrSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level){
  CH_TIME("CdrSolver::readCheckpointLevel");
  if(m_verbosity > 5){
    pout() << m_name + "::readCheckpointLevel" << endl;
  }

  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name, m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);
  read<EBCellFAB>(a_handle, *m_source[a_level], m_name+"_src", m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);
}

Real CdrSolver::computeAdvectionDt(){
  CH_TIME("CdrSolver::computeAdvectionDt");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAdvectionDt" << endl;
  }

  Real min_dt = std::numeric_limits<Real>::max();

  if(m_isMobile){
    const int comp = 0;

    DataOps::setValue(m_scratch, min_dt);
    
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real dx                = m_amr->getDx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	EBCellFAB& dt          = (*m_scratch[lvl])[dit()];
	const EBCellFAB& velo  = (*m_cellVelocity[lvl])[dit()];
	const Box box          = dbl.get(dit());

	BaseFab<Real>& dt_fab         = dt.getSingleValuedFAB();
	const BaseFab<Real>& velo_fab = velo.getSingleValuedFAB();
	FORT_ADVECTION_DT(CHF_FRA1(dt_fab, comp),
			  CHF_CONST_FRA(velo_fab),
			  CHF_CONST_REAL(dx),
			  CHF_BOX(box));


	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	
	  Real vel = 0.0;
	  for (int dir = 0; dir < SpaceDim; dir++){
	    vel += Abs(velo(vof, dir));
	  }

	  dt(vof, comp) = dx/vel;
	}
      }
    }

    Real maxVal = std::numeric_limits<Real>::max();
    
    DataOps::setCoveredValue(m_scratch, comp, maxVal);
    DataOps::getMaxMin(maxVal, min_dt, m_scratch, comp);
  }


  return min_dt;
}

Real CdrSolver::computeDiffusionDt(){
  CH_TIME("CdrSolver::computeDiffusionDt");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDiffusionDt" << endl;
  }

  Real min_dt = std::numeric_limits<Real>::max();

  if(m_isDiffusive){
    const int comp  = 0;

    DataOps::setValue(m_scratch, min_dt);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real dx                = m_amr->getDx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	EBCellFAB& dt                   = (*m_scratch[lvl])[dit()];
	const Box box                   = dbl.get(dit());
	const EBISBox& ebisbox          = ebisl[dit()];
	const BaseIVFAB<Real>& diffcoEB = (*m_ebCenteredDiffusionCoefficient[lvl])[dit()];

	// Regular faces
	for (int dir = 0; dir < SpaceDim; dir++){
	  BaseFab<Real>& dt_fab           = dt.getSingleValuedFAB();
	  const EBFaceFAB& diffco         = (*m_faceCenteredDiffusionCoefficient[lvl])[dit()][dir];
	  const BaseFab<Real>& diffco_fab = diffco.getSingleValuedFAB();
	  
	  FORT_DIFFUSION_DT(CHF_FRA1(dt_fab, comp),
			    CHF_CONST_FRA1(diffco_fab, comp),
			    CHF_CONST_REAL(dx),
			    CHF_CONST_INT(dir),
			    CHF_BOX(box));
	}

	// Irregular faces. 
	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex vof = vofit();

	  Real irregD  = diffcoEB(vof, comp);
	  
	  for (int dir = 0; dir < SpaceDim; dir++){
	    const EBFaceFAB& diffcoFace = (*m_faceCenteredDiffusionCoefficient[lvl])[dit()][dir];
	    for (SideIterator sit; sit.ok(); ++sit){
	      const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());

	      for (int iface = 0; iface < faces.size(); iface++){
		const FaceIndex& face = faces[iface];
		
		irregD = std::max(irregD, diffcoFace(face, comp));
	      }
	    }
	  }

	  dt(vof, comp) = dx*dx/(2.0*SpaceDim*irregD);
	}
      }

      Real maxVal = std::numeric_limits<Real>::max();
    
      DataOps::setCoveredValue(m_scratch, comp, maxVal);
      DataOps::getMaxMin(maxVal, min_dt, m_scratch, comp);
    }
  }

  return min_dt;
}

Real CdrSolver::computeAdvectionDiffusionDt(){
  CH_TIME("CdrSolver::computeAdvectionDiffusionDt");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAdvectionDiffusionDt" << endl;
  }

  Real min_dt = std::numeric_limits<Real>::max();

  if(m_isMobile && !m_isDiffusive){
    min_dt = this->computeAdvectionDt();
  }
  else if(!m_isMobile && m_isDiffusive){
    min_dt = this->computeDiffusionDt();
  }
  else if(m_isMobile && m_isDiffusive){
    const int comp  = 0;

    DataOps::setValue(m_scratch, 0.0);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real dx                = m_amr->getDx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	EBCellFAB& dt                   = (*m_scratch[lvl])[dit()];
	const EBISBox& ebisbox          = ebisl[dit()];
	const EBCellFAB& velo           = (*m_cellVelocity[lvl])[dit()];
	const BaseIVFAB<Real>& diffcoEB = (*m_ebCenteredDiffusionCoefficient[lvl])[dit()];
	const Box box                   = dbl.get(dit());

	dt.setVal(0.0);

	// Regular cells. First compute dt = 2*D*d/(dx*dx)
	BaseFab<Real>& dt_fab         = dt.getSingleValuedFAB();

	// Regular faces
	for (int dir = 0; dir < SpaceDim; dir++){
	  const EBFaceFAB& diffco         = (*m_faceCenteredDiffusionCoefficient[lvl])[dit()][dir];
	  const BaseFab<Real>& diffco_fab = diffco.getSingleValuedFAB();
	  FORT_ADVECTION_DIFFUSION_DT_ONE(CHF_FRA1(dt_fab, comp),
					  CHF_CONST_FRA1(diffco_fab, comp),
					  CHF_CONST_REAL(dx),
					  CHF_CONST_INT(dir),
					  CHF_BOX(box));
	}

	// Add the advective contribution so that dt_fab = (|Vx|+|Vy|+|Vz|)/dx + 2*d*D/(dx*dx)
	const BaseFab<Real>& velo_fab = velo.getSingleValuedFAB();
	FORT_ADVECTION_DIFFUSION_DT_TWO(CHF_FRA1(dt_fab, comp),  
					CHF_CONST_FRA(velo_fab),
					CHF_CONST_REAL(dx),
					CHF_BOX(box));

	// Invert the result. 
	FORT_ADVECTION_DIFFUSION_DT_INVERT(CHF_FRA1(dt_fab, comp),  
					   CHF_BOX(box));


	// Irregular cells
	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex vof = vofit();

	  // Compute advective velocity. 
	  Real vel = 0.0;
	  for (int dir = 0; dir < SpaceDim; dir++){
	    vel += std::abs(velo(vof, dir));
	  }

	  // Find largest diffusion coefficient. 
	  Real irregD  = diffcoEB(vof, comp);
	  for (int dir = 0; dir < SpaceDim; dir++){
	    const EBFaceFAB& diffcoFace = (*m_faceCenteredDiffusionCoefficient[lvl])[dit()][dir];
	    for (SideIterator sit; sit.ok(); ++sit){
	      const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());

	      for (int iface = 0; iface < faces.size(); iface++){
		const FaceIndex& face = faces[iface];
		
		irregD = std::max(irregD, diffcoFace(face, comp));
	      }
	    }
	  }

	  const Real idtA  = vel/dx;
	  const Real idtD  = (2*SpaceDim*irregD)/(dx*dx);

	  dt(vof, comp) = 1.0/(idtA + idtD);
	}
      }
    }

    Real maxVal = std::numeric_limits<Real>::max();
    
    DataOps::setCoveredValue(m_scratch, comp, maxVal);  // Covered cells are bogus. 
    DataOps::getMaxMin(maxVal, min_dt, m_scratch, comp);
  }


  return min_dt;
}

Real CdrSolver::computeSourceDt(const Real a_max, const Real a_tolerance){
  CH_TIME("CdrSolver::computeSourceDt");
  if(m_verbosity > 5){
    pout() << m_name + "::computeSourceDt" << endl;
  }

  const int comp = 0;
  
  Real min_dt = 1.E99;

  if(a_max > 0.0){
    const int finest_level = m_amr->getFinestLevel();
    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real dx                = m_amr->getDx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const EBCellFAB& state  = (*m_phi[lvl])[dit()];
	const EBCellFAB& source = (*m_source[lvl])[dit()];
	const Box box           = dbl.get(dit());

	const BaseFab<Real>& state_fab  = state.getSingleValuedFAB();
	const BaseFab<Real>& source_fab = source.getSingleValuedFAB();


	FORT_SOURCE_DT(CHF_REAL(min_dt),
		       CHF_CONST_FRA1(state_fab, comp),
		       CHF_CONST_FRA1(source_fab, comp),
		       CHF_CONST_REAL(a_tolerance),
		       CHF_CONST_REAL(a_max),
		       CHF_BOX(box));

	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
	for (vofit.reset(); vofit.ok(); ++vofit){
	  const VolIndex vof = vofit();

	  const Real phi = state(vof, comp);
	  const Real src = source(vof, comp);

	  Real thisdt = 1.E99;
	  if(Abs(phi) > a_tolerance*a_max && src > 0.0){
	    thisdt = Abs(phi/src);
	  }
	  min_dt = Min(min_dt, thisdt);
	}
      }
    }
  }

#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&min_dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("CdrSolver::computeSourceDt() - communication error on norm");
  }
  min_dt = tmp;
#endif
  
  return min_dt;

}

Real CdrSolver::computeMass(){
  CH_TIME("CdrSolver::computeMass");
  if(m_verbosity > 5){
    pout() << m_name + "::computeMass" << endl;
  }

  Real mass = 0.;
  const int base = 0;
  const Real dx = m_amr->getDx()[base];
  m_amr->averageDown(m_phi, m_realm, m_phase);
  
  DataOps::kappaSum(mass, *m_phi[base]);
  mass *= pow(dx, SpaceDim);
  
  return mass;
}

Real CdrSolver::computeCharge(){
  CH_TIME("CdrSolver::computeCharge");
  if(m_verbosity > 5){
    pout() << m_name + "::computeCharge" << endl;
  }

  const Real Q = this->computeMass()*m_species->getChargeNumber();

  return Q;
}

bool CdrSolver::extrapolateSourceTerm() const {
  return m_extrapolateSourceTerm;
}

bool CdrSolver::isDiffusive(){
  return m_isDiffusive;
}

bool CdrSolver::isMobile(){
  return m_isMobile;
}

EBAMRCellData& CdrSolver::getPhi(){
  return m_phi;
}

EBAMRCellData& CdrSolver::getSource(){
  return m_source;
}

EBAMRCellData& CdrSolver::getCellCenteredVelocity(){
  return m_cellVelocity;
}

EBAMRFluxData& CdrSolver::getFaceCenteredVelocity(){
  return m_faceVelocity;
}

EBAMRIVData& CdrSolver::getEbCenteredVelocity(){
  return m_ebVelocity;
}

EBAMRFluxData& CdrSolver::getFaceCenteredDiffusionCoefficient(){
  return m_faceCenteredDiffusionCoefficient;
}

EBAMRIVData& CdrSolver::getEbCenteredDiffusionCoefficient(){
  return m_ebCenteredDiffusionCoefficient;
}

EBAMRIVData& CdrSolver::getEbFlux(){
  return m_ebFlux;
}

EBAMRIFData& CdrSolver::getDomainFlux(){
  return m_domainFlux;
}

void CdrSolver::parseDomainBc(){
  ParmParse pp(m_className.c_str());

}

void CdrSolver::parseExtrapolateSourceTerm(){
  CH_TIME("CdrSolver::parseExtrapolateSourceTerm");
  if(m_verbosity > 5){
    pout() << m_name + "::parseExtrapolateSourceTerm" << endl;
  }
  
  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("extrapolateSourceTerm", str);
  m_extrapolateSourceTerm = (str == "true") ? true : false;
}

void CdrSolver::parseDivergenceComputation(){
  CH_TIME("CdrSolver::parseDivergenceComputation");
  if(m_verbosity > 5){
    pout() << m_name + "::parseDivergenceComputation" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("redist_mass_weighted", m_useMassWeightedRedistribution);
  pp.get("blend_conservation",   m_blendConservation);
}

void CdrSolver::parsePlotVariables(){
  ParmParse pp(m_className.c_str());
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  m_plotPhi                  = false;
  m_plotVelocity             = false;
  m_plotDiffusionCoefficient = false;
  m_plotSource               = false;
  m_plotEbFlux               = false;
  
  for (int i = 0; i < num; i++){
    if(     str[i] == "phi")    m_plotPhi                  = true;
    else if(str[i] == "vel")    m_plotVelocity             = true;
    else if(str[i] == "dco")    m_plotDiffusionCoefficient = true; 
    else if(str[i] == "src")    m_plotSource               = true;
    else if(str[i] == "ebflux") m_plotEbFlux               = true;
  }
}

void CdrSolver::gwnDiffusionSource(EBAMRCellData& a_noiseSource, const EBAMRCellData& a_cellPhi){
  CH_TIME("CdrSolver::gwnDiffusionSource");
  if(m_verbosity > 5){
    pout() << m_name + "::gwnDiffusionSource" << endl;
  }

  if(m_isDiffusive){
    const int comp  = 0;
    const int ncomp = 1;

    this->fillGwn(m_scratchFluxTwo, 1.0);                         // Gaussian White Noise = W/sqrt(dV)
    this->smoothHeavisideFaces(m_scratchFluxOne, a_cellPhi); // m_scratchFluxOne = phis
    DataOps::multiply(m_scratchFluxOne, m_faceCenteredDiffusionCoefficient);                // m_scratchFluxOne = D*phis
    DataOps::scale(m_scratchFluxOne, 2.0);                        // m_scratchFluxOne = 2*D*phis
    DataOps::squareRoot(m_scratchFluxOne);                       // m_scratchFluxOne = sqrt(2*D*phis)

#if 0 // Debug
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
  
    DataOps::multiply(m_scratchFluxOne, m_scratchFluxTwo);                     // Holds random, cell-centered flux
    this->conservativeDivergenceNoKappaDivision(a_noiseSource, m_scratchFluxOne, m_ebZero); 
  
#if 0 // Debug
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
  CH_TIME("CdrSolver::smoothHeavisideFaces");
  if(m_verbosity > 5){
    pout() << m_name + "::smoothHeavisideFaces" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    const Real vol               = pow(dx, SpaceDim);
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();

      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& face_states       = (*a_facePhi[lvl])[dit()][dir];
	const EBCellFAB& cell_states = (*a_cellPhi[lvl])[dit()];

	BaseFab<Real>& reg_face       = face_states.getSingleValuedFAB();
	const BaseFab<Real>& reg_cell = cell_states.getSingleValuedFAB();

	Box facebox = box;
	//	facebox.grow(dir,1);
	facebox.surroundingNodes(dir);

	// This will also do irregular cells and boundary faces
	FORT_HEAVISIDE_MEAN(CHF_FRA1(reg_face, comp),  
			    CHF_CONST_FRA1(reg_cell, comp),
			    CHF_CONST_INT(dir),
			    CHF_CONST_REAL(dx),
			    CHF_BOX(facebox));

	// Fix up irregular cell faces
	const IntVectSet& irreg = ebisbox.getIrregIVS(box);
	const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingNoBoundary;
	for (FaceIterator faceit(irreg, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();


	  const VolIndex lovof = face.getVoF(Side::Lo);
	  const VolIndex hivof = face.getVoF(Side::Hi);

	  const Real loval = Max(0.0, cell_states(lovof, comp));
	  const Real hival = Max(0.0, cell_states(hivof, comp));


	  Real Hlo;
	  if(loval*vol <= 0.0){
	    Hlo = 0.0;
	  }
	  else if(loval*vol >= 1.0){
	    Hlo = 1.0;
	  }
	  else{
	    Hlo = loval*vol;
	  }

	  Real Hhi;
	  if(hival*vol <= 0.0){
	    Hhi = 0.0;
	  }
	  else if(hival*vol >= 1.0){
	    Hhi = 1.0;
	  }
	  else{
	    Hhi = hival*vol;
	  }


	  face_states(face, comp) = 0.5*(hival + loval)*Hlo*Hhi;
	}

	// No random flux on domain faces. Reset those. 
	for (SideIterator sit; sit.ok(); ++sit){
	  Box sidebox;
	  if(sit() == Side::Lo){
	    sidebox = bdryLo(domain, dir, 1);
	  }
	  else if(sit() == Side::Hi){
	    sidebox = bdryHi(domain, dir, 1);
	  }
	  
	  sidebox &= facebox;

	  Box cellbox = sidebox.enclosedCells(dir);

	  const IntVectSet ivs(cellbox);
	  const FaceStop::WhichFaces stopcrit = FaceStop::AllBoundaryOnly;
	  for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	    face_states(faceit(), comp) = 0.0;
	  }
	}
      }
    }

    // Covered is bogus
    EBLevelDataOps::setCoveredVal(*a_facePhi[lvl], 0.0);
  }
}

void CdrSolver::fillGwn(EBAMRFluxData& a_noise, const Real a_sigma){
  CH_TIME("CdrSolver::fillGwn");
  if(m_verbosity > 5){
    pout() << m_name + "::fillGwn" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();

  std::normal_distribution<double> GWN(0.0, a_sigma);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    const Real vol               = pow(dx, SpaceDim);
    const Real ivol              = sqrt(1./vol);
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();

      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& noise = (*a_noise[lvl])[dit()][dir];

	Box facebox = box;
	facebox.surroundingNodes(dir);

	noise.setVal(0.0);

	// Regular faces
	BaseFab<Real>& noise_reg = noise.getSingleValuedFAB();
	for (BoxIterator bit(facebox); bit.ok(); ++bit){
	  const IntVect iv = bit();
	  noise_reg(iv, comp) = GWN(m_rng)*ivol;
	}

	// // Irregular faces
	const IntVectSet& irreg = ebisbox.getIrregIVS(box);
	const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingNoBoundary;
	for (FaceIterator faceit(irreg, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  noise(face, comp) = GWN(m_rng)*ivol;
	}
      }
    }
  }
}

void CdrSolver::parseRngSeed(){
  CH_TIME("CdrSolver::parseRngSeed");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRngSeed" << endl;
  }  
  ParmParse pp(m_className.c_str());
  pp.get("seed", m_seed);
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();

  m_rng = std::mt19937_64(m_seed);
}

void CdrSolver::parsePlotMode(){
  CH_TIME("CdrSolver::parsePlotMode");
  if(m_verbosity > 5){
    pout() << m_name + "::parsePlotMode" << endl;
  }  
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
