/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrSolver.cpp
  @brief  Implementation of CD_CdrSolver.H
  @author Robert Marskar
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
  CH_TIME("CdrSolver::setDefaultDomainBC");
  if(m_verbosity > 5){
    pout() << m_name + "::setDefaultDomainBC" << endl;
  }

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
  CH_TIME("CdrSolver::setDomainBcType");
  if(m_verbosity > 5){
    pout() << m_name + "::setDomainBcType" << endl;
  }

  m_domainBC.setBcType(a_domainSide, a_bcType);
}


void CdrSolver::setDomainBcFunction(const CdrDomainBC::DomainSide a_domainSide, const CdrDomainBC::FluxFunction a_fluxFunction){
  CH_TIME("CdrSolver::setDomainBcFunction");
  if(m_verbosity > 5){
    pout() << m_name + "::setDomainBcFunction" << endl;
  }

  m_domainBC.setBcFunction(a_domainSide, a_fluxFunction);
}

std::string CdrSolver::getName() const {
  return m_name;
}

std::string CdrSolver::getRealm() const {
  return m_realm;
}

void CdrSolver::setRealm(const std::string a_realm) {
  m_realm = a_realm;
}

Vector<std::string> CdrSolver::getPlotVariableNames() const {
  CH_TIME("CdrSolver::getPlotVariableNames");
  if(m_verbosity > 5){
    pout() << m_name + "::getPlotVariableNames" << endl;
  }
  
  Vector<std::string> names(0);
  
  if(m_plotPhi)                                   names.push_back(m_name + " phi");
  if(m_plotDiffusionCoefficient && m_isDiffusive) names.push_back(m_name + " diffusion_coefficient");
  if(m_plotSource)                                names.push_back(m_name + " source");
  if(m_plotVelocity && m_isMobile){
    names.push_back("x-Velocity " + m_name);
    names.push_back("y-Velocity " + m_name);
    if(SpaceDim == 3){
      names.push_back("z-Velocity " + m_name);
    }
  }
  if(m_plotEbFlux && m_isMobile) names.push_back(m_name + " eb_flux");
  
  return names;
}

int CdrSolver::getNumberOfPlotVariables() const {
  CH_TIME("CdrSolver::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int num_output = 0;

  if(m_plotPhi)                                   num_output = num_output + 1;
  if(m_plotDiffusionCoefficient && m_isDiffusive) num_output = num_output + 1;
  if(m_plotSource)                                num_output = num_output + 1;
  if(m_plotVelocity && m_isMobile)                num_output = num_output + SpaceDim;
  if(m_plotEbFlux && m_isMobile)                  num_output = num_output + 1;

  return num_output;
}

void CdrSolver::advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt){
  CH_TIME("CdrSolver::advanceEuler(no source)");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceEuler(no source)" << endl;
  }
  
  if(m_isDiffusive){
    
    // Create a source term = S = 0.0 and then call the other version. 
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
  CH_TIME("CdrSolver::advanceTGA(no source)");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceTGA(no source)" << endl;
  }

  if(m_isDiffusive){

    // Create a source term = S = 0.0 and then call the other version. 
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
  CH_TIME("CdrSolver::allocateInternals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocateInternals" << endl;
  }

  const int sca = 1;
  const int vec = SpaceDim;

  // This is allocated no matter what. 
  m_amr->allocate(m_phi,     m_realm, m_phase, sca);
  m_amr->allocate(m_source,  m_realm, m_phase, sca);
  m_amr->allocate(m_scratch, m_realm, m_phase, sca);
  
  DataOps::setValue(m_phi,      0.0);
  DataOps::setValue(m_source,   0.0);
  DataOps::setValue(m_scratch,  0.0);

  // Only allocate memory for cell-centered and face-centered velocities if the solver is mobile. Otherwise, allocate
  // a NULL pointer that we can pass around in TimeStepper in order to handle special cases
  if(m_isMobile){
    m_amr->allocate(m_faceVelocity, m_realm, m_phase, sca);
    m_amr->allocate(m_cellVelocity, m_realm, m_phase, vec);
    m_amr->allocate(m_faceStates,   m_realm, m_phase, sca);
    
    DataOps::setValue(m_faceVelocity,  0.0);
    DataOps::setValue(m_cellVelocity,  0.0);
  }
  else{
    m_amr->allocatePointer(m_faceVelocity);
    m_amr->allocatePointer(m_cellVelocity);
  }

  // Only allocate memory for diffusion coefficients if we need it. Otherwise, allocate a NULL pointer that we can
  // pass around in TimeStepper in order to handle special cases
  if(m_isDiffusive){
    m_amr->allocate(m_faceCenteredDiffusionCoefficient, m_realm, m_phase, sca);
    m_amr->allocate(m_ebCenteredDiffusionCoefficient,   m_realm, m_phase, sca);
    
    DataOps::setValue(m_faceCenteredDiffusionCoefficient, 0.0);
    DataOps::setValue(m_ebCenteredDiffusionCoefficient,   0.0);
  }
  else{
    m_amr->allocatePointer(m_faceCenteredDiffusionCoefficient);
    m_amr->allocatePointer(m_ebCenteredDiffusionCoefficient);
  }

  // Allocate stuff for holding fluxes
  if(m_isDiffusive || m_isMobile){
    m_amr->allocate(m_scratchFluxOne, m_realm, m_phase, sca);
    m_amr->allocate(m_scratchFluxTwo, m_realm, m_phase, sca);
  }

  // These don't consume (much) memory so just allocate them 
  m_amr->allocate(m_ebFlux,              m_realm, m_phase, sca);
  m_amr->allocate(m_ebZero,              m_realm, m_phase, sca);
  m_amr->allocate(m_domainFlux,          m_realm, m_phase, sca);
  m_amr->allocate(m_massDifference,      m_realm, m_phase, sca);
  m_amr->allocate(m_nonConservativeDivG, m_realm, m_phase, sca);
  
  DataOps::setValue(m_ebFlux,     0.0);
  DataOps::setValue(m_ebZero,     0.0);
  DataOps::setValue(m_domainFlux, 0.0);

  // This defines interpolation stencils and space for interpolants
  this->defineInterpolationStencils();
  this->defineInterpolant();
}

void CdrSolver::deallocateInternals(){
  CH_TIME("CdrSolver::deallocateInternals");
  if(m_verbosity > 5){
    pout() << m_name + "::deallocateInternals" << endl;
  }

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
  CH_TIME("CdrSolver::averageVelocityToFaces(public, full)");
  if(m_verbosity > 5){
    pout() << m_name + "::averageVelocityToFaces(public, full)" << endl;
  }

  this->averageVelocityToFaces(m_faceVelocity, m_cellVelocity); // Average velocities to face centers for all levels
}

void CdrSolver::averageVelocityToFaces(EBAMRFluxData& a_faceVelocity, const EBAMRCellData& a_cellVelocity){
  CH_TIME("CdrSolver::averageVelocityToFaces");
  if(m_verbosity > 5){
    pout() << m_name + "::averageVelocityToFaces" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    DataOps::averageCellToFace(*a_faceVelocity[lvl], *a_cellVelocity[lvl], m_amr->getDomains()[lvl]);
    a_faceVelocity[lvl]->exchange();
  }
}

void CdrSolver::preRegrid(const int a_lmin, const int a_oldFinestLevel){
  CH_TIME("CdrSolver::preRegrid");
  if(m_verbosity > 5){
    pout() << m_name + "::preRegrid" << endl;
  }

  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();
  
  m_amr->allocate(m_cachePhi,    m_realm, m_phase, ncomp);
  m_amr->allocate(m_cacheSource, m_realm, m_phase, ncomp);
  
  for (int lvl = 0; lvl <= a_oldFinestLevel; lvl++){
    m_phi[lvl]->localCopyTo(*m_cachePhi[lvl]);
    m_source[lvl]->localCopyTo(*m_cacheSource[lvl]);
  }
}

void CdrSolver::coarseFineIncrement(const EBAMRIVData& a_massDifference){
  CH_TIME("CdrSolver::coarseFineIncrement");
  if(m_verbosity > 5){
    pout() << m_name + "::coarseFineIncrement" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(0,0);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < 0;

    if(has_coar){
      fine2coar_redist->setToZero();

    }
    if(has_fine){
      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      if(has_coar){
	fine2coar_redist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
      }

      if(has_fine){
	coar2fine_redist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
	coar2coar_redist->increment((*a_massDifference[lvl])[dit()], dit(), interv);
      }
    }
  }
}

void CdrSolver::coarseFineRedistribution(EBAMRCellData& a_phi){
  CH_TIME("CdrSolver::coarseFineRedistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::coarseFineRedistribution" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->getDx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];
    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(m_realm, m_phase)[lvl];
    
    if(has_coar){
      fine2coar_redist->redistribute(*a_phi[lvl-1], interv);
      fine2coar_redist->setToZero();
    }

    if(has_fine){
      coar2fine_redist->redistribute(*a_phi[lvl+1], interv);
      coar2coar_redist->redistribute(*a_phi[lvl],   interv);

      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }
  }
}

void CdrSolver::computeDivG(EBAMRCellData& a_divG, EBAMRFluxData& a_G, const EBAMRIVData& a_ebFlux){
  CH_TIME("CdrSolver::computeDivG");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDivG" << endl;
  }

  DataOps::setValue(a_divG, 0.0);
  
  this->conservativeDivergenceNoKappaDivision(a_divG, a_G, a_ebFlux);      // Make the conservative divergence.
  this->nonConservativeDivergence(m_nonConservativeDivG, a_divG);          // Non-conservative divergence
  this->hybridDivergence(a_divG, m_massDifference, m_nonConservativeDivG); // a_divG becomes hybrid divergence. Mass diff computed. 
  this->incrementFluxRegister(a_G);                                        // Increment flux register
  this->incrementRedist(m_massDifference);                                 // Increment level redistribution register

  // Redistribution and reflux magic. 
  this->coarseFineIncrement(m_massDifference);              // Compute C2F, F2C, and C2C mass transfers
  this->incrementRedistFlux();                              // Tell flux register about whats going on
  this->hyperbolicRedistribution(a_divG, m_massDifference); // Level redistribution. Weights is a dummy parameter
  this->coarseFineRedistribution(a_divG);                   // Do the coarse-fine redistribution
  this->reflux(a_divG);                                     // Reflux
}

void CdrSolver::injectEbFlux(EBAMRCellData& a_phi, const EBAMRIVData& a_ebFlux, const Real a_dt){
  CH_TIME("CdrSolver::injectEbFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::injectEbFlux" << endl;
  }

  MayDay::Warning("CdrSolver::injectEbFlux - routine has not been wetted!");

  if(m_useMassWeightedRedistribution){
    this->setRedistWeights(a_phi);
  }

  this->conservativeDivergenceNoKappaDivisionOnlyEbFlux(m_scratch, a_ebFlux);  // Compute conservative divergence, but only EB
  this->nonConservativeDivergence(m_nonConservativeDivG, m_scratch);           // Blend with volume fraction
  this->hybridDivergence(m_scratch, m_massDifference, m_nonConservativeDivG);  // Hybrid divergence
  this->incrementRedist(m_massDifference);                                     // Increment redistribution register

  // Redistribution and reflux magic. 
  this->coarseFineIncrement(m_massDifference);                 // Compute C2F, F2C, and C2C mass transfers
  this->incrementRedistFlux();                                 // Tell flux register about whats going on
  this->hyperbolicRedistribution(m_scratch, m_massDifference); // Level redistribution. 
  this->coarseFineRedistribution(m_scratch);                   // Do the coarse-fine redistribution
  this->reflux(m_scratch);                                     // Reflux

  // Now do the increment. 
  DataOps::incr(a_phi, m_scratch, -a_dt);
}

void CdrSolver::conservativeDivergenceNoKappaDivisionOnlyEbFlux(EBAMRCellData& a_conservativeDivergence, const EBAMRIVData& a_ebFlux){
  CH_TIME("CdrSolver::injectEbFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::injectEbFlux" << endl;
  }

  // TLDR: This sets a_conservativeDivergence = a_ebFlux*area/dx

  const int comp = 0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real dx                = m_amr->getDx()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& divG            = (*a_conservativeDivergence[lvl])[dit()];
      const EBISBox& ebisbox     = ebisl[dit()];
      const BaseIVFAB<Real>& flx = (*a_ebFlux[lvl])[dit()];

      divG.setVal(0.0);
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real area     = ebisbox.bndryArea(vof);

	divG(vof, comp) = flx(vof, comp)*area/dx;
      }
    }
  }
}

void CdrSolver::computeDivergenceIrregular(LevelData<EBCellFAB>&              a_divG,
					   const LevelData<BaseIVFAB<Real> >& a_ebFlux,
					   const int                          a_lvl){
  CH_TIME("CdrSolver::computeDivergenceIrregular");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDivergenceIrregular" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];
  const Real dx                = m_amr->getDx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];

    EBCellFAB& divG               = a_divG[dit()];
    const BaseIVFAB<Real>& ebflux = a_ebFlux[dit()];

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real area     = ebisbox.bndryArea(vof);

      // EB flux
      divG(vof, comp) = ebflux(vof,comp)*area;

      // Face fluxes
      for (int dir = 0; dir < SpaceDim; dir++){

	const BaseIFFAB<Real>& flux = (*(m_interpolant[dir])[a_lvl])[dit()];

	for (SideIterator sit; sit.ok(); ++sit){
	  const int isign = sign(sit());
	  const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());

	  for (int iface = 0; iface < faces.size(); iface++){
	    const FaceIndex face = faces[iface];
	    const Real face_area = ebisbox.areaFrac(face);
	    divG(vof, comp) += isign*face_area*flux(face, comp);
	  }
	}
      }

      // Scale divG by dx but not by kappa.
      divG(vof, comp) *= 1./dx;
    }
  }
}

void CdrSolver::computeAdvectionFlux(EBAMRFluxData&       a_flux,
				     const EBAMRFluxData& a_facePhi,
				     const EBAMRFluxData& a_faceVelocity,
				     const bool           a_addDomainFlux){
  CH_TIME("CdrSolver::computeAdvectionFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::computeAdvectionFlux" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    this->computeAdvectionFlux(*a_flux[lvl], *a_facePhi[lvl], *a_faceVelocity[lvl], lvl);
  }

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

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(box);

    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flx       = a_flux[dit()][dir];
      const EBFaceFAB& phi = a_facePhi[dit()][dir];
      const EBFaceFAB& vel = a_faceVelocity[dit()][dir];

      flx.setVal(0.0, comp);
      flx += phi;
      flx *= vel;

      // Irregular faces
      const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;
      for (FaceIterator faceit(ebisbox.getIrregIVS(box), ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();
	flx(face, comp) = vel(face, comp)*phi(face, comp);
      }
    }
  }
}

void CdrSolver::computeDiffusionFlux(EBAMRFluxData& a_flux, const EBAMRCellData& a_phi, const bool a_addDomainFlux){
  CH_TIME("CdrSolver::computeDiffusionFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDiffusionFlux" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    this->computeDiffusionFlux(*a_flux[lvl], *a_phi[lvl], lvl);
  }

  if(a_addDomainFlux){
    this->fillDomainFlux(a_flux);
  }
  else{
    this->resetDomainFlux(a_flux);
  }
}

void CdrSolver::computeDiffusionFlux(LevelData<EBFluxFAB>& a_flux, const LevelData<EBCellFAB>& a_phi, const int a_lvl){
  CH_TIME("CdrSolver::computeDiffusionFlux(level)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDiffusionFlux(level)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const Real dx                = m_amr->getDx()[a_lvl];
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& cellbox = dbl.get(dit());

    const EBCellFAB& state = a_phi[dit()];
    const EBISBox& ebisbox = state.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();

    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux = a_flux[dit()][dir];
      const EBFaceFAB& dco = (*m_faceCenteredDiffusionCoefficient[a_lvl])[dit()][dir];

      // Do regular cells
      Box facebox = cellbox;
      facebox.surroundingNodes(dir);

      BaseFab<Real>& regFlux        = flux.getSingleValuedFAB();
      const BaseFab<Real>& regState = state.getSingleValuedFAB();
      const BaseFab<Real>& regDco   = dco.getSingleValuedFAB();
      FORT_DFLUX_REG(CHF_FRA1(regFlux, comp),
		     CHF_CONST_FRA1(regState, comp),
		     CHF_CONST_FRA1(regDco, comp),
		     CHF_CONST_INT(dir),
		     CHF_CONST_REAL(dx),
		     CHF_BOX(facebox));


      // Now redo the irregular faces
      const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;
      for (FaceIterator faceit(ebisbox.getIrregIVS(cellbox), ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();

	if(!face.isBoundary()){ 
	  const VolIndex hiVoF = face.getVoF(Side::Hi);
	  const VolIndex loVoF = face.getVoF(Side::Lo);

	  flux(face, comp) = dco(face,comp)*(state(hiVoF,comp) - state(loVoF, comp))/dx;
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

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& flux = (*a_flux[lvl])[dit()][dir];

	for (SideIterator sit; sit.ok(); ++sit){
	  const BaseIFFAB<Real>& domFlux = (*m_domainFlux[lvl])[dit()](dir, sit());

	  const IntVectSet& ivs  = domFlux.getIVS();
	  const EBGraph& ebgraph = domFlux.getEBGraph();

	  for (FaceIterator faceIt(ivs, ebgraph, dir, FaceStop::AllBoundaryOnly); faceIt.ok(); ++faceIt){
	    flux(faceIt(), m_comp) = 0.0;
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
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    this->fillDomainFlux(*a_flux[lvl], lvl);
  }
}

void CdrSolver::fillDomainFlux(LevelData<EBFluxFAB>& a_flux, const int a_level) {
  CH_TIME("CdrSolver::fillDomainFlux(LevelData<EBFluxFAB>, int)");
  if(m_verbosity > 5){
    pout() << m_name + "::fillDomainFlux(LevelData<EBFluxFAB>, int)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const EBISBox& ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_level][dit()];
    
    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux = a_flux[dit()][dir];

      // Iterate over domain cells. I'm REALLY not sure about if this is performant...
      for (SideIterator sit; sit.ok(); ++sit){
	const CdrDomainBC::DomainSide    domainSide = std::make_pair(dir, sit());
	const CdrDomainBC::BcType&       bcType     = m_domainBC.getBcType(domainSide);
	const CdrDomainBC::FluxFunction& bcFunction = m_domainBC.getBcFunction(domainSide);

	const BaseIFFAB<Real>& domflux = (*m_domainFlux[a_level])[dit()](dir, sit());
	const IntVectSet& ivs          = domflux.getIVS();
	const EBGraph& ebgraph         = domflux.getEBGraph();

	// Iterate through faces on the boundary and set the BC condition. 
	for (FaceIterator faceit(ivs, ebgraph, dir, FaceStop::AllBoundaryOnly); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  switch(bcType){
	  case CdrDomainBC::BcType::DataBased:{
	    flux(face, m_comp) = domflux(face, m_comp);
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
	  case CdrDomainBC::BcType::External: // Don't do anything, the solver has responsibility. 
	    break;
	  default:
	    MayDay::Abort("CdrSolver::fillDomainFlux - logic bust");
	  }
	}
      }
    }
  }
}

void CdrSolver::conservativeDivergenceNoKappaDivision(EBAMRCellData& a_conservativeDivergence, EBAMRFluxData& a_flux, const EBAMRIVData& a_ebFlux){
  CH_TIME("CdrSolver::conservativeDivergenceNoKappaDivision");
  if(m_verbosity > 5){
    pout() << m_name + "::conservativeDivergenceNoKappaDivision" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_flux[lvl]->exchange();
    
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    this->conservativeDivergenceRegular(*a_conservativeDivergence[lvl], *a_flux[lvl], lvl);
    this->setupFluxInterpolant(*a_flux[lvl], lvl);                                           // Copy face-centered fluxes in a_flux to m_interpolant
    this->interpolateFluxToFaceCentroids(*a_flux[lvl], lvl);                                 // Interpolate fluxes w m_interpolant. Copy 2 a_flux
    this->computeDivergenceIrregular(*a_conservativeDivergence[lvl], *a_ebFlux[lvl], lvl);   // Recompute divergence on irregular cells

    a_conservativeDivergence[lvl]->exchange();
  }

  m_amr->averageDown(a_conservativeDivergence, m_realm, m_phase);
  m_amr->interpGhost(a_conservativeDivergence, m_realm, m_phase);
}

void CdrSolver::conservativeDivergenceRegular(LevelData<EBCellFAB>& a_divJ, const LevelData<EBFluxFAB>& a_flux, const int a_lvl){
  CH_TIME("CdrSolver::conservativeDivergenceRegular");
  if(m_verbosity > 5){
    pout() << m_name + "::conservativeDivergenceRegular" << endl;
  }

  CH_assert(a_divJ.nComp() == 1);
  CH_assert(a_flux.nComp() == 1);

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain domain  = m_amr->getDomains()[a_lvl];
  const Real dx               = m_amr->getDx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& divJ         = a_divJ[dit()];
    BaseFab<Real>& divJ_fab = divJ.getSingleValuedFAB();
    const Box box = dbl.get(dit());

    divJ.setVal(0.0);
    for (int dir = 0; dir < SpaceDim; dir++){
      const EBFaceFAB& flx         = a_flux[dit()][dir];
      const BaseFab<Real>& flx_fab = flx.getSingleValuedFAB();

      FORT_CONSDIV_REG(CHF_FRA1(divJ_fab, comp),
		       CHF_CONST_FRA1(flx_fab, comp),
		       CHF_CONST_INT(dir),
		       CHF_CONST_REAL(dx),
		       CHF_BOX(box));
    }

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      divJ(vof, comp) = 0.0;
    }
  }

  a_divJ.exchange();
}

void CdrSolver::defineInterpolationStencils(){
  CH_TIME("CdrSolver::defineInterpolationStencils");
  if(m_verbosity > 5){
    pout() << m_name + "::defineInterpolationStencils" << endl;
  }

  const int comp            = 0;
  const int ncomp           = 1;
  const int finest_level    = m_amr->getFinestLevel();
  FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;
  
  for (int dir = 0; dir < SpaceDim; dir++){
    (m_interpStencils[dir]).resize(1 + finest_level);

    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const ProblemDomain& domain  = m_amr->getDomains()[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

      m_interpStencils[dir][lvl] = RefCountedPtr<LayoutData<BaseIFFAB<FaceStencil> > >
	(new LayoutData<BaseIFFAB<FaceStencil> >(dbl));

      LayoutData<IntVectSet> cfivs(dbl);
      EBArith::defineCFIVS(cfivs, dbl, domain);

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	BaseIFFAB<FaceStencil>& sten = (*m_interpStencils[dir][lvl])[dit()];
	const Box box          = dbl.get(dit());
	const EBISBox& ebisbox = ebisl[dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();
	const IntVectSet ivs   = ebisbox.getIrregIVS(box);

	sten.define(ivs, ebisbox.getEBGraph(), dir, ncomp);
	
	for (FaceIterator faceit(ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();
	  FaceStencil& facesten = sten(face,comp);

	  facesten = EBArith::getInterpStencil(face, cfivs[dit()], ebisbox, domain.domainBox());
	}
      }
    }
  }
}

void CdrSolver::incrementRedistFlux(){
  CH_TIME("CdrSolver::incrementRedistFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::incrementRedistFlux" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->getDx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;
    
    if(has_fine){
      RefCountedPtr<EBFluxRegister>& fluxreg = m_amr->getFluxRegister(m_realm, m_phase)[lvl];
      RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];
      RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];
      
      const Real scale = -dx;

      fluxreg->incrementRedistRegister(*coar2fine_redist, interv, scale);
      fluxreg->incrementRedistRegister(*coar2coar_redist, interv, scale);

    }
  }
}

void CdrSolver::initialData(){
  CH_TIME("CdrSolver::initialData");
  if(m_verbosity > 5){
    pout() << m_name + "::initialData" << endl;
  }

  const bool deposit_function  = m_species->initializeWithFunction();
  const bool depositParticles = m_species->initializeWithParticles();

  DataOps::setValue(m_phi, 0.0);
  
  if(depositParticles){
    initialDataParticles();
  }

  // Increment with function values if this is also called for
  if(deposit_function){
    initialDataDistribution();
  }

  m_amr->averageDown(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);
}

void CdrSolver::initialDataDistribution(){
  CH_TIME("CdrSolver::initialDataDistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::initialDataDistribution" << endl;
  }

  const RealVect origin  = m_amr->getProbLo();
  const int finest_level = m_amr->getFinestLevel();

  // Copy this
  DataOps::copy(m_scratch, m_phi);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx = m_amr->getDx()[lvl];

    for (DataIterator dit = m_phi[lvl]->dataIterator(); dit.ok(); ++dit){
      EBCellFAB& state         = (*m_phi[lvl])[dit()];
      const EBCellFAB& scratch = (*m_scratch[lvl])[dit()];
      const Box box            = m_phi[lvl]->disjointBoxLayout().get(dit());
      const EBISBox& ebisbox   = state.getEBISBox();

      BaseFab<Real>& reg_state   = state.getSingleValuedFAB();
      const BaseFab<Real>& reg_scratch = scratch.getSingleValuedFAB();

      // Regular cells
      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();
	const RealVect pos = origin + RealVect(iv)*dx + 0.5*dx*RealVect::Unit;

	for (int comp = 0; comp < state.nComp(); comp++){
	  reg_state(iv, comp) = reg_scratch(iv, comp) + m_species->initialData(pos, m_time);
	}
      }

      // Irreg and multicells
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real kappa    = ebisbox.volFrac(vof);
	const RealVect pos  = EBArith::getVofLocation(vof, m_amr->getDx()[lvl]*RealVect::Unit, origin);
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = scratch(vof, comp) + m_species->initialData(pos, m_time);
	}
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

  const int halo_buffer = 0;
  const int pvr_buffer  = 0;

  const List<Particle>& init_particles = m_species->getInitialParticles();

  if(init_particles.length() > 0){

    ParticleContainer<Particle> particles;
    m_amr->allocate(particles, pvr_buffer, m_realm);
    particles.addParticles(m_species->getInitialParticles());

    // We will deposit onto m_phi, using m_scratch as a scratch holder for interpolation stuff
    DataOps::setValue(m_phi, 0.0);
  
    // Deposit onto mseh
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const RealVect dx            = m_amr->getDx()[lvl]*RealVect::Unit;
      const RealVect origin        = m_amr->getProbLo();
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const ProblemDomain& dom     = m_amr->getDomains()[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    
      // 2. Deposit this levels particles and exchange ghost cells
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const Box box          = dbl.get(dit());
	const EBISBox& ebisbox = ebisl[dit()];
	EbParticleInterp interp(box, ebisbox, dx, origin, true);
	interp.deposit(particles[lvl][dit()].listItems(), (*m_phi[lvl])[dit()].getFArrayBox(), DepositionType::NGP);
      }
    }

#if CH_SPACEDIM==2 // Only do this scaling for planar cartesian
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const Real dx = m_amr->getDx()[lvl];
      DataOps::scale(*m_phi[lvl], 1./dx);
    }
#endif
  }
}

void CdrSolver::hybridDivergence(EBAMRCellData&     a_hybrid_div,
				 EBAMRIVData&       a_massDifference,
				 const EBAMRIVData& a_nonConservativeDivergence){
  CH_TIME("CdrSolver::hybridDivergence(AMR)");
  if(m_verbosity > 5){
    pout() << m_name + "::hybridDivergence(AMR)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    hybridDivergence(*a_hybrid_div[lvl], *a_massDifference[lvl], *a_nonConservativeDivergence[lvl], lvl);
  }
}

void CdrSolver::hybridDivergence(LevelData<EBCellFAB>&              a_hybridDivergence,
				 LevelData<BaseIVFAB<Real> >&       a_massDifference,
				 const LevelData<BaseIVFAB<Real> >& a_nonConservativeDivergence,
				 const int                          a_lvl){
  CH_TIME("CdrSolver::hybridDivergence(level)");
  if(m_verbosity > 5){
    pout() << m_name + "::hybridDivergence(level)" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];
    
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& divH               = a_hybridDivergence[dit()];  // On input, this contains kappa*div(F)
    BaseIVFAB<Real>& deltaM       = a_massDifference[dit()];
    const BaseIVFAB<Real>& divNC  = a_nonConservativeDivergence[dit()]; 

    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real kappa    = ebisbox.volFrac(vof);
      const Real dc       = divH(vof, comp);
      const Real dnc      = divNC(vof, comp);

      // Note to self: deltaM = (1-kappa)*(dc - kappa*dnc) because dc was not divided by kappa,
      // which it would be otherwise. 
      divH(vof, comp)   = dc + (1-kappa)*dnc;          // On output, contains hybrid divergence
      deltaM(vof, comp) = (1-kappa)*(dc - kappa*dnc);
    }
  }
}

void CdrSolver::setRedistWeights(const EBAMRCellData& a_phi){
  CH_TIME("CdrSolver::setRedistWeights");
  if(m_verbosity > 5){
    pout() << m_name + "::setRedistWeights" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    EBLevelRedist& redist = *m_amr->getLevelRedist(m_realm, m_phase)[lvl];
    redist.resetWeights(*a_phi[lvl], comp);

    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_coar){
      EBFineToCoarRedist& fine2coar = *m_amr->getFineToCoarRedist(m_realm, m_phase)[lvl];
      fine2coar.resetWeights(*a_phi[lvl-1], comp);
    }
    if(has_fine){
      EBCoarToCoarRedist& coar2coar = *m_amr->getCoarToCoarRedist(m_realm, m_phase)[lvl];
      EBCoarToFineRedist& coar2fine = *m_amr->getCoarToFineRedist(m_realm, m_phase)[lvl];

      coar2coar.resetWeights(*a_phi[lvl], comp);
      coar2fine.resetWeights(*a_phi[lvl], comp);
    }
  }
}

void CdrSolver::hyperbolicRedistribution(EBAMRCellData& a_divF, const EBAMRIVData&   a_massDifference) {
  CH_TIME("CdrSolver::hyberbolic_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::hyperbolicRedistribution" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    EBLevelRedist& level_redist = *(m_amr->getLevelRedist(m_realm, m_phase)[lvl]);
    level_redist.redistribute(*a_divF[lvl], interv);
    level_redist.setToZero();
  }
}

void CdrSolver::interpolateFluxToFaceCentroids(LevelData<EBFluxFAB>& a_flux, const int a_lvl){
  CH_TIME("CdrSolver::interpolateFluxToFaceCentroids");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolateFluxToFaceCentroids" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];

  const FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);
    
    for (int dir = 0; dir < SpaceDim; dir++){
      BaseIFFAB<Real>& interpolant = (*(m_interpolant[dir])[a_lvl])[dit()]; // Holds the face-centered flux
      BaseIFFAB<Real> centroid_flux;
      EBFaceFAB& faceFlux = a_flux[dit()][dir];

      // Compute face centroid flux
      centroid_flux.define(ivs, ebgraph, dir, ncomp);
      for (FaceIterator faceit(ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	const FaceIndex& face   = faceit();
	const FaceStencil& sten = (*m_interpStencils[dir][a_lvl])[dit()](face, comp);

	centroid_flux(face, comp) = 0.;
	Real sum = 0.0;
	for (int i = 0; i < sten.size(); i++){
	  const FaceIndex& iface = sten.face(i);
	  const Real iweight     = sten.weight(i);
	  sum += iweight;
	  
	  centroid_flux(face, comp) += iweight*interpolant(iface, comp);
	}

	faceFlux(face, comp) = centroid_flux(face, comp);
      }

      // Copy centroid flux into a_flux
      interpolant.setVal(0.0);
      interpolant.copy(box, interv, box, centroid_flux, interv);
    }
  }
}

void CdrSolver::resetFluxRegister(){
  CH_TIME("CdrSolver::resetFluxRegister");
  if(m_verbosity > 5){
    pout() << m_name + "::resetFluxRegister" << endl;
  }
  
  const int finest_level = m_amr->getFinestLevel();
  
  Vector<RefCountedPtr<EBFluxRegister> >& fluxreg = m_amr->getFluxRegister(m_realm, m_phase);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const bool has_fine = lvl < finest_level;
    if(has_fine){
      fluxreg[lvl]->setToZero();
    }
  }
  
}

void CdrSolver::incrementFluxRegister(const EBAMRFluxData& a_flux){
  CH_TIME("CdrSolver::incrementFluxRegister(flux)");
  if(m_verbosity > 5){
    pout() << m_name + "::incrementFluxRegister(flux)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFluxRegister> >& fluxreg = m_amr->getFluxRegister(m_realm, m_phase);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_fine){
      fluxreg[lvl]->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      const Box box          = dbl.get(dit());
      
      for (int dir = 0; dir < SpaceDim; dir++){
	const Real scale      = 1.0;
	const EBFaceFAB& flux = (*a_flux[lvl])[dit()][dir];

	// Increment flux register for irregular/regular. Add both from coarse to fine and from fine to coarse
	if(has_fine){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxreg[lvl]->incrementCoarseBoth(flux, scale, dit(), interv, dir, sit());
	  }
	}
	if(has_coar){
	  for (SideIterator sit; sit.ok(); ++sit){
	    fluxreg[lvl-1]->incrementFineBoth(flux, scale, dit(), interv, dir, sit());
	  }
	}
      }
    }
  }
}

void CdrSolver::incrementRedist(const EBAMRIVData& a_massDifference){
  CH_TIME("CdrSolver::incrementRedist");
  if(m_verbosity > 5){
    pout() << m_name + "::incrementRedist" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    
    EBLevelRedist& level_redist = *(m_amr->getLevelRedist(m_realm, m_phase)[lvl]);
    level_redist.setToZero();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      level_redist.increment((*a_massDifference[lvl])[dit()], dit(), interv);
    }
  }
}

void CdrSolver::nonConservativeDivergence(EBAMRIVData& a_nonConservativeDivergence, const EBAMRCellData& a_divG){
  CH_TIME("CdrSolver::nonConservativeDivergence");
  if(m_verbosity > 5){
    pout() << m_name + "::nonConservativeDivergence" << endl;
  }

  if(m_blendConservation){
    const IrregAmrStencil<NonConservativeDivergenceStencil>& stencils = m_amr->getNonConservativeDivergenceStencils(m_realm, m_phase);
    stencils.apply(a_nonConservativeDivergence, a_divG);
  }
  else{
    DataOps::setValue(a_nonConservativeDivergence, 0.0);
  }
}

void CdrSolver::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("CdrSolver::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  this->allocateInternals();

  Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->getPwlInterpolator(m_realm, m_phase);

  // These levels have not changed
  for (int lvl = 0; lvl <= Max(0,a_lmin-1); lvl++){
    m_cachePhi[lvl]->copyTo(*m_phi[lvl]); // Base level should never change, but ownership can.
    m_cacheSource[lvl]->copyTo(*m_source[lvl]); // Base level should never change, but ownership can.
  }

  // These levels have changed
  for (int lvl = Max(1,a_lmin); lvl <= a_newFinestLevel; lvl++){
    interpolator[lvl]->interpolate(*m_phi[lvl], *m_phi[lvl-1], interv);
    interpolator[lvl]->interpolate(*m_source[lvl], *m_source[lvl-1], interv);
    if(lvl <= Min(a_oldFinestLevel, a_newFinestLevel)){
      m_cachePhi[lvl]->copyTo(*m_phi[lvl]);
      m_cacheSource[lvl]->copyTo(*m_source[lvl]);
    }
  }

  //  DataOps::floor(m_phi, 0.0);
  m_amr->averageDown(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);

  m_amr->averageDown(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void CdrSolver::reflux(EBAMRCellData& a_phi){
  CH_TIME("CdrSolver::reflux");
  if(m_verbosity > 5){
    pout() << m_name + "::reflux" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFluxRegister > >& fluxreg = m_amr->getFluxRegister(m_realm, m_phase);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->getDx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_fine){
      const Real scale = 1.0/dx;
      
      fluxreg[lvl]->reflux(*a_phi[lvl], interv, scale);
      fluxreg[lvl]->setToZero();
    }
  }
}

void CdrSolver::sanityCheck(){
  CH_TIME("CdrSolver::sanityCheck");
  if(m_verbosity > 5){
    pout() << m_name + "::sanityCheck" << endl;
  }

  CH_assert(!m_computationalGeometry.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_species.isNull());
  CH_assert(!m_ebis.isNull());
}

void CdrSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("CdrSolver::setAmr");
  if(m_verbosity > 5){
    pout() << m_name + "::setAmr" << endl;
  }

  m_amr = a_amr;

}

void CdrSolver::registerOperators(){
  CH_TIME("CdrSolver::registerOperators");
  if(m_verbosity > 5){
    pout() << m_name + "::registerOperators" << endl;
  }

  if(m_amr.isNull()){
    MayDay::Abort("CdrSolver::registerOperators - need to set AmrMesh!");
  }
  else{
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, m_phase);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, m_phase);
    m_amr->registerOperator(s_eb_pwl_interp,   m_realm, m_phase);
    m_amr->registerOperator(s_eb_flux_reg,     m_realm, m_phase);
    m_amr->registerOperator(s_eb_redist,       m_realm, m_phase);
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, m_phase);
    m_amr->registerOperator(s_noncons_div,  m_realm, m_phase);
  }
}

void CdrSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry> a_computationalGeometry){
  CH_TIME("CdrSolver::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << m_name + "::setComputationalGeometry" << endl;
  }
  
  m_computationalGeometry = a_computationalGeometry;

  const RefCountedPtr<MultiFluidIndexSpace> MultiFluidIndexSpace = m_computationalGeometry->getMfIndexSpace();
  
  this->setEbIndexSpace(MultiFluidIndexSpace->getEBIndexSpace(m_phase));
}

void CdrSolver::setDiffusionCoefficient(const EBAMRFluxData& a_diffusionCoefficient, const EBAMRIVData& a_ebDiffusionCoefficient){
  CH_TIME("CdrSolver::setDiffusionCoefficient(ebamrflux, ebamriv)");
  if(m_verbosity > 5){
    pout() << m_name + "::setDiffusionCoefficient(ebamrflux, ebamriv)" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_diffusionCoefficient[lvl]->localCopyTo(*m_faceCenteredDiffusionCoefficient[lvl]);
    a_ebDiffusionCoefficient[lvl]->localCopyTo(*m_ebCenteredDiffusionCoefficient[lvl]);
  }

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
  DataOps::setValue(m_faceCenteredDiffusionCoefficient, a_diffCo, m_amr->getProbLo(), m_amr->getDx(), m_comp);
  DataOps::setValue(m_ebCenteredDiffusionCoefficient,   a_diffCo, m_amr->getProbLo(), m_amr->getDx(), m_comp);
}

void CdrSolver::setEbFlux(const EBAMRIVData& a_ebFlux){
  CH_TIME("CdrSolver::setEbFlux(variable)");
  if(m_verbosity > 5){
    pout() << m_name + "::setEbFlux(variable)" << endl;
  }

  const int comp         = 0;
  const int finest_level = m_amr->getFinestLevel();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_ebFlux[lvl]->localCopyTo(interv, *m_ebFlux[lvl], interv);
  }
}

void CdrSolver::setEbFlux(const Real a_ebFlux){
  CH_TIME("CdrSolver::setEbFlux(constant)");
  if(m_verbosity > 5){
    pout() << m_name + "::setEbFlux(constant)" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    DataOps::setValue(*m_ebFlux[lvl], a_ebFlux);
  }
}

void CdrSolver::setEbIndexSpace(const RefCountedPtr<EBIndexSpace>& a_ebis){
  CH_TIME("CdrSolver::setEbIndexSpace");
  if(m_verbosity > 5){
    pout() << m_name + "::setEbIndexSpace" << endl;
  }

  m_ebis = a_ebis;
}

void CdrSolver::setSpecies(const RefCountedPtr<CdrSpecies> a_species){
  CH_TIME("CdrSolver::setSpecies");
  if(m_verbosity > 5){
    pout() << m_name + "::setSpecies" << endl;
  }

  m_species   = a_species;
  m_name      = m_species->getName();
  m_isDiffusive = m_species->isDiffusive();
  m_isMobile    = m_species->isMobile();
}

void CdrSolver::setSource(const EBAMRCellData& a_source){
  CH_TIME("CdrSolver::setSource");
  if(m_verbosity > 5){
    pout() << m_name + "::setSource" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_source[lvl]->localCopyTo(*m_source[lvl]);
  }

  m_amr->averageDown(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void CdrSolver::setSource(const Real a_source){
  CH_TIME("CdrSolver::setSource");
  if(m_verbosity > 5){
    pout() << m_name + "::setSource" << endl;
  }

  const int comp = 0;
  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    DataOps::setValue(*m_source[lvl], a_source, comp);
  }

  m_amr->averageDown(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void CdrSolver::setSource(const std::function<Real(const RealVect a_position)> a_source){
  DataOps::setValue(m_source, a_source, m_amr->getProbLo(), m_amr->getDx(), m_comp);
}

void CdrSolver::setTime(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("CdrSolver::setTime");
  if(m_verbosity > 5){
    pout() << m_name + "::setTime" << endl;
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

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_velo[lvl]->localCopyTo(*m_cellVelocity[lvl]);
  }

  m_amr->averageDown(m_cellVelocity, m_realm, m_phase);
  m_amr->interpGhost(m_cellVelocity, m_realm, m_phase);
}

void CdrSolver::setVelocity(const RealVect a_velo){
  CH_TIME("CdrSolver::setVelocity");
  if(m_verbosity > 5){
    pout() << m_name + "::setVelocity" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (int dir = 0; dir < SpaceDim; dir++){
      DataOps::setValue(*m_cellVelocity[lvl], a_velo[dir], dir);
    }

    m_cellVelocity[lvl]->exchange();
  }

  m_amr->averageDown(m_cellVelocity, m_realm, m_phase);
  m_amr->interpGhost(m_cellVelocity, m_realm, m_phase);
}

void CdrSolver::setVelocity(const std::function<RealVect(const RealVect a_pos)>& a_velo){
  DataOps::setValue(m_cellVelocity, a_velo, m_amr->getProbLo(), m_amr->getDx());
}

void CdrSolver::setPhase(const phase::which_phase a_phase){
  CH_TIME("CdrSolver::setPhase");
  if(m_verbosity > 5){
    pout() << m_name + "::setPhase" << endl;
  }

  m_phase = a_phase;
}

void CdrSolver::setVerbosity(const int a_verbosity){
  CH_TIME("CdrSolver::setVerbosity");
  m_verbosity = a_verbosity;
  
  if(m_verbosity > 5){
    pout() << m_name + "::setVerbosity" << endl;
  }
}

void CdrSolver::setupFluxInterpolant(const LevelData<EBFluxFAB>& a_flux, const int a_lvl){
  CH_TIME("CdrSolver::setupFluxInterpolant");
  if(m_verbosity > 5){
    pout() << m_name + "::setupFluxInterpolant" << endl;
  }

  const int ncomp = 1;
  const int comp  = 0;
  
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];

  for (int dir = 0; dir < SpaceDim; dir++){
    
    const LayoutData<IntVectSet>& grown_set = (*(m_interpSets[dir])[a_lvl]);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      BaseIFFAB<Real>& interpol   = (*(m_interpolant[dir])[a_lvl])[dit()];
      const Box& box              = dbl.get(dit());
      const EBISBox& ebisbox      = ebisl[dit()];
      const EBGraph& ebgraph      = ebisbox.getEBGraph();
      const EBFaceFAB& flux       = a_flux[dit()][dir];
      const IntVectSet ivs        = grown_set[dit()] & box;
      const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;

      for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();
	interpol(face, comp) = flux(face, comp);
      }
    }

    (*(m_interpolant[dir])[a_lvl]).exchange();
  }
}

void CdrSolver::writePlotFile(){
  CH_TIME("CdrSolver::writePlotFile");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotFile" << endl;
  }


  // Number of output components and their names
  const int ncomps = getNumberOfPlotVariables();
  const Vector<std::string> names = getPlotVariableNames();

  // Allocate storage
  EBAMRCellData output;
  m_amr->allocate(output, m_realm, m_phase, ncomps);
  DataOps::setValue(output, 0.0);

  // Copy internal data to be plotted over to 'output'
  int icomp = 0;
  this->writePlotData(output, icomp);

  // Filename
  char file_char[100];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_timeStep, SpaceDim);

  // Alias
  Vector<LevelData<EBCellFAB>* > output_ptr;
  m_amr->alias(output_ptr, output);
  
  Vector<Real> covered_values(ncomps, 0.0);
  string fname(file_char);
  writeEBHDF5(fname,
	      m_amr->getGrids(m_realm),
	      output_ptr,
	      names,
	      m_amr->getDomains()[0].domainBox(),
	      m_amr->getDx()[0],
	      m_dt,
	      m_time,
	      m_amr->getRefinementRatios(),
	      m_amr->getFinestLevel() + 1,
	      true,
	      covered_values,
	      IntVect::Unit);
}

void CdrSolver::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("CdrSolver::writePlotData");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotData" << endl;
  }

  // Plot state
  if(m_plotPhi) {
    this->writeData(a_output, a_comp, m_phi, true);
  }

  // Plot diffusion coefficients
  if(m_plotDiffusionCoefficient && m_isDiffusive) { // Need to compute the cell-centerd stuff first
    DataOps::setValue(m_scratch, 0.0);
    DataOps::averageFaceToCell(m_scratch, m_faceCenteredDiffusionCoefficient, m_amr->getDomains());
    this->writeData(a_output, a_comp, m_scratch,   false);
  }

  // Plot source terms
  if(m_plotSource) {
    this->writeData(a_output, a_comp, m_source, false);
  }

  // Plot velocities
  if(m_plotVelocity && m_isMobile) {
    this->writeData(a_output, a_comp, m_cellVelocity, false);
  }

  // Plot EB fluxes
  if(m_plotEbFlux && m_isMobile){
    DataOps::setValue(m_scratch, 0.0);
    DataOps::incr(m_scratch, m_ebFlux, 1.0);
    this->writeData(a_output, a_comp, m_scratch, false);
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

void CdrSolver::defineInterpolant(){
  CH_TIME("CdrSolver::defineInterpolant");
  if(m_verbosity > 5){
    pout() << m_name + "::defineInterpolant" << endl;
  }

  const int ncomp = 1;


  for (int dir = 0; dir < SpaceDim; dir++){
    //    Vector<RefCountedPtr<LevelData<BaseIFFAB<Real> > > > interpolant = m_interpolant[SpaceDim];
    //    Vector<RefCountedPtr<LayoutData<BaseIFFAB<Real> > > > grown_sets = m_interpSets[SpaceDim];

    m_interpolant[dir].resize(1 + m_amr->getFinestLevel());
    m_interpSets[dir].resize(1 + m_amr->getFinestLevel());
    
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const ProblemDomain& domain  = m_amr->getDomains()[lvl];
      
      (m_interpolant[dir])[lvl] = RefCountedPtr<LevelData<BaseIFFAB<Real> > > (new LevelData<BaseIFFAB<Real> >());
      (m_interpSets[dir])[lvl] = RefCountedPtr<LayoutData<IntVectSet> > (new LayoutData<IntVectSet>());

      EBArith::defineFluxInterpolant(*(m_interpolant[dir])[lvl],
				     *(m_interpSets[dir])[lvl], 
				     dbl,
				     ebisl,
				     domain,
				     ncomp,
				     dir);
    }
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
