/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
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

ExtrapAdvectBCFactory s_physibc;

#include "CD_NamespaceHeader.H"

CdrGodunov::CdrGodunov() : CdrTGA() {
  m_className = "CdrGodunov";
  m_name       = "CdrGodunov";
}

CdrGodunov::~CdrGodunov(){

}

void CdrGodunov::parseOptions(){
  CH_TIME("CdrGodunov::parseOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseOptions" << endl;
  }
  
  parseDomainBc();               // Parses domain BC options
  parseSlopeLimiter();           // Parses slope limiter settings
  parsePlotVariables();          // Parses plot variables
  parsePlotMode();               // Parse plot mdoe
  parseMultigridSettings();          // Parses solver parameters for geometric multigrid
  parseExtrapolateSourceTerm();  // Parse source term extrapolation for time-centering advective comps
  parseRngSeed();                // Get a seed
  parseDivergenceComputation();  // Nonlinear divergence blending
}

void CdrGodunov::parseRuntimeOptions(){
  CH_TIME("CdrGodunov::parseRuntimeOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  parseSlopeLimiter();
  parsePlotVariables();
  parsePlotMode();
  parseMultigridSettings();
  parseDomainBc();
  parseExtrapolateSourceTerm();
  parseDivergenceComputation();
}

void CdrGodunov::parseSlopeLimiter(){
  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("limit_slopes", str);
  m_slopelim = (str == "true") ? true : false;
}

int CdrGodunov::queryGhost() const {
  return 3;
}

void CdrGodunov::averageVelocityToFaces(){
  CH_TIME("CdrGodunov::averageVelocityToFaces(public, full)");
  if(m_verbosity > 5){
    pout() << m_name + "::averageVelocityToFaces(public, full)" << endl;
  }
  this->averageVelocityToFaces(m_faceVelocity, m_cellVelocity); // Average velocities to face centers for all levels
}

void CdrGodunov::averageVelocityToFaces(EBAMRFluxData& a_faceVelocity, const EBAMRCellData& a_cellVelocity){
  CH_TIME("CdrGodunov::averageVelocityToFaces");
  if(m_verbosity > 5){
    pout() << m_name + "::averageVelocityToFaces" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    DataOps::averageCellToFace(*a_faceVelocity[lvl], *a_cellVelocity[lvl], m_amr->getDomains()[lvl]);
    a_faceVelocity[lvl]->exchange();
  }

  // Fix up boundary velocities to ensure no influx. This is (probably) the easiest way to handle this for CdrGodunov
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

      EBFluxFAB& velo = (*a_faceVelocity[lvl])[dit()];
      for (int dir = 0; dir < SpaceDim; dir++){
	for (SideIterator sit; sit.ok(); ++sit){
	  Box box   = dbl.get(dit());
	  box &= domain;
	  box.surroundingNodes(dir); // Convert to facebox

	  const int isign = sign(sit());

	  Box cell_box = box;
	  cell_box.shiftHalf(dir, isign);
	  
	  if(!domain.contains(cell_box)){
	    cell_box &= domain;

	    Box bndry_box = adjCellBox(cell_box, dir, sit(), 1);
	    bndry_box.shift(dir, -isign);

	    IntVectSet ivs(bndry_box);
	    for (VoFIterator vofit(ivs, ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit){
	      const VolIndex& vof = vofit();
	      Vector<FaceIndex> faces = ebisl[dit()].getFaces(vof, dir, sit());

	      for (int iface = 0; iface < faces.size(); iface++){
		const FaceIndex& face = faces[iface];

		if(velo[dir](face, 0)*isign < 0.0){
		  //		  velo[dir](face, 0) = 0.0;
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

void CdrGodunov::allocateInternals(){
  CH_TIME("CdrSolver::allocateInternals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocateInternals" << endl;
  }

  CdrSolver::allocateInternals();

  if(m_isDiffusive){
    this->setupMultigrid();
  }

  if(m_isMobile){
    const Vector<RefCountedPtr<EBLevelGrid> >& eblgs = m_amr->getEBLevelGrid(m_realm, m_phase);
    const Vector<DisjointBoxLayout>& grids           = m_amr->getGrids(m_realm);
    const Vector<int>& ref_ratios                    = m_amr->getRefinementRatios();
    const Vector<Real>& dx                           = m_amr->getDx();
    const int finest_level                           = m_amr->getFinestLevel();
    const bool has_ebcf                              = m_amr->getEbCf();
    m_level_advect.resize(1 + finest_level);
  
    for (int lvl = 0; lvl <= finest_level; lvl++){

      const bool has_coar = lvl > 0;
      const bool has_fine = lvl < finest_level;

      int ref_rat = 1;
      EBLevelGrid eblg_coar;
      if(has_coar){
	eblg_coar = *eblgs[lvl-1];
	ref_rat   = ref_ratios[lvl-1];
      }

      const bool forceNoEBCF = !has_ebcf;
    
      m_level_advect[lvl] = RefCountedPtr<EBAdvectLevelIntegrator> (new EBAdvectLevelIntegrator(*eblgs[lvl],
												eblg_coar,
												ref_rat,
												dx[lvl]*RealVect::Unit,
												has_coar,
												has_fine,
												forceNoEBCF,
												m_slopelim,
												m_ebis));
    }
  }
}
  
void CdrGodunov::advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt){
  CH_TIME("CdrGodunov::advectToFaces");
  if(m_verbosity > 5){
    pout() << m_name + "::advectToFaces" << endl;
  }

  // Compute source for extrapolation
  if(m_extrapolateSourceTerm && a_extrapDt > 0.0){
#if 0 // R.M. April 2020 - disabling this for now. 
    if(m_isDiffusive){
      const int finest_level = m_amr->getFinestLevel();
      Vector<LevelData<EBCellFAB>* > scratchAlias, stateAlias;
      m_amr->alias(scratchAlias, m_scratch);
      m_amr->alias(stateAlias,   a_phi);
      m_multigridSolver->computeAMROperator(scratchAlias, stateAlias, finest_level, 0, false);

      // computeAMROperator fucks my ghost cells. 
      m_amr->interpGhostPwl(const_cast<EBAMRCellData&> (a_phi), m_realm, m_phase);
    }
#endif

    DataOps::copy(m_scratch, m_source);
    m_amr->averageDown(m_scratch,     m_realm, m_phase);
    m_amr->interpGhostPwl(m_scratch, m_realm, m_phase);
  }
  else{
    DataOps::setValue(m_scratch, 0.0);
  }

  // Extrapolate face-centered state on every level
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];
    const Real dx                = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator();dit.ok(); ++dit){

      EBFluxFAB& extrap_state   = (*a_facePhi[lvl])[dit()];
      const EBCellFAB& state    = (*a_phi[lvl])[dit()];
      const EBCellFAB& cell_vel = (*m_cellVelocity[lvl])[dit()];
      const EBFluxFAB& face_vel = (*m_faceVelocity[lvl])[dit()];
      const EBCellFAB& source   = (*m_scratch[lvl])[dit()];
      const EBISBox& ebisbox    = ebisl[dit()];
      const Box& box            = dbl.get(dit());
      const Real time           = 0.0;

      EBAdvectPatchIntegrator& ebpatchad = m_level_advect[lvl]->getPatchAdvect(dit());

      ebpatchad.setVelocities(cell_vel, face_vel);
      ebpatchad.setDoingVel(0);
      ebpatchad.setEBPhysIBC(s_physibc);
      ebpatchad.setCurComp(0);

      ebpatchad.extrapolateBCG(extrap_state, state, source, dit(), time, a_extrapDt);
    }
  }
}

#include <CD_NamespaceFooter.H>
