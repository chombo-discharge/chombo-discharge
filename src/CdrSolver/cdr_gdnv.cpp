/*!
  @file cdr_gdnv.cpp
  @brief Implementation of cdr_gdnv.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "cdr_gdnv.H"
#include "data_ops.H"

#include <ExtrapAdvectBC.H>

#include <EBArith.H>
#include <ParmParse.H>

ExtrapAdvectBCFactory s_physibc;

#include "CD_NamespaceHeader.H"

cdr_gdnv::cdr_gdnv() : cdr_tga() {
  m_className = "cdr_gdnv";
  m_name       = "cdr_gdnv";
}

cdr_gdnv::~cdr_gdnv(){

}

void cdr_gdnv::parseOptions(){
  CH_TIME("cdr_gdnv::parseOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseOptions" << endl;
  }
  
  parseDomainBc();     // Parses domain BC options
  parse_slopelim();      // Parses slope limiter settings
  parsePlotVariables();     // Parses plot variables
  parsePlotMode();      // Parse plot mdoe
  parse_gmg_settings();  // Parses solver parameters for geometric multigrid
  parseExtrapolateSourceTerm(); // Parse source term extrapolation for time-centering advective comps
  parseRngSeed();      // Get a seed
  parseDivergenceComputation();  // Nonlinear divergence blending
}

void cdr_gdnv::parseRuntimeOptions(){
  CH_TIME("cdr_gdnv::parseRuntimeOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  parse_slopelim();
  parsePlotVariables();
  parsePlotMode();
  parse_gmg_settings();
  parseDomainBc();
  parseExtrapolateSourceTerm();
  parseDivergenceComputation();
}

void cdr_gdnv::parse_slopelim(){
  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("limit_slopes", str);
  m_slopelim = (str == "true") ? true : false;
}

int cdr_gdnv::queryGhost() const {
  return 3;
}

void cdr_gdnv::averageVelocityToFaces(){
  CH_TIME("cdr_gdnv::averageVelocityToFaces(public, full)");
  if(m_verbosity > 5){
    pout() << m_name + "::averageVelocityToFaces(public, full)" << endl;
  }
  this->averageVelocityToFaces(m_faceVelocity, m_cellVelocity); // Average velocities to face centers for all levels
}

void cdr_gdnv::averageVelocityToFaces(EBAMRFluxData& a_faceVelocity, const EBAMRCellData& a_cellVelocity){
  CH_TIME("cdr_gdnv::averageVelocityToFaces");
  if(m_verbosity > 5){
    pout() << m_name + "::averageVelocityToFaces" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::average_cell_to_face(*a_faceVelocity[lvl], *a_cellVelocity[lvl], m_amr->getDomains()[lvl]);
    a_faceVelocity[lvl]->exchange();
  }

  // Fix up boundary velocities to ensure no influx. This is (probably) the easiest way to handle this for cdr_gdnv
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

void cdr_gdnv::allocateInternals(){
  CH_TIME("CdrSolver::allocateInternals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocateInternals" << endl;
  }

  CdrSolver::allocateInternals();

  if(m_isDiffusive){
    this->setup_gmg();
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
  
void cdr_gdnv::advect_to_faces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt){
  CH_TIME("cdr_gdnv::advect_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::advect_to_faces" << endl;
  }

  // Compute source for extrapolation
  if(m_extrapolateSourceTerm && a_extrapDt > 0.0){
#if 0 // R.M. April 2020 - disabling this for now. 
    if(m_isDiffusive){
      const int finest_level = m_amr->getFinestLevel();
      Vector<LevelData<EBCellFAB>* > scratchAlias, stateAlias;
      m_amr->alias(scratchAlias, m_scratch);
      m_amr->alias(stateAlias,   a_phi);
      m_gmg_solver->computeAMROperator(scratchAlias, stateAlias, finest_level, 0, false);

      // computeAMROperator fucks my ghost cells. 
      m_amr->interpGhostPwl(const_cast<EBAMRCellData&> (a_phi), m_realm, m_phase);
    }
#endif

    data_ops::copy(m_scratch, m_source);
    m_amr->averageDown(m_scratch,     m_realm, m_phase);
    m_amr->interpGhostPwl(m_scratch, m_realm, m_phase);
  }
  else{
    data_ops::set_value(m_scratch, 0.0);
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
#include "CD_NamespaceFooter.H"
