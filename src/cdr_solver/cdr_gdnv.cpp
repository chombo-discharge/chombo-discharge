/*!
   @file cdr_gdnv.cpp
   @brief Implementation of cdr_gdnv.H
   @author Robert Marskar
   @date Jan. 2018
 */

#include "cdr_gdnv.H"
#include "gdnv_outflow_bc.H"
#include "data_ops.H"

#include <ExtrapAdvectBC.H>

#include <EBArith.H>
#include <ParmParse.H>

ExtrapAdvectBCFactory s_physibc;

cdr_gdnv::cdr_gdnv() : cdr_tga() {
  m_class_name = "cdr_gdnv";
  m_name       = "cdr_gdnv";
}

cdr_gdnv::~cdr_gdnv(){

}

void cdr_gdnv::parse_options(){
  CH_TIME("cdr_gdnv::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }
  
  parse_domain_bc();     // Parses domain BC options
  parse_slopelim();      // Parses slope limiter settings
  parse_plot_vars();     // Parses plot variables
  parse_plotmode();      // Parse plot mdoe
  parse_gmg_settings();  // Parses solver parameters for geometric multigrid
  parse_extrap_source(); // Parse source term extrapolation for time-centering advective comps
  parse_rng_seed();      // Get a seed
  parse_conservation();  // Nonlinear divergence blending
}

void cdr_gdnv::parse_slopelim(){
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("limit_slopes", str);
  m_slopelim = (str == "true") ? true : false;
}

int cdr_gdnv::query_ghost() const {
  return 3;
}

void cdr_gdnv::average_velo_to_faces(){
  CH_TIME("cdr_gdnv::average_velo_to_faces(public, full)");
  if(m_verbosity > 5){
    pout() << m_name + "::average_velo_to_faces(public, full)" << endl;
  }
  this->average_velo_to_faces(m_velo_face, m_velo_cell); // Average velocities to face centers for all levels
}

void cdr_gdnv::average_velo_to_faces(EBAMRFluxData& a_velo_face, const EBAMRCellData& a_velo_cell){
  CH_TIME("cdr_gdnv::average_velo_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::average_velo_to_faces" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::average_cell_to_face(*a_velo_face[lvl], *a_velo_cell[lvl], m_amr->get_domains()[lvl]);
    a_velo_face[lvl]->exchange();
  }

  // Fix up boundary velocities to ensure no influx. This is (probably) the easiest way to handle this for cdr_gdnv
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

      EBFluxFAB& velo = (*a_velo_face[lvl])[dit()];
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

void cdr_gdnv::allocate_internals(){
  CH_TIME("cdr_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  cdr_solver::allocate_internals();

  if(m_diffusive){
    this->setup_gmg();
  }

  if(m_mobile){
    const Vector<RefCountedPtr<EBLevelGrid> >& eblgs = m_amr->get_eblg(m_realm, m_phase);
    const Vector<DisjointBoxLayout>& grids           = m_amr->get_grids(m_realm);
    const Vector<int>& ref_ratios                    = m_amr->get_ref_rat();
    const Vector<Real>& dx                           = m_amr->get_dx();
    const int finest_level                           = m_amr->get_finest_level();
    const bool ebcf                                  = m_amr->get_ebcf();
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
    
      m_level_advect[lvl] = RefCountedPtr<EBAdvectLevelIntegrator> (new EBAdvectLevelIntegrator(*eblgs[lvl],
												eblg_coar,
												ref_rat,
												dx[lvl]*RealVect::Unit,
												has_coar,
												has_fine,
												ebcf,
												m_slopelim,
												m_ebis));
    }
  }
}
  
void cdr_gdnv::advect_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state, const Real a_extrap_dt){
  CH_TIME("cdr_gdnv::advect_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::advect_to_faces" << endl;
  }

  // Compute source for extrapolation
  if(m_extrap_source && a_extrap_dt > 0.0){
#if 0 // R.M. April 2020 - disabling this for now. 
    if(m_diffusive){
      const int finest_level = m_amr->get_finest_level();
      Vector<LevelData<EBCellFAB>* > scratchAlias, stateAlias;
      m_amr->alias(scratchAlias, m_scratch);
      m_amr->alias(stateAlias,   a_state);
      m_gmg_solver->computeAMROperator(scratchAlias, stateAlias, finest_level, 0, false);

      // computeAMROperator fucks my ghost cells. 
      m_amr->interp_ghost_pwl(const_cast<EBAMRCellData&> (a_state), m_realm, m_phase);
    }
#endif

    data_ops::copy(m_scratch, m_source);
    m_amr->average_down(m_scratch,     m_realm, m_phase);
    m_amr->interp_ghost_pwl(m_scratch, m_realm, m_phase);
  }
  else{
    data_ops::set_value(m_scratch, 0.0);
  }

  // Extrapolate face-centered state on every level
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_realm, m_phase)[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator();dit.ok(); ++dit){

      EBFluxFAB& extrap_state   = (*a_face_state[lvl])[dit()];
      const EBCellFAB& state    = (*a_state[lvl])[dit()];
      const EBCellFAB& cell_vel = (*m_velo_cell[lvl])[dit()];
      const EBFluxFAB& face_vel = (*m_velo_face[lvl])[dit()];
      const EBCellFAB& source   = (*m_scratch[lvl])[dit()];
      const EBISBox& ebisbox    = ebisl[dit()];
      const Box& box            = dbl.get(dit());
      const Real time           = 0.0;

      EBAdvectPatchIntegrator& ebpatchad = m_level_advect[lvl]->getPatchAdvect(dit());

      ebpatchad.setVelocities(cell_vel, face_vel);
      ebpatchad.setDoingVel(0);
      ebpatchad.setEBPhysIBC(s_physibc);
      ebpatchad.setCurComp(0);

      ebpatchad.extrapolateBCG(extrap_state, state, source, dit(), time, a_extrap_dt);
    }
  }
}
