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

cdr_gdnv::cdr_gdnv() : cdr_tga() {
  m_class_name = "cdr_gdnv";
  m_name       = "cdr_gdnv";
}

cdr_gdnv::~cdr_gdnv(){
  this->delete_covered();
}

void cdr_gdnv::parse_options(){
  CH_TIME("cdr_gdnv::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }
  
  parse_domain_bc();     // Parses domain BC options
  parse_mass_redist();   // Parses mass redistribution
  parse_hybrid_div();    // Parses options for hybrid divergence
  parse_slopelim();      // Parses slope limiter settings
  parse_plot_vars();     // Parses plot variables
  parse_plotmode();      // Parse plot mdoe
  parse_gmg_settings();  // Parses solver parameters for geometric multigrid
  parse_extrap_source(); // Parse source term extrapolation for time-centering advective comps
  parse_diffusion();     // Parse diffusion
  parse_rng_seed();
}

void cdr_gdnv::parse_hybrid_div(){
  ParmParse pp(m_class_name.c_str());
  std::string str;
  pp.get("divF_nc", str);
  if(str == "conservative_average"){
    m_nonConsDiv = nonConsDiv::conservative_average;
  }
  else if(str == "covered_face"){
    m_nonConsDiv = nonConsDiv::covered_face;
  }
  else {
    MayDay::Abort("cdr_gdnv::cdr_gdnv - unknown non-conservative divergence type requested");
  }
}

void cdr_gdnv::parse_slopelim(){
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("limit_slopes", str);
  m_slopelim = (str == "true") ? true : false;
}

void cdr_gdnv::parse_extrap_source(){
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("extrap_source", str);
  m_extrap_source = (str == "true") ? true : false;
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

  if(m_nonConsDiv == nonConsDiv::covered_face){
    this->extrapolate_vel_to_covered_faces();              // Extrapolate velocities to covered face centers
  }
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
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
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

  if(m_nonConsDiv == nonConsDiv::covered_face && m_mobile){
    this->delete_covered();
    this->allocate_covered();
  }
  if(m_diffusive){
    this->setup_gmg();
  }

  if(m_mobile){
    const Vector<RefCountedPtr<EBLevelGrid> >& eblgs = m_amr->get_eblg(m_phase);
    const Vector<DisjointBoxLayout>& grids           = m_amr->get_grids();
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

void cdr_gdnv::allocate_covered(){
  CH_TIME("cdr_solver::allocate_covered");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_covered" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const int finest_level = m_amr->get_finest_level();

  //  this->delete_covered();

  m_covered_sets_lo.resize(1 + finest_level);
  m_covered_sets_hi.resize(1 + finest_level);
  m_covered_face_lo.resize(1 + finest_level);
  m_covered_face_hi.resize(1 + finest_level);
  m_covered_velo_lo.resize(1 + finest_level);
  m_covered_velo_hi.resize(1 + finest_level);
  m_covered_phi_lo.resize(1 + finest_level);
  m_covered_phi_hi.resize(1 + finest_level);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];

    m_covered_sets_lo[lvl] = RefCountedPtr<LayoutData<Vector<IntVectSet> > > (new LayoutData<Vector<IntVectSet> >());
    m_covered_sets_hi[lvl] = RefCountedPtr<LayoutData<Vector<IntVectSet> > > (new LayoutData<Vector<IntVectSet> >());
    m_covered_face_lo[lvl] = RefCountedPtr<LayoutData<Vector<Vector<VolIndex> > > > (new LayoutData<Vector<Vector<VolIndex> > >());
    m_covered_face_hi[lvl] = RefCountedPtr<LayoutData<Vector<Vector<VolIndex> > > > (new LayoutData<Vector<Vector<VolIndex> > >());
    m_covered_velo_lo[lvl] = RefCountedPtr<LayoutData<Vector<BaseIVFAB<Real>* > > > (new LayoutData<Vector<BaseIVFAB<Real>* > >());
    m_covered_velo_hi[lvl] = RefCountedPtr<LayoutData<Vector<BaseIVFAB<Real>* > > > (new LayoutData<Vector<BaseIVFAB<Real>* > >());
    m_covered_phi_lo[lvl]  = RefCountedPtr<LayoutData<Vector<BaseIVFAB<Real>* > > > (new LayoutData<Vector<BaseIVFAB<Real>* > >());
    m_covered_phi_hi[lvl]  = RefCountedPtr<LayoutData<Vector<BaseIVFAB<Real>* > > > (new LayoutData<Vector<BaseIVFAB<Real>* > >());

    m_covered_sets_lo[lvl]->define(dbl);
    m_covered_sets_hi[lvl]->define(dbl);
    m_covered_face_lo[lvl]->define(dbl);
    m_covered_face_hi[lvl]->define(dbl);
    m_covered_velo_lo[lvl]->define(dbl);
    m_covered_velo_hi[lvl]->define(dbl);
    m_covered_phi_lo[lvl]->define(dbl);
    m_covered_phi_hi[lvl]->define(dbl);

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      
      Box box = dbl.get(dit());
      box.grow(1);
      box &= domain;

      (*m_covered_sets_lo[lvl])[dit()].resize(SpaceDim);
      (*m_covered_sets_hi[lvl])[dit()].resize(SpaceDim);
      (*m_covered_face_lo[lvl])[dit()].resize(SpaceDim);
      (*m_covered_face_hi[lvl])[dit()].resize(SpaceDim);
      (*m_covered_velo_lo[lvl])[dit()].resize(SpaceDim, NULL);
      (*m_covered_velo_hi[lvl])[dit()].resize(SpaceDim, NULL);
      (*m_covered_phi_lo[lvl])[dit()].resize(SpaceDim, NULL);
      (*m_covered_phi_hi[lvl])[dit()].resize(SpaceDim, NULL);


      // Compute m_covered_sets and m_covered_face
      for (int dir = 0; dir < SpaceDim; dir++){
	IntVectSet irreg_lo;
	IntVectSet irreg_hi;
	EBArith::computeCoveredFaces((*m_covered_face_lo[lvl])[dit()][dir],
				     (*m_covered_sets_lo[lvl])[dit()][dir],
				     irreg_lo,
				     dir,
				     Side::Lo,
				     ebisbox,
				     box);
	EBArith::computeCoveredFaces((*m_covered_face_hi[lvl])[dit()][dir],
				     (*m_covered_sets_hi[lvl])[dit()][dir],
				     irreg_hi,
				     dir,
				     Side::Hi,
				     ebisbox,
				     box);


							   
	(*m_covered_velo_lo[lvl])[dit()][dir] = new BaseIVFAB<Real>((*m_covered_sets_lo[lvl])[dit()][dir], ebgraph, ncomp);
	(*m_covered_velo_hi[lvl])[dit()][dir] = new BaseIVFAB<Real>((*m_covered_sets_hi[lvl])[dit()][dir], ebgraph, ncomp);
	(*m_covered_phi_lo[lvl])[dit()][dir]  = new BaseIVFAB<Real>((*m_covered_sets_lo[lvl])[dit()][dir], ebgraph, ncomp);
	(*m_covered_phi_hi[lvl])[dit()][dir]  = new BaseIVFAB<Real>((*m_covered_sets_hi[lvl])[dit()][dir], ebgraph, ncomp);

	(*m_covered_velo_lo[lvl])[dit()][dir]->setVal(0.0);
	(*m_covered_velo_hi[lvl])[dit()][dir]->setVal(0.0);
	(*m_covered_phi_lo[lvl])[dit()][dir]->setVal(0.0);
	(*m_covered_phi_hi[lvl])[dit()][dir]->setVal(0.0);
      }
    }
  }
}

void cdr_gdnv::delete_covered(){
  CH_TIME("cdr_solver::delete_covered");
  if(m_verbosity > 5){
    pout() << m_name + "::delete_covered" << endl;
  }

  for (int lvl = 0; lvl < m_covered_velo_lo.size(); lvl++){
    const BoxLayout& dbl = m_covered_velo_lo[lvl]->boxLayout();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      for (int dir = 0; dir < SpaceDim; dir++){
	if((*m_covered_velo_lo[lvl])[dit()][dir] != NULL){
	  delete (*m_covered_velo_lo[lvl])[dit()][dir];
	  (*m_covered_velo_lo[lvl])[dit()][dir] = NULL;
	}
	if((*m_covered_velo_hi[lvl])[dit()][dir] != NULL){
	  delete (*m_covered_velo_hi[lvl])[dit()][dir];
	  (*m_covered_velo_hi[lvl])[dit()][dir] = NULL;
	}
	if((*m_covered_phi_lo[lvl])[dit()][dir] != NULL){
	  delete (*m_covered_phi_lo[lvl])[dit()][dir];
	  (*m_covered_phi_lo[lvl])[dit()][dir] = NULL;
	}
	if((*m_covered_phi_hi[lvl])[dit()][dir] != NULL){
	  delete (*m_covered_phi_hi[lvl])[dit()][dir];
	  (*m_covered_phi_hi[lvl])[dit()][dir] = NULL;
	}
	// delete (*m_covered_velo_lo[lvl])[dit()][dir];
	// delete (*m_covered_velo_hi[lvl])[dit()][dir];
	// delete (*m_covered_phi_lo[lvl])[dit()][dir];
	// delete (*m_covered_phi_hi[lvl])[dit()][dir];
      }
    }
  }
}

void cdr_gdnv::extrapolate_vel_to_covered_faces(){
  CH_TIME("cdr_gdnv::extrapolate_vel_to_covered_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::extrapolate_vel_to_covered_faces" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  // Need some temporary storage below
  EBAMRCellData velo;
  m_amr->allocate(velo, m_phase, ncomp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl         = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl              = m_amr->get_ebisl(m_phase)[lvl];
    EBAdvectLevelIntegrator& leveladvect = *m_level_advect[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];

      // Local storage for cell-centered velocity component in direction dir
      //      EBCellFAB veldir(ebisbox, box, ncomp);
      EBCellFAB veldir(ebisbox, (*m_velo_cell[lvl])[dit()].getRegion(), ncomp);
      for (int dir = 0; dir < SpaceDim; dir++){
	veldir.setVal(0.0);
	veldir.plus((*m_velo_cell[lvl])[dit()], dir, comp, ncomp);

	EBAdvectPatchIntegrator& patcher = leveladvect.getPatchAdvect(dit());

	patcher.extrapToCoveredFaces(*(*m_covered_velo_lo[lvl])[dit()][dir],
				     (*m_velo_face[lvl])[dit()][dir],
				     veldir,
				     (*m_covered_face_lo[lvl])[dit()][dir],
				     dir,
				     Side::Lo,
				     box);

	patcher.extrapToCoveredFaces(*(*m_covered_velo_hi[lvl])[dit()][dir],
				     (*m_velo_face[lvl])[dit()][dir],
				     veldir,
				     (*m_covered_face_hi[lvl])[dit()][dir],
				     dir,
				     Side::Hi,
				     box);

      }
    }
  }
}
  
void cdr_gdnv::advect_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state, const Real a_extrap_dt){
  CH_TIME("cdr_gdnv::advect_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::advect_to_faces" << endl;
  }

  // Extrapolate phi and vel to covered faces first
  if(m_nonConsDiv == nonConsDiv::covered_face){
    this->extrapolate_vel_to_covered_faces();
  }

  const int finest_level = m_amr->get_finest_level();

  const int comp = 0;
  
  EBAdvectPatchIntegrator::setCurComp(0);
  EBAdvectPatchIntegrator::setDoingVel(0);

  RefCountedPtr<ExtrapAdvectBCFactory> bcfact = RefCountedPtr<ExtrapAdvectBCFactory>
    (new ExtrapAdvectBCFactory());

  // This data is used as a source term for time-extrapolation to edges 
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_phase, 1);
  data_ops::set_value(scratch, 0.0);
  
  if(a_extrap_dt > 0.0 && m_extrap_source){
    //    compute_divD(scratch, a_state);
    if(m_diffusive){
      Vector<LevelData<EBCellFAB>* > scratchAlias, stateAlias;
      m_amr->alias(scratchAlias, scratch);
      m_amr->alias(stateAlias,   a_state);
      m_gmg_solver->computeAMROperator(scratchAlias, stateAlias, finest_level, 0, false);
    }
    data_ops::incr(scratch, m_source, 1.0);

    m_amr->average_down(scratch, m_phase);
    m_amr->interp_ghost(scratch, m_phase);
  }

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const RefCountedPtr<EBAdvectLevelIntegrator>& leveladvect = m_level_advect[lvl];
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

    CH_assert(!leveladvect.isNull());

    const bool has_coar = lvl > 0;
    LevelData<EBCellFAB>* coarstate_old = NULL;
    LevelData<EBCellFAB>* coarstate_new = NULL;
    LevelData<EBCellFAB>* coarvel_old   = NULL;
    LevelData<EBCellFAB>* coarvel_new   = NULL;
    LevelData<EBCellFAB>* coarsrc_old   = NULL;
    LevelData<EBCellFAB>* coarsrc_new   = NULL;

    if(has_coar){
      coarstate_old = a_state[lvl-1];
      coarstate_new = a_state[lvl-1];
      coarvel_old   = m_velo_cell[lvl-1];
      coarvel_new   = m_velo_cell[lvl-1];
      coarsrc_old   = scratch[lvl-1];
      coarsrc_new   = scratch[lvl-1];
    }

    leveladvect->resetBCs(bcfact);

    if(m_nonConsDiv == nonConsDiv::conservative_average){
      leveladvect->advectToFacesBCG(*a_face_state[lvl],
				    *a_state[lvl],
				    *m_velo_cell[lvl],
				    *m_velo_face[lvl],
				    coarstate_old,
				    coarstate_new,
				    coarvel_old,
				    coarvel_new,
				    m_time,
				    m_time,
				    m_time,
				    a_extrap_dt,
				    scratch[lvl],
				    coarsrc_old,
				    coarsrc_new);
    }
    else if(m_nonConsDiv == nonConsDiv::covered_face){
      leveladvect->advectToFacesCol(*a_face_state[lvl],
				    *m_covered_phi_lo[lvl],
				    *m_covered_phi_hi[lvl],
				    *m_covered_face_lo[lvl],
				    *m_covered_face_hi[lvl],
				    *m_covered_sets_lo[lvl],
				    *m_covered_sets_hi[lvl],
				    *a_state[lvl],
				    *m_velo_cell[lvl],
				    *m_velo_face[lvl],
				    coarstate_old,
				    coarstate_new,
				    coarvel_old,
				    coarvel_new,
				    m_time,
				    m_time,
				    m_time,
				    a_extrap_dt,
				    scratch[lvl],
				    coarsrc_old,
				    coarsrc_new);
    }
    else{
      MayDay::Abort("cdr_gdnv::extrapolate_to_faces - unknown method requested");
    }
  }
}

void cdr_gdnv::nonconservative_divergence(EBAMRIVData&         a_div_nc,
					  const EBAMRCellData& a_divF,
					  const EBAMRFluxData& a_face_state){
  CH_TIME("cdr_gdnv::nonconservative_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::nonconservative_divergence" << endl;
  }

  if(m_nonConsDiv = nonConsDiv::conservative_average){
    cdr_solver::nonconservative_divergence(a_div_nc, a_divF);
  }
  else if(m_nonConsDiv = nonConsDiv::covered_face){
    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl         = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl              = m_amr->get_ebisl(m_phase)[lvl];
      EBAdvectLevelIntegrator& leveladvect = *(m_level_advect[lvl]);

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const Box box          = dbl.get(dit());
	const EBISBox& ebisbox = ebisl[dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();
	const IntVectSet& ivs  = ebisbox.getIrregIVS(box);
      
	EBAdvectPatchIntegrator& patcher = leveladvect.getPatchAdvect(dit());
	EBCellFAB divF(ebisbox, box, ncomp); // Temporary storage

	divF.setVal(0.0);
	patcher.advectiveDerivative(divF,
				    (*a_face_state[lvl])[dit()],
				    (*m_velo_face[lvl])[dit()],
				    (*m_covered_phi_lo[lvl])[dit()],
				    (*m_covered_phi_hi[lvl])[dit()],
				    (*m_covered_velo_lo[lvl])[dit()],
				    (*m_covered_velo_hi[lvl])[dit()],
				    (*m_covered_face_lo[lvl])[dit()],
				    (*m_covered_face_hi[lvl])[dit()],
				    box);

	// Copy results for irregular cells over to a_div_nc
	for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	  (*a_div_nc[lvl])[dit()](vof, comp) = divF(vof, comp);
	}
      }
    }
  }
  else{
    MayDay::Abort("cdr_gdnv::nonconservative_divergence - unknown type for div(F)_nc");
  }
}

void cdr_gdnv::compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt){
  CH_TIME("cdr_gdnv::compute_divJ(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divJ(divF, state)" << endl;
  }

  if(m_nonConsDiv == nonConsDiv::conservative_average){
    cdr_tga::compute_divJ(a_divJ, a_state, a_extrap_dt);
  }
  else if(m_nonConsDiv == nonConsDiv::covered_face){

    const int comp  = 0;
    const int ncomp = 1;

    // Have to do advection and diffusion separately because the use differnet nonconservative divergences
    if(m_mobile || m_diffusive){
      EBAMRCellData divG;
      if(m_mobile){
	m_amr->allocate(divG, m_phase, ncomp);
	compute_divF(divG, a_state, a_extrap_dt);
	data_ops::incr(a_divJ, divG, 1.0);
      }
      if(m_diffusive){
	compute_divD(divG, a_state);
	data_ops::scale(divG, -1.0);
	data_ops::incr(a_divJ, divG, 1.0);
      }
    }
    else{
      data_ops::set_value(a_divJ, 0.0);
    }
  }

  return;
}

void cdr_gdnv::compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt){
  CH_TIME("cdr_gdnv::compute_divF");
  if(m_verbosity > 5){
    pout() << "cdr_gdnv::compute_divF" << endl;
  }


  if(m_mobile){
    if(m_nonConsDiv == nonConsDiv::conservative_average){
      cdr_tga::compute_divF(a_divF, a_state, a_extrap_dt);
    }
    else if(m_nonConsDiv == nonConsDiv::covered_face){
      // TLDR: This is almost a complete overwwrite of cdr_tga::compute_divF, with the exception that
      //       the nonconservative divergence is computed differentenly

      const int comp       = 0;
      const int ncomp      = 1;
      const int redist_rad = m_amr->get_redist_rad();

      EBAMRFluxData face_state;
      EBAMRFluxData flux;
      EBAMRIVData   div_nc;
      EBAMRIVData   mass_diff;
      EBAMRCellData weights;

      m_amr->allocate(face_state, m_phase, ncomp);
      m_amr->allocate(flux, m_phase, ncomp);
      m_amr->allocate(div_nc,     m_phase, ncomp);
      m_amr->allocate(mass_diff,  m_phase, ncomp);
      m_amr->allocate(weights,    m_phase, ncomp, 2*redist_rad);

      data_ops::set_value(a_divF,     0.0); 
      data_ops::set_value(face_state, 0.0);
      data_ops::set_value(div_nc,     0.0);
      data_ops::set_value(mass_diff,  0.0);
      data_ops::set_value(weights,    0.0);

      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	a_state[lvl]->localCopyTo(*weights[lvl]);
      }

      // Compute the advective derivative
      this->average_velo_to_faces();
      this->advect_to_faces(face_state, a_state, a_extrap_dt);          // Face extrapolation to cell-centered faces
      this->compute_flux(flux, m_velo_face, face_state, m_domainflux);  // Compute face-centered fluxes
      this->conservative_divergence(a_divF, flux, m_ebflux);            // Compute conservative divergence
      this->nonconservative_divergence(div_nc, a_divF, face_state);     // Compute alternative non-conservative divergence
      this->hybrid_divergence(a_divF, mass_diff, div_nc);               // Make divF = hybrid divergence. Compute mass diff.
      this->increment_flux_register(flux);                              // Increment flux registers
      this->increment_redist(mass_diff);                                // Increment redistribution objects

      const bool ebcf = m_amr->get_ebcf();
      if(ebcf){ 
	this->coarse_fine_increment(mass_diff);                         // Increment the coarse-fine redistribution objects
	this->increment_redist_flux();                                  // Increment flux registers with the redistribution stuff
	this->hyperbolic_redistribution(a_divF, mass_diff, weights);    // Redistribute mass into hybrid divergence
	this->coarse_fine_redistribution(a_divF);                       // Redistribute
	this->reflux(a_divF);
      }
      else{
	this->hyperbolic_redistribution(a_divF, mass_diff, weights);    // Redistribute that was left out
	this->reflux(a_divF);                                           // Reflux at coarse-fine interfaces
      }
    }
    else{
      MayDay::Abort("cdr_gdnv::compute_divF - unknown nonconsdiv");
    }
  }
  else{
    data_ops::set_value(a_divF, 0.0);
  }
}

void cdr_gdnv::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("cdr_gdnv::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  this->delete_covered();

  cdr_solver::regrid(a_lmin, a_old_finest_level, a_new_finest_level);
}
