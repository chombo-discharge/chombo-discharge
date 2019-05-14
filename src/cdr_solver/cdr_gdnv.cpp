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
  m_name = "cdr_gdnv";
  m_slopelim = false;
  std::string str = "covered_face";
  int which = 0;

  // Get options from input script
  {
    ParmParse pp("cdr_gdnv");
    pp.query("divF_nc", str);
    if(str == "conservative_average"){
      which = 0;
      MayDay::Warning("cdr_gdnv::cdr_gdnv - we are decommissioning the conservative average stuff for this solver");
    }
    else if(str == "covered_face"){
      which = 1;
    }
    else {
      MayDay::Abort("cdr_gdnv::cdr_gdnv - unknown non-conservative divergence type requested");
    }

    if(pp.contains("limit_slopes")){
      pp.get("limit_slopes", str);
      if(str == "true"){
	m_slopelim = true;
      }
      else if(str == "false"){
	m_slopelim = false;
      }
    }
  }
  
  this->set_divF_nc(which);
}

cdr_gdnv::~cdr_gdnv(){
  this->delete_covered();
}

int cdr_gdnv::query_ghost() const {
  return 3;
}

void cdr_gdnv::advance_advect(EBAMRCellData& a_state, const Real a_dt){
  CH_TIME("cdr_gdnv::advance_advect");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_advect" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();
  
  EBAMRCellData divF;
  m_amr->allocate(divF, m_phase, ncomp);

  
  this->compute_divF(divF, a_state, 0.5*a_dt, true); // Compute div(n*v). Use semi-implicit extrapolation to half time steps
  data_ops::incr(a_state, divF, a_dt);               // phi_new = phi_old + dt*divF

  m_amr->average_down(a_state, m_phase);
  m_amr->interp_ghost(a_state, m_phase);

  data_ops::floor(a_state, 0.0);
}

void cdr_gdnv::average_velo_to_faces(){
  CH_TIME("cdr_gdnv::average_velo_to_faces(public, full)");
  if(m_verbosity > 5){
    pout() << m_name + "::average_velo_to_faces(public, full)" << endl;
  }

  this->average_velo_to_faces(m_velo_face, m_velo_cell); // Average velocities to face centers for all levels
  this->extrapolate_vel_to_covered_faces();              // Extrapolate velocities to covered face centers
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

void cdr_gdnv::set_divF_nc(const int a_which_divFnc){

  CH_assert(a_which_divFnc == 0 || a_which_divFnc == 1);
  m_which_divFnc = a_which_divFnc;
}

void cdr_gdnv::compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt){
  cdr_tga::compute_divJ(a_divJ, a_state, a_extrap_dt);
}

void cdr_gdnv::allocate_internals(){
  CH_TIME("cdr_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  cdr_solver::allocate_internals();

  if(m_which_divFnc == 1){
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
  if(m_which_divFnc == 1){
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
#if 0 // Development feature
  if(a_extrap_dt > 0.0){
    compute_divD(scratch, a_state);
    data_ops::incr(scratch, m_source, 1.0);

    m_amr->average_down(scratch, m_phase);
    m_amr->interp_ghost(scratch, m_phase);
  }
#endif

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

    if(m_which_divFnc == 0){
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
    else if(m_which_divFnc == 1){
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

  if(m_which_divFnc == 0){
    cdr_solver::nonconservative_divergence(a_div_nc, a_divF, a_face_state);
    return;
  }
  else if(m_which_divFnc == 1){
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

void cdr_gdnv::nonconservative_divergence(LevelData<BaseIVFAB<Real> >& a_divF_nc,
					  const LevelData<EBFluxFAB>&  a_face_state,
					  const int                    a_lvl){

  const int comp         = 0;
  const int ncomp        = 1;

  if(m_which_divFnc == 0){
    MayDay::Abort("cdr_gdnv::nonconservative_divergence(per level) - this only works for covered cells stuff");
  }

  const DisjointBoxLayout& dbl         = m_amr->get_grids()[a_lvl];
  const EBISLayout& ebisl              = m_amr->get_ebisl(m_phase)[a_lvl];
  EBAdvectLevelIntegrator& leveladvect = *(m_level_advect[a_lvl]);

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet& ivs  = ebisbox.getIrregIVS(box);
      
    EBAdvectPatchIntegrator& patcher = leveladvect.getPatchAdvect(dit());
    EBCellFAB divF(ebisbox, box, ncomp); // Temporary storage

    patcher.advectiveDerivative(divF,
				a_face_state[dit()],
				(*m_velo_face[a_lvl])[dit()],
				(*m_covered_phi_lo[a_lvl])[dit()],
				(*m_covered_phi_hi[a_lvl])[dit()],
				(*m_covered_velo_lo[a_lvl])[dit()],
				(*m_covered_velo_hi[a_lvl])[dit()],
				(*m_covered_face_lo[a_lvl])[dit()],
				(*m_covered_face_hi[a_lvl])[dit()],
				box);

    // Copy results for irregular cells over to a_div_nc
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      a_divF_nc[dit()](vof, comp) = divF(vof, comp);
    }
  }
}

void cdr_gdnv::regrid(const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("cdr_gdnv::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  this->delete_covered();

  cdr_solver::regrid(a_old_finest_level, a_new_finest_level);
}


void cdr_gdnv::eulerF_subcycle(EBAMRCellData& a_state, const Real a_dt, const bool a_redist){
  CH_TIME("cdr_gdnv::eulerF_subcycle(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::eulerF_subcycle(divF, state)" << endl;
  }

  if(m_mobile){
    const int comp       = 0;
    const int ncomp      = 1;
    const int redist_rad = m_amr->get_redist_rad();

    EBAMRFluxData flux;
    EBAMRFluxData face_state;
    EBAMRCellData divF_c;
    EBAMRCellData weights;
    EBAMRIVData   divF_nc;
    EBAMRIVData   mass_diff;

    EBAMRCellData coar_old;
    EBAMRCellData coar_new;
    
    m_amr->allocate(face_state, m_phase, ncomp);
    m_amr->allocate(flux,       m_phase, ncomp);
    m_amr->allocate(divF_nc,    m_phase, ncomp);
    m_amr->allocate(mass_diff,  m_phase, ncomp);
    m_amr->allocate(divF_c,     m_phase, ncomp);
    m_amr->allocate(coar_old,   m_phase, ncomp);
    m_amr->allocate(coar_new,   m_phase, ncomp);
    m_amr->allocate(weights,    m_phase, ncomp, 2*redist_rad);

    data_ops::set_value(face_state, 0.0);
    data_ops::set_value(flux,       0.0);
    data_ops::set_value(divF_nc,    0.0);
    data_ops::set_value(divF_c,     0.0);
    data_ops::set_value(coar_old,   0.0);
    data_ops::set_value(mass_diff,  0.0);
    data_ops::set_value(weights,    0.0);

    this->average_velo_to_faces(m_velo_face, m_velo_cell); // Average velocities to face centers for all levels
    this->extrapolate_vel_to_covered_faces();              // Extrapolate velocities to covered face centers

    // We now have all the things that we need. Cycle through levels
    const int coar_lvl = 0;
    const int fine_lvl = m_amr->get_finest_level();

    // Do subcycle advance between coarsest_level and finest_level. Use scratch storage.
    Vector<Real> times_new(1+fine_lvl, 0.0); // These are the CURRENT grid times 
    Vector<Real> times_old(1+fine_lvl, 0.0); // These are the PREVIOUS grid times
    advance_subcycle_amr(a_state,
			 flux,
			 face_state,
			 coar_old,
			 divF_c,
			 weights,
			 divF_nc,
			 mass_diff,
			 coar_lvl,
			 times_new,
			 times_old,
			 coar_lvl,
			 fine_lvl,
			 a_dt);

  }
}

void cdr_gdnv::advect_to_faces(LevelData<EBFluxFAB>&        a_face_state,
			       const LevelData<EBCellFAB>&  a_cell_state,
			       const LevelData<EBCellFAB>*  a_state_coarse_old,
			       const LevelData<EBCellFAB>*  a_state_coarse_new,
			       const Real                   a_time,
			       const Real                   a_coarse_time_old,
			       const Real                   a_coarse_time_new,
			       const int                    a_lvl,
			       const Real                   a_extrap_dt){
  CH_TIME("cdr_gdnv::advect_to_faces(subcycle)");
  if(m_verbosity > 5){
    pout() << m_name + "::advect_to_faces(subcycle)" << endl;
  }


  const int comp = 0;
  const bool has_coar = a_lvl > 0;
  
  EBAdvectPatchIntegrator::setCurComp(0);
  EBAdvectPatchIntegrator::setDoingVel(0);

  RefCountedPtr<ExtrapAdvectBCFactory> bcfact = RefCountedPtr<ExtrapAdvectBCFactory> (new ExtrapAdvectBCFactory());

  const RefCountedPtr<EBAdvectLevelIntegrator>& leveladvect = m_level_advect[a_lvl];
  
  CH_assert(!leveladvect.isNull());
  
  leveladvect->resetBCs(bcfact);


  // Level source on this level and coarser level
  data_ops::set_value(m_scratch, 0.0);
  if(a_extrap_dt > 0.0){
#if 0 // Development feature....
    // Compute DivD on this level
    if(m_diffusive){
      AMRLevelOp<LevelData<EBCellFAB> >* oper = m_gmg_solver->getAMROperators()[a_lvl];
      oper->applyOp(*m_scratch[a_lvl], a_cell_state);


      if(has_coar){
	AMRLevelOp<LevelData<EBCellFAB> >* coar_oper = m_gmg_solver->getAMROperators()[a_lvl-1];
	oper->applyOp(*m_resid[a_lvl-1], a_state_coar_old);
	data_ops::incr(*m_scratch[a_lvl-1], *m_resid[a_lvl-1], 1.0);

      }

    }
    data_ops::incr(*m_scratch[a_lvl], *m_source[a_lvl], 1.0);
    if(has_coar){
      data_ops::incr(*m_scratch[a_lvl-1], *m_source[a_lvl-1], 1.0);
    }
#endif
  }

  const LevelData<EBCellFAB>* source_ptr   = m_scratch[a_lvl];
  const LevelData<EBCellFAB>* coar_src_old = has_coar ? &(*m_scratch[a_lvl-1])   : NULL;
  const LevelData<EBCellFAB>* coar_src_new = has_coar ? &(*m_scratch[a_lvl-1])   : NULL;
  const LevelData<EBCellFAB>* coar_vel_old = has_coar ? &(*m_velo_cell[a_lvl-1]) : NULL;
  const LevelData<EBCellFAB>* coar_vel_new = has_coar ? &(*m_velo_cell[a_lvl-1]) : NULL;

  // Extrapolate to faces. Interpolate to a_time between a_coarse_time_old and a_coarse_time new
  if(m_which_divFnc == 0){
    leveladvect->advectToFacesBCG(a_face_state,
				  a_cell_state,
				  *m_velo_cell[a_lvl],
				  *m_velo_face[a_lvl],
				  a_state_coarse_old,
				  a_state_coarse_new,
				  coar_vel_old,
				  coar_vel_new,
				  a_coarse_time_old,
				  a_coarse_time_new,
				  a_time,
				  a_extrap_dt,
				  source_ptr,
				  coar_src_old,
				  coar_src_new);
  }
  else if(m_which_divFnc == 1){
    leveladvect->advectToFacesCol(a_face_state,
				  *m_covered_phi_lo[a_lvl],
				  *m_covered_phi_hi[a_lvl],
				  *m_covered_face_lo[a_lvl],
				  *m_covered_face_hi[a_lvl],
				  *m_covered_sets_lo[a_lvl],
				  *m_covered_sets_hi[a_lvl],
				  a_cell_state,
				  *m_velo_cell[a_lvl],
				  *m_velo_face[a_lvl],
				  a_state_coarse_old,
				  a_state_coarse_new,
				  coar_vel_old,
				  coar_vel_new,
				  a_coarse_time_old,
				  a_coarse_time_new,
				  a_time,
				  a_extrap_dt,
				  source_ptr,
				  coar_src_old,
				  coar_src_new);
  }
  else{
    MayDay::Abort("cdr_gdnv::extrapolate_to_faces - unknown method requested");
  }
}

void cdr_gdnv::advance_subcycle_amr(EBAMRCellData& a_state,
				    EBAMRFluxData& a_flux,
				    EBAMRFluxData& a_face_states,
				    EBAMRCellData& a_coar_old_state,
				    EBAMRCellData& a_divF_c,
				    EBAMRCellData& a_weights,
				    EBAMRIVData&   a_divF_nc,
				    EBAMRIVData&   a_mass_diff,
				    int            a_lvl,
				    Vector<Real>&  a_tnew,
				    Vector<Real>&  a_told,
				    const int      a_coarsest_level,
				    const int      a_finest_level,
				    const Real     a_dt){

  // TLDR: When we come in here we want to advance a_dt using subcycling over finer grid levels. Normal and advective velocities
  //       have been filled, and boundary conditions are fixed. 

  // Grid times for the coarser grid
  Real coar_time_old = 0.0;
  Real coar_time_new = 0.0;

  const bool has_fine = a_lvl < a_finest_level;
  const bool has_coar = a_lvl > a_coarsest_level;

  LevelData<EBCellFAB>* coar_state_old = NULL;
  LevelData<EBCellFAB>* coar_state_new = NULL;
  if(has_coar){
    
    coar_time_old = a_told[a_lvl-1];
    coar_time_new = a_tnew[a_lvl-1];

    coar_state_old = &(*a_coar_old_state[a_lvl-1]);
    coar_state_new = &(*a_state[a_lvl-1]);

  }

  // This advances this level from a_tnew[a_lvl] to a_tnew[a_lvl] + a_dt. Copy the old state to a_coar_old_state first
  a_state[a_lvl]->copyTo(*a_coar_old_state[a_lvl]); // Need to back up old states
  advance_subcycle_level(*a_state[a_lvl],
			 *a_flux[a_lvl],
			 *a_face_states[a_lvl],
  			 *a_divF_c[a_lvl],
  			 *a_weights[a_lvl],
  			 *a_divF_nc[a_lvl],
  			 *a_mass_diff[a_lvl],
  			 coar_state_old,
  			 coar_state_new,
  			 a_lvl,
  			 a_coarsest_level,
  			 a_finest_level,
  			 coar_time_old,
  			 coar_time_new,
			 a_tnew[a_lvl],
  			 a_dt);

  // We have advanced this level. Update new times
  a_told[a_lvl] = a_tnew[a_lvl];
  a_tnew[a_lvl] = a_tnew[a_lvl] + a_dt;

  // If there is a coarser level, advance it nref times so that we synchronize. a_state[lvl]

  if(has_fine){ 
    const int nref    = m_amr->get_ref_rat()[a_lvl];
    const Real dt_ref = a_dt/(nref);
    for (int i=0; i < nref; i++){
      advance_subcycle_amr(a_state, a_flux, a_face_states, a_coar_old_state, a_divF_c,a_weights, a_divF_nc,
			   a_mass_diff, a_lvl+1, a_tnew, a_told, a_coarsest_level, a_finest_level, dt_ref);
    }

    // Finer level has reached this level, average down finer level onto this one and reflux
    m_amr->average_down(a_state, m_phase, a_lvl);
    reflux_level(a_state, a_lvl, a_coarsest_level, a_finest_level, 1.0);
  }
}

void cdr_gdnv::advance_subcycle_level(LevelData<EBCellFAB>&        a_state,
				      LevelData<EBFluxFAB>&        a_flux,
				      LevelData<EBFluxFAB>&        a_face_states,
				      LevelData<EBCellFAB>&        a_divF_c,
				      LevelData<EBCellFAB>&        a_weights,
				      LevelData<BaseIVFAB<Real> >& a_mass_diff,
				      LevelData<BaseIVFAB<Real> >& a_divF_nc,
				      const LevelData<EBCellFAB>*  a_coar_old_state,
				      const LevelData<EBCellFAB>*  a_coar_new_state,
				      const int                    a_lvl,
				      const int                    a_coarsest_level,
				      const int                    a_finest_level,
				      const Real                   a_coar_time_old,
				      const Real                   a_coar_time_new,
				      const Real                   a_time,
				      const Real                   a_dt){

  // This is the level advance. It advances with hyperbolic redistribution and increments flux registers on the way. 
  //
  // The actual advance is state = state - dt*div(F). Reflux registers are incremented 
  //
  // 

  const bool ebcf = m_amr->get_ebcf();
  if(ebcf){
    MayDay::Abort("cdr_gdnv::advance_subcycle_level - ebcf not (yet) supported");
  }

  //  a_state.exchange(); // Necessary???

  // Advect to faces and compute fluxes on face centers, and compute the conservative divergence on regular cells
  advect_to_faces(a_face_states, a_state, a_coar_old_state, a_coar_new_state, a_time, a_coar_time_old, a_coar_time_new, a_lvl);
  new_compute_flux(a_flux, a_face_states, *m_velo_face[a_lvl], *m_domainflux[a_lvl], a_lvl);
  consdiv_regular(a_divF_c, a_flux, a_lvl);

  // Set up flux interpolant and compute conservative flux on irregular cells. 
  LevelData<BaseIFFAB<Real> > flux_interp[SpaceDim];
  setup_flux_interpolant(flux_interp, a_flux, a_lvl); // Compute interpolant
  interpolate_flux_to_centroids(flux_interp, a_lvl);  // Interpolant now holds face centroid-centered fluxes
  compute_divF_irreg(a_divF_c, flux_interp, *m_ebflux[a_lvl], a_lvl); // Update the conservative divergence in the cut cells

  // Recess: So far the conservative divergence is scaled by 1/dx but not yet divided by the volume fraction. This
  // means that the actual advance without the hybrid stuff would be
  //
  // new_state -= dt*a_divF_c/kappa
  //
  // The stuff below was originally written for d(phi)/dt = Div(F) rather than d(phi)/dt = -Div(F) so that's why
  // there's a (-a_dt) in all the stuff below. This design choice was made because I am, in fact, an ass. 

  // Compute the nonconservative and hybrid divergences (hybrid put on storage for conservative, which is lost)
  nonconservative_divergence(a_divF_nc, a_face_states, a_lvl);
  hybrid_divergence(a_divF_c, a_mass_diff, a_divF_nc, a_lvl); // Puts hybrid in a_divF_c. mass_diff as usual without dt, but it
  data_ops::scale(a_mass_diff, -a_dt);                        // Was given the opposite sign of the Chombo convention

  // Update flux registers and redistribution registers
  update_flux_registers(a_flux, a_lvl, a_coarsest_level, a_finest_level, a_dt);
  update_redist_register(a_mass_diff, a_lvl);
  if(m_mass_redist){
    data_ops::incr(a_weights, a_state, 1.0);
    a_weights.exchange();
  }

  // Euler advance and then redistribute mass on this level
  data_ops::incr(a_state, a_divF_c, -a_dt);
  a_state.exchange();
  redist_level(a_state, a_weights, a_lvl);

}

void cdr_gdnv::update_flux_registers(LevelData<EBFluxFAB>& a_flux,
				     const int             a_lvl,
				     const int             a_coarsest_level,
				     const int             a_finest_level,
				     const Real            a_dt){

#if 0 // Debug
  if(procID() == 0){
    std::cout << "update fluxreg on a_lvl = " << a_lvl << "\t dt = " << a_dt << std::endl;
  }
#endif

    // Increment the coarser flux register and initialize the finer flux register. a_flux holds phi*vel which we can use
  const bool has_fine = a_lvl < a_finest_level;
  const bool has_coar = a_lvl > a_coarsest_level;

  const int comp = 0;
  const Interval interv(comp, comp);
  EBFluxRegister* fluxreg_fine = NULL;
  EBFluxRegister* fluxreg_coar = NULL;

  // Remember, register on a_lvl holds flux between level a_lvl and a_lvl+1
  if(has_fine) fluxreg_fine = m_amr->get_flux_reg(m_phase)[a_lvl];   
  if(has_coar) fluxreg_coar = m_amr->get_flux_reg(m_phase)[a_lvl-1]; 

  // Initialize flux register
  if(has_fine) {

    fluxreg_fine->setToZero();
  }

#if 0 // Debug
  if(procID() == 0){
    std::cout << "------------------------------ " << std::endl;
    std::cout << "Updating flux registers from " << a_lvl << std::endl;
    if(has_fine){
      std::cout << "Resetting flux reg on level = " << a_lvl << std::endl;
      std::cout << "Incrementing flux reg on level = " << a_lvl << " with coarse flux" << std::endl;
    }
    if(has_coar){
      std::cout << "incrementing flux reg on coarser level = " << a_lvl-1 << " with fine flux" << std::endl;
    }
    std::cout << "------------------------------ " << std::endl;
    std::cout << std::endl;
  }
#endif

  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const EBISBox& ebisbox = ebisl[dit()];
    const Box box          = dbl.get(dit());
      
    for (int dir = 0; dir < SpaceDim; dir++){

      // Initialize finer flux register with flux leaving this level and into the finer level (or invalid region of that level)
      if(has_fine) {
	for (SideIterator sit; sit.ok(); ++sit){
	  fluxreg_fine->incrementCoarseBoth(a_flux[dit()][dir], a_dt, dit(), interv, dir, sit());
	}
      }

      // Increment coarser flux register with flux entering that level from this level 
      if(has_coar){ // The coarser level has already been initialized with the coarse side flux
	for (SideIterator sit; sit.ok(); ++sit){
	  fluxreg_coar->incrementFineBoth(a_flux[dit()][dir], a_dt, dit(), interv, dir, sit());
	}
      }
    }
  }
}

void cdr_gdnv::update_redist_register(const LevelData<BaseIVFAB<Real> >& a_mass_diff, const int a_lvl){
  CH_TIME("cdr_solver::update_redist_register");
  if(m_verbosity > 5){
    pout() << m_name + "::update_redist_register" << endl;
  }

  const int comp  = 0;
  const Interval interv(comp, comp);

  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
    
  EBLevelRedist& level_redist = *(m_amr->get_level_redist(m_phase)[a_lvl]);
  level_redist.setToZero();

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    level_redist.increment(a_mass_diff[dit()], dit(), interv);
  }
}

void cdr_gdnv::reflux_level(EBAMRCellData& a_state,
			    const int      a_lvl,
			    const int      a_coarsest_level,
			    const int      a_finest_level,
			    const Real     a_scale){
  CH_TIME("cdr_gdnv::reflux");
  if(m_verbosity > 5){
    pout() << m_name + "::reflux(cdr_gdnv, level stuff)" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  // Remember, the flux register for smooshing mass from a_level+1 to a_level lives on a_level
  RefCountedPtr<EBFluxRegister >& fluxreg = m_amr->get_flux_reg(m_phase)[a_lvl];

  // Reflux this level, then empty the register.
  const Real dx = m_amr->get_dx()[a_lvl];
  //  fluxreg->reflux(*a_state[a_lvl], interv, a_scale);
  fluxreg->reflux(*a_state[a_lvl], interv, 1./dx);
  fluxreg->setToZero();
}

void cdr_gdnv::redist_level(LevelData<EBCellFAB>& a_state, const LevelData<EBCellFAB>& a_weights, const int a_lvl){
  CH_TIME("cdr_gdnv::redist_level");
  if(m_verbosity > 5){
    pout() << m_name + "::redist_level" << endl;
  }

  const int comp = 0;
  const Interval interv(comp, comp);

  EBLevelRedist& level_redist = *(m_amr->get_level_redist(m_phase)[a_lvl]);
  if(m_mass_redist){
    level_redist.resetWeights(a_weights, comp);
  }
  level_redist.redistribute(a_state, interv);
  level_redist.setToZero();
}
