/*!
  @file cdr_gdnv.cpp
  @brief Implementation of cdr_gdnv.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "cdr_gdnv.H"

#include <ExtrapAdvectBC.H>
#include <EBArith.H>

cdr_gdnv::cdr_gdnv() : cdr_tga() {
  m_name = "cdr_gdnv";

  this->set_divF_nc(0);
}


cdr_gdnv::~cdr_gdnv(){

}

int cdr_gdnv::query_ghost() const {
  return 3;
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
    this->allocate_covered();
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
  m_amr->allocate(velo, m_phase, ncomp, 2);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl         = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl              = m_amr->get_ebisl(m_phase)[lvl];
    EBAdvectLevelIntegrator& leveladvect = *(m_amr->get_level_advect(m_phase)[lvl]);
    
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

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const RefCountedPtr<EBAdvectLevelIntegrator>& leveladvect = m_amr->get_level_advect(m_phase)[lvl];
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
      coarsrc_old   = m_source[lvl-1];
      coarsrc_new   = m_source[lvl-1];
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
				    m_source[lvl], 
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
				    m_source[lvl],
				    coarsrc_old,
				    coarsrc_new);
    }
    else{
      MayDay::Abort("cdr_gdnv::extrapolate_to_faces - unknown method requested");
    }


#if 0    // Do not fall below 0
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      for (int dir = 0; dir < SpaceDim; dir++){
	BaseIVFAB<Real>& phi_lo = *(*m_covered_phi_lo[lvl])[dit()][dir];
	BaseIVFAB<Real>& phi_hi = *(*m_covered_phi_hi[lvl])[dit()][dir];


	for(VoFIterator vofit(phi_lo.getIVS(), ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	  phi_lo(vof, comp) = Max(phi_lo(vof, comp), 0.0);
	}
	for(VoFIterator vofit(phi_hi.getIVS(), ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	  phi_hi(vof, comp) = Max(phi_hi(vof, comp), 0.0);
	}


	const Box box = dbl.get(dit());
	const IntVectSet ivs(box);
	FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;

	EBFaceFAB& phi = (*a_face_state[lvl])[dit()][dir];
	for (FaceIterator faceit(ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();
	  if(phi(face, comp) < 0.){
	    phi(face, comp) = 0.;
	  }
	}
      
      }
    }
#endif
  }
}

void cdr_gdnv::nonconservative_divergence(EBAMRIVData& a_div_nc, const EBAMRCellData& a_divF, const EBAMRFluxData& a_face_state){
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
      EBAdvectLevelIntegrator& leveladvect = *(m_amr->get_level_advect(m_phase)[lvl]);

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
