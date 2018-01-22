/*!
  @file cdr_solver.cpp
  @brief Implementation of cdr_solver.H
  @author Robert Marskar
  @date Nov. 2017
  @todo See if we can make EBAdvectLevelIntegrator a static in cdr_gdnv
*/

#include "cdr_solver.H"
#include "cdr_solverF_F.H"
#include "data_ops.H"

#include <EBArith.H>
#include <EBAMRIO.H>

#if 1 // Move to cdr_gdnv
#include <ExtrapAdvectBC.H>
#endif

cdr_solver::cdr_solver(){

  this->set_phase(phase::gas);
  this->set_time(0, 0., 0.);
  
  m_name = "cdr_solver";
}

cdr_solver::~cdr_solver(){

}

int cdr_solver::query_ghost() const {
  CH_TIME("cdr_solver::query_ghost");
  if(m_verbosity > 5){
    pout() << m_name + "::query_ghost" << endl;
  }

  return 3;
}

void cdr_solver::sanity_check(){
  CH_TIME("cdr_solver::sanity_check");
  if(m_verbosity > 5){
    pout() << m_name + "::sanity_check" << endl;
  }

  CH_assert(!m_compgeom.isNull());
  CH_assert(!m_physdom.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_species.isNull());
  CH_assert(!m_ebis.isNull());
}

void cdr_solver::set_species(const RefCountedPtr<species> a_species){
  CH_TIME("cdr_solver::set_species");
  if(m_verbosity > 5){
    pout() << m_name + "::set_species" << endl;
  }

  m_species   = a_species;
  m_name      = m_species->get_name();
  m_diffusive = m_species->is_diffusive();
}

void cdr_solver::set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom){
  CH_TIME("cdr_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << m_name + "::set_computational_geometry" << endl;
  }
  m_compgeom = a_compgeom;

  const RefCountedPtr<mfis> mfis = m_compgeom->get_mfis();
  
  this->set_ebis(mfis->get_ebis(m_phase));
}

void cdr_solver::set_ebis(const RefCountedPtr<EBIndexSpace>& a_ebis){
  CH_TIME("cdr_solver::set_ebis");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebis" << endl;
  }

  m_ebis = a_ebis;
}

void cdr_solver::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("cdr_solver::set_physical_domain");
  if(m_verbosity > 5){
    pout() << m_name + "::set_physical_domain" << endl;
  }

  m_physdom = a_physdom;
}

void cdr_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("cdr_solver::set_amr");
  if(m_verbosity > 5){
    pout() << m_name + "::set_amr" << endl;
  }

  m_amr = a_amr;
}

void cdr_solver::set_time(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("cdr_solver::set_time");
  if(m_verbosity > 5){
    pout() << m_name + "::set_time" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void cdr_solver::set_velocity(const EBAMRCellData& a_velo){
  CH_TIME("cdr_solver::set_velocity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_velocity" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_velo[lvl]->copyTo(*m_velo_cell[lvl]);
  }

  m_amr->average_down(m_velo_cell, m_phase);
  m_amr->interp_ghost(m_velo_cell, m_phase);
}

void cdr_solver::set_velocity(const RealVect a_velo){
  CH_TIME("cdr_solver::set_velocity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_velocity" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (int dir = 0; dir < SpaceDim; dir++){
      data_ops::set_value(*m_velo_cell[lvl], a_velo[dir], dir);
    }
  }

  m_amr->average_down(m_velo_cell, m_phase);
  m_amr->interp_ghost(m_velo_cell, m_phase);
}

void cdr_solver::set_diffco(const EBAMRFluxData& a_diffco, const EBAMRIVData& a_diffco_eb){
  CH_TIME("cdr_solver::set_diffco");
  if(m_verbosity > 5){
    pout() << m_name + "::set_diffco" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_diffco[lvl]->copyTo(*m_diffco[lvl]);
    a_diffco_eb[lvl]->copyTo(*m_diffco_eb[lvl]);
  }

  m_amr->average_down(m_diffco, m_phase);
  m_amr->average_down(m_diffco_eb, m_phase);
}

void cdr_solver::set_diffco(const Real a_diffco){
  CH_TIME("cdr_solver::set_diffco");
  if(m_verbosity > 5){
    pout() << m_name + "::set_diffco" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(*m_diffco[lvl],    a_diffco);
    data_ops::set_value(*m_diffco_eb[lvl], a_diffco);
  }

  m_amr->average_down(m_diffco, m_phase);
  m_amr->average_down(m_diffco_eb, m_phase);
}

void cdr_solver::set_source(const EBAMRCellData& a_source){
  CH_TIME("cdr_solver::set_source");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_source[lvl]->copyTo(*m_source[lvl]);
  }

  m_amr->average_down(m_source, m_phase);
  m_amr->interp_ghost(m_source, m_phase);
}

void cdr_solver::set_source(const Real a_source){
  CH_TIME("cdr_solver::set_source");
  if(m_verbosity > 5){
    pout() << m_name + "::set_source" << endl;
  }

  const int comp = 0;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(*m_source[lvl], a_source, comp);
  }

  m_amr->average_down(m_source, m_phase);
  m_amr->interp_ghost(m_source, m_phase);
}

void cdr_solver::set_ebflux(const EBAMRIVData& a_ebflux){
  CH_TIME("cdr_solver::set_ebflux(variable)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebflux(variable)" << endl;
  }

  const int comp         = 0;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_ebflux[lvl]->copyTo(interv, *m_ebflux[lvl], interv);
  }
}

void cdr_solver::set_ebflux(const Real a_ebflux){
  CH_TIME("cdr_solver::set_ebflux(constant)");
  if(m_verbosity > 5){
    pout() << m_name + "::set_ebflux(constant)" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::set_value(*m_ebflux[lvl], a_ebflux);
  }
}

void cdr_solver::set_phase(const phase::which_phase a_phase){
  CH_TIME("cdr_solver::set_phase");
  if(m_verbosity > 5){
    pout() << m_name + "::set_phase" << endl;
  }

  m_phase = a_phase;
}

void cdr_solver::set_verbosity(const int a_verbosity){
  CH_TIME("cdr_solver::set_verbosity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_verbosity" << endl;
  }
  
  m_verbosity = a_verbosity;
}

void cdr_solver::initial_data(){
  CH_TIME("cdr_solver::initial_data");
  if(m_verbosity > 5){
    pout() << m_name + "::initial_data" << endl;
  }

  const RealVect origin  = m_physdom->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (DataIterator dit = m_state[lvl]->dataIterator(); dit.ok(); ++dit){
      EBCellFAB& state       = (*m_state[lvl])[dit()];
      const Box box          = m_state[lvl]->disjointBoxLayout().get(dit());
      const EBISBox& ebisbox = state.getEBISBox();
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = EBArith::getVofLocation(vof, m_amr->get_dx()[lvl]*RealVect::Unit, origin);
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = m_species->initial_data(pos, m_time);
	}
      }
    }
  }

  m_amr->average_down(m_state, m_phase);
}

void cdr_solver::allocate_internals(){
  CH_TIME("cdr_solver::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  const int sca = 1;
  const int vec = SpaceDim;

  m_amr->allocate(m_state,      m_phase, sca);
  m_amr->allocate(m_source,     m_phase, sca);
  m_amr->allocate(m_velo_face,  m_phase, sca);
  m_amr->allocate(m_velo_cell,  m_phase, vec); 
  m_amr->allocate(m_ebflux,     m_phase, sca);
  m_amr->allocate(m_diffco,     m_phase, sca);
  m_amr->allocate(m_diffco_eb,  m_phase, sca);

  data_ops::set_value(m_state,      0.0);
  data_ops::set_value(m_source,     0.0);
  data_ops::set_value(m_velo_face,  0.0);
  data_ops::set_value(m_velo_cell,  0.0);
  data_ops::set_value(m_ebflux,     0.0);
  data_ops::set_value(m_diffco,     0.0);
  data_ops::set_value(m_diffco_eb,  0.0);

  this->define_interp_stencils();
  this->define_divFnc_stencils();
}

void cdr_solver::define_interp_stencils(){
  CH_TIME("cdr_solver::define_interp_stencils");
  if(m_verbosity > 5){
    pout() << m_name + "::define_interp_stencils" << endl;
  }

  const int comp            = 0;
  const int ncomp           = 1;
  const int finest_level    = m_amr->get_finest_level();
  FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;
  
  for (int dir = 0; dir < SpaceDim; dir++){
    (m_interp_stencils[dir]).resize(1 + finest_level);

    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const ProblemDomain& domain  = m_amr->get_domains()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

      m_interp_stencils[dir][lvl] = RefCountedPtr<LayoutData<BaseIFFAB<FaceStencil> > >
	(new LayoutData<BaseIFFAB<FaceStencil> >(dbl));

      LayoutData<IntVectSet> cfivs(dbl);
      EBArith::defineCFIVS(cfivs, dbl, domain);

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	BaseIFFAB<FaceStencil>& sten = (*m_interp_stencils[dir][lvl])[dit()];
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

void cdr_solver::define_divFnc_stencils(){
  CH_TIME("cdr_solver::define_divFnc_stencils");
  if(m_verbosity > 5){
    pout() << m_name + "::define_divFnc_stencils" << endl;
  }

  const int rad   = 1;
  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  m_stencils_nc.resize(1 + finest_level);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

    m_stencils_nc[lvl] = RefCountedPtr<LayoutData<BaseIVFAB<VoFStencil> > > (new LayoutData<BaseIVFAB<VoFStencil> >(dbl));
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = ebisbox.getIrregIVS(box);

      BaseIVFAB<VoFStencil>& stens = (*m_stencils_nc[lvl])[dit()];
      stens.define(ivs, ebgraph, ncomp);
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	VoFStencil& sten = stens(vof, comp);
	sten.clear();
	
	Real norm = 0.;
	Vector<VolIndex> vofs;
	EBArith::getAllVoFsInMonotonePath(vofs, vof, ebisbox, rad);
	for (int i = 0; i < vofs.size(); i++){
	  const VolIndex& ivof = vofs[i];
	  const Real iweight   = ebisbox.volFrac(ivof);

	  norm += iweight;
	  sten.add(ivof, 1.0);
	}

	sten *= 1./norm;
      }
    }
  }
}

void cdr_solver::advance(const Real& a_dt){
  CH_TIME("cdr_solver::advance(dt)");
  if(m_verbosity > 1){
    pout() << m_name + "::advance(dt)" << endl;
  }
  
  this->advance(m_state, a_dt);
}

void cdr_solver::advance(EBAMRCellData& a_state, const Real& a_dt){
  CH_TIME("cdr_solver::state(dt)");
  if(m_verbosity > 1){
    pout() << m_name + "::state(dt)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();


  // Compute right-hand-side
  EBAMRCellData k1, k2, phi;
  m_amr->allocate(k1,  m_phase, ncomp);
  m_amr->allocate(k2,  m_phase, ncomp);
  m_amr->allocate(phi, m_phase, ncomp);

#if 1
  // Compute f(yn, tn);
  this->compute_rhs(k1, a_state, a_dt);

  // Make phi(n+1) = phi(n) + dt*f(n)
  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_state[lvl]->copyTo(*phi[lvl]);
    data_ops::incr(*phi[lvl], *k1[lvl], a_dt);
  }

  // Compute f(n+1)
  this->compute_rhs(k2, phi, a_dt);

  // Make phi = phi(n) + 0.5*(f(n) + f(n+1))
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::incr(*a_state[lvl], *k1[lvl], 0.5*a_dt);
    data_ops::incr(*a_state[lvl], *k2[lvl], 0.5*a_dt);
  }
#else
  MayDay::Abort("cdr_solver::advance - review exponenents and stuff like that");
  EBAMRCellData phi_new, phi_old, rhs;
  m_amr->allocate(rhs    ,  m_phase, ncomp);
  m_amr->allocate(phi_new,  m_phase, ncomp);
  m_amr->allocate(phi_old,  m_phase, ncomp);

  this->compute_rhs(rhs, a_state, a_dt);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_state[lvl]->copyTo(*phi_new[lvl]);
    a_state[lvl]->copyTo(*phi_old[lvl]);

    data_ops::ln(*phi_old[lvl]);
    data_ops::ln(*phi_new[lvl]);

    data_ops::multiply(*rhs[lvl], *phi_old[lvl]);
    data_ops::incr(*phi_new[lvl], *rhs[lvl], 1.0*a_dt);
    data_ops::exponentiate(*phi_new[lvl], 1.0);

    phi_new[lvl]->copyTo(*a_state[lvl]);
    
  }
#endif

  m_amr->average_down(a_state, m_phase);
  m_amr->interp_ghost(a_state, m_phase);

  m_time += a_dt;
  m_step++;
}

void cdr_solver::compute_rhs(EBAMRCellData& a_rhs, const Real& a_dt){
  CH_TIME("cdr_solver::compute_rhs(rhs, dt)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_rhs(rhs, dt)" << endl;
  }
  
  this->compute_rhs(a_rhs, m_state, a_dt);
}

void cdr_solver::compute_rhs(EBAMRCellData& a_rhs, const EBAMRCellData& a_state, const Real& a_dt){
  CH_TIME("cdr_solver::compute_rhs(rhs, state, dt)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_rhs(rhs, state, dt)" << endl;
  }

  data_ops::set_value(a_rhs, 0.0);  

  const int comp  = 0;
  const int ncomp = 1;

  // Advective derivative
  EBAMRCellData advective_term;
  EBAMRCellData diffusion_term;

  m_amr->allocate(advective_term, m_phase, ncomp);
  if(this->is_diffusive()){
    m_amr->allocate(diffusion_term, m_phase, ncomp);
  }

  // Advective term
  this->compute_divF(advective_term, a_state);
  data_ops::incr(a_rhs, advective_term, -1.0);

  // Diffusion term
  if(this->is_diffusive()){
    this->compute_diffusion_term(diffusion_term, a_state);
    data_ops::incr(a_rhs, diffusion_term, 1.0);
  }

  // Source term
  data_ops::incr(a_rhs, m_source, 1.0);

  m_amr->average_down(a_rhs, m_phase);
  m_amr->interp_ghost(a_rhs, m_phase);
}

void cdr_solver::compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state){
  CH_TIME("cdr_solver::compute_divF(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divF(divF, state)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  EBAMRFluxData face_state;
  EBAMRIVData   div_nc;
  EBAMRIVData   mass_diff;

  m_amr->allocate(face_state, m_phase, ncomp);
  m_amr->allocate(div_nc,     m_phase, ncomp);
  m_amr->allocate(mass_diff,  m_phase, ncomp);

  data_ops::set_value(a_divF,     0.0); 
  data_ops::set_value(face_state, 0.0);
  data_ops::set_value(div_nc,     0.0);
  data_ops::set_value(mass_diff,  0.0);

  // Compute the advective derivative
  this->average_velo_to_faces(m_velo_face, m_velo_cell);            // Average cell-centered velocities to face centers
  this->extrapolate_to_faces(face_state, a_state);                  // Face extrapolation to cell-centered faces
  this->conservative_divergence(a_divF, face_state, m_velo_face);   // a_divF holds the conservative divergence
  this->nonconservative_divergence(div_nc, a_divF);                 // Compute non-conservative divergence
  this->hybrid_divergence(a_divF, mass_diff, div_nc);               // Make divF = hybrid divergence. Compute mass diff.
  this->increment_flux_register(face_state, m_velo_face);           // Increment flux registers
  this->increment_redist(mass_diff);                                // Increment redistribution objects

  const bool ebcf = m_amr->get_ebcf();
  if(ebcf){ 
    this->coarse_fine_increment(mass_diff);     // Increment the coarse-fine redistribution objects
    this->increment_redist_flux();              // Increment flux registers with the redistribution stuff
    this->coarse_fine_redistribution(a_divF);   // Redistribute
  }
  else{
    this->hyperbolic_redistribution(a_divF, mass_diff, a_state);  // Redistribute mass into hybrid divergence
    //    this->reflux(a_divF);                                         // Reflux at coarse-fine interfaces
  }
}

void cdr_solver::compute_diffusion_term(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state){
  CH_TIME("cdr_solver::compute_diffusion_term");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_term" << endl;
  }

  MayDay::Abort("cdr_solver::compute_diffusion_term - not implemented");
}

void cdr_solver::average_velo_to_faces(EBAMRFluxData& a_velo_face, const EBAMRCellData& a_velo_cell){
  CH_TIME("cdr_solver::average_velo_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::average_velo_to_faces" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::average_cell_to_face(*a_velo_face[lvl], *a_velo_cell[lvl], m_amr->get_domains()[lvl]);
  }
}

#if 1 // This should be moved to cdr_gdnv
void cdr_solver::extrapolate_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state){
  CH_TIME("cdr_solver::extrapolate_to_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::extrapolate_to_faces" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  const int comp = 0;
  
  EBAdvectPatchIntegrator::setCurComp(0);
  EBAdvectPatchIntegrator::setDoingVel(0);

  RefCountedPtr<ExtrapAdvectBCFactory> bcfact = RefCountedPtr<ExtrapAdvectBCFactory>
    (new ExtrapAdvectBCFactory());


  for (int lvl = 0; lvl <= finest_level; lvl++){
    const RefCountedPtr<EBAdvectLevelIntegrator>& leveladvect = m_amr->get_level_advect(m_phase)[lvl];
    leveladvect->resetBCs(bcfact);
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
				  m_dt,
				  m_source[lvl], 
				  coarsrc_old,
				  coarsrc_new);
  }
}
#endif

void cdr_solver::conservative_divergence(EBAMRCellData&       a_cons_div,
					 const EBAMRFluxData& a_face_vel,
					 const EBAMRFluxData& a_face_state){
  CH_TIME("cdr_solver::conservative_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::conservative_divergence" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];

    // Compute div(F) on regular cells
    this->advective_derivative(*a_cons_div[lvl], *a_face_vel[lvl], *a_face_state[lvl], lvl); 

    // Interpolate fluxes to face centroids on irregular cells
    LevelData<BaseIFFAB<Real> > flux[SpaceDim];
    this->compute_flux_interpolant(flux, *a_face_state[lvl], *a_face_vel[lvl], lvl);
    this->interpolate_flux_to_centroids(flux, lvl);
    this->compute_divF_irreg(*a_cons_div[lvl], flux, *m_ebflux[lvl], lvl);
  }
}

void cdr_solver::advective_derivative(LevelData<EBCellFAB>&       a_divF,
				      const LevelData<EBFluxFAB>& a_vel,
				      const LevelData<EBFluxFAB>& a_phi,
				      const int                   a_lvl){

  CH_TIME("cdr_solver::advective_derivative");
  if(m_verbosity > 5){
    pout() << m_name + "::advective_derivative" << endl;
  }

  CH_assert(a_divF.nComp() == 1);
  CH_assert(a_vel.nComp()  == 1);
  CH_assert(a_phi.nComp()  == 1);

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain domain  = m_amr->get_domains()[a_lvl];
  const Real dx               = m_amr->get_dx()[a_lvl];

  for (DataIterator dit = a_divF.dataIterator(); dit.ok(); ++dit){
    EBCellFAB& divF         = a_divF[dit()];
    BaseFab<Real>& divF_fab = divF.getSingleValuedFAB();
    const Box box = dbl.get(dit());

    for (int dir = 0; dir < SpaceDim; dir++){
      const EBFaceFAB& vel = a_vel[dit()][dir];
      const EBFaceFAB& phi = a_phi[dit()][dir];

      const BaseFab<Real>& vel_fab = vel.getSingleValuedFAB();
      const BaseFab<Real>& phi_fab = phi.getSingleValuedFAB();

      FORT_ADVECTIVEDERIV(CHF_FRA1(divF_fab, comp),
			  CHF_CONST_FRA1(vel_fab, comp),
			  CHF_CONST_FRA1(phi_fab, comp),
			  CHF_CONST_INT(dir),
			  CHF_CONST_INT(ncomp),
			  CHF_CONST_REAL(dx),
			  CHF_BOX(box));;
    }


    // Reset irregular cells - these contain bogus values
    const EBISBox& ebisbox = divF.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      divF(vof, comp) = 0.;
    }
  }
}

void cdr_solver::compute_flux_interpolant(LevelData<BaseIFFAB<Real> >   a_interpolant[SpaceDim],
					  const LevelData<EBFluxFAB>&   a_vel,
					  const LevelData<EBFluxFAB>&   a_state,
					  const int                     a_lvl){

  const int ncomp = 1;
  const int comp  = 0;
  
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];

  for (int dir = 0; dir < SpaceDim; dir++){
    // Define interpolant
    LayoutData<IntVectSet> grown_set;
    EBArith::defineFluxInterpolant(a_interpolant[dir],
				   grown_set,
				   dbl,
				   ebisl,
				   domain,
				   ncomp,
				   dir);

    // Compute interpolant
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      BaseIFFAB<Real>& interpol   = a_interpolant[dir][dit()];
      const Box& box              = dbl.get(dit());
      const EBISBox& ebisbox      = ebisl[dit()];
      const EBGraph& ebgraph      = ebisbox.getEBGraph();
      const EBFaceFAB& vel        = a_vel[dit()][dir];
      const EBFaceFAB& phi        = a_state[dit()][dir];
      const IntVectSet ivs        = grown_set[dit()] & box;
      const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingWithBoundary;

      for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();
	interpol(face, comp) = vel(face, comp)*phi(face, comp);
      }
    }

    a_interpolant[dir].exchange();
  }
}

void cdr_solver::interpolate_flux_to_centroids(LevelData<BaseIFFAB<Real> >       a_flux[SpaceDim],
					       const int                         a_lvl){
  CH_TIME("cdr_solver::interpolate_flux_to_centroids");
  if(m_verbosity > 5){
    pout() << m_name + "::interpolate_flux_to_centroids" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);
  
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];

  const FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);
    
    for (int dir = 0; dir < SpaceDim; dir++){
      BaseIFFAB<Real>& flux = a_flux[dir][dit()];
      BaseIFFAB<Real> centroid_flux;

      // Compute centroid flux
      centroid_flux.define(ivs, ebgraph, dir, ncomp);
      for (FaceIterator faceit(ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	const FaceIndex& face  = faceit();
	const FaceStencil& sten = (*m_interp_stencils[dir][a_lvl])[dit()](face, comp);

#if 1 // This can probably be removed
	if(sten.size() == 0){
	  centroid_flux(face, comp) = flux(face, comp);
	}
	else{
#endif
	  centroid_flux(face, comp) = 0.;
	  for (int i = 0; i < sten.size(); i++){
	    const FaceIndex& iface = sten.face(i);
	    const Real iweight     = sten.weight(i);
	  
	    centroid_flux(face, comp) += iweight*flux(iface, comp);
	  }
#if 1
	}
#endif
      }

      // Copy centroid flux into a_flux
      flux.setVal(0.0);
      flux.copy(box, interv, box, centroid_flux, interv);
    }
  }
}

void cdr_solver::compute_divF_irreg(LevelData<EBCellFAB>&              a_divF,
				    const LevelData<BaseIFFAB<Real> >  a_flux[SpaceDim],
				    const LevelData<BaseIVFAB<Real> >& a_ebflux,
				    const int                          a_lvl){
  CH_TIME("cdr_solver::compute_divF_irreg");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divF_irreg" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];
  const Real dx                = m_amr->get_dx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);

    EBCellFAB& divF               = a_divF[dit()];
    const BaseIVFAB<Real>& ebflux = a_ebflux[dit()];

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const Real area     = ebisbox.bndryArea(vof);

      // EB flux
      divF(vof, comp) = ebflux(vof,comp)*area;

      // Face fluxes
      for (int dir = 0; dir < SpaceDim; dir++){
	const BaseIFFAB<Real>& flux = a_flux[dir][dit()];

	for (SideIterator sit; sit.ok(); ++sit){
	  const int isign = sign(sit());
	  const Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());

	  for (int iface = 0; iface < faces.size(); iface++){
	    const FaceIndex face = faces[iface];
	    const Real face_area = ebisbox.areaFrac(face);
	    divF(vof, comp) += isign*face_area*flux(face, comp);
	  }
	}
      }

      // Scale divF by dx but not by kappa.
      divF(vof, comp) *= 1./dx;
    }
  }
}

void cdr_solver::nonconservative_divergence(EBAMRIVData& a_div_nc, const EBAMRCellData& a_divF){
  CH_TIME("cdr_solver::nonconservative_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::nonconservative_divergence" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    
    for (DataIterator dit = a_div_nc[lvl]->dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& div_nc = (*a_div_nc[lvl])[dit()];
      EBCellFAB& divF         = (*a_divF[lvl])[dit()];

      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = ebisbox.getIrregIVS(box);

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof    = vofit();
	const VoFStencil& sten = (*m_stencils_nc[lvl])[dit()](vof, comp);

	div_nc(vof, comp) = 0.;
	for (int i = 0; i < sten.size(); i++){
	  const VolIndex& ivof = sten.vof(i);
	  const Real& iweight  = sten.weight(i);
	  div_nc(vof,comp) += iweight*divF(ivof, comp);
	}
      }
    }
  }
}

void cdr_solver::hybrid_divergence(EBAMRCellData&     a_hybrid_div,
				   EBAMRIVData&       a_mass_diff,
				   const EBAMRIVData& a_noncons_div){
  CH_TIME("cdr_solver::hybrid_divergence");
  if(m_verbosity > 5){
    pout() << m_name + "::hybrid_divergence" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBCellFAB& divH               = (*a_hybrid_div[lvl])[dit()];  // On input, this contains kappa*div(F)
      BaseIVFAB<Real>& deltaM       = (*a_mass_diff[lvl])[dit()];
      const BaseIVFAB<Real>& divNC  = (*a_noncons_div[lvl])[dit()]; 

      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = ebisbox.getIrregIVS(box);

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const Real kappa    = ebisbox.volFrac(vof);
	const Real dc       = divH(vof, comp);
	const Real dnc      = divNC(vof, comp);

	divH(vof, comp)   = (dc + (1-kappa)*dnc);          // On output, contains hybrid divergence
	deltaM(vof, comp) = (1-kappa)*(dc - kappa*dnc);
      }
    }
  }
}

void cdr_solver::increment_flux_register(const EBAMRFluxData& a_face_state, const EBAMRFluxData& a_velo_face){
  CH_TIME("cdr_solver::increment_flux_register");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_flux_register" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFastFR> >& fluxreg = m_amr->get_flux_reg(m_phase);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_fine){
      fluxreg[lvl]->setToZero();
    }

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      const Box box          = dbl.get(dit());
      
      for (int dir = 0; dir < SpaceDim; dir++){
	const Real scale     = 1.;
	const EBFaceFAB& phi = (*a_face_state[lvl])[dit()][dir];
	const EBFaceFAB& vel = (*a_velo_face[lvl])[dit()][dir];

	// Compute flux
	EBFaceFAB flux(ebisbox, box, dir, ncomp);
	flux.setVal(0.0);
	flux += phi;
	flux *= vel;

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

void cdr_solver::increment_redist(const EBAMRIVData& a_mass_diff){
  CH_TIME("cdr_solver::increment_redist");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_redist" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    
    EBLevelRedist& level_redist = *(m_amr->get_level_redist(m_phase)[lvl]);
    level_redist.setToZero();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      level_redist.increment((*a_mass_diff[lvl])[dit()], dit(), interv);
    }
  }
}

void cdr_solver::hyperbolic_redistribution(EBAMRCellData&       a_divF,
					   const EBAMRIVData&   a_mass_diff,
					   const EBAMRCellData& a_redist_weights){
  CH_TIME("cdr_solver::hyberbolic_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::hyperbolic_redistribution" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);


  for (int lvl = 0; lvl <= finest_level; lvl++){
    EBLevelRedist& level_redist = *(m_amr->get_level_redist(m_phase)[lvl]);
    level_redist.redistribute(*a_divF[lvl], interv);
    level_redist.setToZero();
  }
}

void cdr_solver::reflux(EBAMRCellData& a_state){
  CH_TIME("cdr_solver::reflux");
  if(m_verbosity > 5){
    pout() << m_name + "::reflux" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFastFR > > fluxreg = m_amr->get_flux_reg(m_phase);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Real dx       = m_amr->get_dx()[lvl];
    const bool has_coar = lvl > 0;
    const bool has_fine = lvl < finest_level;

    if(has_fine){
      const Real scale = -1.0/dx;
      fluxreg[lvl]->reflux(*a_state[lvl], interv, scale);
      fluxreg[lvl]->setToZero();
    }
  }
}

void cdr_solver::increment_redist_flux(){
  CH_TIME("cdr_solver::increment_redist_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::increment_redist_flux" << endl;
  }

  MayDay::Abort("cdr_solver::increment_redist_flux - not implemented");
}

void cdr_solver::coarse_fine_increment(const EBAMRIVData& m_mass_diff){
  CH_TIME("cdr_solver::coarse_fine_increment");
  if(m_verbosity > 5){
    pout() << m_name + "::coarse_fine_increment" << endl;
  }

  MayDay::Abort("cdr_solver::coarse_fine_increment - not implemented");
}

void cdr_solver::coarse_fine_redistribution(EBAMRCellData& a_state){
  CH_TIME("cdr_solver::coarse_fine_redistribution");
  if(m_verbosity > 5){
    pout() << m_name + "::coarse_fine_redistribution" << endl;
  }

  MayDay::Abort("cdr_solver::coarse_fine_redistribution - not implemented");
}

#ifdef CH_USE_HDF5
void cdr_solver::write_plot_file(){
  CH_TIME("cdr_solver::write_plot_file");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_file" << endl;
  }

  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_step, SpaceDim);

  const int ncomps = 2 + SpaceDim;
  Vector<string> names(ncomps);
  names[0] = "density";
  names[1] = "source";
  names[2] = "x-velocity";
  names[3] = "y-velocity";
  if(SpaceDim == 3){
    names[4] = "y-velocity";
  }


  EBAMRCellData output;
  m_amr->allocate(output, m_phase, ncomps);

  for (int lvl = 0; lvl < output.size(); lvl++){
    LevelData<EBCellFAB>& state  = *m_state[lvl];
    LevelData<EBCellFAB>& source = *m_source[lvl];
    LevelData<EBCellFAB>& velo   = *m_velo_cell[lvl];

    state.copyTo(Interval(0,0),           *output[lvl],  Interval(0,0));
    source.copyTo(Interval(0,0),          *output[lvl],  Interval(1,1));
    velo.copyTo(Interval(0,SpaceDim - 1), *output[lvl],  Interval(2, 2 + (SpaceDim-1)));
  }

#if 0 // Possibly removed. 
  // Transform to centroids
  irreg_amr_stencil<centroid_interp>& sten = m_amr->get_centroid_interp_stencils(phase::gas);
  sten.apply(output, true);
#endif


  Vector<LevelData<EBCellFAB>* > output_ptr;
  m_amr->alias(output_ptr, output);

#if 0 // Should be removed
  data_ops::floor(output, 0.0);
#endif
  
  Vector<Real> covered_values(ncomps, 0.0);
  string fname(file_char);
  writeEBHDF5(fname,
	      m_amr->get_grids(),
	      output_ptr,
	      names,
	      m_amr->get_domains()[0].domainBox(),
	      m_amr->get_dx()[0],
	      m_dt,
	      m_time,
	      m_amr->get_ref_rat(),
	      m_amr->get_finest_level() + 1,
	      true,
	      covered_values);
}
#endif

Real cdr_solver::compute_dt(){
  CH_TIME("cdr_solver::compute_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_dt" << endl;
  }

  Real min_dt = 1.E99;

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  Vector<RefCountedPtr<EBFastFR > > fluxreg = m_amr->get_flux_reg(m_phase);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& velo  = (*m_velo_cell[lvl])[dit()];
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs   = ebisbox.getIrregIVS(box);

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex vof = vofit();
	const RealVect u  = RealVect(D_DECL(velo(vof, 0), velo(vof, 1), velo(vof, 2)));
	const Real thisdt = dx/u.vectorLength();

	min_dt = Min(min_dt, thisdt);
      }
    }
  }
  
#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&min_dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("EBAMRAdvectDiffuse::computeAdvectiveDt() - communication error on norm");
  }
  min_dt = tmp;
#endif


  return min_dt;
}

Real cdr_solver::compute_mass(){
  CH_TIME("cdr_solver::compute_mass");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_mass" << endl;
  }

  Real mass = 0.;
  const int base = 0;
  const Real dx = m_amr->get_dx()[base];
  m_amr->average_down(m_state, m_phase);
  
  data_ops::kappa_sum(mass, *m_state[base]);
  mass *= pow(dx, SpaceDim);
  
  return mass;
}

bool cdr_solver::is_diffusive(){
  CH_TIME("cdr_solver::is_diffusive");
  if(m_verbosity > 5){
    pout() << m_name + "::is_diffusive" << endl;
  }
  
  return m_diffusive;
}

EBAMRCellData& cdr_solver::get_state(){
  return m_state;
}

EBAMRCellData& cdr_solver::get_source(){
  return m_source;
}

EBAMRFluxData& cdr_solver::get_velo_face(){
  return m_velo_face;
}

EBAMRIVData& cdr_solver::get_velo_eb(){
  return m_velo_eb;
}

EBAMRFluxData& cdr_solver::get_diffco_face(){
  return m_diffco;
}

EBAMRIVData& cdr_solver::get_diffco_eb(){
  return m_diffco_eb;
}
