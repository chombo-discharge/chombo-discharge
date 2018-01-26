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
#include <NeumannConductivityEBBC.H>
#include <NeumannConductivityDomainBC.H>

cdr_solver::cdr_solver(){

  this->set_phase(phase::gas);
  this->set_time(0, 0., 0.);
  this->set_mass_redist(false);
  this->set_gmg_solver_parameters();
  this->set_bottom_solver(1);
  this->set_bottom_drop(2);
  this->set_tga(true);
  
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

void cdr_solver::set_bottom_solver(const int a_whichsolver){
  CH_TIME("cdr_solver::set_bottom_solver");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bottom_solver" << endl;
  }
  if(a_whichsolver == 0 || a_whichsolver == 1){
    m_bottomsolver = a_whichsolver;
  }
  else{
    MayDay::Abort("cdr_solver::set_bottom_solver - Unsupported solver type requested");
  }
}

void cdr_solver::set_botsolver_smooth(const int a_numsmooth){
  CH_TIME("cdr_solver::set_botsolver_smooth");
  if(m_verbosity > 5){
    pout() << m_name + "::set_botsolver_smooth" << endl;
  }
  CH_assert(a_numsmooth > 0);
  
  m_numsmooth = a_numsmooth;
}

void cdr_solver::set_bottom_drop(const int a_bottom_drop){
  CH_TIME("cdr_solver::set_bottom_drop");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bottom_drop" << endl;
  }
  
  m_bottom_drop = a_bottom_drop;
}

void cdr_solver::set_tga(const bool a_use_tga){
  CH_TIME("cdr_solver::set_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::set_tga" << endl;
  }
  
  m_use_tga = a_use_tga;
}

void cdr_solver::set_gmg_solver_parameters(relax::which_relax a_relax_type,
					   amrmg::which_mg    a_gmg_type,      
					   const int          a_verbosity,          
					   const int          a_pre_smooth,         
					   const int          a_post_smooth,       
					   const int          a_bot_smooth,         
					   const int          a_max_iter,
					   const int          a_min_iter,
					   const Real         a_eps,               
					   const Real         a_hang){
  CH_TIME("cdr_solver::set_gmg_solver_parameters");
  if(m_verbosity > 5){
    pout() << m_name + "::set_gmg_solver_parameters" << endl;
  }

  m_gmg_relax_type  = a_relax_type;
  m_gmg_type        = a_gmg_type;
  m_gmg_verbosity   = a_verbosity;
  m_gmg_pre_smooth  = a_pre_smooth;
  m_gmg_post_smooth = a_post_smooth;
  m_gmg_bot_smooth  = a_bot_smooth;
  m_gmg_max_iter    = a_max_iter;
  m_gmg_min_iter    = a_min_iter;
  m_gmg_eps         = a_eps;
  m_gmg_hang        = a_hang;
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
  m_amr->allocate(m_scratch,       m_phase, sca);

  data_ops::set_value(m_state,      0.0);
  data_ops::set_value(m_source,     0.0);
  data_ops::set_value(m_velo_face,  0.0);
  data_ops::set_value(m_velo_cell,  0.0);
  data_ops::set_value(m_ebflux,     0.0);
  data_ops::set_value(m_diffco,     0.0);
  data_ops::set_value(m_diffco_eb,  0.0);
  data_ops::set_value(m_scratch,       0.0);

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

void cdr_solver::setup_gmg(){
  CH_TIME("cdr_solver::setup_gmg");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_gmg" << endl;
  }

  this->setup_operator_factory();
  this->setup_multigrid();

  if(m_use_tga){
    this->setup_tga();
  }
  else{
    this->setup_euler();
  }
}

void cdr_solver::setup_operator_factory(){
  CH_TIME("cdr_solver::setup_operator_factory");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_operator_factory" << endl;
  }

  const int finest_level                 = m_amr->get_finest_level();
  const int ghost                        = m_amr->get_num_ghost();
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids();
  const Vector<int>& refinement_ratios   = m_amr->get_ref_rat();
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real>& dx                 = m_amr->get_dx();
  const RealVect& origin                 = m_physdom->get_prob_lo();
  const Vector<EBISLayout>& ebisl        = m_amr->get_ebisl(m_phase);
  const Vector<RefCountedPtr<EBQuadCFInterp> >& quadcfi  = m_amr->get_old_quadcfi(m_phase);

  Vector<EBLevelGrid> levelgrids;
  for (int lvl = 0; lvl <= finest_level; lvl++){ 
    levelgrids.push_back(*(m_amr->get_eblg(m_phase)[lvl])); // amr_mesh uses RefCounted levelgrids. EBConductivityOp does not. 
  }

  // Appropriate coefficients. 
  const Real alpha =  0.0;
  const Real beta  =  1.0;

  // Default is Neumann. This might change in the future. 
  RefCountedPtr<NeumannConductivityDomainBCFactory> domfact = RefCountedPtr<NeumannConductivityDomainBCFactory>
    (new NeumannConductivityDomainBCFactory());
  RefCountedPtr<NeumannConductivityEBBCFactory> ebfact  = RefCountedPtr<NeumannConductivityEBBCFactory>
    (new NeumannConductivityEBBCFactory());

  domfact->setValue(0.0);
  ebfact->setValue(0.0);
  
  // Create operator factory.
  data_ops::set_value(m_scratch, 1.0);
  m_opfact = RefCountedPtr<ebconductivityopfactory> (new ebconductivityopfactory(levelgrids,
										 quadcfi,
										 alpha,
										 beta,
										 m_scratch,
										 m_diffco,
										 m_diffco_eb,
										 dx[0],
										 refinement_ratios,
										 domfact,
										 ebfact,
										 ghost*IntVect::Unit,
										 ghost*IntVect::Unit,
										 m_gmg_relax_type,
										 m_bottom_drop));
}

void cdr_solver::setup_multigrid(){
  CH_TIME("cdr_solver::setup_multigrid");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_multigrid" << endl;
  }

  const int finest_level       = m_amr->get_finest_level();
  const ProblemDomain coar_dom = m_amr->get_domains()[0];

  // Select bottom solver
  LinearSolver<LevelData<EBCellFAB> >* botsolver = NULL;
  if(m_bottomsolver == 0){
    m_simple_solver.setNumSmooths(m_numsmooth);
    botsolver = &m_simple_solver;
  }
  else{
    botsolver = &m_bicgstab;
  }

  m_gmg_solver = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >());
  m_gmg_solver->m_imin = m_gmg_min_iter;
  m_gmg_solver->m_verbosity = m_gmg_verbosity;
  m_gmg_solver->define(coar_dom, *m_opfact, botsolver, 1 + finest_level);
  m_gmg_solver->setSolverParameters(m_gmg_pre_smooth,
				    m_gmg_post_smooth,
				    m_gmg_bot_smooth,
				    m_gmg_type,
				    m_gmg_max_iter,
				    m_gmg_eps,
				    m_gmg_hang,
				    1.E-99); // Residue set through other means

}

void cdr_solver::setup_tga(){
  CH_TIME("cdr_solver::setup_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_tga" << endl;
  }
  
  const int finest_level       = m_amr->get_finest_level();
  const ProblemDomain coar_dom = m_amr->get_domains()[0];
  const Vector<int> ref_rat    = m_amr->get_ref_rat();

  m_tgasolver = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
    (new AMRTGA<LevelData<EBCellFAB> > (m_gmg_solver, *m_opfact, coar_dom, ref_rat, 1 + finest_level, m_gmg_solver->m_verbosity));

  // Must init gmg for TGA
  Vector<LevelData<EBCellFAB>* > phi, rhs;
  m_amr->alias(phi, m_state);
  m_amr->alias(rhs, m_source);
  m_gmg_solver->init(phi, rhs, finest_level, 0);
}

void cdr_solver::setup_euler(){
  CH_TIME("cdr_solver::setup_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_euler" << endl;
  }

  const int finest_level       = m_amr->get_finest_level();
  const ProblemDomain coar_dom = m_amr->get_domains()[0];
  const Vector<int> ref_rat    = m_amr->get_ref_rat();

  m_eulersolver = RefCountedPtr<EBBackwardEuler> 
    (new EBBackwardEuler (m_gmg_solver, *m_opfact, coar_dom, ref_rat, 1 + finest_level, m_gmg_solver->m_verbosity));

  // Note: If this crashes, try to init gmg first
}

void cdr_solver::advance(const Real& a_dt){
  CH_TIME("cdr_solver::advance(dt)");
  if(m_verbosity > 1){
    pout() << m_name + "::advance(dt)" << endl;
  }
  
  this->advance(m_state, a_dt);
}

void cdr_solver::advance(EBAMRCellData& a_state, const Real& a_dt){
  CH_TIME("cdr_solver::advance(state, dt)");
  if(m_verbosity > 1){
    pout() << m_name + "::advance(state, dt)" << endl;
  }

  if(this->is_diffusive()){
    bool m_needs_setup = true;
    if(m_needs_setup){
      this->setup_gmg();
    }
    if(m_use_tga){
      this->advance_tga(a_state, a_dt);
    }
    else{
      MayDay::Abort("cdr_solver::advance - advance_euler is not implemented (yet)");
    }
  }
  else{
    const Real midpoint = 0.5;
    this->advance_rk2(a_state, a_dt, midpoint);
  }

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::floor(*a_state[lvl], 0.0); // This removes mass!
  }
  
  m_amr->average_down(a_state, m_phase);
  m_amr->interp_ghost(a_state, m_phase);

  m_time += a_dt;
  m_step++;
}

void cdr_solver::advance_rk2(EBAMRCellData& a_state, const Real a_dt, const Real a_alpha){
  CH_TIME("cdr_solver::advance_rk2");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_rk2" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  // Compute right-hand-side
  EBAMRCellData k1, k2, phi;
  m_amr->allocate(k1,  m_phase, ncomp);
  m_amr->allocate(k2,  m_phase, ncomp);
  m_amr->allocate(phi, m_phase, ncomp);

  // Compute f(yn, tn);
  this->compute_rhs(k1, a_state, a_dt);

  //  phi = yn + a_alpha*dt*f(tn, yn)
  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_state[lvl]->copyTo(*phi[lvl]);
    data_ops::incr(*phi[lvl], *k1[lvl], a_alpha*a_dt);
  }

  // Compute f(tn + alpha*dt, yn + alpha*dt*f(yn, tn)) = f(tn + alpha*dt, phi)
  this->compute_rhs(k2, phi, a_dt);

  // Make y(n+1) = y(n) 
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::incr(*a_state[lvl], *k1[lvl], (1. - 1./(2*a_alpha))*a_dt);
    data_ops::incr(*a_state[lvl], *k2[lvl], 1./(2*a_alpha)*a_dt);
  }
}

void cdr_solver::advance_tga(EBAMRCellData& a_state, const Real a_dt){
  CH_TIME("cdr_solver::advance_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga" << endl;
  }

  bool converged = false;

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  // Create a source term = S - div(n*v)
  EBAMRCellData src, phi, divF;
  m_amr->allocate(src,  m_phase, ncomp);
  m_amr->allocate(phi,  m_phase, ncomp);
  m_amr->allocate(divF, m_phase, ncomp);

  data_ops::set_value(src,  0.0);
  data_ops::set_value(phi,  0.0);
  data_ops::set_value(divF, 0.0);

  this->compute_advection_term(divF, a_state, 0.5*a_dt); // extrapolate div(n*v) to 0.5*a_dt and use it as a source term for TGA

  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::incr(*phi[lvl], *a_state[lvl],  1.0);
    data_ops::incr(*src[lvl], *m_source[lvl], 1.0);
    data_ops::incr(*src[lvl], *divF[lvl],    -1.0);
  }

  // Do the aliasing stuff
  Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
  m_amr->alias(new_state, phi);
  m_amr->alias(old_state, a_state);
  m_amr->alias(source,    src);
  

  // Advance
  const Real alpha = 0.0;
  const Real beta  = 1.0;
  
  if(m_use_tga){
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, 0.0);
  }
  else{
    m_eulersolver->resetAlphaAndBeta(alpha, beta);
    m_eulersolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level);
  }

  const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
  if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
    converged = true;
  }

  data_ops::copy(a_state, phi);
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
  
  this->compute_divJ(a_rhs, a_state, a_dt); // a_rhs = div(J)
  data_ops::incr(a_rhs, m_source, 1.0);     // a_rhs = div(J) + S

  m_amr->average_down(a_rhs, m_phase);
  m_amr->interp_ghost(a_rhs, m_phase);
}

void cdr_solver::compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt){
  CH_TIME("cdr_solver::compute_divJ(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divJ(divF, state)" << endl;
  }


  const int comp  = 0;
  const int ncomp = 1;

  data_ops::set_value(a_divJ, 0.0);

  EBAMRCellData advective_term;
  EBAMRCellData diffusion_term;
  m_amr->allocate(advective_term, m_phase, ncomp);
  if(this->is_diffusive()){
    m_amr->allocate(diffusion_term, m_phase, ncomp);
  }

  // Compute advective term
  this->compute_advection_term(advective_term, a_state, 0.0);
  data_ops::incr(a_divJ, advective_term, -1.0);

  // Add in diffusion term
  if(this->is_diffusive()){
    this->compute_diffusion_term(diffusion_term, a_state);
    data_ops::incr(a_divJ, diffusion_term, 1.0);
  }
}

void cdr_solver::compute_advection_term(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt){
  CH_TIME("cdr_solver::compute_advection_term(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_advection_term(divF, state)" << endl;
  }

  const int comp       = 0;
  const int ncomp      = 1;
  const int redist_rad = m_amr->get_redist_rad();

  EBAMRFluxData face_state;
  EBAMRIVData   div_nc;
  EBAMRIVData   mass_diff;
  EBAMRCellData weights;

  m_amr->allocate(face_state, m_phase, ncomp);
  m_amr->allocate(div_nc,     m_phase, ncomp);
  m_amr->allocate(mass_diff,  m_phase, ncomp);
  m_amr->allocate(weights,    m_phase, ncomp, 2*redist_rad);

  data_ops::set_value(a_divF,     0.0); 
  data_ops::set_value(face_state, 0.0);
  data_ops::set_value(div_nc,     0.0);
  data_ops::set_value(mass_diff,  0.0);
  data_ops::set_value(weights,    0.0);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    a_state[lvl]->copyTo(*weights[lvl]);
  }

  // Compute the advective derivative
  this->average_velo_to_faces(m_velo_face, m_velo_cell);            // Average cell-centered velocities to face centers
  this->extrapolate_to_faces(face_state, a_state, a_extrap_dt);     // Face extrapolation to cell-centered faces
  this->conservative_divergence(a_divF, face_state, m_velo_face);   // a_divF holds the conservative divergence
  this->nonconservative_divergence(div_nc, a_divF, face_state);     // Compute non-conservative divergence
  this->hybrid_divergence(a_divF, mass_diff, div_nc);               // Make divF = hybrid divergence. Compute mass diff.
  this->increment_flux_register(face_state, m_velo_face);           // Increment flux registers
  this->increment_redist(mass_diff);                                // Increment redistribution objects

  // Mass weights.
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    data_ops::incr(*weights[lvl], *a_state[lvl], 1.0);
  }

  const bool ebcf = m_amr->get_ebcf();
  if(ebcf){ 
    this->coarse_fine_increment(mass_diff);     // Increment the coarse-fine redistribution objects
    this->increment_redist_flux();              // Increment flux registers with the redistribution stuff
    this->coarse_fine_redistribution(a_divF);   // Redistribute
  }
  else{
    this->hyperbolic_redistribution(a_divF, mass_diff, weights);  // Redistribute mass into hybrid divergence
    this->reflux(a_divF);                                         // Reflux at coarse-fine interfaces
  }
}

void cdr_solver::compute_diffusion_term(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state){
  CH_TIME("cdr_solver::compute_diffusion_term");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_diffusion_term" << endl;
  }

  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  // Copy a_state because AMRMultiGrid may do exchange operations
  EBAMRCellData clone;
  m_amr->allocate(clone, m_phase, ncomp);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    a_state[lvl]->copyTo(*clone[lvl]);
  }

  // Multigrid doesn't want smart pointers, so do aliasing. 
  Vector<LevelData<EBCellFAB>* > res, phi, zero;
  m_amr->alias(res,  a_diffusive_term);
  m_amr->alias(phi,  clone);
  m_amr->alias(zero, m_scratch);

  // TGA can mess with alpha and beta so I need to reset them to the appropriate values
  const Real alpha =  0.0;
  const Real beta  = -1.0; // Minus one because the AMRMultiGrid computes the residual as resid = rhs - L(phi)
  
  if(m_use_tga){
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
  }
  else{
    m_eulersolver->resetAlphaAndBeta(alpha, beta);
  }
  
  data_ops::set_value(m_scratch, 0.0);
  m_gmg_solver->computeAMRResidual(res, phi, zero, finest_level, 0); // Computes res = L(phi) - zero
  data_ops::set_value(m_scratch, 0.0);

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
			  CHF_BOX(box));
    }


    // Reset irregular cells - these contain bogus values
    const EBISBox& ebisbox = divF.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs   = ebisbox.getIrregIVS(box);

    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      divF(vof, comp) = 0.12345678E89;
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

      // Compute face centroid flux
      centroid_flux.define(ivs, ebgraph, dir, ncomp);
      for (FaceIterator faceit(ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	const FaceIndex& face   = faceit();
	const FaceStencil& sten = (*m_interp_stencils[dir][a_lvl])[dit()](face, comp);

	centroid_flux(face, comp) = 0.;
	for (int i = 0; i < sten.size(); i++){
	  const FaceIndex& iface = sten.face(i);
	  const Real iweight     = sten.weight(i);
	  
	  centroid_flux(face, comp) += iweight*flux(iface, comp);
	}
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

void cdr_solver::nonconservative_divergence(EBAMRIVData&         a_div_nc,
					    const EBAMRCellData& a_divF,
					    const EBAMRFluxData& a_face_state){
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

	divH(vof, comp)   = dc + (1-kappa)*dnc;          // On output, contains hybrid divergence
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
    if(m_mass_redist){
      level_redist.resetWeights(*a_redist_weights[lvl], comp);
    }
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

void cdr_solver::set_mass_redist(const bool a_mass_redist){
  CH_TIME("cdr_solver::set_mass_redist");
  if(m_verbosity > 5){
    pout() << m_name + "::set_mass_redist" << endl;
  }

  m_mass_redist = a_mass_redist;
}

Real cdr_solver::compute_cfl_dt(){
  CH_TIME("cdr_solver::compute_cfl_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_cfl_dt" << endl;
  }

  Real min_dt = 1.E99;

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& velo  = (*m_velo_cell[lvl])[dit()];
      const Box box          = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet ivs(box);

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

Real cdr_solver::compute_diffusive_dt(){
  CH_TIME("cdr_solver::compute_dt");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_dt" << endl;
  }

  Real min_dt = 1.E99;

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();
  const Interval interv(comp, comp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& velo     = (*m_velo_cell[lvl])[dit()];
      const Box box             = dbl.get(dit());
      const EBISBox& ebisbox    = ebisl[dit()];
      const EBGraph& ebgraph    = ebisbox.getEBGraph();
      const IntVectSet& ivs     = ebisbox.getIrregIVS(box);
      FaceStop::WhichFaces stop = FaceStop::SurroundingWithBoundary;
      const IntVectSet norm_ivs(box);

      // Regular faces
      for (int dir = 0; dir < SpaceDim; dir++){
	const EBFaceFAB& diffco = (*m_diffco[lvl])[dit()][dir];

	for (FaceIterator faceit(norm_ivs, ebgraph, dir, stop); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();
	  const Real thisdt = dx*dx/(2*SpaceDim*diffco(face, comp));

	  min_dt = Min(thisdt, min_dt);
	}
      }

      // Irregular faces
      const BaseIVFAB<Real>& diffco = (*m_diffco_eb[lvl])[dit()];
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
      	const VolIndex vof = vofit();
	const Real thisdt = dx*dx/(2*SpaceDim*diffco(vof, comp));

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
