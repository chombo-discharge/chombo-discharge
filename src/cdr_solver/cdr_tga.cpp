/*!
  @file   cdr_tga.cpp
  @brief  Implementation of cdr_tga.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "cdr_tga.H"
#include "data_ops.H"

#include <ParmParse.H>
#include <EBArith.H>
#include <EBAMRIO.H>
#include <NeumannConductivityEBBC.H>
#include <NeumannConductivityDomainBC.H>

#define CDR_TGA_DEBUG_TIMER 1

cdr_tga::cdr_tga() : cdr_solver() {
  m_name       = "cdr_tga";
  m_class_name = "cdr_tga";
}

cdr_tga::~cdr_tga(){

}

void cdr_tga::advance(EBAMRCellData& a_state, const Real& a_dt){
  CH_TIME("cdr_tga::advance(state, dt)");
  if(m_verbosity > 1){
    pout() << m_name + "::advance(state, dt)" << endl;
  }
  
  if(this->is_diffusive()){
    bool m_needs_setup = false;
    if(m_needs_setup){
      this->setup_gmg();
    }
    if(m_use_tga){
      this->advance_tga(a_state, a_dt);
    }
  }
  else{
    const Real midpoint = 0.5;
    cdr_solver::advance_rk2(a_state, a_dt, midpoint);
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

void cdr_tga::advance_diffusion(EBAMRCellData& a_state, EBAMRCellData& a_error, const Real a_dt){
  CH_TIME("cdr_tga::advance_diffusion");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_diffusion" << endl;
  }

  if(m_diffusive){
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    // Create a source term = S = 0.0;
    EBAMRCellData src, phi_euler, phi_tga;
    m_amr->allocate(src,       m_phase, ncomp);
    m_amr->allocate(phi_euler, m_phase, ncomp);
    m_amr->allocate(phi_tga,   m_phase, ncomp); 

    data_ops::set_value(src, 0.0);
    data_ops::copy(phi_tga, a_state);

    // Do the aliasing stuff
    Vector<LevelData<EBCellFAB>* > euler_state, tga_state, old_state, source;
    m_amr->alias(euler_state, phi_euler);
    m_amr->alias(tga_state,   phi_tga);
    m_amr->alias(old_state, a_state);
    m_amr->alias(source,    src);

    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);

    // TGA solve
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(tga_state, old_state, source, a_dt, 0, finest_level, false);

    // Euler solve
    data_ops::copy(phi_euler, phi_tga); // I assume that phi_tga is a very good initial guess
    m_eulersolver->resetAlphaAndBeta(alpha, beta);
    m_eulersolver->oneStep(euler_state, old_state, source, a_dt, 0, finest_level, false);

    // Compute error
    data_ops::copy(a_error, phi_euler);
    data_ops::scale(a_error, -1.0);
    data_ops::incr(a_error, phi_tga, 1.0);

    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }

    data_ops::copy(a_state, phi_tga);
  }
}


void cdr_tga::advance_tga(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt){
  CH_TIME("cdr_tga::advance_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga" << endl;
  }

  if(m_diffusive){
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    // Create a source term = S = 0.0;
    EBAMRCellData src;
    m_amr->allocate(src, m_phase, ncomp);
    data_ops::set_value(src, 0.0);

    // Do the aliasing stuff
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    src);

    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);

    // TGA solve
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);

    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
}

void cdr_tga::advance_euler(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt){
  CH_TIME("cdr_tga::advance_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_euler" << endl;
  }
  if(m_diffusive){
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    // Create a source term = S = 0.0;
    EBAMRCellData src;
    m_amr->allocate(src, m_phase, ncomp);
    data_ops::set_value(src, 0.0);

    // Do the aliasing stuff
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    src);

    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);

    // TGA solve
    m_eulersolver->resetAlphaAndBeta(alpha, beta);
    m_eulersolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);

    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
}

void cdr_tga::advance_euler(EBAMRCellData&       a_new_state,
			    const EBAMRCellData& a_old_state,
			    const EBAMRCellData& a_source,
			    const Real           a_dt){
  CH_TIME("cdr_tga::advance_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_euler" << endl;
  }
  
  if(m_diffusive){
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    // Do the aliasing stuff
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    a_source);

#if 0 // Debug
    setup_gmg();
#endif
    const Real alpha = 0.0;
    const Real beta  = 1.0;

    m_eulersolver->resetAlphaAndBeta(alpha, beta);
    data_ops::set_value(m_diffco_eb, 0.0);
    
    // Euler solve
    m_eulersolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);


    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
}

void cdr_tga::advance_tga(EBAMRCellData& a_state, const Real a_dt){
  CH_TIME("cdr_tga::advance_tga");
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

  this->compute_divF(divF, a_state, 0.5*a_dt, true); // extrapolate div(n*v) to 0.5*a_dt and use it as a source term for TGA

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

  data_ops::set_value(m_diffco_eb, 0.0);
  
  if(m_use_tga){
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
  }
  else{
    m_eulersolver->resetAlphaAndBeta(alpha, beta);
    m_eulersolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
  }

  const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
  if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
    converged = true;
  }

  data_ops::copy(a_state, phi);
}

void cdr_tga::advance_tga(EBAMRCellData&       a_new_state,
			    const EBAMRCellData& a_old_state,
			    const EBAMRCellData& a_source,
			    const Real           a_dt){
  CH_TIME("cdr_tga::advance_tga(full)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga(full)" << endl;
  }
  
  if(m_diffusive){
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    // Do the aliasing stuff
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    a_source);

    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);

    // TGA solve
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);

    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{
    data_ops::copy(a_new_state, a_old_state);
  }
}

void cdr_tga::set_bottom_solver(const int a_whichsolver){
  CH_TIME("cdr_tga::set_bottom_solver");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bottom_solver" << endl;
  }

  if(a_whichsolver == 0 || a_whichsolver == 1){
    m_bottomsolver = a_whichsolver;

    std::string str;
    ParmParse pp("cdr_tga");
    pp.get("gmg_bottom_solver", str);
    if(str == "simple"){
      m_bottomsolver = 0;
    }
    else if(str == "bicgstab"){
      m_bottomsolver = 1;
    }
  }
  else{
    MayDay::Abort("cdr_tga::set_bottom_solver - Unsupported solver type requested");
  }
}

void cdr_tga::set_botsolver_smooth(const int a_numsmooth){
  CH_TIME("cdr_tga::set_botsolver_smooth");
  if(m_verbosity > 5){
    pout() << m_name + "::set_botsolver_smooth" << endl;
  }
  CH_assert(a_numsmooth > 0);


  ParmParse pp("cdr_tga");
  pp.query("gmg_bottom_relax", m_numsmooth);

  if(m_numsmooth < 0){
    m_numsmooth = 0;
  }
}

void cdr_tga::set_bottom_drop(const int a_bottom_drop){
  CH_TIME("cdr_tga::set_bottom_drop");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bottom_drop" << endl;
  }

  ParmParse pp("cdr_tga");
  pp.get("gmg_bottom_drop", m_bottom_drop);

  if(m_bottom_drop < 2){
    m_bottom_drop = 2;
  }
}

void cdr_tga::set_tga(const bool a_use_tga){
  CH_TIME("cdr_tga::set_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::set_tga" << endl;
  }
  
  m_use_tga = a_use_tga;

  ParmParse pp("cdr_tga");
  std::string str;
  if(pp.contains("use_tga")){
    pp.get("use_tga", str);
    if(str == "true"){
      m_use_tga = true;
    }
    else if(str == "false"){
      m_use_tga = false;
    }
  }
}

void cdr_tga::set_gmg_solver_parameters(relax::which_relax a_relax_type,
					amrmg::which_mg    a_gmg_type,      
					const int          a_verbosity,          
					const int          a_pre_smooth,         
					const int          a_post_smooth,       
					const int          a_bot_smooth,         
					const int          a_max_iter,
					const int          a_min_iter,
					const Real         a_eps,               
					const Real         a_hang){
  CH_TIME("cdr_tga::set_gmg_solver_parameters");
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

  ParmParse pp("cdr_tga");

  pp.get("gmg_verbosity",   m_gmg_verbosity);
  pp.get("gmg_pre_smooth",  m_gmg_pre_smooth);
  pp.get("gmg_post_smooth", m_gmg_post_smooth);
  pp.get("gmg_bott_smooth", m_gmg_bot_smooth);
  pp.get("gmg_max_iter",    m_gmg_max_iter);
  pp.get("gmg_min_iter",    m_gmg_min_iter);
  pp.get("gmg_tolerance",   m_gmg_eps);
  pp.get("gmg_hang",        m_gmg_hang);
}

void cdr_tga::setup_gmg(){
  CH_TIME("cdr_tga::setup_gmg");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_gmg" << endl;
  }

  this->setup_operator_factory();
  this->setup_multigrid();

  this->setup_tga();
  this->setup_euler();
}

void cdr_tga::setup_operator_factory(){
  CH_TIME("cdr_tga::setup_operator_factory");
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

  Vector<EBLevelGrid> mg_levelgrids;
  Vector<RefCountedPtr<EBLevelGrid> >& mg_eblg = m_amr->get_mg_eblg(m_phase);
  for (int lvl = 0; lvl < mg_eblg.size(); lvl++){
    mg_levelgrids.push_back(*mg_eblg[lvl]);
  }

  // Appropriate coefficients. These don't matter right now. 
  const Real alpha =  1.0;
  const Real beta  =  1.0;

  // Default is Neumann. This might change in the future. 
  RefCountedPtr<NeumannConductivityDomainBCFactory> domfact = RefCountedPtr<NeumannConductivityDomainBCFactory>
    (new NeumannConductivityDomainBCFactory());
  RefCountedPtr<NeumannConductivityEBBCFactory> ebfact  = RefCountedPtr<NeumannConductivityEBBCFactory>
    (new NeumannConductivityEBBCFactory());

  domfact->setValue(0.0);
  ebfact->setValue(0.0);
  
  // Create operator factory.
  data_ops::set_value(m_aco, 1.0); // We're usually solving (1 - dt*nabla^2)*phi^(k+1) = phi^k + dt*S^k so aco=1
  data_ops::set_value(m_diffco_eb, 0.0);
  m_opfact = RefCountedPtr<ebconductivityopfactory> (new ebconductivityopfactory(levelgrids,
										 quadcfi,
										 alpha,
										 beta,
										 m_aco,
										 m_diffco,
										 m_diffco_eb,
										 dx[0],
										 refinement_ratios,
										 domfact,
										 ebfact,
										 ghost*IntVect::Unit,
										 ghost*IntVect::Unit,
										 m_gmg_relax_type,
										 m_bottom_drop,
										 -1,
										 mg_levelgrids));
}

void cdr_tga::setup_multigrid(){
  CH_TIME("cdr_tga::setup_multigrid");
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

void cdr_tga::setup_tga(){
  CH_TIME("cdr_tga::setup_tga");
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

void cdr_tga::setup_euler(){
  CH_TIME("cdr_tga::setup_euler");
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

void cdr_tga::compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt){
  CH_TIME("cdr_tga::compute_divJ(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divJ(divF, state)" << endl;
  }

  const Real t0 = MPI_Wtime();
  const int comp  = 0;
  const int ncomp = 1;

  data_ops::set_value(a_divJ, 0.0);

  EBAMRCellData advective_term;
  EBAMRCellData diffusion_term;
  m_amr->allocate(advective_term, m_phase, ncomp);
  if(this->is_diffusive()){
    m_amr->allocate(diffusion_term, m_phase, ncomp);
  }

  const Real t1 = MPI_Wtime();

  // Compute advective term
  if(this->is_mobile()){
    this->compute_divF(advective_term, a_state, 0.0, true);
    data_ops::incr(a_divJ, advective_term, 1.0);
  }

  const Real t2 = MPI_Wtime();

  // Add in diffusion term
  if(this->is_diffusive()){
    this->compute_divD(diffusion_term, a_state); // This already does refluxing. 
    data_ops::incr(a_divJ, diffusion_term, -1.0);
  }

  const Real t3 = MPI_Wtime();
#if CDR_TGA_DEBUG_TIMER
  const Real T = t3 - t0;
  pout() << endl;
  pout() << "cdr_tga::compute_divJ" << endl;
  pout() << "Allocations:  " << 100.*(t1 - t0)/T << "%" << endl;
  pout() << "Compute divF: " << 100.*(t2 - t1)/T << "%" << endl;
  pout() << "compute divD: " << 100.*(t3 - t2)/T << "%" << endl;
  pout() << "Total:        " << T << endl;
#endif
}

void cdr_tga::compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_redist){
  CH_TIME("cdr_tga::compute_divF(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divF(divF, state)" << endl;
  }

  if(m_mobile){
    const Real t0 = MPI_Wtime();
  
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
      a_state[lvl]->localCopyTo(*weights[lvl]);
    }

    const Real t1 = MPI_Wtime();

    // Compute the advective derivative
    this->average_velo_to_faces(m_velo_face, m_velo_cell);            // Average cell-centered velocities to face centers
    const Real t10 = MPI_Wtime();
    this->advect_to_faces(face_state, a_state, a_extrap_dt);          // Face extrapolation to cell-centered faces
    const Real t11 = MPI_Wtime();
    this->conservative_divergence(a_divF, face_state, m_velo_face);   // a_divF holds the conservative divergence
    const Real t12 = MPI_Wtime();
    this->nonconservative_divergence(div_nc, a_divF, face_state);     // Compute non-conservative divergence
    const Real t13 = MPI_Wtime();
    this->hybrid_divergence(a_divF, mass_diff, div_nc);               // Make divF = hybrid divergence. Compute mass diff.
    const Real t14 = MPI_Wtime();
    this->increment_flux_register(face_state, m_velo_face);           // Increment flux registers
    const Real t15 = MPI_Wtime();
    this->increment_redist(mass_diff);                                // Increment redistribution objects
    const Real t16 = MPI_Wtime();
    const Real t2 = MPI_Wtime();

    // Mass weights.
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      //      data_ops::incr(*weights[lvl], *a_state[lvl], 1.0);
    }

    const bool ebcf = m_amr->get_ebcf();
    if(ebcf){ 
      this->coarse_fine_increment(mass_diff);     // Increment the coarse-fine redistribution objects
      this->increment_redist_flux();              // Increment flux registers with the redistribution stuff
      this->coarse_fine_redistribution(a_divF);   // Redistribute
      this->reflux(a_divF);
    }
    else{
      this->hyperbolic_redistribution(a_divF, mass_diff, weights);  // Redistribute mass into hybrid divergence
      this->reflux(a_divF);                                         // Reflux at coarse-fine interfaces
    }

    const Real t3 = MPI_Wtime();

#if CDR_TGA_DEBUG_TIMER
    const Real T = t3 - t0;
    pout() << endl;
    pout() << "cdr_tga::compute_divF breakdown" << endl;
    pout() << "Allocations:      " << 100.*(t1 - t0)/T   << "%" << endl;
    pout() << "Average to faces: " << 100.*(t10 - t1)/T  << "%" << endl;
    pout() << "Advect to faces:  " << 100.*(t11 - t10)/T << "%" << endl;
    pout() << "Cons. Div.:       " << 100.*(t12 - t11)/T << "%" << endl;
    pout() << "NonCons. Div.:    " << 100.*(t13 - t12)/T << "%" << endl;
    pout() << "Hybrid. Div.:     " << 100.*(t14 - t13)/T << "%" << endl;
    pout() << "Flux reg:         " << 100.*(t15 - t14)/T << "%" << endl;
    pout() << "Incr redist:      " << 100.*(t16 - t15)/T << "%" << endl;
    pout() << "Reflux-redist:    " << 100.*(t3 - t2)/T   << "%" << endl;
    pout() << "Total time:       " << T << endl;
#endif


  }
  else{
    data_ops::set_value(a_divF, 0.0);
  }

#if 0 // Debug
  Real volume;
  const Real sum = EBLevelDataOps::noKappaSumLevel(volume, *a_divF[0], 0, m_amr->get_domains()[0]);

  std::cout << sum << std::endl;
#endif

}

void cdr_tga::compute_divD(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state){
  CH_TIME("cdr_tga::compute_divD");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divD" << endl;
  }

  if(m_diffusive){
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    // Copy a_state because AMRMultiGrid may do exchange operations
    EBAMRCellData clone;
    m_amr->allocate(clone, m_phase, ncomp);
    for (int lvl = 0; lvl <= finest_level; lvl++){
      a_state[lvl]->localCopyTo(*clone[lvl]);
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

    m_amr->average_down(a_diffusive_term, m_phase);
    m_amr->interp_ghost(a_diffusive_term, m_phase);
  }
  else{
    data_ops::set_value(a_diffusive_term, 0.0);
  }
}

void cdr_tga::parse_gmg_settings(){
  ParmParse pp(m_class_name.c_str());

  std::string str;

  pp.get("gmg_verbosity",   m_gmg_verbosity);
  pp.get("gmg_pre_smooth",  m_gmg_pre_smooth);
  pp.get("gmg_post_smooth", m_gmg_post_smooth);
  pp.get("gmg_bott_smooth", m_gmg_bot_smooth);
  pp.get("gmg_max_iter",    m_gmg_max_iter);
  pp.get("gmg_min_iter",    m_gmg_min_iter);
  pp.get("gmg_tolerance",   m_gmg_eps);
  pp.get("gmg_hang",        m_gmg_hang);
  pp.get("gmg_bottom_drop", m_bottom_drop);

  // Bottom solver
  pp.get("gmg_bottom_solver", str);
  if(str == "simple"){
    m_bottomsolver = 0;
  }
  else if(str == "bicgstab"){
    m_bottomsolver = 1;
  }
  else{
    MayDay::Abort("cdr_tga::parse_gmg_settings - unknown bottom solver requested");
  }

  // Relaxation type
  pp.get("gmg_relax_type", str);
  if(str == "gsrb"){
    m_gmg_relax_type = relax::gsrb_fast;
  }
  else if(str == "jacobi"){
    m_gmg_relax_type = relax::jacobi;
  }
  else if(str == "gauss_seidel"){
    m_gmg_relax_type = relax::gauss_seidel;
  }
  else{
    MayDay::Abort("cdr_tga::parse_gmg_settings - unknown relaxation method requested");
  }

  // Cycle type
  pp.get("gmg_cycle", str);
  if(str == "vcycle"){
    m_gmg_type = amrmg::vcycle;
  }
  else{
    MayDay::Abort("cdr_tga::parse_gmg_settings - unknown cycle type requested");
  }

  // No lower than 2. 
  if(m_bottom_drop < 2){
    m_bottom_drop = 2;
  }
}

void cdr_tga::parse_diffusion(){
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("stochastic_diffusion", str);
  
  m_stochastic_diffusion = (str == "true") ? true : false;
  m_stochastic_advection = false;
}

void cdr_tga::parse_rng_seed(){
  ParmParse pp(m_class_name.c_str());
  pp.get("seed", m_seed);
  if(m_seed < 0) m_seed = std::chrono::system_clock::now().time_since_epoch().count();

  m_rng = new std::mt19937_64(m_seed);
}

void cdr_tga::parse_plotmode(){
  ParmParse pp(m_class_name.c_str());

  m_plot_numbers = false;
  
  std::string str;
  pp.get("plot_mode", str);
  if(str == "density"){
    m_plot_numbers = false;
  }
  else if(str == "numbers"){
    m_plot_numbers = true;
  }
}

void cdr_tga::GWN_diffusion_source(EBAMRCellData& a_ransource, const EBAMRCellData& a_cell_states){
  CH_TIME("cdr_tga::GWN_diffusion_source");
  if(m_verbosity > 5){
    pout() << m_name + "::GWN_diffusion_source" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  
  EBAMRFluxData ranflux;
  EBAMRFluxData GWN;

  // Putting this in here because it is really, really important that we don't send in any stupid ghost cells or negative
  // values for the diffusion advance
  EBAMRCellData& states = const_cast<EBAMRCellData&> (a_cell_states);
  m_amr->average_down(states, m_phase);
  m_amr->interp_ghost(states, m_phase);
  data_ops::floor(states, 0.0);
  
  m_amr->allocate(ranflux,   m_phase, ncomp);
  m_amr->allocate(GWN,       m_phase, ncomp);
  
  data_ops::set_value(a_ransource, 0.0);
  data_ops::set_value(ranflux,     0.0);
  data_ops::set_value(GWN,         0.0);

  this->fill_GWN(GWN, 1.0);                             // Gaussian White Noise
  this->smooth_heaviside_faces(ranflux, a_cell_states); // ranflux = phis
   data_ops::multiply(ranflux, m_diffco);                // ranflux = D*phis
   data_ops::scale(ranflux, 2.0);                        // ranflux = 2*D*phis
   data_ops::square_root(ranflux);                       // ranflux = sqrt(2*D*phis)

#if 1 // Debug
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    for (int dir = 0; dir <SpaceDim; dir++){
      Real max, min;
      EBLevelDataOps::getMaxMin(max, min, *ranflux[lvl], 0, dir);
      if(min < 0.0 || max < 0.0){
	MayDay::Abort("stop - negative face value");
      }
    }
  }
#endif
  data_ops::multiply(ranflux, GWN);                     // Holds random, cell-centered flux

  // Source term. 
  // I want to re-use conservative_divergence(), but that also computes with the EB fluxes. Back that up
  // first and use the already written and well-tested routine. Then copy back
  EBAMRIVData backup;
  m_amr->allocate(backup, m_phase, 1);
  data_ops::copy(backup, m_ebflux);
  data_ops::set_value(m_ebflux, 0.0);
  data_ops::set_value(a_ransource, 0.0);
  conservative_divergence(a_ransource, ranflux); // Compute the conservative divergence. This also refluxes. 
  data_ops::copy(m_ebflux, backup);

  //  data_ops::set_value(a_ransource, 0.0);

  
#if 1 // Debug
  m_amr->average_down(a_ransource, m_phase);
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    if(EBLevelDataOps::checkNANINF(*a_ransource[lvl])){
      MayDay::Abort("something is wrong");
    }
  }
#endif
}

void cdr_tga::smooth_heaviside_faces(EBAMRFluxData& a_face_states, const EBAMRCellData& a_cell_states){
  CH_TIME("cdr_tga::smooth_heaviside_faces");
  if(m_verbosity > 5){
    pout() << m_name + "::smooth_heaviside_faces" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    const Real vol               = pow(dx, SpaceDim);
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box& box = dbl.get(dit());
      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();

      for (int dir = 0; dir < SpaceDim; dir++){
	EBFaceFAB& face_states       = (*a_face_states[lvl])[dit()][dir];
	const EBCellFAB& cell_states = (*a_cell_states[lvl])[dit()];

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
    EBLevelDataOps::setCoveredVal(*a_face_states[lvl], 0.0);
  }
}

void cdr_tga::fill_GWN(EBAMRFluxData& a_noise, const Real a_sigma){
  CH_TIME("cdr_tga::fill_GWN");
  if(m_verbosity > 5){
    pout() << m_name + "::fill_GWN" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  std::normal_distribution<double> GWN(0.0, a_sigma);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
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
	  noise_reg(iv, comp) = GWN(*m_rng)*ivol;
	}

	// Irregular faces
	const IntVectSet& irreg = ebisbox.getIrregIVS(box);
	const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingNoBoundary;
	for (FaceIterator faceit(irreg, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  noise(face, comp) = GWN(*m_rng)*ivol;
	}
      }
    }
  }
}

void cdr_tga::GWN_advection_source(EBAMRCellData& a_ransource, const EBAMRCellData& a_cell_states){
  CH_TIME("cdr_tga::GWN_advection_source");
  if(m_verbosity > 5){
    pout() << m_name + "::GWN_advection_source" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  
  EBAMRFluxData ranflux;
  EBAMRFluxData GWN;
  
  m_amr->allocate(ranflux,   m_phase, ncomp);
  m_amr->allocate(GWN,       m_phase, ncomp);
  
  data_ops::set_value(a_ransource, 0.0);
  data_ops::set_value(ranflux,   0.0);
  data_ops::set_value(GWN,       0.0);
      
  this->fill_GWN(GWN, 1.0);                             // Gaussian White Noise
  this->smooth_heaviside_faces(ranflux, a_cell_states); // ranflux = phis
  data_ops::square_root(ranflux);                       // ranflux = sqrt(phis)
  data_ops::multiply(ranflux, GWN);                     // ranflux = sqrt(phis)*W/sqrt(dV)
  data_ops::multiply(ranflux, m_velo_face);             // ranflux = v*sqrt(phis)*W/sqrt(dV)

  // Source term. 
  // I want to re-use conservative_divergence(), but that also computes with the EB fluxes. Back that up
  // first and use the already written and well-tested routine. Then copy back.
  EBAMRIVData backup;
  m_amr->allocate(backup, m_phase, 1);
  data_ops::copy(backup, m_ebflux);
  data_ops::set_value(m_ebflux, 0.0);
  conservative_divergence(a_ransource, ranflux); // Compute the conservative divergence. This also refluxes. 
  data_ops::copy(m_ebflux, backup);
}

void cdr_tga::write_plot_data(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("cdr_tga::write_plot_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_data" << endl;
  }

  if(m_plot_phi) {
    if(!m_plot_numbers){ // Regular write
      write_data(a_output, a_comp, m_state, true);
    }
    else{ // Scale, write, and scale back
     
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_state[lvl], (pow(m_amr->get_dx()[lvl], 3)));
      }
      write_data(a_output, a_comp, m_state,    false);
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_state[lvl], 1./(pow(m_amr->get_dx()[lvl], 3)));
      }
    }
  }

  if(m_plot_dco && m_diffusive) { // Need to compute the cell-centerd stuff first
    data_ops::set_value(m_scratch, 0.0);
    data_ops::average_face_to_cell(m_scratch, m_diffco, m_amr->get_domains());
    write_data(a_output, a_comp, m_scratch,   false);
  }

  if(m_plot_src) {
    if(!m_plot_numbers){
      write_data(a_output, a_comp, m_source,    false);
    }
    else { // Scale, write, and scale back
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_source[lvl], (pow(m_amr->get_dx()[lvl], 3)));
      }
      write_data(a_output, a_comp, m_source,    false);
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_source[lvl], 1./(pow(m_amr->get_dx()[lvl], 3)));
      }
    }
  }

  if(m_plot_vel && m_mobile) {
    write_data(a_output, a_comp, m_velo_cell, false);
  }

  // Plot EB fluxes
  if(m_plot_ebf && m_mobile){
    data_ops::set_value(m_scratch, 0.0);
    data_ops::incr(m_scratch, m_ebflux, 1.0);
    write_data(a_output, a_comp, m_scratch, false);
  }
}
