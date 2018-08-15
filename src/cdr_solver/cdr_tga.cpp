
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

#define CDR_TGA_DEBUG_TIMER 0

cdr_tga::cdr_tga() : cdr_solver() {
  this->set_gmg_solver_parameters();
  this->set_bottom_solver(1);
  this->set_bottom_drop(2);
  this->set_tga(true);

  m_name = "cdr_tga";
}

cdr_tga::~cdr_tga(){

}

void cdr_tga::advance(EBAMRCellData& a_state, const Real& a_dt){
  CH_TIME("cdr_tga::advance(state, dt)");
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

void cdr_tga::advance_diffusion(EBAMRCellData& a_state, const Real a_dt){
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
    EBAMRCellData src, phi, divF;
    m_amr->allocate(src,  m_phase, ncomp);
    m_amr->allocate(phi,  m_phase, ncomp); 

    for (int lvl = 0; lvl <= finest_level; lvl++){
      data_ops::set_value(*src[lvl],  0.0);
      data_ops::set_value(*phi[lvl],  0.0);
      data_ops::incr(*phi[lvl], *a_state[lvl], 1.0);
    }

    // Do the aliasing stuff
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, phi);
    m_amr->alias(old_state, a_state);
    m_amr->alias(source,    src);

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

void cdr_tga::set_bottom_solver(const int a_whichsolver){
  CH_TIME("cdr_tga::set_bottom_solver");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bottom_solver" << endl;
  }
  
  if(a_whichsolver == 0 || a_whichsolver == 1){
    m_bottomsolver = a_whichsolver;
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
  
  m_numsmooth = a_numsmooth;
}

void cdr_tga::set_bottom_drop(const int a_bottom_drop){
  CH_TIME("cdr_tga::set_bottom_drop");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bottom_drop" << endl;
  }
  
  m_bottom_drop = a_bottom_drop;
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
}

void cdr_tga::setup_gmg(){
  CH_TIME("cdr_tga::setup_gmg");
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
  pout() << endl;
  pout() << "cdr_tga::compute_divJ" << endl;
  pout() << "t1 - t0 = " << t1 - t0 << endl;
  pout() << "t2 - t1 = " << t2 - t1 << endl;
  pout() << "t3 - t2 = " << t3 - t2 << endl;
  pout() << "Total = " << t3 - t0 << endl;
#endif
}

void cdr_tga::compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_redist){
  CH_TIME("cdr_tga::compute_divF(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divF(divF, state)" << endl;
  }

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
    a_state[lvl]->copyTo(*weights[lvl]);
  }

  const Real t1 = MPI_Wtime();

  // Compute the advective derivative
  this->average_velo_to_faces(m_velo_face, m_velo_cell);            // Average cell-centered velocities to face centers
  this->advect_to_faces(face_state, a_state, a_extrap_dt);          // Face extrapolation to cell-centered faces
  this->conservative_divergence(a_divF, face_state, m_velo_face);   // a_divF holds the conservative divergence
  this->nonconservative_divergence(div_nc, a_divF, face_state);     // Compute non-conservative divergence
  this->hybrid_divergence(a_divF, mass_diff, div_nc);               // Make divF = hybrid divergence. Compute mass diff.
  this->increment_flux_register(face_state, m_velo_face);           // Increment flux registers
  this->increment_redist(mass_diff);                                // Increment redistribution objects

  const Real t2 = MPI_Wtime();

  // Mass weights.
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    data_ops::incr(*weights[lvl], *a_state[lvl], 1.0);
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
  pout() << endl;
  pout() << "Calilng cdr_tga::compute_divF" << endl;
  pout() << "t1 - t0 = " << t1 - t0 << endl;
  pout() << "t2 - t1 = " << t2 - t1 << endl;
  pout() << "t3 - t2 = " << t3 - t2 << endl;
  pout() << "Total = " << t3 - t0 << endl;
  pout() << endl;
#endif
}

void cdr_tga::compute_divD(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state){
  CH_TIME("cdr_tga::compute_divD");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divD" << endl;
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
