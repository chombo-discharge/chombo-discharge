/*!
  @file   eddington_sp1.cpp
  @brief  Implementation of eddington_sp1.H
  @author Robert Marskar
  @date   Jan. 2018
  @todo   Also plot kappa on output
*/

#include "eddington_sp1.H"
#include "data_ops.H"
#include "units.H"

#include <ParmParse.H>
#include <EBAMRIO.H>

#define eddington_sp1_feature 1 // Comment Feb. 14 2018: I think we can keep this - it appears to produce the correct physics. 

eddington_sp1::eddington_sp1() : rte_solver() {
  
  this->set_verbosity(-1);
  this->set_stationary(true);
  this->set_gmg_solver_parameters();
  this->set_bottom_solver(1);
  this->set_bottom_drop(16);
  this->set_tga(true);
  this->set_reflectivity(0.0);
  this->set_time(0, 0., 0.);

  m_needs_setup = true;
}

eddington_sp1::~eddington_sp1(){
}

int eddington_sp1::query_ghost() const{
  return 3;
}

void eddington_sp1::cache_state(){
  CH_TIME("eddington_sp1::cache_state");
  if(m_verbosity > 5){
    pout() << m_name + "::cache_state" << endl;
  }

  const int ncomp = 1;
  const int finest_level = m_amr->get_finest_level();

  m_amr->allocate(m_cache, m_phase, ncomp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    m_state[lvl]->copyTo(*m_cache[lvl]);
  }
}

void eddington_sp1::set_reflectivity(const Real a_reflectivity){
  CH_TIME("eddington_sp1::set_reflectivity");
  if(m_verbosity > 5){
    pout() << m_name + "::set_reflectivity" << endl;
  }

  Real r = a_reflectivity;

  { // Get from input script
    ParmParse pp("eddington_sp1");
    pp.query("reflectivity", r);
  }

  m_r1 = r/(2.0);
  m_r2 = r/(3.0);
  
}

void eddington_sp1::set_reflection_coefficients(const Real a_r1, const Real a_r2){
  CH_TIME("eddington_sp1::set_reflection_coefficients");
  if(m_verbosity > 5){
    pout() << m_name + "::set_reflection_coefficients" << endl;
  }
  
  m_r1 = a_r1;
  m_r2 = a_r2;
}

void eddington_sp1::set_bottom_solver(const int a_whichsolver){
  CH_TIME("eddington_sp1::set_bottom_solver");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bottom_solver" << endl;
  }
  if(a_whichsolver == 0 || a_whichsolver == 1){
    m_bottomsolver = a_whichsolver;

    std::string str;
    ParmParse pp("eddington_sp1");
    pp.query("gmg_bottom_solver", str);
    if(str == "simple"){
      m_bottomsolver = 0;
    }
    else if(str == "bicgstab"){
      m_bottomsolver = 1;
    }
    
  }
  else{
    MayDay::Abort("eddington_sp1::set_bottom_solver - Unsupported solver type requested");
  }
}

void eddington_sp1::set_botsolver_smooth(const int a_numsmooth){
  CH_TIME("eddington_sp1::set_botsolver_smooth");
  if(m_verbosity > 5){
    pout() << m_name + "::set_botsolver_smooth" << endl;
  }
  CH_assert(a_numsmooth > 0);
  
  m_numsmooth = a_numsmooth;

  ParmParse pp("eddington_sp1");
  pp.query("gmg_bottom_relax", m_numsmooth);
}

void eddington_sp1::set_bottom_drop(const int a_bottom_drop){
  CH_TIME("eddington_sp1::set_bottom_drop");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bottom_drop" << endl;
  }
  
  m_bottom_drop = a_bottom_drop;

  ParmParse pp("eddington_sp1");
  pp.query("gmg_bottom_drop", m_bottom_drop);
  if(m_bottom_drop < 2){
    m_bottom_drop = 2;
  }
}

void eddington_sp1::set_tga(const bool a_use_tga){
  CH_TIME("eddington_sp1::set_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::set_tga" << endl;
  }
  
  m_use_tga = a_use_tga;
}

void eddington_sp1::set_gmg_solver_parameters(relax::which_relax a_relax_type,
					      amrmg::which_mg    a_gmg_type,      
					      const int          a_verbosity,          
					      const int          a_pre_smooth,         
					      const int          a_post_smooth,       
					      const int          a_bot_smooth,         
					      const int          a_max_iter,
					      const int          a_min_iter,
					      const Real         a_eps,               
					      const Real         a_hang){
  CH_TIME("eddington_sp1::set_gmg_solver_parameters");
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

  ParmParse pp("eddington_sp1");

  pp.query("gmg_verbosity",   m_gmg_verbosity);
  pp.query("gmg_pre_smooth",  m_gmg_pre_smooth);
  pp.query("gmg_post_smooth", m_gmg_post_smooth);
  pp.query("gmg_bott_smooth", m_gmg_bot_smooth);
  pp.query("gmg_max_iter",    m_gmg_max_iter);
  pp.query("gmg_min_iter",    m_gmg_min_iter);
  pp.query("gmg_tolerance",   m_gmg_eps);
  pp.query("gmg_hang",        m_gmg_hang);
}

void eddington_sp1::allocate_internals(){
  CH_TIME("eddington_sp1::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }
  
  const int ncomp = 1;

  m_amr->allocate(m_aco,        m_phase, ncomp);
  m_amr->allocate(m_bco,        m_phase, ncomp);
  m_amr->allocate(m_bco_irreg,  m_phase, ncomp);
  m_amr->allocate(m_state,  m_phase, ncomp);
  m_amr->allocate(m_source, m_phase, ncomp);
  m_amr->allocate(m_resid,  m_phase, ncomp);

  data_ops::set_value(m_resid,  0.0);
  data_ops::set_value(m_state,  0.0);
  data_ops::set_value(m_source, 0.0);

  this->set_aco_and_bco();
}

void eddington_sp1::deallocate_internals(){
  m_amr->deallocate(m_aco);
  m_amr->deallocate(m_bco);
  m_amr->deallocate(m_bco_irreg);
  m_amr->deallocate(m_state);
  m_amr->deallocate(m_source);
  m_amr->deallocate(m_resid);
}

void eddington_sp1::regrid(const int a_old_finest_level, const int a_new_finest_level) {
  CH_TIME("eddington_sp1::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  this->allocate_internals();

  Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->get_eb_pwl_interp(m_phase);

  m_cache[0]->copyTo(*m_state[0]); // Base level should never change. 
  for (int lvl = 1; lvl <= a_new_finest_level; lvl++){
    interpolator[lvl]->interpolate(*m_state[lvl], *m_state[lvl-1], interv);

    if(lvl <= a_old_finest_level){
      m_cache[lvl]->copyTo(*m_state[lvl]);
    }
  }

  m_needs_setup = true;
}

bool eddington_sp1::advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source, const bool a_zerophi){
  CH_TIME("eddington_sp1::advance(ebamrcell, ebamrcell)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance(ebamrcell, ebamrcell)" << endl;
  }

  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();

  if(m_needs_setup){
    this->setup_gmg();
  }

  bool converged;

  // Must have a dummy for chekcing intiial residual
  EBAMRCellData dummy;
  EBAMRCellData source;
  m_amr->allocate(dummy, m_phase, ncomp);
  m_amr->allocate(source, m_phase, ncomp);
  data_ops::set_value(dummy, 0.0);

  // Various source term manipulations. 
  data_ops::set_value(source, 0.0);
  data_ops::incr(source, a_source, 1.0);
#if eddington_sp1_feature
  data_ops::scale(source, 1./units::s_c0); // Source should be scaled by 1./c0
#endif
  if(m_stationary){ // Should kappa-scale for transient solvres
    data_ops::kappa_scale(source);
  }



  Vector<LevelData<EBCellFAB>* > phi, rhs, res, zero;
  m_amr->alias(phi,  a_state);
  m_amr->alias(rhs,  source);
  m_amr->alias(res,  m_resid);
  m_amr->alias(zero, dummy);

  if(m_stationary){
    m_gmg_solver->init(phi, rhs, finest_level, 0);

    const Real phi_resid  = m_gmg_solver->computeAMRResidual(phi,  rhs, finest_level, 0); // Incoming residual
    const Real zero_resid = m_gmg_solver->computeAMRResidual(zero, rhs, finest_level, 0); // Zero residual

    if(phi_resid > zero_resid*m_gmg_eps){ // Residual is too large
      m_gmg_solver->m_convergenceMetric = zero_resid;
      m_gmg_solver->solveNoInitResid(phi, res, rhs, finest_level, 0, a_zerophi);
      
      const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
      if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
	converged = true;
      }
    }
    else{ // Solution is already good enough
      converged = true;
    }
    m_gmg_solver->revert(phi, rhs, finest_level, 0);
  }
  else{
    if(m_use_tga){
#if eddington_sp1_feature
      m_tgasolver->oneStep(res, phi, rhs, units::s_c0*a_dt, 0, finest_level, m_time);
#else
      m_tgasolver->oneStep(res, phi, rhs, a_dt, 0, finest_level, m_time);
#endif
    }
    else{
#if eddington_sp1_feature
      m_eulersolver->oneStep(res, phi, rhs, units::s_c0*a_dt, 0, finest_level, false);
#else
      m_eulersolver->oneStep(res, phi, rhs, a_dt, 0, finest_level, false);
#endif
    }
    
    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }

    // We solve onto res, copy back to state
    data_ops::copy(a_state, m_resid);
  }

  m_amr->average_down(a_state, m_phase);
  m_amr->interp_ghost(a_state, m_phase);

  data_ops::floor(a_state, 0.0);

  return converged;
}

void eddington_sp1::setup_gmg(){
  CH_TIME("eddington_sp1::setup_gmg");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_gmg" << endl;
  }
  
  this->set_coefficients();       // Set coefficients, kappa, aco, bco
  this->setup_operator_factory(); // Set the operator factory
  this->setup_multigrid();        // Set up the AMR multigrid solver

  if(!m_stationary){
    if(m_use_tga){
      this->setup_tga();
    }
    else{
      this->setup_euler();
    }
  }

  m_needs_setup = false;
}

void eddington_sp1::set_coefficients(){
  CH_TIME("eddington_sp1::set_coefficients");
  if(m_verbosity > 5){
    pout() << m_name + "::set_coefficients" << endl;
  }

  const int ncomp = 1;
  const int ghost = 1;

  m_amr->allocate(m_aco,        m_phase, ncomp, ghost);
  m_amr->allocate(m_bco,        m_phase, ncomp, ghost);
  m_amr->allocate(m_bco_irreg,  m_phase, ncomp, ghost);

  this->set_aco_and_bco();
}

void eddington_sp1::set_aco_and_bco(){
  CH_TIME("eddington_sp1::set_aco_and_bco");
  if(m_verbosity > 5){
    pout() << m_name + "::set_aco_and_bco" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  for (int lvl = 0; lvl <= finest_level; lvl++){

    const RealVect origin = m_physdom->get_prob_lo();
    const Real dx         = m_amr->get_dx()[lvl];

    LevelData<EBCellFAB>& aco            = *m_aco[lvl];
    LevelData<EBFluxFAB>& bco            = *m_bco[lvl];
    LevelData<BaseIVFAB<Real> >& bco_irr = *m_bco_irreg[lvl];

    for (DataIterator dit = aco.dataIterator(); dit.ok(); ++dit){
      this->set_aco(aco[dit()],        origin, dx);   // Set aco = kappa
      this->set_bco_face(bco[dit()],   origin, dx);   // Set bco = 1./kappa
      this->set_bco_eb(bco_irr[dit()], origin, dx);   // Set bco = 1./kappa
    }
  }

#if eddington_sp1_feature // Different scaling for the RTE
  data_ops::scale(m_aco,       1.0);       // aco = c*kappa
  data_ops::scale(m_bco,       1.0/(3.0)); // bco = c/(3*kappa)
  data_ops::scale(m_bco_irreg, 1.0/(3.0)); // bco = c/(3*kappa)
#else // Original code before different scaling
  data_ops::scale(m_aco,       units::s_c0);       // aco = c*kappa
  data_ops::scale(m_bco,       units::s_c0/(3.0)); // bco = c/(3*kappa)
  data_ops::scale(m_bco_irreg, units::s_c0/(3.0)); // bco = c/(3*kappa)
#endif
}

void eddington_sp1::set_aco(EBCellFAB& a_aco, const RealVect a_origin, const Real a_dx){
  CH_TIME("eddington_sp1::set_aco");
  if(m_verbosity > 10){
    pout() << m_name + "::set_aco" << endl;
  }

  const int comp = 0;

  const IntVectSet ivs(a_aco.getRegion());
  const EBISBox& ebisbox = a_aco.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  
  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, a_origin);

    a_aco(vof, comp) = m_photon_group->get_kappa(pos);
  }
}

void eddington_sp1::set_bco_face(EBFluxFAB& a_bco, const RealVect a_origin, const Real a_dx){
  CH_TIME("eddington_sp1::set_bco_face");
  if(m_verbosity > 10){
    pout() << m_name + "::set_bco_face" << endl;
  }

  const int comp         = 0;
  const IntVectSet ivs   = IntVectSet(a_bco.getRegion());
  const EBISBox& ebisbox = a_bco.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  FaceStop::WhichFaces stop_crit = FaceStop::SurroundingWithBoundary;
  
  for (int dir = 0; dir < SpaceDim; dir++){
    for (FaceIterator faceit(ivs, ebgraph, dir, stop_crit); faceit.ok(); ++faceit){
      const FaceIndex& face  = faceit();
      const IntVect iv       = face.gridIndex(Side::Lo);
      const RealVect pos     = a_origin + a_dx*RealVect(iv) + 0.5*a_dx*RealVect(BASISV(dir));

      const Real kappa = m_photon_group->get_kappa(pos);
      
      a_bco[dir](face, comp) = 1./kappa;
    }
  }
}

void eddington_sp1::set_bco_eb(BaseIVFAB<Real>&          a_bco,
			       const RealVect           a_origin,
			       const Real               a_dx){
  CH_TIME("eddington_sp1::set_bco_eb");
  if(m_verbosity > 10){
    pout() << m_name + "::set_bco_eb" << endl;
  }
  
  const int comp         = 0;
  const IntVectSet ivs   = a_bco.getIVS();
  const EBGraph& ebgraph = a_bco.getEBGraph();

  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx, a_origin); // This is strictly speaking not on the boundary...

    const Real kappa = m_photon_group->get_kappa(pos);
    a_bco(vof, comp) = 1./kappa;
  }
}

void eddington_sp1::setup_operator_factory(){
  CH_TIME("eddington_sp1::setup_operator_factory");
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
  const Real alpha =  1.0;
  const Real beta  = -1.0;

  m_domfact = RefCountedPtr<robinconductivitydomainbcfactory> (new robinconductivitydomainbcfactory());
  m_ebfact  = RefCountedPtr<robinconductivityebbcfactory> (new robinconductivityebbcfactory(origin));
  m_robinco = RefCountedPtr<larsen_coefs> (new larsen_coefs(m_photon_group, m_r1, m_r2));

  m_domfact->set_coefs(m_robinco);
  m_ebfact->set_coefs(m_robinco);
  m_ebfact->set_type(stencil_type::lsq);

  // Create operator factory. 
  m_opfact = RefCountedPtr<ebconductivityopfactory> (new ebconductivityopfactory(levelgrids,
										 quadcfi,
										 alpha,
										 beta,
										 m_aco,
										 m_bco,
										 m_bco_irreg,
										 dx[0],
										 refinement_ratios,
										 m_domfact,
										 m_ebfact,
										 ghost*IntVect::Unit,
										 ghost*IntVect::Unit,
										 m_gmg_relax_type,
										 m_bottom_drop));
}

void eddington_sp1::setup_multigrid(){
  CH_TIME("eddington_sp1::setup_multigrid");
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
  m_gmg_solver->define(coar_dom, *m_opfact, botsolver, 1 + finest_level);
  m_gmg_solver->setSolverParameters(m_gmg_pre_smooth,
				    m_gmg_post_smooth,
				    m_gmg_bot_smooth,
				    m_gmg_type,
				    m_gmg_max_iter,
				    m_gmg_eps,
				    m_gmg_hang,
				    1.E-99); // Residue set through other means
  m_gmg_solver->m_imin    = m_gmg_min_iter;
  m_gmg_solver->m_verbosity = m_gmg_verbosity;
}

void eddington_sp1::setup_tga(){
  CH_TIME("eddington_sp1::setup_tga");
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

void eddington_sp1::setup_euler(){
  CH_TIME("eddington_sp1::setup_euler");
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

void eddington_sp1::compute_boundary_flux(EBAMRIVData& a_ebflux, const EBAMRCellData& a_state){
  CH_TIME("eddington_sp1::compute_boundary_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_boundary_flux" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  
  irreg_amr_stencil<eb_centroid_interp>& sten = m_amr->get_eb_centroid_interp_stencils(phase::gas);
  for(int lvl = 0; lvl <= finest_level; lvl++){
    sten.apply(*a_ebflux[lvl], *a_state[lvl], lvl, true);
  }

  m_amr->average_down(a_ebflux, m_phase);

  data_ops::scale(a_ebflux, 0.5*units::s_c0);
}

void eddington_sp1::compute_flux(EBAMRCellData& a_flux, const EBAMRCellData& a_state){
  CH_TIME("eddington_sp1::compute_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_flux" << endl;
  }

  const int finest_level = m_amr->get_finest_level();

  m_amr->compute_gradient(a_flux, a_state); // flux = grad(phi)
  for (int lvl = 0; lvl <= finest_level; lvl++){
    data_ops::divide_scalar(*a_flux[lvl], *m_aco[lvl]);   // flux = grad(phi)/(c*kappa)
    data_ops::scale(*a_flux[lvl], -units::s_c0*units::s_c0/3.0);  // flux = -c*grad(phi)/3.
  }

  m_amr->average_down(a_flux, m_phase);
  m_amr->interp_ghost(a_flux, m_phase);
}


void eddington_sp1::compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_state){
  CH_TIME("eddington_sp1::compute_density");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_density" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Interval interv(0,0);
    a_state[lvl]->copyTo(interv, *a_isotropic[lvl], interv);
  }
}

#ifdef CH_USE_HDF5
void eddington_sp1::write_plot_file(){
  CH_TIME("eddington_sp1::write_plot_file");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_file" << endl;
  }

  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_step, SpaceDim);

  const int ncomps = 3 + SpaceDim;
  Vector<string> names(ncomps);
  names[0] = "density";
  names[1] = "x-flux";
  names[2] = "y-flux";
  if(SpaceDim == 3){
    names[3] = "z-flux";
    names[4] = "isotropic source";
    names[5] = "residue";
  }
  else{
    names[3] = "isotropic source";
    names[4] = "residue";
  }

  // Compute the flux
  EBAMRCellData flux;
  m_amr->allocate(flux, m_phase, SpaceDim);
  this->compute_flux(flux, m_state);

  // Allocate output storage
  Vector<RefCountedPtr<LevelData<EBCellFAB> > > output;
  m_amr->allocate(output, m_phase, ncomps, 1);


  for (int lvl = 0; lvl < output.size(); lvl++){
    LevelData<EBCellFAB>& state  = *m_state[lvl];
    LevelData<EBCellFAB>& source = *m_source[lvl];
    LevelData<EBCellFAB>& flx    = *flux[lvl];
    LevelData<EBCellFAB>& res    = *m_resid[lvl];

    state.copyTo(Interval(0,0),          *output[lvl], Interval(0,0));
    flx.copyTo(Interval(0,SpaceDim - 1), *output[lvl], Interval(1, SpaceDim));
    source.copyTo(Interval(0,0),         *output[lvl], Interval(1 + SpaceDim, 1 + SpaceDim));
    res.copyTo(Interval(0,0),            *output[lvl], Interval(2+SpaceDim, 2+SpaceDim));
  }

  // Transform to centroid-centered
  irreg_amr_stencil<centroid_interp>& sten = m_amr->get_centroid_interp_stencils(phase::gas);
  sten.apply(output);

  // Alias this stuff
  Vector<LevelData<EBCellFAB>* > output_ptr;
  m_amr->alias(output_ptr, output);

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
