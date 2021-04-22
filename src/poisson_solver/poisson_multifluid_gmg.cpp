/*!
  @file poisson_multifluid_gmg.cpp
  @brief Implementation of poisson_multifluid_gmg.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "poisson_multifluid_gmg.H"
#include "dirichlet_func.H"
#include "data_ops.H"
#include "MFQuadCFInterp.H"
#include "MFInterfaceFAB.H"
#include "jump_bc.H"
#include "amr_mesh.H"
#include "conductivitydomainbc_wrapper_factory.H"
#include "units.H"

#include <Stencils.H>
#include <MFCellFAB.H>
#include <LayoutData.H>
#include <MFLevelDataOps.H>
#include <DirichletConductivityDomainBC.H>
#include <DirichletConductivityEBBC.H>
#include <ParmParse.H>
#include <BRMeshRefine.H>

#define POISSON_MF_GMG_TIMER 0

poisson_multifluid_gmg::poisson_multifluid_gmg(){
  m_needs_setup = true;
  m_has_mg_stuff = false;
  m_class_name  = "poisson_multifluid_gmg";
}

poisson_multifluid_gmg::~poisson_multifluid_gmg(){

}

Real poisson_multifluid_gmg::s_constant_one(const RealVect a_pos){
  return 1.0;
}

void poisson_multifluid_gmg::parse_options(){

  parse_autotune();
  parse_domain_bc();
  parse_plot_vars();
  parse_gmg_settings();
  parse_kappa_source();
}

void poisson_multifluid_gmg::parse_runtime_options(){
  parse_domain_bc();
  parse_plot_vars();
  parse_gmg_settings();
  parse_kappa_source();
}

void poisson_multifluid_gmg::parse_kappa_source(){
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("kappa_source", m_kappa_source);

}

void poisson_multifluid_gmg::parse_autotune(){
  ParmParse pp(m_class_name.c_str());

  std::string str;
  pp.get("auto_tune", str);
  m_autotune = (str == "true") ? true : false;
}

void poisson_multifluid_gmg::parse_domain_bc(){
  ParmParse pp(m_class_name.c_str());

  this->allocate_wall_bc();

  // Check each side in each direction
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Side::LoHiSide side = sit();
	
      std::string str_dir;
      if(dir == 0){
	str_dir = "x";
      }
      else if(dir == 1){
	str_dir = "y";
      }
      else if(dir == 2){
	str_dir = "z";
      }

      // Get dir/side and set accordingly
      if(side == Side::Lo){
	std::string type;
	std::string bc_string = "bc_" + str_dir + "_low";
	pp.get(bc_string.c_str(), type);
	if(type == "dirichlet_ground"){
	  this->set_dirichlet_wall_bc(dir, Side::Lo, potential::ground);
	}
	else if(type == "dirichlet_live"){
	  this->set_dirichlet_wall_bc(dir, Side::Lo, potential::live);
	}
	else if(type == "neumann"){
	  this->set_neumann_wall_bc(dir, Side::Lo, 0.0);
	}
	else if(type == "robin"){
	  this->set_robin_wall_bc(dir, Side::Lo, 0.0);
	}
	else {
	  std::string error = "poisson_multifluid_gmg::poisson_multifluid_gmg - unknown bc requested for " + bc_string;
	  MayDay::Abort(error.c_str());
	}
      }
      else if(side == Side::Hi){
	std::string type;
	std::string bc_string = "bc_" + str_dir + "_high";
	pp.get(bc_string.c_str(), type);
	if(type == "dirichlet_ground"){
	  this->set_dirichlet_wall_bc(dir, Side::Hi, potential::ground);
	}
	else if(type == "dirichlet_live"){
	  this->set_dirichlet_wall_bc(dir, Side::Hi, potential::live);
	}
	else if(type == "neumann"){
	  this->set_neumann_wall_bc(dir, Side::Hi, 0.0);
	}
	else if(type == "robin"){
	  this->set_robin_wall_bc(dir, Side::Hi, 0.0);
	}
	else {
	  std::string error = "poisson_multifluid_gmg::poisson_multifluid_gmg - unknown bc requested for " + bc_string;
	  MayDay::Abort(error.c_str());
	}
      }
    }
  }

  // Set default distribution on domain edges
  m_wall_func_x_lo = poisson_solver::s_constant_one;
  m_wall_func_x_hi = poisson_solver::s_constant_one;
  m_wall_func_y_lo = poisson_solver::s_constant_one;
  m_wall_func_y_hi = poisson_solver::s_constant_one;
#if CH_SPACEDIM==3
  m_wall_func_z_lo = poisson_solver::s_constant_one;
  m_wall_func_z_hi = poisson_solver::s_constant_one;
#endif
}

void poisson_multifluid_gmg::parse_plot_vars(){
  ParmParse pp(m_class_name.c_str());
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);
  
  m_plot_phi = false;
  m_plot_rho = false;
  m_plot_E   = false;
  m_plot_res = false;

  for (int i = 0; i < num; i++){
    if(     str[i] == "phi") m_plot_phi = true;
    else if(str[i] == "rho") m_plot_rho = true;
    else if(str[i] == "res") m_plot_res = true; 
    else if(str[i] == "E")   m_plot_E   = true;
  }
}

void poisson_multifluid_gmg::parse_gmg_settings(){
  ParmParse pp(m_class_name.c_str());

  std::string str;

  pp.get("gmg_coarsen",     m_mg_coarsen);
  pp.get("gmg_verbosity",   m_gmg_verbosity);
  pp.get("gmg_pre_smooth",  m_gmg_pre_smooth);
  pp.get("gmg_post_smooth", m_gmg_post_smooth);
  pp.get("gmg_bott_smooth", m_gmg_bot_smooth);
  pp.get("gmg_max_iter",    m_gmg_max_iter);
  pp.get("gmg_min_iter",    m_gmg_min_iter);
  pp.get("gmg_tolerance",   m_gmg_eps);
  pp.get("gmg_hang",        m_gmg_hang);
  pp.get("gmg_bottom_drop", m_bottom_drop);
  pp.get("gmg_bc_order",    m_bc_order);

  if(!(m_bc_order == 1 || m_bc_order == 2)){
    MayDay::Abort("poisson_multifluid_gmg::parse_gmg_settings - boundary condition order must be 1 or 2");
  }

  // Bottom solver
  pp.get("gmg_bottom_solver", str);
  if(str == "simple"){
    m_bottomsolver = 0;
  }
  else if(str == "bicgstab"){
    m_bottomsolver = 1;
  }
  else if(str == "gmres"){
    m_bottomsolver = 2;
    }
  else{
    MayDay::Abort("poisson_multifluid_gmg::parse_gmg_settings - unknown bottom solver requested");
  }

  // Relaxation type
  pp.get("gmg_relax_type", str);
  if(str == "gsrb"){
    m_gmg_relax_type = relax::gsrb_fast;
  }
  else if( str == "gauss_seidel"){
    m_gmg_relax_type = relax::gauss_seidel;
  }
  else{
    MayDay::Abort("poisson_multifluid_gmgcdr_gdnv::parse_gmg_settings - unsupported relaxation method requested");
  }

  // Cycle type
  pp.get("gmg_cycle", str);
  if(str == "vcycle"){
    m_gmg_type = amrmg::vcycle;
  }
  else{
    MayDay::Abort("cdr_gdnv::parse_gmg_settings - unknown cycle type requested");
  }

  // No lower than 2. 
  if(m_bottom_drop < 2){
    m_bottom_drop = 2;
  }
}

void poisson_multifluid_gmg::post_checkpoint(){
  CH_TIME("poisson_multifluid_gmg::post_checkpoint");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::post_checkpoint" << endl;
  }

  // Define the MG levels. 
  this->define_mg_levels();
}

void poisson_multifluid_gmg::allocate_internals(){
  CH_TIME("poisson_multifluid_gmg::allocate_internals");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::allocate_internals" << endl;
  }

  poisson_solver::allocate_internals();

  const int ncomp = 1;

  m_amr->allocate(m_zero,          m_realm, ncomp);
  m_amr->allocate(m_scaled_source, m_realm, ncomp);
  m_amr->allocate(m_scaled_sigma,  m_realm, phase::gas, ncomp);

  data_ops::set_value(m_zero, 0.0);
}

bool poisson_multifluid_gmg::solve(const bool a_zerophi){
  CH_TIME("poisson_multifluid_gmg::solve");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::solve" << endl;
  }
  
  const bool converged = this->solve(m_state, m_source, m_sigma, a_zerophi);


  return converged;
}

bool poisson_multifluid_gmg::solve(MFAMRCellData&       a_state,
				   const MFAMRCellData& a_source,
				   const EBAMRIVData&   a_sigma,
				   const bool           a_zerophi){
  CH_TIME("poisson_multifluid_gmg::solve(mfamrcell, mfamrcell");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::solve(mfamrcell, mfamrcell)" << endl;
  }

  bool converged = false;

  if(m_needs_setup){
    this->setup_gmg(); // This does everything, allocates coefficients, gets bc stuff and so on
  }

  const Real t0 = MPI_Wtime();

  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();


  const Real t1 = MPI_Wtime();

  // Do the scaled space charge density
  data_ops::copy(m_scaled_source, a_source);
  data_ops::scale(m_scaled_source, 1./(units::s_eps0));
  data_ops::scale(m_scaled_source, 1./(m_length_scale*m_length_scale));

  if(m_kappa_source){ // Scale source by kappa
    data_ops::kappa_scale(m_scaled_source);
  }

  // Do the scaled surface charge
  data_ops::copy(m_scaled_sigma, a_sigma);
  data_ops::scale(m_scaled_sigma, 1./(m_length_scale*m_length_scale));
  m_opfact->set_jump(m_scaled_sigma, 1.0/units::s_eps0);

  const Real t2 = MPI_Wtime();
  
#if 0 // Debug
  MayDay::Warning("poisson_multifluid_gmg::solve - debug mode");
  m_opfact->set_jump(0.0, 1.0);
  data_ops::set_value(source, 0.0);
#endif

#if 0 // Check NaN/Inf input
  EBAMRCellData phiGas = m_amr->alias(phase::gas,   a_state);
  EBAMRCellData phiSol = m_amr->alias(phase::solid, a_state);

  EBAMRCellData srcGas = m_amr->alias(phase::gas,   m_scaled_source);
  EBAMRCellData srcSol = m_amr->alias(phase::solid, m_scaled_source);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    EBLevelDataOps::checkNANINF(*phiGas[lvl]);
    EBLevelDataOps::checkNANINF(*phiSol[lvl]);

    EBLevelDataOps::checkNANINF(*srcGas[lvl]);
    EBLevelDataOps::checkNANINF(*srcSol[lvl]);
  }
#endif

  // Aliasing
  Vector<LevelData<MFCellFAB>* > phi, cpy, rhs, res, zero;
  m_amr->alias(phi,     a_state);
  m_amr->alias(rhs,     m_scaled_source);
  m_amr->alias(res,     m_resid);
  m_amr->alias(zero,    m_zero);

  // GMG solve. Use phi = zero as initial metric. Want to reduce this by m_gmg_eps
  //  m_gmg_solver.init(phi, rhs, finest_level, 0);
  const Real phi_resid  = m_gmg_solver.computeAMRResidual(phi,  rhs, finest_level, 0);
  const Real zero_resid = m_gmg_solver.computeAMRResidual(zero, rhs, finest_level, 0);

  m_converged_resid = zero_resid*m_gmg_eps;

  const Real t3 = MPI_Wtime();

  if(phi_resid > m_converged_resid){ // Residual is too large, recompute solution
    m_gmg_solver.m_convergenceMetric = zero_resid;
    m_gmg_solver.solveNoInitResid(phi, res, rhs, finest_level, 0, a_zerophi);

    const int status = m_gmg_solver.m_exitStatus;   // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{ // Solution is already converged
    converged = true;

  }

  const Real t4 = MPI_Wtime();

#if 0 // Why is this required??? Is it because of op->zeroCovered()????
  Real new_resid = m_gmg_solver.computeAMRResidual(phi, rhs, finest_level, 0);
  new_resid = m_gmg_solver.computeAMRResidual(phi, rhs, finest_level, 0);
#endif

  
#if 0 // Debug. Solve again
  m_gmg_solver.solveNoInitResid(phi, res, rhs, finest_level, 0, a_zerophi);
#endif


  m_gmg_solver.revert(phi, rhs, finest_level, 0);

  m_amr->average_down(a_state, m_realm);
  m_amr->interp_ghost(a_state, m_realm);

  const Real t5 = MPI_Wtime();

#if POISSON_MF_GMG_TIMER
  const Real T = t5 - t0;
  pout() << endl;
  pout() << "poisson_multifluid_gmg::solve breakdown" << endl;
  pout() << "alloc :     " << 100.*(t1-t0)/T << "%" << endl;
  pout() << "set jump:   " << 100.*(t2-t1)/T << "%" << endl;
  pout() << "resid:      " << 100.*(t3-t2)/T << "%" << endl;
  pout() << "solve:      " << 100.*(t4-t3)/T << "%" << " = " << t4 - t3 << endl;
  pout() << "revert/avg: " << 100.*(t5-t4)/T << "%" << endl;
  pout() << "Total time: " << T << endl;
#endif

  return converged;
}

int poisson_multifluid_gmg::query_ghost() const {
  return 3; // Need this many cells
}

void poisson_multifluid_gmg::set_potential(Real (*a_potential)(const Real a_time)){
  poisson_solver::set_potential(a_potential);

  const RealVect origin  = m_amr->get_prob_lo();

  m_bcfunc = RefCountedPtr<dirichlet_func>(new dirichlet_func(m_potential, s_constant_one, origin));

  const int ixlo = wall_bc::map_bc(0, Side::Lo);
  const int ixhi = wall_bc::map_bc(0, Side::Hi);
  const int iylo = wall_bc::map_bc(1, Side::Lo);
  const int iyhi = wall_bc::map_bc(1, Side::Hi);
#if CH_SPACEDIM==3
  const int izlo = wall_bc::map_bc(2, Side::Lo);
  const int izhi = wall_bc::map_bc(2, Side::Hi);
#endif

  m_wall_bcfunc.resize(2*SpaceDim);
  m_wall_bcfunc[ixlo] = RefCountedPtr<BaseBCFuncEval>(new dirichlet_func(m_potential, m_wall_func_x_lo, origin));
  m_wall_bcfunc[ixhi] = RefCountedPtr<BaseBCFuncEval>(new dirichlet_func(m_potential, m_wall_func_x_hi, origin));
  m_wall_bcfunc[iylo] = RefCountedPtr<BaseBCFuncEval>(new dirichlet_func(m_potential, m_wall_func_y_lo, origin));
  m_wall_bcfunc[iyhi] = RefCountedPtr<BaseBCFuncEval>(new dirichlet_func(m_potential, m_wall_func_y_hi, origin));
#if CH_SPACEDIM==3
  m_wall_bcfunc[izlo] = RefCountedPtr<BaseBCFuncEval>(new dirichlet_func(m_potential, m_wall_func_z_lo, origin));
  m_wall_bcfunc[izhi] = RefCountedPtr<BaseBCFuncEval>(new dirichlet_func(m_potential, m_wall_func_z_hi, origin));
#endif
}

void poisson_multifluid_gmg::set_time(const int a_step, const Real a_time, const Real a_dt){
  poisson_solver::set_time(a_step, a_time, a_dt);
  m_bcfunc->set_time(a_time);
  for (int i = 0; i < m_wall_bcfunc.size(); i++){
    dirichlet_func* func = static_cast<dirichlet_func*> (&(*m_wall_bcfunc[i]));
    func->set_time(a_time);
  }
  
  if(!m_opfact.isNull()){
    m_opfact->set_time(&m_time);
  }
}

void poisson_multifluid_gmg::auto_tune(){
  CH_TIME("poisson_multifluid_gmg::auto_tune");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::auto_tune" << endl;
  }

  if(m_autotune){
    const ProblemDomain coarsest = m_amr->get_domains()[0];

    if(m_verbosity > 2){
      pout() << "poisson_multifluid_gmg - autotuning solver parameters (bottom drop only)" << endl;
    }
      
      
    Real best_time;
    int best_drop;
    int best_relax;
    
    best_time = 1.E99;
    for (int drop = 2; drop <= coarsest.size(0)/2; drop *= 2){
      this->set_bottom_drop(drop);
      const Real t1 = MPI_Wtime();
      data_ops::set_value(m_state, 0.0);
      this->solve(true);
      const Real t2 = MPI_Wtime();

      if(t2 - t1 < best_time){
	best_time = t2 - t1;
	best_drop = drop;
      }
    }

    this->set_bottom_drop(best_drop);

#if 0 // This is very dangerous because it might overestimate the necessary number of relaxations
    best_time = 1.E99;
    for (int relax = 8; relax <= 32; relax += 4){
      this->set_gmg_solver_parameters(m_gmg_relax_type,
				      m_gmg_type,
				      m_gmg_verbosity,
				      relax,
				      relax,
				      relax,
				      m_gmg_max_iter,
				      m_gmg_min_iter,
				      m_gmg_eps,
				      m_gmg_hang);

      data_ops::set_value(m_state, 0.0);
      const Real t1 = MPI_Wtime();
      this->solve(true);
      const Real t2 = MPI_Wtime();
      pout() << t2 - t1 << endl;

      if(t2 - t1 < best_time){
	best_time = t2 - t1;
	best_relax = relax;
      }
    }
#endif
  }
}

void poisson_multifluid_gmg::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("poisson_multifluid_gmg::regrid");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::regrid" << endl;
  }
  poisson_solver::regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  m_needs_setup = true;
}

void poisson_multifluid_gmg::register_operators(){
  CH_TIME("poisson_multifluid_gmg::register_operators");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::register_operators" << endl;
  }

  if(m_amr.isNull()){
    MayDay::Abort("poisson_multifluid_gmg::register_operators - need to set amr_mesh!");
  }
  else{
    m_amr->register_operator(s_eb_coar_ave,     m_realm, phase::gas);
    m_amr->register_operator(s_eb_coar_ave,     m_realm, phase::solid);
    m_amr->register_operator(s_eb_fill_patch,   m_realm, phase::gas);
    m_amr->register_operator(s_eb_fill_patch,   m_realm, phase::solid);
    m_amr->register_operator(s_eb_pwl_interp,   m_realm, phase::gas);
    m_amr->register_operator(s_eb_pwl_interp,   m_realm, phase::solid);
    m_amr->register_operator(s_eb_quad_cfi,     m_realm, phase::gas);
    m_amr->register_operator(s_eb_quad_cfi,     m_realm, phase::solid);
    m_amr->register_operator(s_eb_irreg_interp, m_realm, phase::gas);
    m_amr->register_operator(s_eb_irreg_interp, m_realm, phase::solid);
    m_amr->register_operator(s_eb_flux_reg,     m_realm, phase::gas);
    m_amr->register_operator(s_eb_flux_reg,     m_realm, phase::solid);
  }
}

void poisson_multifluid_gmg::set_bottom_solver(const int a_whichsolver){
  CH_TIME("poisson_multifluid_gmg::set_bottom_solver");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::set_bottom_solver" << endl;
  }
  
  if(a_whichsolver == 0 || a_whichsolver == 1 || a_whichsolver == 2){
    m_bottomsolver = a_whichsolver;

    std::string str;
    ParmParse pp("poisson_multifluid");
    pp.query("gmg_bottom_solver", str);
    if(str == "simple"){
      m_bottomsolver = 0;
    }
    else if(str == "bicgstab"){
      m_bottomsolver = 1;
    }
    else if(str == "gmres"){
      m_bottomsolver=2;
    }
  }
  else{
    MayDay::Abort("poisson_multifluid_gmg::set_bottom_solver - Unsupported solver type requested");
  }
}

void poisson_multifluid_gmg::set_botsolver_smooth(const int a_numsmooth){
  CH_TIME("poisson_multifluid_gmg::set_botsolver_smooth");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::set_botsolver_smooth" << endl;
  }
  CH_assert(a_numsmooth > 0);
  
  m_numsmooth = a_numsmooth;

  ParmParse pp("poisson_multifluid");
  pp.query("gmg_bottom_relax", m_numsmooth);
}

void poisson_multifluid_gmg::set_bottom_drop(const int a_bottom_drop){
  CH_TIME("poisson_multifluid_gmg::set_bottom_drop");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::set_bottom_drop" << endl;
  }
  
  m_bottom_drop = a_bottom_drop;

  ParmParse pp("poisson_multifluid");
  pp.query("gmg_bottom_drop", m_bottom_drop);
  if(m_bottom_drop < 2){
    m_bottom_drop = 2;
  }
}

void poisson_multifluid_gmg::set_gmg_solver_parameters(relax      a_relax_type,
						       amrmg      a_gmg_type,      
						       const int  a_verbosity,          
						       const int  a_pre_smooth,         
						       const int  a_post_smooth,       
						       const int  a_bot_smooth,         
						       const int  a_max_iter,
						       const int  a_min_iter,
						       const Real a_eps,               
						       const Real a_hang){
  CH_TIME("poisson_multifluid_gmg::set_gmg_solver_parameters");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::set_gmg_solver_parameters" << endl;
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

  ParmParse pp("poisson_multifluid");

  pp.query("gmg_verbosity",   m_gmg_verbosity);
  pp.query("gmg_pre_smooth",  m_gmg_pre_smooth);
  pp.query("gmg_post_smooth", m_gmg_post_smooth);
  pp.query("gmg_bott_smooth", m_gmg_bot_smooth);
  pp.query("gmg_max_iter",    m_gmg_max_iter);
  pp.query("gmg_min_iter",    m_gmg_min_iter);
  pp.query("gmg_tolerance",   m_gmg_eps);
  pp.query("gmg_hang",        m_gmg_hang);

#if 1 // Overriding these because this appears to be the only things that works with mfconductivityop
  m_gmg_relax_type = relax::gsrb_fast;
  m_gmg_type       = amrmg::vcycle;
#endif
}



void poisson_multifluid_gmg::set_coefficients(){
  CH_TIME("poisson_multifluid_gmg::set_coefficients");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::set_coefficients" << endl;
  }

  const int ncomps = 1;
  const int ghosts = 1;
  const Real eps0  = m_compgeom->get_eps0();
  
  m_amr->allocate(m_aco,       m_realm, ncomps);
  m_amr->allocate(m_bco,       m_realm, ncomps);
  m_amr->allocate(m_bco_irreg, m_realm, ncomps);

  data_ops::set_value(m_aco,       0.0);  // Always zero for poisson equation, but that is done from alpha. 
  data_ops::set_value(m_bco,       eps0); // Will override this later
  data_ops::set_value(m_bco_irreg, eps0); // Will override this later

  this->set_permittivities(m_compgeom->get_dielectrics());
}

void poisson_multifluid_gmg::set_permittivities(const Vector<dielectric>& a_dielectrics){
  CH_TIME("poisson_multifluid_gmg::set_permittivities");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::set_permittivities" << endl;
  }

  if(a_dielectrics.size() > 0){
    const RealVect origin  = m_amr->get_prob_lo();
    const Vector<Real> dx  = m_amr->get_dx();
    const int finest_level = m_amr->get_finest_level();

    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[lvl];
      
      LevelData<EBFluxFAB> bco;
      LevelData<BaseIVFAB<Real> > bco_irreg;

      mfalias::aliasMF(bco,       phase::solid, *m_bco[lvl]);
      mfalias::aliasMF(bco_irreg, phase::solid, *m_bco_irreg[lvl]);

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	EBFluxFAB& perm          = bco[dit()];
	BaseIVFAB<Real>& perm_eb = bco_irreg[dit()];
	const Box box            = dbl.get(dit());

	this->set_face_perm(perm,  box, origin, dx[lvl], a_dielectrics);
	this->set_eb_perm(perm_eb, box, origin, dx[lvl], a_dielectrics);
      }
    }
  }
}

void poisson_multifluid_gmg::set_face_perm(EBFluxFAB&                a_perm,
					   const Box&                a_box,
					   const RealVect&           a_origin,
					   const Real&               a_dx,
					   const Vector<dielectric>& a_dielectrics){
  CH_TIME("poisson_multifluid_gmg::set_face_perm");
  if(m_verbosity > 10){
    pout() << "poisson_multifluid_gmg::set_face_perm" << endl;
  }

  const int comp         = 0;
  const IntVectSet ivs   = IntVectSet(a_perm.getRegion());
  const EBISBox& ebisbox = a_perm.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  FaceStop::WhichFaces stop_crit = FaceStop::SurroundingWithBoundary;
  
  for (int dir = 0; dir < SpaceDim; dir++){

    // Regular clels
    Box facebox = a_box;
    facebox.surroundingNodes(dir);
    BaseFab<Real>& perm_fab = a_perm[dir].getSingleValuedFAB();
    for (BoxIterator bit(facebox); bit.ok(); ++bit){
      const IntVect iv = bit();
      const RealVect pos     = a_origin + a_dx*RealVect(iv) + 0.5*a_dx*RealVect(BASISV(dir));

	Real dist   = 1.E99;
	int closest = 0;
	for (int i = 0; i < a_dielectrics.size(); i++){
	  const RefCountedPtr<BaseIF> func = a_dielectrics[i].get_function();

	  const Real cur_dist = func->value(pos);
	
	  if(cur_dist <= dist){
	    dist = cur_dist;
	    closest = i;
	  }
	}
	perm_fab(iv, comp) = a_dielectrics[closest].get_permittivity(pos);
    }
    

    // Irregular cells
    const IntVectSet irreg = ebisbox.getIrregIVS(a_box);
    for (FaceIterator faceit(irreg, ebgraph, dir, stop_crit); faceit.ok(); ++faceit){
      const FaceIndex& face  = faceit();
      const IntVect iv       = face.gridIndex(Side::Lo);
      const RealVect pos     = a_origin + a_dx*RealVect(iv) + 0.5*a_dx*RealVect(BASISV(dir));
      
      Real dist   = 1.E99;
      int closest = 0;
      for (int i = 0; i < a_dielectrics.size(); i++){
	const RefCountedPtr<BaseIF> func = a_dielectrics[i].get_function();

	const Real cur_dist = func->value(pos);
	
	if(cur_dist <= dist){
	  dist = cur_dist;
	  closest = i;
	}
      }
      a_perm[dir](face, comp) = a_dielectrics[closest].get_permittivity(pos);
    }
  }
}


void poisson_multifluid_gmg::set_eb_perm(BaseIVFAB<Real>&          a_perm,
					 const Box&                a_box,
					 const RealVect&           a_origin,
					 const Real&               a_dx,
					 const Vector<dielectric>& a_dielectrics){
  CH_TIME("poisson_multifluid_gmg::set_eb_perm");
  if(m_verbosity > 10){
    pout() << "poisson_multifluid_gmg::set_eb_perm" << endl;
  }
  
  const int comp         = 0;
  const IntVectSet ivs   = a_perm.getIVS();
  const EBGraph& ebgraph = a_perm.getEBGraph();

  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx, a_origin); // This is strictly speaking not on the boundary...
      
    Real dist   = 1.E99;
    int closest = 0;
    for (int i = 0; i < a_dielectrics.size(); i++){
      const RefCountedPtr<BaseIF> func = a_dielectrics[i].get_function();

      const Real cur_dist = func->value(pos);
	
      if(cur_dist <= dist){
	dist = cur_dist;
	closest = i;
      }
    }

    a_perm(vof, comp) = a_dielectrics[closest].get_permittivity(pos);
  }
}

void poisson_multifluid_gmg::define_mg_levels(){
  CH_TIME("poisson_multifluid_gmg::define_mg_level");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::define_mg_levels" << endl;
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  const int coar_ref = 2;
  
  m_mg_domains.resize(0);
  m_mg_grids.resize(0);
  m_mg_mflg.resize(0);
  
  m_mg_eblg.resize(phase::num_phases);

  m_mg_eblg[phase::gas].resize(0);
  m_mg_eblg[phase::solid].resize(0);

  // Get some stuff from amr_mesh on how to decompose the levels
  const Vector<ProblemDomain>& domains = m_amr->get_domains();
  const int max_box_size               = m_amr->get_max_box_size();
  const int blocking_factor            = m_amr->get_blocking_factor();
  const int num_ebghost                = m_amr->get_eb_ghost();

  int num_coar       = 0;
  bool has_coar      = true;
  ProblemDomain fine = domains[0];

  // Coarsen problem domains and create grids
  while(num_coar < m_mg_coarsen || !has_coar){

    // Check if we can coarsen
    const ProblemDomain coar = fine.coarsen(coar_ref);
    const Box coar_box       = coar.domainBox();
    for (int dir = 0; dir < SpaceDim; dir++){
      if(coar_box.size()[dir] < max_box_size || coar_box.size()[dir]%max_box_size != 0){
	has_coar = false;
      }
    }

    if(has_coar){
      // Split the domain into pieces, then order and load balance them
      Vector<Box> boxes;
      Vector<int> proc_assign;
      domainSplit(coar, boxes, max_box_size, blocking_factor);
      mortonOrdering(boxes);
      load_balance::make_balance(proc_assign, boxes);

      // Add problem domain and grid
      m_mg_domains.push_back(coar);
      m_mg_grids.push_back(DisjointBoxLayout(boxes, proc_assign, coar));

      // Define the EBLevelGrids
      const int idx = m_mg_grids.size() - 1; // Last element added
      if(!ebis_gas.isNull()){
	m_mg_eblg[phase::gas].push_back(RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_mg_grids[idx],
										    m_mg_domains[idx],
										    num_ebghost,
										    ebis_gas)));
      }
      if(!ebis_sol.isNull()){
	m_mg_eblg[phase::solid].push_back(RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_mg_grids[idx],
										      m_mg_domains[idx],
										      num_ebghost,
										      ebis_sol)));
      }

      // Define the MFLevelGrid object
      Vector<EBLevelGrid> eblgs;
      if(!ebis_gas.isNull()) eblgs.push_back(*m_mg_eblg[phase::gas][idx]);
      if(!ebis_sol.isNull()) eblgs.push_back(*m_mg_eblg[phase::solid][idx]);
      
      m_mg_mflg.push_back(RefCountedPtr<MFLevelGrid> (new MFLevelGrid(m_mfis, eblgs)));

      // Next iterate
      fine = coar;
      num_coar++;
    }
    else{
      break;
    }
  }
  
  m_has_mg_stuff = true;
}

void poisson_multifluid_gmg::setup_gmg(){
  CH_TIME("poisson_multifluid_gmg::setup_gmg");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::setup_gmg" << endl;
  }

  if(!m_has_mg_stuff){
    this->define_mg_levels();     // Define MG levels. These don't change during regrids so we only need to set them once. 
  }
  this->set_coefficients();       // Set coefficients
  this->setup_operator_factory(); // Set the operator factory
  this->setup_solver();           // Set up the AMR multigrid solver

  m_needs_setup = false;
}

void poisson_multifluid_gmg::setup_operator_factory(){
  CH_TIME("poisson_multifluid_gmg::setup_operator_factory");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::setup_operator_factory" << endl;
  }

  const int nphases                      = m_mfis->num_phases();
  const int finest_level                 = m_amr->get_finest_level();
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids(m_realm);
  const Vector<int>& refinement_ratios   = m_amr->get_ref_rat();
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real>& dx                 = m_amr->get_dx();
  const RealVect& origin                 = m_amr->get_prob_lo();

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_mfis->get_ebis(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_mfis->get_ebis(phase::solid);

  // This stuff is needed for the operator factory
  Vector<MFLevelGrid>    mflg(1 + finest_level);
  Vector<MFQuadCFInterp> mfquadcfi(1 + finest_level);
  Vector<MFFastFluxReg>  mffluxreg(1 + finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<EBLevelGrid>                    eblg_phases(nphases);
    Vector<RefCountedPtr<EBQuadCFInterp> > quadcfi_phases(nphases);
    Vector<RefCountedPtr<EBFluxRegister> > fluxreg_phases(nphases);

    if(!ebis_gas.isNull()) eblg_phases[phase::gas]   = *(m_amr->get_eblg(m_realm, phase::gas)[lvl]);
    if(!ebis_sol.isNull()) eblg_phases[phase::solid] = *(m_amr->get_eblg(m_realm, phase::solid)[lvl]);

    if(!ebis_gas.isNull()) quadcfi_phases[phase::gas]   = (m_amr->get_old_quadcfi(m_realm, phase::gas)[lvl]);
    if(!ebis_sol.isNull()) quadcfi_phases[phase::solid] = (m_amr->get_old_quadcfi(m_realm, phase::solid)[lvl]);

    if(!ebis_gas.isNull()) fluxreg_phases[phase::gas]   = (m_amr->get_flux_reg(m_realm, phase::gas)[lvl]);
    if(!ebis_sol.isNull()) fluxreg_phases[phase::solid] = (m_amr->get_flux_reg(m_realm, phase::solid)[lvl]);
    
    mflg[lvl].define(m_mfis, eblg_phases);
    mfquadcfi[lvl].define(quadcfi_phases);
    mffluxreg[lvl].define(fluxreg_phases);
  }

  // Appropriate coefficients for poisson equation
  const Real alpha =  1.0; // Aco is zero. 
  const Real beta  = -1.0;

  RefCountedPtr<BaseDomainBCFactory> domfact = RefCountedPtr<BaseDomainBCFactory> (NULL);

  const IntVect ghost_phi = this->query_ghost()*IntVect::Unit;
  const IntVect ghost_rhs = this->query_ghost()*IntVect::Unit;

  conductivitydomainbc_wrapper_factory* bcfact = new conductivitydomainbc_wrapper_factory();
  RefCountedPtr<dirichlet_func> pot = RefCountedPtr<dirichlet_func> (new dirichlet_func(m_potential,
											s_constant_one,
											RealVect::Zero));
  
  bcfact->set_wallbc(m_wallbc);
  bcfact->set_potentials(m_wall_bcfunc);
  domfact = RefCountedPtr<BaseDomainBCFactory> (bcfact);

  Vector<MFLevelGrid> mg_levelgrids;
  for (int lvl = 0; lvl < m_mg_mflg.size(); lvl++){
    mg_levelgrids.push_back(*m_mg_mflg[lvl]);
  }

  // Set the length scale for the Poisson equation. This is equivalent to solving the Poisson equation
  // on the domain [-1,1] in the x-direction
  m_length_scale = 2./(domains[0].size(0)*dx[0]);

  int relax_type;
  if(m_gmg_relax_type == relax::jacobi){
    relax_type = 0;
  }
  else if(m_gmg_relax_type == relax::gauss_seidel){
    relax_type = 1;
  }
  else if(m_gmg_relax_type == relax::gsrb_fast){
    relax_type = 2;
  }

  // Create factory and set potential
  m_opfact = RefCountedPtr<mfconductivityopfactory> (new mfconductivityopfactory(m_mfis,
										 mflg,
										 mfquadcfi,
										 mffluxreg,
										 refinement_ratios,
										 grids,
										 m_aco,
										 m_bco,
										 m_bco_irreg,
										 alpha,
										 beta,
										 m_length_scale,
										 dx[0]*m_length_scale,
										 domains[0],
										 domfact,
										 origin,
										 ghost_phi,
										 ghost_rhs,
										 m_bc_order,
										 relax_type,
										 m_bottom_drop,
										 1 + finest_level,
										 mg_levelgrids));
  CH_assert(!m_bcfunc.isNull());
  m_opfact->set_electrodes(m_compgeom->get_electrodes(), m_bcfunc);
}

void poisson_multifluid_gmg::setup_solver(){
  CH_TIME("poisson_multifluid_gmg::setup_solver");
  if(m_verbosity > 5){
    pout() << "poisson_multifluid_gmg::setup_solver" << endl;
  }

  const int finest_level       = m_amr->get_finest_level();
  const ProblemDomain coar_dom = m_amr->get_domains()[0];

  // Select bottom solver
  LinearSolver<LevelData<MFCellFAB> >* botsolver = NULL;
  if(m_bottomsolver == 0){
    m_mfsolver.setNumSmooths(m_numsmooth);
    botsolver = &m_mfsolver;
  }
  else if(m_bottomsolver==1){
    botsolver = &m_bicgstab;
#if 1
    if(m_mfis->num_phases() == 2){ // BiCGStab doesn't work with multifluid (yet)
      botsolver = &m_mfsolver;
      
      if(m_verbosity > 0){
	pout() << "poisson_multifluid_gmg::poisson_multifluid_gmg - BiCGStab not supported for multifluid" << endl;
      }
    }
#endif
  }
  else if(m_bottomsolver==2){
    botsolver = &m_gmres;
  }

  int gmg_type;
  if(m_gmg_type == amrmg::full){
    gmg_type = 0;
  }
  else if(m_gmg_type == amrmg::vcycle){
    gmg_type = 1;
  }
  else if(m_gmg_type == amrmg::fcycle){
    gmg_type = 2;
  }


  m_gmg_solver.define(coar_dom, *m_opfact, botsolver, 1 + finest_level);
  m_gmg_solver.setSolverParameters(m_gmg_pre_smooth,
				   m_gmg_post_smooth,
				   m_gmg_bot_smooth,
				   gmg_type,
				   m_gmg_max_iter,
				   m_gmg_eps,
				   m_gmg_hang,
				   1.E-99); // Norm thresh will be set via eps
  m_gmg_solver.m_imin = m_gmg_min_iter;
  m_gmg_solver.m_verbosity = m_gmg_verbosity;


  // Dummies for init
  const int ncomp = 1;
  MFAMRCellData dummy1, dummy2;
  m_amr->allocate(dummy1, m_realm, ncomp);
  m_amr->allocate(dummy2, m_realm, ncomp);
  data_ops::set_value(dummy1, 0.0);
  data_ops::set_value(dummy2, 0.0);

  // Aliasing
  Vector<LevelData<MFCellFAB>* > phi, rhs;
  m_amr->alias(phi, dummy1);
  m_amr->alias(rhs, dummy2);

  // Init solver
  m_gmg_solver.init(phi, rhs, finest_level, 0);
}

MFAMRCellData& poisson_multifluid_gmg::get_aco(){
  return m_aco;
}

MFAMRFluxData& poisson_multifluid_gmg::get_bco(){
  return m_bco;
}

MFAMRIVData& poisson_multifluid_gmg::get_bco_irreg(){
  return m_bco_irreg;
}

void poisson_multifluid_gmg::set_needs_setup(const bool& a_needs_setup){
  m_needs_setup = a_needs_setup;
}
