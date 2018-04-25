/*!
  @file   splitstep_rk2_tga_trapz.cpp
  @brief  Implementation of splitstep_rk2_tga_trapz.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "splitstep_rk2_tga_trapz.H"
#include "splitstep_rk2_tga_trapz_storage.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "data_ops.H"
#include "units.H"

#include <EBArith.H>
#include <ParmParse.H>

extern "C" void newton_point_trapz_fort(int *i);

typedef splitstep_rk2_tga_trapz::cdr_storage     cdr_storage;
typedef splitstep_rk2_tga_trapz::poisson_storage poisson_storage;
typedef splitstep_rk2_tga_trapz::rte_storage     rte_storage;
typedef splitstep_rk2_tga_trapz::sigma_storage   sigma_storage;

splitstep_rk2_tga_trapz::splitstep_rk2_tga_trapz(){
  m_alpha    = 1.0;
  m_tol_x    = 1.E-6;
  m_tol_f    = 1.E-6;
  m_max_iter = 10;
  m_simpi    = true;

  // Basically only for debugging
  m_do_advection = true;
  m_do_diffusion = true;
  m_do_source    = true;
  m_do_rte       = true;
  m_do_poisson   = true;
  
  {
    ParmParse pp("splitstep_rk2_tga_trapz");

    std::string str;
    
    pp.query("alpha",                 m_alpha);
    pp.query("solution_tolerance",    m_tol_x);
    pp.query("function_tolerance",    m_tol_f);
    pp.query("max_newton_iter",       m_max_iter);
    
    if(pp.contains("coupling")){
      pp.query("coupling", str);
      if(str == "semi_implicit"){
	m_simpi = true; 
      }
      else if(str == "implicit"){
	m_simpi = false;
      }
    }

    if(pp.contains("turn_off_advection")){
      pp.get("turn_off_advection", str);
      if(str == "true"){
	m_do_advection = false;
	if(m_verbosity > 2){
	  pout() << "splitstep_rk2_tga_trapz::splitstep_rk2_tga_trapz - Turning off advection" << endl;
	}
      }
    }
    
    if(pp.contains("turn_off_diffusion")){
      pp.get("turn_off_diffusion", str);
      if(str == "true"){
	m_do_diffusion = false;
	if(m_verbosity > 2){
	  pout() << "splitstep_rk2_tga_trapz::splitstep_rk2_tga_trapz - Turning off diffusion" << endl;
	}
      }
    }
    
    if(pp.contains("turn_off_source")){
      pp.get("turn_off_source", str);
      if(str == "true"){
	m_do_source = false;

	if(m_verbosity > 2){
	  pout() << "splitstep_rk2_tga_trapz::splitstep_rk2_tga_trapz - Turning off source" << endl;
	}
      }
    }

    if(pp.contains("turn_off_rte")){
      pp.get("turn_off_rte", str);
      if(str == "true"){
	m_do_rte = false;

	if(m_verbosity > 2){
	  pout() << "splitstep_rk2_tga_trapz::splitstep_rk2_tga_trapz - Turning off rte" << endl;
	}
      }
    }

    if(pp.contains("turn_off_poisson")){
      pp.get("turn_off_poisson", str);
      if(str == "true"){
	m_do_poisson = false;

	if(m_verbosity > 2){
	  pout() << "splitstep_rk2_tga_trapz::splitstep_rk2_tga_trapz - Turning off poisson" << endl;
	}
      }
    }
  }
}

splitstep_rk2_tga_trapz::~splitstep_rk2_tga_trapz(){

}

RefCountedPtr<cdr_storage>& splitstep_rk2_tga_trapz::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& splitstep_rk2_tga_trapz::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real splitstep_rk2_tga_trapz::restrict_dt(){
  return 1.E99;
}

Real splitstep_rk2_tga_trapz::advance(const Real a_dt){
  CH_TIME("splitstep_rk2_tga_trapz::advance");
  if(m_verbosity > 2){
    pout() << "splitstep_rk2_tga_trapz::advance" << endl;
  }

  this->cache_solutions(); // Cache old solutions. Used in case the time step is rejected. 

  if(m_do_advection){              // Rules for advective advance: cdr solvers contain the advected states, the 
    this->advance_advection(a_dt); // poisson solver contains the potential after advection. Sigma solver is updated.
  }                                // Scratch storage contains nothing but junk (except for cached stuff)
  
  if(m_do_diffusion){              // Rules for diffusion advance: cdr solvers contain the advected AND diffused states,
    this->advance_diffusion(a_dt); // the poisson solver contains the potential after this. Nothing happens to sigma
  }                                // and RTE. Scratch storage contains nothing but junk.

  
  this->compute_E_at_start_of_time_step(); // Compute the electric field using the advected/diffused states. This is
                                           // needed for both the RTE update and the first source term iteration

#if 0 // Todo
  if(m_do_rte){ // Must be implemented. 
    if(m_rte->is_stationary()){              // Rules for RTE update: After this advance, the RTE solvers contain
      this->solve_rte_using_solver_states(); // the solution by using the 
    }
    else{
      MayDay::Abort("splitstep_rk2_tga_trapz::advance - transient RTE is not supported for this time stepper");
    }
  }
#endif

  if(m_do_source){
    this->advance_sources(a_dt);   // Source term advance. 
  }

  // Put solver back in useable state so that we can reliably compute the next time step. 
  this->compute_cdr_velocities();
  this->compute_cdr_diffusion();
  this->compute_cdr_sources();
  
  return a_dt;
}

void splitstep_rk2_tga_trapz::advance_advection(const Real a_dt){
  CH_TIME("splitstep_rk2_tga_trapz::advance_advection");
  if(m_verbosity > 2){
    pout() << "splitstep_rk2_tga_trapz::advance_advection" << endl;
  }

  // Not implemented...
}

void splitstep_rk2_tga_trapz::advance_diffusion(const Real a_dt){
  CH_TIME("splitstep_rk2_tga_trapz::advance_diffusion");
  if(m_verbosity > 2){
    pout() << "splitstep_rk2_tga_trapz::advance_diffusion" << endl;
  }

  // Not implemented...
}

void splitstep_rk2_tga_trapz::advance_sources(const Real a_dt){
  CH_TIME("splitstep_rk2_tga_trapz::advance_sources");
  if(m_verbosity > 2){
    pout() << "splitstep_rk2_tga_trapz::advance_sources" << endl;
  }

  const Real time       = m_time + a_dt;
  const int num_photons = m_plaskin->get_num_photons();
  const int num_species = m_plaskin->get_num_species();


  // Extra storage for grad(|E|) and grad(n) which we will need later
  EBAMRCellData grad_E;
  EBAMRCellData E_norm;
  Vector<EBAMRCellData*> grad_cdr(num_species); 
  m_amr->allocate(grad_E, m_cdr->get_phase(), SpaceDim);
  m_amr->allocate(E_norm, m_cdr->get_phase(), 1);
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    grad_cdr[idx] = new EBAMRCellData();                            
    m_amr->allocate(*grad_cdr[idx], m_cdr->get_phase(), SpaceDim);  
  }


  // Preparation for Newton iteration
  this->compute_epsj();                      // Get tolerances for the finite difference evaluation 
  this->advance_semi_implicit_newton(a_dt);  // Advance cdr semi-implicitly. Initial guess in cdr_storage->phi. Updates source term
  this->compute_trapz_rhs(a_dt);             // Compute the right-hand side for the trapezoidal discretization. Uses source term. 
  if(m_do_poisson){
    this->compute_semi_implicit_potential(); // Poisson solve using semi-implicitly advanced states
    this->compute_E_at_start_of_time_step(); // Recompute field by using the semi-implicit Poisson advance
  }
  if(m_do_rte){
    this->compute_semi_implicit_rte(a_dt);   // RTE solve by using semi-implictly advances states
  }


  // Here, we have advanced the cdr states explictly recomputed the electric field and RTE solutions in a semi-implicit way
  // so that we have estimates for all variables at time (k+1). We can now start Newton iteration using.
  // 
  // RULES: Here is where storage should go
  // cdr_storage->m_phi:       Newton iterates
  // cdr_storage->m_k1:        RHS of the trapezoidal discretization (i.e. n_k - 0.5*dt_S_k)
  // poisson_solver->m_phi:    Electric field updates (implicit or semi-implicit)
  // poisson_scratch->E_cell:  Electric field (implicit or semi-implicit)
  // rte_solver->m_state:      RTE solutions (implicit or semi-implicit)




  int iter = 0;           // Number of iterations
  bool converged = false; // Convergence track

  
  while(!converged && iter < m_max_iter){

    Vector<Real> x(num_species); // Newton solution
    Vector<Real> p(num_species); // Increment
    const EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();

    // Compute grad(|E|)
    data_ops::vector_length(E_norm, E_cell);          // Compute |E|
    m_amr->average_down(E_norm, m_cdr->get_phase());  // Average down
    m_amr->interp_ghost(E_norm, m_cdr->get_phase());  // Interpolate ghost cells
    m_amr->compute_gradient(grad_E, E_norm);          // Compute gradient
    m_amr->average_down(grad_E, m_cdr->get_phase());  // Average down gradient
    m_amr->interp_ghost(grad_E, m_cdr->get_phase());  // Interpolate gradient ghost cells

    // Compute grad(n)
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      const RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
      const EBAMRCellData& iter_state = storage->get_phi();
      
      m_amr->compute_gradient(*grad_cdr[idx], iter_state);     // Compute grad()
      m_amr->average_down(*grad_cdr[idx], m_cdr->get_phase()); // Average down
      m_amr->interp_ghost(*grad_cdr[idx], m_cdr->get_phase()); // Interpolate ghost cells
    }

    // Compute cell-centered grad(E), grad(n)

    const int finest_level = m_amr->get_finest_level();
    for (int lvl = 0; lvl <= finest_level; lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids()[lvl];
      const EBISLayout& ebisl      = m_amr->get_ebisl(phase::gas)[lvl];
      const Real dx                = m_amr->get_dx()[lvl];
      
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const Box box            = dbl.get(dit());
	const EBISBox& ebisbox   = ebisl[dit()];
	const EBGraph& ebgraph   = ebisbox.getEBGraph();
	const IntVectSet ivs_irr = ebisbox.getIrregIVS(box);
	const IntVectSet ivs_reg = IntVectSet(box) - ivs_irr;

	
	RealVect pos;
	RealVect E;
	RealVect Egrad;
	Vector<Real> x(num_species); // Iterates
	Vector<Real> p(num_species); // Corrections
	Vector<Real> rhs(num_species);
	Vector<Real> rte_densities(num_photons);
	Vector<RealVect> cdr_gradients(num_species);

	// Regular cells. Straightforward stuff
	for (VoFIterator vofit(ivs_reg, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();

	  pos   = EBArith::getVofLocation(vof, dx*RealVect::Unit, m_physdom->get_prob_lo());
	  E     = RealVect(D_DECL((*E_cell[lvl])[dit()](vof, 0),
				  (*E_cell[lvl])[dit()](vof, 1),
				  (*E_cell[lvl])[dit()](vof, 2)));
	  Egrad = RealVect(D_DECL((*grad_E[lvl])[dit()](vof, 0),
				  (*grad_E[lvl])[dit()](vof, 1),
				  (*grad_E[lvl])[dit()](vof, 2)));


	  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.get_solver();
	    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
	    EBAMRCellData& phi = storage->get_phi();
	    EBAMRCellData& RHS = storage->get_k1();
	    x[idx] = Max(0.0, (*phi[lvl])[dit()](vof, 0));
	    cdr_gradients[idx] = RealVect(D_DECL((*(*grad_cdr[idx])[lvl])[dit()](vof, 0),
						 (*(*grad_cdr[idx])[lvl])[dit()](vof, 1),
						 (*(*grad_cdr[idx])[lvl])[dit()](vof, 2)));
	    rhs[idx] = (*RHS[lvl])[dit()](vof, 0);
	  }

	  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.get_solver();
	    RefCountedPtr<rte_solver>& solver = solver_it();
	    const EBAMRCellData& state = solver->get_state();

	    rte_densities[idx] = Max(0.0, (*state[lvl])[dit()](vof, 0));
	  }

	  this->newton_point_trapz(p, rhs, x, cdr_gradients, E, Egrad, rte_densities, pos, time, a_dt);
	}

	// Irregular cells. These require more since the source term must be extrapolated
	for (VoFIterator vofit(ivs_irr, ebgraph); vofit.ok(); ++vofit){

	}
      }
    }

    if(!m_simpi){ // Update E after Newton iteration
      this->recompute_newton_E();   // Update the electric field
      this->recompute_newton_rte(); // Update the RTE equations
    }
    

    iter++;
    converged = true;
  }

  // Delete grad_cdr. It didn't use smart pointers. 
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    delete grad_cdr[idx];
  }
}

void splitstep_rk2_tga_trapz::allocate_cdr_storage(){
  CH_TIME("splitstep_rk2_tga_trapz::allocate_cdr_storage");
  if(m_verbosity > 2){
    pout() << "splitstep_rk2_tga_trapz::allocate_cdr_storage" << endl;
  }

  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void splitstep_rk2_tga_trapz::allocate_poisson_storage(){
  CH_TIME("splitstep_rk2_tga_trapz::allocate_poisson_storage");
  if(m_verbosity > 2){
    pout() << "splitstep_rk2_tga_trapz::allocate_poisson_storage" << endl;
  }

  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage();
}

void splitstep_rk2_tga_trapz::allocate_rte_storage(){
  CH_TIME("splitstep_rk2_tga_trapz::allocate_rte_storage");
  if(m_verbosity > 2){
    pout() << "splitstep_rk2_tga_trapz::allocate_rte_storage" << endl;
  }

  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void splitstep_rk2_tga_trapz::allocate_sigma_storage(){
  CH_TIME("splitstep_rk2_tga_trapz::allocate_sigma_storage");
  if(m_verbosity > 2){
    pout() << "splitstep_rk2_tga_trapz::allocate_sigma_storage" << endl;
  }

  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void splitstep_rk2_tga_trapz::cache_solutions(){
  CH_TIME("splitstep_rk2_tga_trapz::cache_solutions");
  if(m_verbosity > 2){
    pout() << "splitstep_rk2_tga_trapz::cache_solutions" << endl;
  }

  // Not implemented...
}

void splitstep_rk2_tga_trapz::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  CH_TIME("splitstep_rk2_tga_trapz::compute_dt");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::compute_dt" << endl;
  }

  Real dt = 1.E99;

  const Real dt_cfl = m_cfl*m_cdr->compute_cfl_dt();
  if(dt_cfl < dt){
    dt = dt_cfl;
    a_timecode = time_code::cfl;
  }

  const Real dt_src = m_src_growth*m_cdr->compute_source_dt();
  if(dt_src < dt){
    dt = dt_src;
    a_timecode = time_code::source;
  }

  const Real dt_relax = m_relax_time*this->compute_relaxation_time();
  if(dt_relax < dt){
    dt = dt_relax;
    a_timecode = time_code::relaxation_time;
  }

  const Real dt_restrict = this->restrict_dt();
  if(dt_restrict < dt){
    dt = dt_restrict;
    a_timecode = time_code::restricted;
  }

  if(dt < m_min_dt){
    dt = m_min_dt;
    a_timecode = time_code::hardcap;
  }

  if(dt > m_max_dt){
    dt = m_max_dt;
    a_timecode = time_code::hardcap;
  }

  a_dt = dt;
}

void splitstep_rk2_tga_trapz::compute_cdr_sources_for_simp_newt(){
  CH_TIME("splitstep_rk2_tga_trapz::compute_cdr_sources_for_simp_newt");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::compute_cdr_sources_for_simp_newt" << endl;
  }
  
  // The solvers contain the advected/diffused states
  Vector<EBAMRCellData*> cdr_sources = m_cdr->get_sources();
  Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();
  Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
  EBAMRCellData& E                   = m_poisson_scratch->get_E_cell();

  this->compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, m_time, centering::cell_center);
}

void splitstep_rk2_tga_trapz::compute_trapz_rhs(const Real a_dt){
  CH_TIME("splitstep_rk2_tga_trapz::compute_trapz_rhs");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::compute_trapz_rhs" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    RefCountedPtr<cdr_solver>& solver = solver_it();

    const EBAMRCellData& state  = solver->get_state();
    const EBAMRCellData& source = solver->get_source();


    EBAMRCellData& rhs = storage->get_k1();
    data_ops::copy(rhs, state);
    data_ops::incr(rhs, source, -0.5*a_dt);
    
  }
}

void splitstep_rk2_tga_trapz::compute_epsj(){
  CH_TIME("splitstep_rk2_tga_trapz::compute_epsj");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::compute_epsj" << endl;
  }

  m_epsj = 1.E6;
}

void splitstep_rk2_tga_trapz::compute_E_at_start_of_time_step(){
  CH_TIME("splitstep_rk2_tga_trapz::compute_E_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::compute_E_at_start_of_time_step" << endl;
  }
  
  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();

  const MFAMRCellData& phi = m_poisson->get_state();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void splitstep_rk2_tga_trapz::advance_semi_implicit_newton(const Real a_dt){
  CH_TIME("splitstep_rk2_tga_trapz::advance_semi_implicit_newton");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::advance_semi_implicit_newton" << endl;
  }

  this->compute_E_at_start_of_time_step();   // Compute the electric field using the available cdr/sigma solver states
  this->compute_cdr_sources_for_simp_newt(); // Compute cdr sources for semi-implicit newton advance

  // Explicit Euler advance
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& adv_state       = storage->get_phi();
    const EBAMRCellData& old_state = solver->get_state();
    const EBAMRCellData& source    = solver->get_source();

    data_ops::copy(adv_state, old_state);
    data_ops::incr(adv_state, source, a_dt);

    m_amr->average_down(adv_state, m_cdr->get_phase());

    data_ops::floor(adv_state, 0.0);
  }
}

void splitstep_rk2_tga_trapz::compute_semi_implicit_potential(){
  CH_TIME("splitstep_rk2_tga_trapz::compute_semi_implicit_potential");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::compute_semi_implicit_potential" << endl;
  }

  // TLDR: This routine solves the Poisson equation by using the semi-implicitly advanced cdr states and the sigma solver
    
  Vector<EBAMRCellData*> cdr_densities;
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    cdr_densities.push_back(&(storage->get_phi()));
  }
  

  bool converged = this->solve_poisson(m_poisson->get_state(),
				       m_poisson->get_source(),
				       cdr_densities,
				       m_sigma->get_state(),
				       centering::cell_center);
  if(!converged){
    pout() << "spltistep_rk2_tga_trapz - compute_semi_implicit_potential - solver did not converge at step " << m_step << endl;
  }
}

void splitstep_rk2_tga_trapz::compute_semi_implicit_rte(const Real a_dt){
  CH_TIME("splitstep_rk2_tga_trapz::compute_semi_implicit_rte");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::compute_semi_implicit_rte" << endl;
  }

  const Real time = m_time + a_dt;

  Vector<EBAMRCellData*> rte_states = m_rte->get_states();
  Vector<EBAMRCellData*> rte_sources = m_rte->get_sources();
  Vector<EBAMRCellData*> cdr_states;

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    cdr_states.push_back(&(storage->get_phi()));
  }

  EBAMRCellData& E = m_poisson_scratch->get_E_cell();

  if((m_step + 1) % m_fast_rte == 0){
    const Real dummy_dt = 0.0;
    this->solve_rte(rte_states, rte_sources, cdr_states, E, time, dummy_dt, centering::cell_center);
  }
}

void splitstep_rk2_tga_trapz::deallocate_internals(){
  CH_TIME("splitstep_rk2_tga_trapz::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::deallocate_internals" << endl;
  }
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx]->deallocate_storage();
  }

  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx]->deallocate_storage();
  }

  m_poisson_scratch->deallocate_storage();
  m_sigma_scratch->deallocate_storage();
}

void splitstep_rk2_tga_trapz::newton_point_trapz(Vector<Real>&           a_p,
						 const Vector<Real>&     a_rhs,
						 const Vector<Real>&     a_x,
						 const Vector<RealVect>& a_gradx,
						 const RealVect&         a_E,
						 const RealVect&         a_grad_E,
						 const Vector<Real>&     a_rte,
						 const RealVect&         a_pos,
						 const Real&             a_time,
						 const Real&             a_dt){

  const int num_species = a_p.size();
  const int num_photons = a_rte.size();

  // Compute Fi
  Vector<Real> F(num_species, 0.0);
  Vector<Real> S = m_plaskin->compute_cdr_source_terms(a_time, a_pos, a_E, a_grad_E, a_x, a_rte, a_gradx);
  for (int i = 0; i < num_species; i++){
    F[i] = a_x[i] - 0.5*a_dt*S[i] - a_rhs[i];
  }
  
  Vector<Vector<Real> > jac(num_species);
  for (int i = 0; i < num_species; i++){
    jac[i].resize(num_species);
  }

  // Compute Jacobian
  for (int i = 0; i < num_species; i++){
    const Real Fi = F[i];
    
    for (int j = 0; j < num_species; j++){
      Vector<Real> xj = a_x;
      xj[j] += m_epsj;

      Vector<Real> Sj = m_plaskin->compute_cdr_source_terms(a_time, a_pos, a_E, a_grad_E, xj, a_rte, a_gradx);
      const Real Fii = xj[i] - 0.5*a_dt*Sj[i] - a_rhs[i];
      jac[i][j] = (Fii - Fi)/m_epsj;
    }
  }

  

  MayDay::Abort("splitstep_rk2_tga_trapz::newton_point_trapz - not implemented");
}

void splitstep_rk2_tga_trapz::recompute_newton_E(){
  CH_TIME("splitstep_rk2_tga_trapz::recompute_newton_E");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::recompute_newton_E" << endl;
  }

  MayDay::Abort("splitstep_rk2_tga_trapz::recompute_newton_E - not implemented");
}

void splitstep_rk2_tga_trapz::recompute_newton_rte(){
  CH_TIME("splitstep_rk2_tga_trapz::recompute_newton_rte");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::recompute_newton_rte" << endl;
  }

    MayDay::Abort("splitstep_rk2_tga_trapz::recompute_newton_rte - not implemented");
}

void splitstep_rk2_tga_trapz::regrid_internals(){
  CH_TIME("splitstep_rk2_tga_trapz::regrid_internals");
  if(m_verbosity > 5){
    pout() << "splitstep_rk2_tga_trapz::regrid_internals" << endl;
  }
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();
}
