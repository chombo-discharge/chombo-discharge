/*!
  @file   rk2_stiff.cpp
  @brief  Implementation of rk2_stiff.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "rk2_stiff.H"
#include "rk2_stiff_storage.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "data_ops.H"
#include "units.H"

#include <EBArith.H>
#include <ParmParse.H>

#define debug_source 0

extern "C" void FORT_SOLVE_LU(int* N, double* J, double* F, int* INFO);

typedef rk2_stiff::cdr_storage     cdr_storage;
typedef rk2_stiff::poisson_storage poisson_storage;
typedef rk2_stiff::rte_storage     rte_storage;
typedef rk2_stiff::sigma_storage   sigma_storage;

rk2_stiff::rk2_stiff(){
  m_alpha    = 1.0;
  m_tol_x    = 1.E-8;
  m_tol_f    = 1.E-8;
  m_EPS      = 1.E-3;
  m_max_iter = 10;
  m_simpi    = true;

  // Basically only for debugging
  m_do_adv_diff  = true;
  m_do_source    = true;
  m_do_rte       = true;
  m_do_poisson   = true;
  
  {
    ParmParse pp("rk2_stiff");

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
	MayDay::Abort("rk2_stiff::rk2_stiff - fully implicit not yet supported");
      }
    }

    if(pp.contains("turn_off_adv_diff")){
      pp.get("turn_off_adv_diff", str);
      if(str == "true"){
	m_do_adv_diff = false;
	if(m_verbosity > 2){
	  pout() << "rk2_stiff::rk2_stiff - Turning off advection-diffusion" << endl;
	}
      }
    }
    
    if(pp.contains("turn_off_source")){
      pp.get("turn_off_source", str);
      if(str == "true"){
	m_do_source = false;

	if(m_verbosity > 2){
	  pout() << "rk2_stiff::rk2_stiff - Turning off source" << endl;
	}
      }
    }

    if(pp.contains("turn_off_rte")){
      pp.get("turn_off_rte", str);
      if(str == "true"){
	m_do_rte = false;

	if(m_verbosity > 2){
	  pout() << "rk2_stiff::rk2_stiff - Turning off rte" << endl;
	}
      }
    }

    if(pp.contains("turn_off_poisson")){
      pp.get("turn_off_poisson", str);
      if(str == "true"){
	m_do_poisson = false;

	if(m_verbosity > 2){
	  pout() << "rk2_stiff::rk2_stiff - Turning off poisson" << endl;
	}
      }
    }
  }
}

rk2_stiff::~rk2_stiff(){

}

RefCountedPtr<cdr_storage>& rk2_stiff::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& rk2_stiff::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real rk2_stiff::restrict_dt(){
  return 1.E99;
}

Real rk2_stiff::advance(const Real a_dt){
  CH_TIME("rk2_stiff::advance");
  if(m_verbosity > 2){
    pout() << "rk2_stiff::advance" << endl;
  }

  this->cache_solutions(); // Cache old solutions. Used in case the time step is rejected. 

  if(m_do_adv_diff) {                         // Rules for advective advance: 
    this->advance_advection_diffusion(a_dt);  // states, and the poisson solver contains the updated potential after advection.
  }                                           // The sigma solver is also updated, while scratch storage contains junk.

  bool converged_source = true;
  if(m_do_source){
    converged_source = this->advance_sources(a_dt);   // Source term advance. 
  }

  // Put solver back in useable state so that we can reliably compute the next time step. 
  this->compute_cdr_velocities();
  this->compute_cdr_diffusion();
  this->compute_cdr_sources();
  
  return a_dt;
}

bool rk2_stiff::advance_sources(const Real a_dt){
  CH_TIME("rk2_stiff::advance_sources");
  if(m_verbosity > 2){
    pout() << "rk2_stiff::advance_sources" << endl;
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

  // Prepare for Newton iteration
  Vector<EBAMRCellData*> iterates;
  this->compute_dnj();                        // Get tolerances for the finite difference evaluation 
  this->setup_newton_iterates(iterates);       // Set up target iterates. This iterates directly in the solver
  this->compute_E_into_scratch();              // Compute electric field at time step k
  this->compute_cdr_sources_for_newton_pred(); // Compute source terms at time step k
  this->compute_trapz_rhs(a_dt);               // Compute the right-hand side for the trapezoid rule = n_k + 0.5*dt*S_k
  this->explicit_euler_predict_newton(a_dt);   // Explicit Euler advance as initial guess

  if(m_do_poisson){
    this->solve_poisson();          // Solve Poisson equation with explicit Euler data. Updated potential lies in solver.
    this->compute_E_into_scratch(); // Recompute field by using the semi-implicit Poisson advance
  }
  if(m_do_rte){
    const Real dummy_dt = 0.0;  
    this->solve_rte(dummy_dt); // Update the RTE equations. 
  }


  int iter = 0;           // Number of iterations
  bool converged = false; // Convergence track
  Real preF = 0.0;

  while(!converged && iter < m_max_iter){

    // Compute grad(|E|). This is needed for source term computation. 
    const EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
    data_ops::vector_length(E_norm, E_cell);          // Compute |E|
    m_amr->average_down(E_norm, m_cdr->get_phase());  // Average down
    m_amr->interp_ghost(E_norm, m_cdr->get_phase());  // Interpolate ghost cells
    m_amr->compute_gradient(grad_E, E_norm);          // Compute gradient
    m_amr->average_down(grad_E, m_cdr->get_phase());  // Average down gradient
    m_amr->interp_ghost(grad_E, m_cdr->get_phase());  // Interpolate gradient ghost cells

    // Compute grad(n). This is needed for source term computation. 
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      const RefCountedPtr<cdr_solver>& solver = solver_it();
      const EBAMRCellData& iter_state = solver->get_state();
      
      m_amr->compute_gradient(*grad_cdr[idx], iter_state);     // Compute grad()
      m_amr->average_down(*grad_cdr[idx], m_cdr->get_phase()); // Average down
      m_amr->interp_ghost(*grad_cdr[idx], m_cdr->get_phase()); // Interpolate ghost cells
    }

    converged = true; // Set to false if point ODEs don't converge.
    Real maxF = 0.0;  // Error bound
    Real maxX = 0.0;  // Error bound

    // Extrapolation stencils
    const irreg_amr_stencil<centroid_interp>& interp_stencils = m_amr->get_centroid_interp_stencils(m_cdr->get_phase());

    // Level & grid loops
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

	// Irregular cells. These require more since the source term must be extrapolated. We MUST do these
	// before regular cells because if we overwrite the regular cells the extrapolations become fuzzy
#if debug_source
	for (VoFIterator vofit(ivs_irr, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	  const VoFStencil& stencil = interp_stencils[lvl][dit()](vof, 0);
	  
	  pos   = EBArith::getVofLocation(vof, dx*RealVect::Unit, m_physdom->get_prob_lo());

	  E     = RealVect::Zero;
	  Egrad = RealVect::Zero;
	  for (int i = 0; i < stencil.size(); i++){
	    const VolIndex& ivof = stencil.vof(i);
	    const Real& iweight  = stencil.weight(i);

	    for (int dir = 0; dir < SpaceDim; dir++){
	      E[dir] += (*E_cell[lvl])[dit()](ivof, dir)*iweight;
	      Egrad[dir] += (*grad_E[lvl])[dit()](ivof, dir)*iweight;
	    }
	  }

	  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.get_solver();
	    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
	    EBAMRCellData& RHS = storage->get_k1();

	    Real phi = 0.0;
	    Real cur_rhs = 0.0;
	    RealVect grad = RealVect::Zero;


	    for (int i = 0; i < stencil.size(); i++){
	      const VolIndex& ivof = stencil.vof(i);
	      const Real& iweight  = stencil.weight(i);

	      phi += (*(*iterates[idx])[lvl])[dit()](ivof, 0)*iweight;

	      for (int dir = 0; dir < SpaceDim; dir++){
		grad[dir] += (*(*grad_cdr[idx])[lvl])[dit()](ivof, 0)*iweight;
	      }
	    }

	    x[idx] = (*(*iterates[idx])[lvl])[dit()](vof, 0);
	    cdr_gradients[idx] = grad;
	    rhs[idx] = (*RHS[lvl])[dit()](vof, 0);
	  }

	  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.get_solver();
	    RefCountedPtr<rte_solver>& solver = solver_it();
	    const EBAMRCellData& state = solver->get_state();

	    Real phi = 0.0;
	    for (int i = 0; i < stencil.size(); i++){
	      const VolIndex& ivof = stencil.vof(i);
	      const Real& iweight  = stencil.weight(i);

	      phi += (*state[lvl])[dit()](ivof, 0)*iweight;
	    }

	    rte_densities[idx] = Max(phi, 0.0);
	  }

	  // Newton solve for correction => p
	  const Real sumF = this->newton_point_trapz(p, rhs, x, cdr_gradients, E, Egrad, rte_densities, pos, time, a_dt);
	  
	  // Increment data. 
	  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.get_solver();
	    (*(*iterates[idx])[lvl])[dit()](vof, 0) += p[idx];
	  }

	  // Errors
	  Real sumX = 0.0;
	  for (int i = 0; i < p.size(); i++){
	    sumX += Abs(p[i]);
	  }
	  maxF = Max(sumF, maxF);
	  maxX = Max(maxX, sumX);


	}

#endif

	// Regular cells. Straightforward stuff
	for (VoFIterator vofit(ivs_reg, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();

	  // Position and electric field
	  pos   = EBArith::getVofLocation(vof, dx*RealVect::Unit, m_physdom->get_prob_lo());
	  E     = RealVect(D_DECL((*E_cell[lvl])[dit()](vof, 0),
				  (*E_cell[lvl])[dit()](vof, 1),
				  (*E_cell[lvl])[dit()](vof, 2)));
	  Egrad = RealVect(D_DECL((*grad_E[lvl])[dit()](vof, 0),
				  (*grad_E[lvl])[dit()](vof, 1),
				  (*grad_E[lvl])[dit()](vof, 2)));


	  // Get previous iterate, the gradients and the right-hand side of the trapezoidal equation in the current cell
	  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.get_solver();
	    RefCountedPtr<cdr_solver>& solver   = solver_it();
	    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
	    EBAMRCellData& phi = solver->get_state();
	    EBAMRCellData& RHS = storage->get_k1();
	    x[idx] = Max(0.0, (*phi[lvl])[dit()](vof, 0));
	    cdr_gradients[idx] = RealVect(D_DECL((*(*grad_cdr[idx])[lvl])[dit()](vof, 0),
						 (*(*grad_cdr[idx])[lvl])[dit()](vof, 1),
						 (*(*grad_cdr[idx])[lvl])[dit()](vof, 2)));
	    rhs[idx] = (*RHS[lvl])[dit()](vof, 0);
	  }

	  // Get isotropic photon densities
	  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.get_solver();
	    RefCountedPtr<rte_solver>& solver = solver_it();
	    const EBAMRCellData& state = solver->get_state();

	    rte_densities[idx] = Max(0.0, (*state[lvl])[dit()](vof, 0));
	  }

	  // Newton solve for correction => p
	  const Real sumF = this->newton_point_trapz(p, rhs, x, cdr_gradients, E, Egrad, rte_densities, pos, time, a_dt);
	  
	  // Increment data. 
	  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	    const int idx = solver_it.get_solver();
	    (*(*iterates[idx])[lvl])[dit()](vof, 0) += p[idx];
	  }

	  // Errors
	  Real sumX = 0.0;
	  for (int i = 0; i < p.size(); i++){
	    sumX += Abs(p[i]);
	  }
	  maxF = Max(sumF, maxF);
	  maxX = Max(maxX, sumX);
	}
      }
    }

    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      m_amr->interp_ghost(*iterates[idx], m_cdr->get_phase());
    }

#ifdef CH_MPI
    Real lmaxF = maxF;
    Real lmaxX = maxX;
    int result  = MPI_Allreduce(&lmaxF, &maxF, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
    int result2 = MPI_Allreduce(&lmaxX, &maxX, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
    if(result != MPI_SUCCESS || result2 != MPI_SUCCESS){
      MayDay::Error("rk2_stiff::advance_sources - communication error in advance_sources");
    }
#endif

    converged = (maxF < m_tol_f*m_nmax) || (maxX < m_tol_x*m_nmax);

    if( iter == 0){
      preF = maxF;
    }
    if(m_verbosity > 0){
      pout() << "rk2_stiff::advance_sources - Newton iteration " << iter
	     << " \t Converged = " << converged
	     << " \t Error_F = " << maxF
	     << " \t Error_x = " << maxX
	     << " \t Rate = " << preF/maxF
	     << endl;
    }

    // Prepare for next iteration
    preF = maxF;
    iter = iter + 1;

    if(!m_simpi && !converged){ // Update E after Newton iteration
      this->recompute_newton_E();   // Update the electric field
      this->recompute_newton_rte(); // Update the RTE equations
    }
  }

  // Delete grad_cdr. It didn't use smart pointers. 
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    delete grad_cdr[idx];
  }

  return converged;
}

void rk2_stiff::advance_advection_diffusion(const Real a_dt){
  CH_TIME("rk2_stiff::advance_advection_diffusion");
  if(m_verbosity > 2){
    pout() << "rk2_stiff::advance_advection_diffusion" << endl;
  }

  this->set_cdr_sources_to_zero();              // Set cdr sources to zero

  // Prepare for k1 advance
  this->compute_E_into_scratch();   // Compute E into scratch storage by using m_poisson->m_state
  this->compute_cdr_velo_using_solver_states(m_time);   // Compute cdr velocities using what is available in solvers and E_scratch
  this->compute_cdr_eb_states_using_solver_states();    // Compute cdr EB states using what is available in solvers and E_scratch
  this->compute_cdr_diffco_using_solver_states(m_time); // Compute cdr diffusion coefficients using what is available in solvers..
  this->compute_cdr_fluxes_using_solver_states(m_time); // Compute cdr fluxes. Call plasma_kinetics and update boundary conditions
  this->compute_sigma_flux_using_solver_states();       // Compute sigma flux. Does the sum of cdr_fluxes

  // Do the k1-advance into scratch storages
  this->advance_cdr_k1(a_dt);     // First Runge Kutta stage. This computes k1 into scratch and a temporary state into the solver
  this->advance_sigma_k1(a_dt);   // First Runge Kutta stage. This computes k1 into scratch and a temporary state into the solver
  this->solve_poisson();          // Solve Poisson equation using the intermediately advanced states (reside in solvers)
  this->compute_E_into_scratch(); // Compute E using the potential that we got. Put the result in scratch. 
  if(m_do_rte){
    if(m_rte->is_stationary()){
      const Real dummy = 0.0;
      this->solve_rte(dummy);      // Compute RTE equations after the k1 advance
    }
    else{
      MayDay::Abort("rk2_stiff::advance_advection_diffusion - transient RTE is not supported for this time_stepper");
    }
  }

  // Do the same shit all over again, but solvers now contain the intermediately advanced states. Also, we
  // need to take the remainder of the time step. 
  const Real time_k1 = m_time + m_alpha*a_dt;
  this->compute_cdr_velo_using_solver_states(time_k1);
  this->compute_cdr_diffco_using_solver_states(time_k1);
  this->compute_cdr_eb_states_using_solver_states();
  this->compute_cdr_fluxes_using_solver_states(time_k1);
  this->compute_sigma_flux_using_solver_states();
  
  // Do the k2 advance back into solver states. 
  this->advance_cdr_k2(a_dt);
  this->advance_sigma_k2(a_dt);
  this->solve_poisson();
  this->compute_E_into_scratch();
  if(m_do_rte){
    if(m_rte->is_stationary()){
      const Real dummy = 0.0;
      this->solve_rte(dummy); // Compute RTE equations after the k2 advance
    }
    else{
      MayDay::Abort("rk2_stiff::advance_advection_diffusion - transient RTE is not supported for this time_stepper");
    }
  }
}

void rk2_stiff::set_cdr_sources_to_zero(){
  CH_TIME("rk2_stiff::set_cdr_sources_to_zero");
  if(m_verbosity > 2){
    pout() << "rk2_stiff::set_cdr_sources_to_zero" << endl;
  }
  
  m_cdr->set_source(0.0);
}

void rk2_stiff::compute_cdr_velo_using_solver_states(const Real a_time){
  CH_TIME("rk2_stiff::compute_cdr_velo_using_solver_states");
  if(m_verbosity > 2){
    pout() << "rk2_stiff::compute_cdr_velo_using_solver_states" << endl;
  }

  Vector<EBAMRCellData*> states     = m_cdr->get_states();
  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  this->compute_cdr_velocities(velocities, states, m_poisson_scratch->get_E_cell(), a_time);
}

void rk2_stiff::compute_cdr_eb_states_using_solver_states(){
  CH_TIME("rk2_stiff::compute_cdr_eb_states_using_solver_states");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::compute_cdr_eb_states_using_solver_states" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_states;
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);


    cdr_states.push_back(&(solver->get_state()));
    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
  }

  this->extrapolate_to_eb(eb_states,          m_cdr->get_phase(), cdr_states);
  this->compute_gradients_at_eb(eb_gradients, m_cdr->get_phase(), cdr_states);

}

void rk2_stiff::compute_cdr_diffco_using_solver_states(const Real a_time){
  CH_TIME("rk2::compute_cdr_diffco_using_solver_states");
  if(m_verbosity > 5){
    pout() << "rk2::compute_cdr_diffco_using_solver_states" << endl;
  }

  const int num_species = m_plaskin->get_num_species();

  Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();
  Vector<EBAMRFluxData*> diffco_face = m_cdr->get_diffco_face();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->get_diffco_eb();

  const EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  const EBAMRIVData& E_eb     = m_poisson_scratch->get_E_eb();

  // Get extrapolated states
  Vector<EBAMRIVData*> eb_states(num_species);
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    eb_states[idx] = &(storage->get_eb_state());
  }
  
  this->compute_cdr_diffco_face(diffco_face, cdr_states, E_cell, a_time);
  this->compute_cdr_diffco_eb(diffco_eb,     eb_states,  E_eb,   a_time);
}

void rk2_stiff::compute_cdr_fluxes_using_solver_states(const Real a_time){
    CH_TIME("rk2::compute_cdr_fluxes_using_solver_states");
  if(m_verbosity > 5){
    pout() << "rk2::compute_cdr_fluxes_using_solver_states" << endl;
  }
  
  Vector<EBAMRIVData*> cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_densities;
  Vector<EBAMRIVData*> extrap_cdr_velocities;
  Vector<EBAMRIVData*> extrap_cdr_gradients;
  Vector<EBAMRIVData*> extrap_rte_fluxes;

  cdr_fluxes = m_cdr->get_ebflux();

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRIVData& dens_eb = storage->get_eb_state();
    EBAMRIVData& velo_eb = storage->get_eb_velo();
    EBAMRIVData& flux_eb = storage->get_eb_flux();
    EBAMRIVData& grad_eb = storage->get_eb_grad();

    extrap_cdr_densities.push_back(&dens_eb);  // Already been computed
    extrap_cdr_velocities.push_back(&velo_eb);
    extrap_cdr_fluxes.push_back(&flux_eb);
    extrap_cdr_gradients.push_back(&grad_eb);  // Already been computed
  }

  // Extrapolate densities, velocities, and fluxes
  Vector<EBAMRCellData*> cdr_densities = m_cdr->get_states();
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->get_velocities();
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, cdr_densities, cdr_velocities, m_cdr->get_phase());
  this->extrapolate_to_eb(extrap_cdr_velocities, m_cdr->get_phase(), cdr_velocities);
  //  this->extrapolate_to_eb(extrap_cdr_densities,  m_cdr->get_phase(), cdr_densities); // Already been done, no?

  // Compute RTE flux on the boundary
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIVData& flux_eb = storage->get_eb_flux();
    solver->compute_boundary_flux(flux_eb, solver->get_state());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_poisson_scratch->get_E_eb();

  this->compute_cdr_fluxes(cdr_fluxes,
			   extrap_cdr_fluxes,
			   extrap_cdr_densities,
			   extrap_cdr_velocities,
			   extrap_cdr_gradients,
			   extrap_rte_fluxes,
			   E,
			   a_time);
}

void rk2_stiff::compute_sigma_flux_using_solver_states(){
  CH_TIME("rk2::compute_sigma_flux_using_solver_states");
  if(m_verbosity > 5){
    pout() << "rk2::compute_sigma_flux_using_solver_states" << endl;
  }

  EBAMRIVData& flux = m_sigma->get_flux();
  data_ops::set_value(flux, 0.0);

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const RefCountedPtr<species>& spec      = solver_it.get_species();
    const EBAMRIVData& solver_flux          = solver->get_ebflux();

    data_ops::incr(flux, solver_flux, spec->get_charge()*units::s_Qe);
  }

  m_sigma->reset_cells(flux);
}

void rk2_stiff::advance_cdr_k1(const Real a_dt){
  CH_TIME("rk2_stiff::advance_cdr_k1");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::advance_cdr_k1" << endl;
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& k1          = storage->get_k1();
    EBAMRCellData& phi         = solver->get_state();
    const EBAMRCellData& state = storage->get_cache();

    data_ops::set_value(k1, 0.0);
    solver->compute_rhs(k1, state, a_dt);

    data_ops::set_value(phi, 0.0);
    data_ops::incr(phi, state, 1.0);
    data_ops::incr(phi, k1,    m_alpha*a_dt);

    m_amr->average_down(phi, m_cdr->get_phase());
    m_amr->interp_ghost(phi, m_cdr->get_phase());

    data_ops::floor(phi, 0.0);
  }
}

void rk2_stiff::advance_sigma_k1(const Real a_dt){
    CH_TIME("rk2_stiff::advance_sigma_k1");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::advance_sigma_k1" << endl;
  }

  const EBAMRIVData& state = m_sigma_scratch->get_cache();
  
  EBAMRIVData& k1  = m_sigma_scratch->get_k1();
  EBAMRIVData& phi = m_sigma->get_state();
  m_sigma->compute_rhs(k1);
  data_ops::set_value(phi,   0.0);
  data_ops::incr(phi, state, 1.0);
  data_ops::incr(phi, k1,    m_alpha*a_dt);

  m_amr->average_down(phi, m_cdr->get_phase());
  
  m_sigma->reset_cells(k1);
  m_sigma->reset_cells(phi);
}

void rk2_stiff::advance_cdr_k2(const Real a_dt){
  CH_TIME("rk2::advance_cdr_k2");
  if(m_verbosity > 5){
    pout() << "rk2::advance_cdr_k2" << endl;
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state   = storage->get_cache();    // Old solution
    EBAMRCellData& k1      = storage->get_k1();       // k1
    EBAMRCellData& k2      = storage->get_k2();       // k2
    EBAMRCellData& phi     = solver->get_state();     // Intermediate state, stored in solver

    solver->compute_rhs(k2, phi, a_dt);               // Compute k2

    // Copy old solution into solver state and do the Runge-Kutta advance
    data_ops::set_value(phi, 1.0);
    data_ops::incr(phi, state, 1.0);                  
    data_ops::incr(state, k1, a_dt*(1 - 1./(2.*m_alpha))); 
    data_ops::incr(state, k2, a_dt*1./(2.*m_alpha));       

    m_amr->average_down(state, m_cdr->get_phase());
    m_amr->interp_ghost(state, m_cdr->get_phase());

    data_ops::floor(phi, 0.0);
  }
}

void rk2_stiff::advance_sigma_k2(const Real a_dt){
  CH_TIME("rk2_stiff::advance_sigma_k2");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::advance_sigma_k2" << endl;
  }
  
  EBAMRIVData& k1  = m_sigma_scratch->get_k1();
  EBAMRIVData& k2  = m_sigma_scratch->get_k2();
  m_sigma->compute_rhs(k2);

  EBAMRIVData& phi   = m_sigma->get_state();
  EBAMRIVData& state = m_sigma_scratch->get_cache();
  data_ops::set_value(phi, 0.0);
  data_ops::incr(phi, state, 1.0);
  data_ops::incr(phi, k1, a_dt*(1 - 1./(2.*m_alpha)));
  data_ops::incr(phi, k2, a_dt*1./(2.*m_alpha));

  m_sigma->reset_cells(phi);
}

void rk2_stiff::allocate_cdr_storage(){
  CH_TIME("rk2_stiff::allocate_cdr_storage");
  if(m_verbosity > 2){
    pout() << "rk2_stiff::allocate_cdr_storage" << endl;
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

void rk2_stiff::allocate_poisson_storage(){
  CH_TIME("rk2_stiff::allocate_poisson_storage");
  if(m_verbosity > 2){
    pout() << "rk2_stiff::allocate_poisson_storage" << endl;
  }

  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage();
}

void rk2_stiff::allocate_rte_storage(){
  CH_TIME("rk2_stiff::allocate_rte_storage");
  if(m_verbosity > 2){
    pout() << "rk2_stiff::allocate_rte_storage" << endl;
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

void rk2_stiff::allocate_sigma_storage(){
  CH_TIME("rk2_stiff::allocate_sigma_storage");
  if(m_verbosity > 2){
    pout() << "rk2_stiff::allocate_sigma_storage" << endl;
  }

  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void rk2_stiff::cache_solutions(){
  CH_TIME("rk2_stiff::cache_solutions");
  if(m_verbosity > 2){
    pout() << "rk2_stiff::cache_solutions" << endl;
  }
  return;
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    const EBAMRCellData& state = solver->get_state();
    EBAMRCellData cache = storage->get_cache();

    data_ops::set_value(cache, 0.0);
    data_ops::incr(cache, state, 1.0);
  }

  //  MayDay::Abort("rk2_stiff::cache solutions - really need to copy everything over into the respective caches...");
}

void rk2_stiff::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  CH_TIME("rk2_stiff::compute_dt");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::compute_dt" << endl;
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

void rk2_stiff::compute_cdr_sources_for_newton_pred(){
  CH_TIME("rk2_stiff::compute_cdr_sources_for_newton_pred");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::compute_cdr_sources_for_newton_pred" << endl;
  }
  
  // The solvers contain the advected/diffused states
  Vector<EBAMRCellData*> cdr_sources = m_cdr->get_sources();
  Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();
  Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
  EBAMRCellData& E                   = m_poisson_scratch->get_E_cell();

  this->compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, m_time, centering::cell_center);
}

void rk2_stiff::compute_trapz_rhs(const Real a_dt){
  CH_TIME("rk2_stiff::compute_trapz_rhs");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::compute_trapz_rhs" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    RefCountedPtr<cdr_solver>& solver = solver_it();

    const EBAMRCellData& state  = solver->get_state();
    const EBAMRCellData& source = solver->get_source();


    EBAMRCellData& rhs = storage->get_k1();
    data_ops::copy(rhs, state);
    data_ops::incr(rhs, source, 0.5*a_dt);

    m_amr->average_down(rhs, m_cdr->get_phase());
    m_amr->interp_ghost(rhs, m_cdr->get_phase());
  }
}

void rk2_stiff::compute_dnj(){
  CH_TIME("rk2_stiff::compute_dnj");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::compute_dnj" << endl;
  }

  int comp = 0;

  m_nmax = 0.0;
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    EBAMRCellData& state = solver_it()->get_state();

    Real cur_max;
    Real cur_min;
    
    data_ops::get_max_min(cur_max, cur_min, state, comp);

    m_nmax = Max(cur_max, m_nmax);
  }

  m_dnj = m_EPS*m_nmax;
  if(m_nmax == 0.0) m_dnj = m_EPS;

}

void rk2_stiff::compute_E_into_scratch(){
  CH_TIME("rk2_stiff::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();

  const MFAMRCellData& phi = m_poisson->get_state();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void rk2_stiff::deallocate_internals(){
  CH_TIME("rk2_stiff::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::deallocate_internals" << endl;
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

void rk2_stiff::explicit_euler_predict_newton(const Real a_dt){
  CH_TIME("rk2_stiff::explicit_euler_predict_newton");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::explicit_euler_predict_newton" << endl;
  }

  // Explicit Euler advance over dt
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();

    EBAMRCellData& state        = solver->get_state();
    const EBAMRCellData& source = solver->get_source();

    data_ops::incr(state, source, a_dt);

    m_amr->average_down(state, m_cdr->get_phase());
    m_amr->interp_ghost(state, m_cdr->get_phase());

    data_ops::floor(state, 0.0);
  }
}

Real rk2_stiff::newton_point_trapz(Vector<Real>&           a_p,
				   const Vector<Real>&     a_rhs,
				   const Vector<Real>&     a_x,
				   const Vector<RealVect>& a_gradx,
				   const RealVect&         a_E,
				   const RealVect&         a_grad_E,
				   const Vector<Real>&     a_rte,
				   const RealVect&         a_pos,
				   const Real&             a_time,
				   const Real&             a_dt){

  int N = a_p.size();

  // Compute Fi
  Real sumF = 0.0;
  double F[N];
  Vector<Real> S = m_plaskin->compute_cdr_source_terms(a_time, a_pos, a_E, a_grad_E, a_x, a_rte, a_gradx);
  for (int i = 0; i < N; i++){
    F[i] = a_x[i] - 0.5*a_dt*S[i] - a_rhs[i];
    sumF += Abs(F[i]);
  }

  // Compute Jacobian. This must be done in Fortran major order
  double J[N*N];
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      const int n = i + j*N;
      
      Vector<Real> xp = a_x;
      Vector<Real> xm = a_x;
      xp[j] += m_dnj;
      xm[j] -= m_dnj;
      Vector<Real> Sp = m_plaskin->compute_cdr_source_terms(a_time, a_pos, a_E, a_grad_E, xp, a_rte, a_gradx);
      //      Vector<Real> Sm = m_plaskin->compute_cdr_source_terms(a_time, a_pos, a_E, a_grad_E, xm, a_rte, a_gradx);

      const Real Fp = xp[i] - 0.5*a_dt*Sp[i] - a_rhs[i];
      //      const Real Fm = xm[i] - 0.5*a_dt*Sm[i] - a_rhs[i];

      
      J[n] = (Fp - F[i])/(m_dnj); // Centered difference Jacobian
    }
  }

  int info;
  FORT_SOLVE_LU(&N, J, F, &info);
  if(info != 0){
    MayDay::Abort("splitstep_rks_tga_trapz::newton_point_trapz - LU decomposition solver failed");
  }

  for (int i = 0; i < N; i++){
    a_p[i] = F[i];
  }

  return sumF;
}

void rk2_stiff::recompute_newton_E(){
  CH_TIME("rk2_stiff::recompute_newton_E");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::recompute_newton_E" << endl;
  }

  MayDay::Abort("rk2_stiff::recompute_newton_E - not implemented");
}

void rk2_stiff::recompute_newton_rte(){
  CH_TIME("rk2_stiff::recompute_newton_rte");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::recompute_newton_rte" << endl;
  }

  MayDay::Abort("rk2_stiff::recompute_newton_rte - not implemented");
}

void rk2_stiff::regrid_internals(){
  CH_TIME("rk2_stiff::regrid_internals");
  if(m_verbosity > 5){
    pout() << "rk2_stiff::regrid_internals" << endl;
  }
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();
}

void rk2_stiff::setup_newton_iterates(Vector<EBAMRCellData*>& a_iterates){
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    a_iterates.push_back(&(solver->get_state()));
  }
}
