/*!
  @file   implicit_trapezoidal.cpp
  @brief  Implementation of implicit_trapezoidal.H
  @author Robert Marskar
  @date   May. 2018
*/


#include "implicit_trapezoidal.H"
#include "implicit_trapezoidal_storage.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "data_ops.H"
#include "units.H"

#include <ParmParse.H>

typedef implicit_trapezoidal::cdr_storage     cdr_storage;
typedef implicit_trapezoidal::poisson_storage poisson_storage;
typedef implicit_trapezoidal::rte_storage     rte_storage;
typedef implicit_trapezoidal::sigma_storage   sigma_storage;

implicit_trapezoidal::implicit_trapezoidal(){
  m_tol_x           = 1.E-6;
  m_tol_f           = 1.E-6;
  m_EPS             = 1.E-3;
  m_cfl_redu        = 0.9;
  m_max_newton_iter = 10;
  m_max_pc_iter     = 10;
  
  m_semi_implicit = true;
  m_do_convection = true; // Basically only for debugging
  m_do_diffusion  = true; // Basically only for debugging
  m_do_reaction   = true; // Basically only for debugging
  m_do_rte        = true; // Basically only for debugging
  m_do_poisson    = true; // Basically only for debugging

  m_sequence.resize(3);
  m_sequence[0] = "convection";
  m_sequence[1] = "diffusion";
  m_sequence[2] = "reaction";

  {
    ParmParse pp("implicit_trapezoidal");

    std::string str;
    
    pp.query("solution_tolerance",    m_tol_x);
    pp.query("function_tolerance",    m_tol_f);
    pp.query("max_pc_iter",           m_max_pc_iter);
    pp.query("max_newton_iter",       m_max_newton_iter);
    pp.query("diff_eps",              m_EPS);
    pp.query("cfl_redu",              m_cfl_redu);
    pp.queryarr("sequence",           m_sequence, 0, 3);
    
    if(pp.contains("coupling")){
      pp.query("coupling", str);
      if(str == "semi_implicit"){
	m_semi_implicit = true; 
      }
      else if(str == "implicit"){
	m_semi_implicit = false;
      }
    }

    if(pp.contains("turn_off_advection")){
      pp.get("turn_off_advection", str);
      if(str == "true"){
	m_do_convection = false;
	if(m_verbosity > 2){
	  pout() << "implicit_trapezoidal::implicit_trapezoidal - Turning off advection" << endl;
	}
      }
    }
    
    if(pp.contains("turn_off_diffusion")){
      pp.get("turn_off_diffusion", str);
      if(str == "true"){
	m_do_diffusion = false;
	if(m_verbosity > 2){
	  pout() << "implicit_trapezoidal::implicit_trapezoidal - Turning off diffusion" << endl;
	}
      }
    }
    
    if(pp.contains("turn_off_source")){
      pp.get("turn_off_source", str);
      if(str == "true"){
	m_do_reaction = false;

	if(m_verbosity > 2){
	  pout() << "implicit_trapezoidal::implicit_trapezoidal - Turning off source" << endl;
	}
      }
    }

    if(pp.contains("turn_off_rte")){
      pp.get("turn_off_rte", str);
      if(str == "true"){
	m_do_rte = false;

	if(m_verbosity > 2){
	  pout() << "implicit_trapezoidal::implicit_trapezoidal - Turning off rte" << endl;
	}
      }
    }

    if(pp.contains("turn_off_poisson")){
      pp.get("turn_off_poisson", str);
      if(str == "true"){
	m_do_poisson = false;

	if(m_verbosity > 2){
	  pout() << "implicit_trapezoidal::implicit_trapezoidal - Turning off poisson" << endl;
	}
      }
    }
  }
}

implicit_trapezoidal::~implicit_trapezoidal(){

}

RefCountedPtr<cdr_storage>& implicit_trapezoidal::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& implicit_trapezoidal::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real implicit_trapezoidal::restrict_dt(){
  return 1.E99;
}

Real implicit_trapezoidal::advance(const Real a_dt){
  CH_TIME("implicit_trapezoidal::advance");
  if(m_verbosity > 2){
    pout() << "implicit_trapezoidal::advance" << endl;
  }

  Real dt = a_dt;

  bool converged_convection = false;
  bool converged_diffusion  = true;
  bool converged_reaction   = true;
  
  this->store_states(); // Cache old solutions. Used in cases where the time step is rejected. 
  
  for (int i = 0; i < m_sequence.size(); i++){
    if(m_sequence[i] == "convection" && m_do_convection){
      converged_convection = this->advance_convection(a_dt);
    }
    else if(m_sequence[i] == "diffusion" && m_do_diffusion){
      converged_diffusion = this->advance_diffusion(a_dt);
    }
    else if(m_sequence[i] == "reaction" && m_do_reaction){
      converged_reaction = this->advance_reaction(a_dt);
    }
    else{
      pout() << "implicit_trapezoidal::advance - unknown sequence requested = " << m_sequence[i] << endl;
      MayDay::Abort("implicit_trapezoidal::advance - unknown sequence requested");
    }
  }

  const bool converged = converged_convection && converged_diffusion && converged_reaction;
  
  if(!converged){ // Reject time step, reduce CFL by a little bit. 
    m_cfl = m_cfl*m_cfl_redu;
    this->restore_states();
    dt = 0.0;
  }


  return dt;
}



void implicit_trapezoidal::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void implicit_trapezoidal::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void implicit_trapezoidal::allocate_poisson_storage(){
  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage();
}

void implicit_trapezoidal::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void implicit_trapezoidal::compute_E_into_scratch(const MFAMRCellData& a_phi){
  CH_TIME("splitstep_tga::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();
  
  this->compute_E(E_cell, m_cdr->get_phase(), a_phi);   // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void implicit_trapezoidal::store_states(){
  CH_TIME("implicit_trapezoidal::store_states");
  if(m_verbosity > 2){
    pout() << "implicit_trapezoidal::store_states" << endl;
  }

  {// Cache cdr solutions
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const RefCountedPtr<cdr_solver>& solver = solver_it();

      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
      EBAMRCellData& cache = storage->get_cache();

      data_ops::copy(cache, solver->get_state());
    }
  }
  

  {// Cache RTE solutions
    for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const RefCountedPtr<rte_solver>& solver = solver_it();

      RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
      EBAMRCellData& cache = storage->get_cache();

      data_ops::copy(cache, solver->get_state());
    }
  }
  

  {// Cache Poisson solution
    MFAMRCellData& cache = m_poisson_scratch->get_cache();
    data_ops::copy(cache, m_poisson->get_state());
  }
  

  { // Cache sigma
    EBAMRIVData& cache = m_sigma_scratch->get_cache();
    data_ops::copy(cache, m_sigma->get_state());
  }
  
}

void implicit_trapezoidal::deallocate_internals(){
  CH_TIME("implicit_trapezoidal::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "implicit_trapezoidal::deallocate_internals" << endl;
  }

  m_poisson_scratch->deallocate_storage();
  m_sigma_scratch->deallocate_storage();

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx]->deallocate_storage();
  }

  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx]->deallocate_storage();
  }
}

void implicit_trapezoidal::regrid_internals(){
  CH_TIME("implicit_trapezoidal::regrid_internals");
  if(m_verbosity > 5){
    pout() << "implicit_trapezoidal::regrid_internals" << endl;
  }
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();
}

void implicit_trapezoidal::restore_states(){
  CH_TIME("implicit_trapezoidal::restore_states");
  if(m_verbosity > 2){
    pout() << "implicit_trapezoidal::restore_states" << endl;
  }


  {// Cache cdr solutions
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      const RefCountedPtr<cdr_solver>& solver = solver_it();

      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
      EBAMRCellData& cache = storage->get_cache();

      data_ops::copy(solver->get_state(), cache);
    }
  }
  

  {// Cache RTE solutions
    for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const RefCountedPtr<rte_solver>& solver = solver_it();

      RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
      EBAMRCellData& cache = storage->get_cache();

      data_ops::copy(solver->get_state(), cache);
    }
  }
  

  {// Cache Poisson solution
    MFAMRCellData& cache = m_poisson_scratch->get_cache();
    data_ops::copy(m_poisson->get_state(), cache);
  }
  

  { // Cache sigma
    EBAMRIVData& cache = m_sigma_scratch->get_cache();
    data_ops::copy(m_sigma->get_state(), cache);
  }
}

void implicit_trapezoidal::predictor_convection(const Real a_dt){
  CH_TIME("implicit_trapezoidal::predictor_convection");
  if(m_verbosity > 2){
    pout() << "implicit_trapezoidal::predictor_convection" << endl;
  }

  m_cdr->set_source(0.0);
  m_cdr->set_diffco(0.0);


  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_velocities;
  Vector<EBAMRIVData*>   cdr_fluxes;
  Vector<EBAMRIVData*>   extrap_cdr_states;
  Vector<EBAMRIVData*>   extrap_cdr_velo;
  Vector<EBAMRIVData*>   extrap_cdr_gradients;
  Vector<EBAMRIVData*>   extrap_cdr_fluxes;
  Vector<EBAMRIVData*>   extrap_rte_fluxes;

  cdr_fluxes = m_cdr->get_ebflux();
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver> solver   = solver_it();
    RefCountedPtr<cdr_storage> storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state = solver->get_state();
    EBAMRCellData& velo  = solver->get_velo_cell();

    EBAMRIVData& dens_eb = storage->get_eb_state();
    EBAMRIVData& velo_eb = storage->get_eb_velo();
    EBAMRIVData& flux_eb = storage->get_eb_flux();
    EBAMRIVData& grad_eb = storage->get_eb_grad();

    cdr_states.push_back(&state);
    cdr_velocities.push_back(&velo);

    extrap_cdr_states.push_back(&dens_eb);
    extrap_cdr_velo.push_back(&velo_eb);
    extrap_cdr_gradients.push_back(&flux_eb);
    extrap_cdr_fluxes.push_back(&grad_eb);
  }

  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIVData& flux_eb = storage->get_eb_flux();
    solver->compute_boundary_flux(flux_eb, solver->get_state());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  // Compute cdr fluxes at the boundary
  this->compute_E_into_scratch(m_poisson->get_state()); 
  this->compute_cdr_velocities(cdr_velocities,         cdr_states,         m_poisson_scratch->get_E_cell(), m_time);
  this->extrapolate_to_eb(extrap_cdr_states,           m_cdr->get_phase(), cdr_states);
  this->extrapolate_to_eb(extrap_cdr_velo,             m_cdr->get_phase(), cdr_velocities);
  this->compute_gradients_at_eb(extrap_cdr_gradients,  m_cdr->get_phase(), cdr_states);
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, cdr_states,         cdr_velocities, m_cdr->get_phase());
  this->compute_cdr_fluxes(cdr_fluxes, 
			   extrap_cdr_fluxes,
			   extrap_cdr_states,
			   extrap_cdr_velo,
			   extrap_cdr_gradients,
			   extrap_rte_fluxes,
			   m_poisson_scratch->get_E_eb(),
			   m_time);
  this->compute_charge_flux(m_sigma->get_flux(), cdr_fluxes);

  // Euler advance. Store -0.5*dt*div(nv)^k in scratch storage. Solver contains explicitly advanced state
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& phi_old = storage->get_cache();
    EBAMRCellData& rhs     = storage->get_scratch1();
    EBAMRCellData& phi_new = solver->get_state();

    solver->compute_divF(rhs, phi_old, 0.0, true);
    data_ops::scale(rhs, -0.5*a_dt);
    data_ops::incr(phi_new, rhs, 2.0);

    m_amr->average_down(phi_new, m_cdr->get_phase());
    m_amr->interp_ghost(phi_new, m_cdr->get_phase());

    data_ops::floor(phi_new, 0.0);
  }


  // Euler advance for sigma
  MayDay::Abort("implicit_trapezoidal::predictor_convection - missing Euler advance for sigma");
}

bool implicit_trapezoidal::advance_convection(const Real a_dt){
  CH_TIME("implicit_trapezoidal::advance_convection");
  if(m_verbosity > 2){
    pout() << "implicit_trapezoidal::advance_convection" << endl;
  }

  bool converged = true;

  this->predictor_convection(a_dt);

  return false;
}

bool implicit_trapezoidal::advance_diffusion(const Real a_dt){
  CH_TIME("implicit_trapezoidal::advance_diffusion");
  if(m_verbosity > 2){
    pout() << "implicit_trapezoidal::advance_diffusion" << endl;
  }

  return false;
}

bool implicit_trapezoidal::advance_reaction(const Real a_dt){
  CH_TIME("implicit_trapezoidal::advance_reaction");
  if(m_verbosity > 2){
    pout() << "implicit_trapezoidal::advance_reaction" << endl;
  }

  return false;
}
