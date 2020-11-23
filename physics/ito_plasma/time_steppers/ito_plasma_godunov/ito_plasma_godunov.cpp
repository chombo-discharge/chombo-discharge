/*!
  @file   ito_plasma_godunov.cpp
  @author Robert Marskar
  @date   June 2020
  @brief  Implementation of ito_plasma_godunov
*/

#include "ito_plasma_godunov.H"
#include "data_ops.H"
#include "units.H"
#include "poisson_multifluid_gmg.H"

#include <ParmParse.H>

using namespace physics::ito_plasma;

ito_plasma_godunov::ito_plasma_godunov(RefCountedPtr<ito_plasma_physics>& a_physics){
  m_name    = "ito_plasma_godunov";
  m_physics = a_physics;

  m_dt_relax = 1.E99;

  ParmParse pp("ito_plasma_godunov");
  pp.get("particle_realm", m_particle_realm);
  pp.get("profile", m_profile);
  pp.get("load_ppc", m_load_ppc);

  m_avg_cfl = 0.0;


}

ito_plasma_godunov::~ito_plasma_godunov(){

}

int ito_plasma_godunov::get_num_plot_vars() const {
  CH_TIME("ito_plasma_godunov::get_num_plot_vars");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::get_num_plot_vars" << endl;
  }

  int ncomp = ito_plasma_stepper::get_num_plot_vars();

  ncomp++; // Add conductivity

  return ncomp;
}

void ito_plasma_godunov::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {
  CH_TIME("ito_plasma_godunov::write_conductivity");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::write_conductivity" << endl;
  }

  ito_plasma_stepper::write_plot_data(a_output, a_plotvar_names, a_icomp);

  // Do conductivity
  this->write_conductivity(a_output, a_icomp);
  a_plotvar_names.push_back("conductivity");
}

void ito_plasma_godunov::write_conductivity(EBAMRCellData& a_output, int& a_icomp) const {
  CH_TIME("ito_plasma_stepper::write_conductivity");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_conductivity" << endl;
  }

  const Interval src_interv(0, 0);
  const Interval dst_interv(a_icomp, a_icomp);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    if(m_conduct_cell.get_realm() == a_output.get_realm()){
      m_conduct_cell[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else {
      m_conduct_cell[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }
  a_icomp += 1;
}

void ito_plasma_godunov::allocate(){
  CH_TIME("ito_plasma_godunov::allocate");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::allocate" << endl;
  }

  m_ito->allocate_internals();
  m_rte->allocate_internals();
  m_poisson->allocate_internals();
  m_sigma->allocate_internals();

  // Now allocate for the conductivity particles and rho^dagger particles
  const int num_ito_species = m_physics->get_num_ito_species();
  
  m_conductivity_particles.resize(num_ito_species);
  m_rho_dagger_particles.resize(num_ito_species);
  
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver = solver_it();

    const int idx         = solver_it.get_solver();
    const int pvr_buffer  = solver->get_pvr_buffer();
    const int halo_buffer = solver->get_halo_buffer();

    m_conductivity_particles[idx] = new particle_container<godunov_particle>();
    m_rho_dagger_particles[idx]   = new particle_container<godunov_particle>();
    
    m_amr->allocate(*m_conductivity_particles[idx], pvr_buffer, halo_buffer, m_particle_realm);
    m_amr->allocate(*m_rho_dagger_particles[idx],   pvr_buffer, halo_buffer, m_particle_realm);
  }
}

void ito_plasma_godunov::parse_options() {
  CH_TIME("ito_plasma_godunov::parse_options");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_options" << endl;
  }

  ParmParse pp(m_name.c_str());
  std::string str;
  
  pp.get("verbosity",      m_verbosity);
  pp.get("ppc",            m_ppc);
  pp.get("max_cells_hop",  m_max_cells_hop);
  pp.get("merge_interval", m_merge_interval);
  pp.get("relax_factor",   m_relax_factor);
  pp.get("regrid_super",   m_regrid_superparticles);
  pp.get("algorithm",      str);
  pp.get("load_balance",   m_load_balance);
  pp.get("min_dt",         m_min_dt);
  pp.get("max_dt",         m_max_dt);
  pp.get("halo_buffer",    m_halo_buffer);
  pp.get("pvr_buffer",     m_pvr_buffer);


  if(str == "euler_maruyama"){
    m_algorithm = which_algorithm::euler_maruyama;
  }
  else if(str ==  "trapezoidal"){
    m_algorithm = which_algorithm::trapezoidal;
  }
  else{
    MayDay::Abort("ito_plasma_godunov::parse_options - unknown algorithm requested");
  }


  pp.get("which_dt", str);
  if(str == "advection"){
    m_whichDt = which_dt::advection;
  }
  else if(str == "diffusion"){
    m_whichDt = which_dt::diffusion;
  }
  else if(str == "advection_diffusion"){
    m_whichDt = which_dt::advection_diffusion;
  }
  else{
    MayDay::Abort("ito_plasma_godunov::parse_options - unknown 'which_dt' requested");
  }

  // Setup runtime storage (requirements change with algorithm)
  this->setup_runtime_storage();
}

void ito_plasma_godunov::allocate_internals(){
  CH_TIME("ito_plasma_godunov::allocate_internals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocate_internals" << endl;
  }

  m_amr->allocate(m_fluid_scratch1,    m_fluid_realm,    m_phase, 1);
  m_amr->allocate(m_fluid_scratchD,    m_fluid_realm,    m_phase, SpaceDim);
  
  m_amr->allocate(m_particle_scratch1,  m_particle_realm, m_phase, 1);
  m_amr->allocate(m_particle_scratchD,  m_particle_realm, m_phase, SpaceDim);
  m_amr->allocate(m_particle_E,         m_particle_realm, m_phase, SpaceDim);

  m_amr->allocate(m_J,            m_fluid_realm, m_phase, SpaceDim);
  m_amr->allocate(m_scratch1,     m_fluid_realm, m_phase, 1);
  m_amr->allocate(m_scratch2,     m_fluid_realm, m_phase, 1);
  m_amr->allocate(m_conduct_cell, m_fluid_realm, m_phase, 1);
  m_amr->allocate(m_conduct_face, m_fluid_realm, m_phase, 1);
  m_amr->allocate(m_conduct_eb,   m_fluid_realm, m_phase, 1);
  m_amr->allocate(m_fluid_E,      m_fluid_realm, m_phase, SpaceDim);

  // Allocate for energy sources
  const int num_ito_species = m_physics->get_num_ito_species();
  m_energy_sources.resize(num_ito_species);
  for (int i = 0; i < m_energy_sources.size(); i++){
    m_amr->allocate(m_energy_sources[i],  m_particle_realm, m_phase, 1);
  }

  // Allocate fluid scratch storage
  m_fscratch1.resize(num_ito_species);
  m_fscratch2.resize(num_ito_species);
  for (int i = 0; i < num_ito_species; i++){
    m_amr->allocate(m_fscratch1[i], m_fluid_realm, m_phase, 1);
    m_amr->allocate(m_fscratch2[i], m_fluid_realm, m_phase, 1);
  }
}

Real ito_plasma_godunov::advance(const Real a_dt) {
  CH_TIME("ito_plasma_godunov::advance");
  if(m_verbosity > 5){
    pout() << m_name + "::advance" << endl;
  }


  Real particle_time = 0.0;
  Real relax_time    = 0.0;
  Real photon_time   = 0.0;
  Real sort_time     = 0.0;
  Real super_time    = 0.0;
  Real reaction_time = 0.0;
  Real clear_time    = 0.0;
  Real deposit_time  = 0.0;
  Real velo_time     = 0.0;
  Real diff_time     = 0.0;

  Real total_time    = -MPI_Wtime();
  
  // Particle algorithms
  MPI_Barrier(Chombo_MPI::comm);
  particle_time -= MPI_Wtime();
  switch(m_algorithm){
  case which_algorithm::euler_maruyama:
    this->advance_particles_euler_maruyama(a_dt);
    break;
  case which_algorithm::trapezoidal:
    this->advance_particles_trapezoidal(a_dt);
    break;
  default:
    MayDay::Abort("ito_plasma_godunov::advance - logic bust");
  }
  particle_time += MPI_Wtime();

  // Compute current and relaxation time.
  MPI_Barrier(Chombo_MPI::comm);
  relax_time = -MPI_Wtime();
  this->compute_J(m_J, a_dt);
  m_dt_relax = this->compute_relaxation_time(); // This is for the restricting the next step.
  relax_time += MPI_Wtime();

  // Move photons
  MPI_Barrier(Chombo_MPI::comm);
  photon_time = -MPI_Wtime();
  this->advance_photons(a_dt);
  photon_time += MPI_Wtime();

  // If we are using the LEA, we must compute the Ohmic heating term. This must be done
  // BEFORE sorting the particles per cell. 
  if(m_physics->get_coupling() == ito_plasma_physics::coupling::LEA){
    this->compute_EdotJ_source();
  }
  
  // Sort the particles and photons per cell so we can call reaction algorithms
  MPI_Barrier(Chombo_MPI::comm);
  sort_time = -MPI_Wtime();
  m_ito->sort_particles_by_cell();
  this->sort_bulk_photons_by_cell();
  this->sort_source_photons_by_cell();
  sort_time += MPI_Wtime();

  // Chemistry kernel.
  MPI_Barrier(Chombo_MPI::comm);
  reaction_time = -MPI_Wtime();
  this->advance_reaction_network(a_dt);
  reaction_time += MPI_Wtime();

  // Make superparticles
  MPI_Barrier(Chombo_MPI::comm);
  super_time = -MPI_Wtime();
  if((m_step+1) % m_merge_interval == 0 && m_merge_interval > 0){
    m_ito->make_superparticles(m_ppc);
  }
  super_time += MPI_Wtime();

  // Sort particles per patch.
  MPI_Barrier(Chombo_MPI::comm);
  sort_time -= MPI_Wtime();
  m_ito->sort_particles_by_patch();
  this->sort_bulk_photons_by_patch();
  this->sort_source_photons_by_patch();
  sort_time += MPI_Wtime();

  // Clear other data holders for now. BC comes later...
  MPI_Barrier(Chombo_MPI::comm);
  clear_time = -MPI_Wtime();
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->clear(solver_it()->get_eb_particles());
    solver_it()->clear(solver_it()->get_domain_particles());
  }
  clear_time += MPI_Wtime();

  //
  MPI_Barrier(Chombo_MPI::comm);
  deposit_time -= MPI_Wtime();
  m_ito->deposit_particles();
  deposit_time += MPI_Wtime();

  // Prepare next step
  MPI_Barrier(Chombo_MPI::comm);
  velo_time -= MPI_Wtime();
  this->compute_ito_velocities();
  velo_time += MPI_Wtime();

  MPI_Barrier(Chombo_MPI::comm);
  diff_time -= MPI_Wtime();
  this->compute_ito_diffusion();
  diff_time += MPI_Wtime();

  total_time += MPI_Wtime();

  // Convert to %
  particle_time *= 100./total_time;
  relax_time    *= 100./total_time;
  photon_time   *= 100./total_time;
  sort_time     *= 100./total_time;
  super_time    *= 100./total_time;
  reaction_time *= 100./total_time;
  clear_time    *= 100./total_time;
  deposit_time  *= 100./total_time;
  velo_time     *= 100./total_time;
  diff_time     *= 100./total_time;



  if(m_profile){

    // Dummy-check total percentage
    Real percentage = 0.0;
    percentage += particle_time;
    percentage += relax_time;
    percentage += photon_time;
    percentage += sort_time;
    percentage += super_time;
    percentage += reaction_time;
    percentage += clear_time;
    percentage += deposit_time;
    percentage += velo_time;
    percentage += diff_time;
    
    pout() << endl
      	   << "ito_plasma_godunov::advance breakdown:" << endl
	   << "======================================" << endl
	   << "particle time = " << particle_time << "%" << endl
	   << "relax time    = " << relax_time << "%" << endl
	   << "photon time   = " << photon_time << "%" << endl
	   << "sort time     = " << sort_time << "%" << endl
	   << "super time    = " << super_time << "%" << endl
	   << "reaction time = " << reaction_time << "%" << endl
	   << "clear time    = " << clear_time << "%" << endl
      	   << "deposit time  = " << deposit_time << "%" << endl
	   << "velo time     = " << velo_time << "%" << endl
      	   << "diff time     = " << diff_time << "%" << endl
	   << "total percent = " << percentage << endl
	   << "imbalance     = " << 100. - percentage << endl
	   << "total time    = " << total_time << " (seconds)" << endl
	   << endl;
  }
  
  return a_dt;
}

void ito_plasma_godunov::compute_dt(Real& a_dt, time_code& a_timecode){
  CH_TIME("ito_plasma_godunov::compute_dt");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::compute_dt" << endl;
  }

  a_dt = 1.E99;
  
  if(m_whichDt == which_dt::advection){
    a_dt = m_ito->compute_advective_dt();
  }
  else if(m_whichDt == which_dt::diffusion){
    a_dt = m_ito->compute_diffusive_dt();
  }
  else if(m_whichDt == which_dt::advection_diffusion){
    a_dt = m_ito->compute_dt();
  }
    
  a_dt = a_dt*m_max_cells_hop;
  a_timecode = time_code::cfl;


  // Physics-based restriction
  const Real physicsDt = this->compute_physics_dt();
  if(physicsDt < a_dt){
    a_dt = physicsDt;
    a_timecode = time_code::physics;
  }

  if(a_dt < m_min_dt){
    a_dt = m_min_dt;
    a_timecode = time_code::hardcap;
  }

  if(a_dt > m_max_dt){
    a_dt = m_max_dt;
    a_timecode = time_code::hardcap;
  }

  m_timecode = a_timecode;

#if 0 // Debug code
  const Real dtCFL = m_ito->compute_dt();
  m_avg_cfl += a_dt/dtCFL;
  if(procID() == 0) std::cout << "dt = " << a_dt
			      << "\t relax dt = " << m_dt_relax
			      << "\t factor = " << a_dt/m_dt_relax
			      << "\t CFL = " << a_dt/dtCFL
			      << "\t avgCFL = " << m_avg_cfl/(1+m_step)
			      << std::endl;
#endif
}

void ito_plasma_godunov::pre_regrid(const int a_lmin, const int a_old_finest_level){
  CH_TIME("ito_plasma_godunov::pre_regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::pre_regrid" << endl;
  }

  ito_plasma_stepper::pre_regrid(a_lmin, a_old_finest_level);


  // Copy conductivity to scratch storage
  const int ncomp        = 1;
  const int finest_level = m_amr->get_finest_level();
  m_amr->allocate(m_cache,  m_fluid_realm, m_phase, ncomp);
  for (int lvl = 0; lvl <= a_old_finest_level; lvl++){
    m_conduct_cell[lvl]->localCopyTo(*m_cache[lvl]);
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
      
    m_conductivity_particles[idx]->pre_regrid(a_lmin);
    m_rho_dagger_particles[idx]->pre_regrid(a_lmin);
  }
}

void ito_plasma_godunov::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) {
  CH_TIME("ito_plasma_godunov::regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::regrid" << endl;
  }

  // Regrid solvers
  m_ito->regrid(a_lmin,     a_old_finest_level, a_new_finest_level);
  m_poisson->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  m_rte->regrid(a_lmin,     a_old_finest_level, a_new_finest_level);
  m_sigma->regrid(a_lmin,   a_old_finest_level, a_new_finest_level);

  // Allocate internal memory for ito_plasma_godunov now....
  this->allocate_internals();

  // We need to remap/regrid the stored particles as well. 
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids(m_particle_realm);
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real>& dx                 = m_amr->get_dx();
  const Vector<int>& ref_rat             = m_amr->get_ref_rat();

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rho_dagger_particles[idx]->regrid(  grids, domains, dx, ref_rat, a_lmin, a_new_finest_level);
    m_conductivity_particles[idx]->regrid(grids, domains, dx, ref_rat, a_lmin, a_new_finest_level);
  }
  
  // Recompute the conductivity and space charge densities. 
  this->compute_regrid_conductivity();
  this->compute_regrid_rho();
  this->setup_semi_implicit_poisson(m_prevDt);

  const bool converged = this->solve_poisson();
  if(!converged){
    MayDay::Abort("ito_plasma_godunov::regrid - Poisson solve did not converge after regrid!!!");
  }

  // Regrid superparticles. 
  if(m_regrid_superparticles){
    m_ito->sort_particles_by_cell();
    m_ito->make_superparticles(m_ppc);
    m_ito->sort_particles_by_patch();
  }

  // Now let the ito solver deposit its actual particles... In the above it deposit m_rho_dagger_particles. 
  m_ito->deposit_particles();

  // Recompute new velocities and diffusion coefficients
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_godunov::setup_runtime_storage(){
    CH_TIME("ito_plasma_godunov::setup_runtime_storage");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_runtime_storage" << endl;
  }

  switch (m_algorithm){
  case which_algorithm::euler_maruyama:
    ito_particle::set_num_runtime_vectors(1);
    break;
  case which_algorithm::trapezoidal:
    ito_particle::set_num_runtime_vectors(2); // For V^k and the diffusion hop. 
    break;
  default:
    MayDay::Abort("ito_plasma_godunov::setup_runtime_storage - logic bust");
  }
}



void ito_plasma_godunov::set_old_positions(){
  CH_TIME("ito_plasma_godunov::set_old_positions()");
  if(m_verbosity > 5){
    pout() << m_name + "::set_old_positions()" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    
    if(solver->is_diffusive()){
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();


	  for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
	    ito_particle& p = particleList[lit];
	    p.oldPosition() = p.position();
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::remap_godunov_particles(Vector<particle_container<godunov_particle>* >& a_particles, const which_particles a_which_particles){
  CH_TIME("ito_plasma_godunov::remap_godunov_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::remap_godunov_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>&   solver = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();
    const bool charged   = species->get_charge() != 0;

    switch(a_which_particles) {
    case which_particles::all:
      a_particles[idx]->remap();
      break;
    case which_particles::all_mobile:
      if(mobile) a_particles[idx]->remap();
      break;
    case which_particles::all_diffusive:
      if(diffusive) a_particles[idx]->remap();
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) a_particles[idx]->remap();
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) a_particles[idx]->remap();
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) a_particles[idx]->remap();
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) a_particles[idx]->remap();
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) a_particles[idx]->remap();
      break;
    default:
      MayDay::Abort("ito_plasma_godunov::remap_godunov_particles - logic bust");
    }
  }
}

void ito_plasma_godunov::deposit_godunov_particles(const Vector<particle_container<godunov_particle>* >& a_particles, const which_particles a_which_particles){
  CH_TIME("ito_plasma_godunov::deposit_godunov_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_godunov_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>&   solver = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();
    const bool charged   = species->get_charge() != 0;

    switch(a_which_particles) {
    case which_particles::all:
      solver->deposit_particles(solver->get_state(), *a_particles[idx]);
      break;
    case which_particles::all_mobile:
      if(mobile) solver->deposit_particles(solver->get_state(), *a_particles[idx]);
      break;
    case which_particles::all_diffusive:
      if(diffusive) solver->deposit_particles(solver->get_state(), *a_particles[idx]);
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) solver->deposit_particles(solver->get_state(), *a_particles[idx]);
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) solver->deposit_particles(solver->get_state(), *a_particles[idx]);
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) solver->deposit_particles(solver->get_state(), *a_particles[idx]);
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) solver->deposit_particles(solver->get_state(), *a_particles[idx]);
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) solver->deposit_particles(solver->get_state(), *a_particles[idx]);
      break;
    default:
      MayDay::Abort("ito_plasma_godunov::deposit_godunov_particles - logic bust");
    }
  }
}

void ito_plasma_godunov::clear_godunov_particles(const Vector<particle_container<godunov_particle>* >& a_particles, const which_particles a_which_particles){
  CH_TIME("ito_plasma_godunov::clear_godunov_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_clear_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>&   solver = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();
    const bool charged   = species->get_charge() != 0;

    switch(a_which_particles) {
    case which_particles::all:
      a_particles[idx]->clear_particles();
      break;
    case which_particles::all_mobile:
      if(mobile) a_particles[idx]->clear_particles();
      break;
    case which_particles::all_diffusive:
      if(diffusive) a_particles[idx]->clear_particles();
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) a_particles[idx]->clear_particles();
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) a_particles[idx]->clear_particles();
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) a_particles[idx]->clear_particles();
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) a_particles[idx]->clear_particles();
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) a_particles[idx]->clear_particles();
      break;
    default:
      MayDay::Abort("ito_plasma_godunov::clear_godunov_particles - logic bust");
    }
  }
}

void ito_plasma_godunov::compute_all_conductivities(const Vector<particle_container<godunov_particle>* >& a_particles){
  CH_TIME("ito_plasma_godunov::compute_all_conductivities");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_all_conductivities" << endl;
  }

  this->compute_cell_conductivity(m_conduct_cell, a_particles);

  // Now do the faces
  this->compute_face_conductivity();
}

void ito_plasma_godunov::compute_cell_conductivity(EBAMRCellData& a_conductivity, const Vector<particle_container<godunov_particle>* >& a_particles){
  CH_TIME("ito_plasma_godunov::compute_cell_conductivity(conductivity, godunov_particle");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_cell_conductivity(conductivity, godunov_particle)" << endl;
  }

  data_ops::set_value(a_conductivity, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>&   solver = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();
    
    const int idx = solver_it.get_solver();
    const int q   = species->get_charge();

    if(Abs(q) > 0 && solver->is_mobile()){
      data_ops::set_value(m_particle_scratch1, 0.0);

      solver->deposit_particles(m_particle_scratch1, *a_particles[idx]); // The particles should have "masses" = m*mu 

      // Copy to fluid realm and add to fluid stuff
      m_fluid_scratch1.copy(m_particle_scratch1);
      data_ops::incr(a_conductivity, m_fluid_scratch1, Abs(q));
    }
  }

  data_ops::scale(a_conductivity, units::s_Qe);

  m_amr->average_down(a_conductivity, m_fluid_realm, m_phase);
  m_amr->interp_ghost_pwl(a_conductivity, m_fluid_realm, m_phase);

  // See if this helps....
  m_amr->interpolate_to_centroids(a_conductivity, m_fluid_realm, m_phase);

}

void ito_plasma_godunov::compute_face_conductivity(){
  CH_TIME("ito_plasma_godunov::compute_face_conductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_face_conductivity" << endl;
  }

  // This code does averaging from cell to face. 
  data_ops::average_cell_to_face_allcomps(m_conduct_face, m_conduct_cell, m_amr->get_domains());

  // This code extrapolates the conductivity to the EB. This should actually be the EB centroid but since the stencils
  // for EB extrapolation can be a bit nasty (e.g. negative weights), we do the centroid instead and take that as an approximation.
#if 0
  const irreg_amr_stencil<centroid_interp>& ebsten = m_amr->get_centroid_interp_stencils(m_fluid_realm, m_phase);
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    ebsten.apply(m_conduct_eb, m_conduct_cell, lvl);
  }
#else
  data_ops::set_value(m_conduct_eb, 0.0);
  data_ops::incr(m_conduct_eb, m_conduct_cell, 1.0);
#endif

}

void ito_plasma_godunov::setup_semi_implicit_poisson(const Real a_dt){
  CH_TIME("ito_plasma_godunov::setup_semi_implicit_poisson");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_semi_implicit_poisson" << endl;
  }

  poisson_multifluid_gmg* poisson = (poisson_multifluid_gmg*) (&(*m_poisson));

  // Set coefficients as usual
  poisson->set_coefficients();

  // Get bco and increment with mobilities
  MFAMRFluxData& bco   = poisson->get_bco();
  MFAMRIVData& bco_irr = poisson->get_bco_irreg();
  
  EBAMRFluxData bco_gas;
  EBAMRIVData   bco_irr_gas;
  
  m_amr->allocate_ptr(bco_gas);
  m_amr->allocate_ptr(bco_irr_gas);
  
  m_amr->alias(bco_gas,     phase::gas, bco);
  m_amr->alias(bco_irr_gas, phase::gas, bco_irr);

  data_ops::scale(m_conduct_face, a_dt/units::s_eps0);
  data_ops::scale(m_conduct_eb,   a_dt/units::s_eps0);

  data_ops::multiply(m_conduct_face, bco_gas);
  data_ops::multiply(m_conduct_eb,   bco_irr_gas);

  data_ops::incr(bco_gas,     m_conduct_face, 1.0);
  data_ops::incr(bco_irr_gas, m_conduct_eb,   1.0);

  // Set up the multigrid solver
  poisson->setup_operator_factory();
  poisson->setup_solver();
  poisson->set_needs_setup(false);
}

void ito_plasma_godunov::setup_standard_poisson(){
  CH_TIME("ito_plasma_godunov::setup_standard_poisson");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_standard_poisson" << endl;
  }

  poisson_multifluid_gmg* poisson = (poisson_multifluid_gmg*) (&(*m_poisson));

  // Set coefficients as usual
  poisson->set_coefficients();

  // Set up the multigrid solver
  poisson->setup_operator_factory();
  poisson->setup_solver();
  poisson->set_needs_setup(false);
}

void ito_plasma_godunov::copy_conductivity_particles(Vector<particle_container<godunov_particle>* >& a_conductivity_particles){
  CH_TIME("ito_plasma_godunov::copy_conductivity_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::copy_conductivity_particles" << endl;
  }

  this->clear_godunov_particles(a_conductivity_particles, which_particles::all);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver   = solver_it();
    const RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();
    const int q   = species->get_charge();

    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const List<ito_particle>& ito_parts = solver->get_particles()[lvl][dit()].listItems();
	List<godunov_particle>& gdnv_parts  = (*a_conductivity_particles[idx])[lvl][dit()].listItems();

	if(q != 0 && solver->is_mobile()){
	  for (ListIterator<ito_particle> lit(ito_parts); lit.ok(); ++lit){
	    const ito_particle& p = lit();
	    const RealVect& pos   = p.position();
	    const Real& mass      = p.mass();
	    const Real& mobility  = p.mobility();

	    gdnv_parts.add(godunov_particle(pos, mass*mobility));
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::copy_rho_dagger_particles(Vector<particle_container<godunov_particle>* >& a_rho_dagger_particles){
  CH_TIME("ito_plasma_godunov::copy_rho_dagger_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::copy_rho_dagger_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver   = solver_it();
    const RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();
    const int q   = species->get_charge();

    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	const List<ito_particle>& ito_parts = solver->get_particles()[lvl][dit()].listItems();
	List<godunov_particle>& gdnv_parts  = (*a_rho_dagger_particles[idx])[lvl][dit()].listItems();

	gdnv_parts.clear();

	if(q != 0){
	  for (ListIterator<ito_particle> lit(ito_parts); lit.ok(); ++lit){
	    const ito_particle& p = lit();
	    const RealVect& pos   = p.position();
	    const Real& mass      = p.mass();

	    gdnv_parts.add(godunov_particle(pos, mass));
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::compute_regrid_conductivity(){
  CH_TIME("ito_plasma_godunov::compute_regrid_conductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_regrid_conductivity" << endl;
  }

  this->compute_all_conductivities(m_conductivity_particles);
}

void ito_plasma_godunov::compute_regrid_rho(){
  CH_TIME("ito_plasma_godunov::compute_regrid_rho");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_regrid_rho" << endl;
  }

  this->deposit_godunov_particles(m_rho_dagger_particles, which_particles::all);
}

void ito_plasma_godunov::advance_particles_euler_maruyama(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advance_particles_euler_maruyama");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_particles_euler_maruyama" << endl;
  }

  Real posTime = 0.0;
  Real diffuseTime = 0.0;
  Real remapGdnvTime = 0.0;
  Real depositGdnvTime = 0.0;
  Real copyCondTime = 0.0;
  Real condTime = 0.0;
  Real setupTime = 0.0;
  Real poissonTime = 0.0;
  Real velocityTime = 0.0;
  Real particleTime = 0.0;
  Real remapTime = 0.0;
  Real isectTime = 0.0;
  Real depositTime = 0.0;
  Real totalTime = 0.0;


  totalTime -= MPI_Wtime();

  m_prevDt = a_dt; // Needed for regrids.

  // 1. Store X^k positions.
  MPI_Barrier(Chombo_MPI::comm);
  posTime -= MPI_Wtime();
  this->set_old_positions();
  posTime += MPI_Wtime();

  // 2. Diffuse the particles. This copies onto m_rho_dagger_particles and stores the hop on the full particles.
  MPI_Barrier(Chombo_MPI::comm);
  diffuseTime -= MPI_Wtime();
  this->diffuse_particles_euler_maruyama(m_rho_dagger_particles, a_dt);
  diffuseTime += MPI_Wtime();

  MPI_Barrier(Chombo_MPI::comm);
  remapGdnvTime -= MPI_Wtime();
  this->remap_godunov_particles(m_rho_dagger_particles,   which_particles::all_diffusive);
  remapGdnvTime += MPI_Wtime();

  MPI_Barrier(Chombo_MPI::comm);
  depositGdnvTime -= MPI_Wtime();
  this->deposit_godunov_particles(m_rho_dagger_particles, which_particles::all_diffusive);
  depositGdnvTime += MPI_Wtime();

  // 3. Solve the semi-implicit Poisson equation. Also, copy the particles used for computing the conductivity to scratch.
  MPI_Barrier(Chombo_MPI::comm);
  copyCondTime -= MPI_Wtime();
  this->copy_conductivity_particles(m_conductivity_particles); // Sets particle "weights" = w*mu
  copyCondTime += MPI_Wtime();

  MPI_Barrier(Chombo_MPI::comm);
  condTime -= MPI_Wtime();
  this->compute_all_conductivities(m_conductivity_particles);  // Deposits q_e*Z*w*mu on the mesh
  condTime += MPI_Wtime();

  MPI_Barrier(Chombo_MPI::comm);
  setupTime -= MPI_Wtime();
  this->setup_semi_implicit_poisson(a_dt);                     // Multigrid setup
  setupTime += MPI_Wtime();

  MPI_Barrier(Chombo_MPI::comm);
  poissonTime -= MPI_Wtime();
  this->solve_poisson();                                       // Solve the stinking equation.
  poissonTime += MPI_Wtime();

  // 4. Recompute velocities with the new electric field, then do the actual semi-implicit Euler-Maruyama update.
  MPI_Barrier(Chombo_MPI::comm);
  velocityTime -= MPI_Wtime();
  this->set_ito_velocity_funcs();
  m_ito->interpolate_velocities();
  velocityTime += MPI_Wtime();

  MPI_Barrier(Chombo_MPI::comm);
  particleTime -= MPI_Wtime();
  this->step_euler_maruyama(a_dt);
  particleTime += MPI_Wtime();

  MPI_Barrier(Chombo_MPI::comm);
  remapTime -= MPI_Wtime();
  this->remap_particles(which_particles::all_mobile_or_diffusive);
  remapTime += MPI_Wtime();

  // 5. Do intersection test and remove EB particles. These particles are NOT allowed to react later.
  MPI_Barrier(Chombo_MPI::comm);
  isectTime -= MPI_Wtime();
  this->intersect_particles(which_particles::all_mobile_or_diffusive, EB_representation::implicit_function);
  this->remove_eb_particles(which_particles::all_mobile_or_diffusive, EB_representation::discrete);         
  isectTime += MPI_Wtime();

  // 6. Deposit particles. This shouldn't be necessary unless we want to compute (E,J)
  MPI_Barrier(Chombo_MPI::comm);
  depositTime -= MPI_Wtime();
  this->deposit_particles(which_particles::all_mobile_or_diffusive);
  depositTime += MPI_Wtime();

  totalTime += MPI_Wtime();

  if(m_profile){

    posTime *= 100./totalTime;
    diffuseTime *= 100./totalTime;
    remapGdnvTime *= 100./totalTime;
    depositGdnvTime *= 100./totalTime;
    copyCondTime *= 100./totalTime;
    condTime *= 100./totalTime;
    setupTime *= 100./totalTime;
    poissonTime *= 100./totalTime;
    velocityTime *= 100./totalTime;
    particleTime *= 100./totalTime;
    remapTime *= 100./totalTime;
    isectTime *= 100./totalTime;
    depositTime *= 100./totalTime;

    Real percent = 0.0;
    percent += posTime;
    percent += diffuseTime;
    percent += remapGdnvTime;
    percent += copyCondTime;
    percent += condTime;
    percent += setupTime;
    percent += poissonTime;
    percent += velocityTime;
    percent += particleTime;
    percent += remapTime;
    percent += isectTime;
    percent += depositTime;
    
    pout() << endl
      	   << "ito_plasma_godunov::euler_maruyama breakdown:" << endl
	   << "======================================" << endl
	   << "posTime         = " << posTime << "%" << endl
      	   << "diffuseTime     = " << diffuseTime << "%" << endl
	   << "remapGdnvTime   = " << remapGdnvTime << "%" << endl
      	   << "depositGdnvTime = " << depositGdnvTime << "%" << endl
	   << "copyCondTime    = " << copyCondTime << "%" << endl
      	   << "condTime        = " << condTime << "%" << endl
	   << "setupTime       = " << setupTime << "%" << endl
      	   << "poissonTime     = " << poissonTime << "%" << endl
	   << "velocityTime    = " << velocityTime << "%" << endl
      	   << "particleTime    = " << particleTime << "%" << endl
	   << "remapTime       = " << remapTime << "%" << endl
      	   << "isectTime       = " << isectTime << "%" << endl
	   << "depositTime     = " << depositTime << "%" << endl
	   << "total percent   = " << percent << "%" << endl
      	   << "imbalance       = " << 100. - percent << "%" << endl
	   << "total time      = " << totalTime << " (seconds)" << endl
	   << endl;
  }
}

void ito_plasma_godunov::diffuse_particles_euler_maruyama(Vector<particle_container<godunov_particle>* >& a_rho_dagger, const Real a_dt){
  CH_TIME("ito_plasma_godunov::diffuse_particles_euler_maruyama");
  if(m_verbosity > 5){
    pout() << m_name + "::diffuse_particles_euler_maruyama" << endl;
  }

  this->clear_godunov_particles(a_rho_dagger, which_particles::all);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver   = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();

    const Real f = mobile    ? 1.0 : 0.0; // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion
    
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
      ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	List<ito_particle>& ito_particles  = particles[dit()].listItems();
	List<godunov_particle>& gdnv_parts = (*a_rho_dagger[idx])[lvl][dit()].listItems();

	// Store the diffusion hop, and add the godunov particles
	for (ListIterator<ito_particle> lit(ito_particles); lit.ok(); ++lit){
	  ito_particle& p     = lit();
	  const Real factor   = g*sqrt(2.0*p.diffusion()*a_dt);
	  const RealVect hop  = factor*solver->random_gaussian();
	  const RealVect& pos = p.position();
	  const Real& mass    = p.mass();

	  // Store the diffusion hop
	  p.runtime_vector(0) = hop;

	  // Add simpler particle
	  gdnv_parts.add(godunov_particle(pos + hop, mass));
	}
      }
    }
  }
}

void ito_plasma_godunov::step_euler_maruyama(const Real a_dt){
  CH_TIME("ito_plasma_godunov::step_euler_maruyama");
  if(m_verbosity > 5){
    pout() << m_name + "::step_euler_maruyama" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();

    const Real f = mobile    ? 1.0 : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;
    
    if(mobile || diffusive){
      
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();

	  for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
	    ito_particle& p = lit();

	    // Add diffusion hop again. The position after the diffusion hop is oldPosition() and X^k is in position()
	    const RealVect& hop = p.runtime_vector(0);
	    p.position()        = p.oldPosition() + f*p.velocity()*a_dt + g*hop;
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::advance_particles_trapezoidal(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advance_particles_trapezoidal");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_particles_trapezoidal" << endl;
  }

  m_prevDt = 0.5*a_dt;  // Needed for regrids. 

  this->set_old_positions();

  // ====== PREDICTOR BEGIN ======
  this->pre_trapezoidal_predictor(m_rho_dagger_particles, a_dt);
  this->remap_godunov_particles(m_rho_dagger_particles,   which_particles::all_diffusive); // Particles that were copied but not moved are in the right box.
  this->deposit_godunov_particles(m_rho_dagger_particles, which_particles::all);           // All copies need to deposit. 

  this->copy_conductivity_particles(m_conductivity_particles); 
  this->compute_all_conductivities(m_conductivity_particles);        
  this->setup_semi_implicit_poisson(a_dt);                     
  this->solve_poisson();                                       

  this->set_ito_velocity_funcs();
  m_ito->interpolate_velocities();
  this->trapezoidal_predictor(a_dt); 
  this->remap_particles(which_particles::all_mobile_or_diffusive);
  // ====== PREDICTOR END ======

  // ====== CORRECTOR BEGIN =====
  this->pre_trapezoidal_corrector(m_rho_dagger_particles, a_dt);                                    // Mobile or diffusive moves to X^dagger = X^k + 0.5*dt*V^k + hop
  this->remap_godunov_particles(m_rho_dagger_particles, which_particles::all_mobile_or_diffusive);  // Only need to remap particles that were mobile or diffusive
  this->deposit_godunov_particles(m_rho_dagger_particles, which_particles::all);                    // Everything needs to deposit...

  this->copy_conductivity_particles(m_conductivity_particles); 
  this->compute_all_conductivities(m_conductivity_particles);        
  this->setup_semi_implicit_poisson(0.5*a_dt);                 
  this->solve_poisson();                                       

  this->set_ito_velocity_funcs();
  m_ito->interpolate_velocities();
  this->trapezoidal_corrector(a_dt); 
  this->remap_particles(which_particles::all_mobile_or_diffusive);
  // ====== CORRECTOR END =====

  // Do particle-boundary intersection. 
  this->intersect_particles(which_particles::all_mobile_or_diffusive, EB_representation::implicit_function);
  this->remove_eb_particles(which_particles::all_mobile_or_diffusive, EB_representation::implicit_function);

  // Finally, deposit particles. 
  this->deposit_particles(which_particles::all_mobile_or_diffusive);
}

void ito_plasma_godunov::pre_trapezoidal_predictor(Vector<particle_container<godunov_particle>* >& a_rho_dagger, const Real a_dt){
  CH_TIME("ito_plasma_godunov::pre_trapezoidal_predictor");
  if(m_verbosity > 5){
    pout() << m_name + "::pre_trapezoidal_predictor" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver   = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();

    const Real f = mobile    ? 1.0 : 0.0; // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion
    
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
      ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	List<ito_particle>& ito_particles  = particles[dit()].listItems();
	List<godunov_particle>& gdnv_parts = (*a_rho_dagger[idx])[lvl][dit()].listItems();

	gdnv_parts.clear();

	// Store the diffusion hop, and add the godunov particles
	for (ListIterator<ito_particle> lit(ito_particles); lit.ok(); ++lit){
	  ito_particle& p     = lit();
	  const Real factor   = sqrt(2.0*p.diffusion()*a_dt);
	  const RealVect hop  = factor*solver->random_gaussian();
	  const RealVect& Xk  = p.oldPosition();
	  const Real& mass    = p.mass();

	  // Store the diffusion hop and the current velocity. 
	  p.runtime_vector(0) = g*hop;
	  p.runtime_vector(1) = f*p.velocity();

	  // Add simpler particle
	  gdnv_parts.add(godunov_particle(Xk + g*hop, mass));
	}
      }
    }
  }
}

void ito_plasma_godunov::trapezoidal_predictor(const Real a_dt){
  CH_TIME("ito_plasma_godunov::trapezoidal_predictor");
  if(m_verbosity > 5){
    pout() << m_name + "::trapezoidal_predictor" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();

    const Real f = mobile    ? 1.0 : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;
    
    if(mobile || diffusive){
      
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();

	  for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
	    ito_particle& p = lit();

	    // Add diffusion hop again. The position after the diffusion hop is oldPosition() and X^k is in position()
	    const RealVect& hop = p.runtime_vector(0);
	    const RealVect& Vk  = p.runtime_vector(1);
	    
	    p.position()        = p.oldPosition() + f*p.velocity()*a_dt + g*hop;
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::pre_trapezoidal_corrector(Vector<particle_container<godunov_particle>* >& a_rho_dagger, const Real a_dt){
  CH_TIME("ito_plasma_godunov::pre_trapezoidal_corrector");
  if(m_verbosity > 5){
    pout() << m_name + "::pre_trapezoidal_corrector" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver   = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();

    const Real f = mobile    ? 1.0 : 0.0; // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion
    
    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
      ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	List<ito_particle>& ito_particles  = particles[dit()].listItems();
	List<godunov_particle>& gdnv_parts = (*a_rho_dagger[idx])[lvl][dit()].listItems();

	gdnv_parts.clear();

	// Store the diffusion hop, and add the godunov particles
	for (ListIterator<ito_particle> lit(ito_particles); lit.ok(); ++lit){
	  ito_particle& p     = lit();
	  
	  const Real& mass    = p.mass();
	  const RealVect& Xk  = p.oldPosition();
	  const RealVect& hop = p.runtime_vector(0);
	  const RealVect& Vk  = p.runtime_vector(1);

	  // Move particle. 
	  const RealVect pos  = Xk + 0.5*a_dt*f*Vk + g*hop;
	  gdnv_parts.add(godunov_particle(pos, mass));
	}
      }
    }
  }
}

void ito_plasma_godunov::trapezoidal_corrector(const Real a_dt){
  CH_TIME("ito_plasma_godunov::trapezoidal_corrector");
  if(m_verbosity > 5){
    pout() << m_name + "::trapezoidal_corrector" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();

    const Real f = mobile    ? 1.0 : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;
    
    if(mobile || diffusive){
      
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();

	  for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
	    ito_particle& p = lit();

	    const RealVect& Xk  = p.oldPosition();
	    const RealVect& hop = p.runtime_vector(0);
	    const RealVect& Vk  = p.runtime_vector(1);
	    const RealVect& Vk1 = p.velocity();
	    
	    p.position() = Xk + 0.5*f*a_dt*(Vk + Vk1) + g*hop;
	  }
	}
      }
    }
  }
}
