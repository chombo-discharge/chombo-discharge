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

  m_avg_cfl = 0.0;
}

ito_plasma_godunov::~ito_plasma_godunov(){

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

  if(str == "euler"){
    m_algorithm = which_algorithm::euler;
  }
  else if(str == "semi_implicit"){
    m_algorithm = which_algorithm::semi_implicit;
  }
  else if(str == "split_semi"){
    m_algorithm = which_algorithm::split_semi_implicit;
  }
  else if(str == "split_pc"){
    m_algorithm = which_algorithm::split_pc;
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

    const int idx        = solver_it.get_solver();
    const int pvr_buffer = solver->get_pvr_buffer();

    m_conductivity_particles[idx] = new particle_container<godunov_particle>();
    m_rho_dagger_particles[idx]   = new particle_container<godunov_particle>();
    
    m_amr->allocate(*m_conductivity_particles[idx], pvr_buffer, m_particle_realm);
    m_amr->allocate(*m_rho_dagger_particles[idx],   pvr_buffer, m_particle_realm);
  }
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
  m_amr->allocate(m_particle_E, m_particle_realm, m_phase, SpaceDim);

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

  // Euler needs to limit by relaxation time
  if(m_algorithm == which_algorithm::euler){
    const Real dtRelax = m_relax_factor*m_dt_relax;
    if(dtRelax < a_dt){
      a_dt = dtRelax;
      a_timecode = time_code::relaxation_time;
    }
  }

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

  if(m_algorithm == which_algorithm::semi_implicit){

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
}

void ito_plasma_godunov::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) {
  CH_TIME("ito_plasma_godunov::regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::regrid" << endl;
  }

  if(m_algorithm == which_algorithm::euler){
    ito_plasma_stepper::regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  }
  else{
    this->regrid_si(a_lmin, a_old_finest_level, a_new_finest_level);
  }
}

void ito_plasma_godunov::regrid_si(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) {
  CH_TIME("ito_plasma_godunov::regrid_si");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::regrid_si" << endl;
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
    m_conductivity_particles[idx]->regrid(grids, domains, dx, ref_rat, a_lmin, a_new_finest_level);
    m_rho_dagger_particles[idx]->regrid(grids, domains, dx, ref_rat, a_lmin, a_new_finest_level);
  }


  // Recompute the conductivity with the other particles
  this->compute_regrid_conductivity();

  // Set up semi-implicit poisson again, with particles after diffusion jump (the rho^\dagger)
  this->setup_semi_implicit_poisson(m_prevDt);

  // Recompute the space charge. This deposits the m_rho_dagger particles onto the states. 
  this->compute_regrid_rho();
  
  // Compute the electric field
  const bool converged = this->solve_poisson();
  if(!converged){
    MayDay::Abort("ito_plasma_stepper::regrid - Poisson solve did not converge after regrid!!!");
  }

  // Now let the ito solver deposit its actual particles...
  m_ito->deposit_particles();

  // Recompute new velocities and diffusion coefficients
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_godunov::compute_regrid_conductivity(){
  CH_TIME("ito_plasma_godunov::compute_regrid_conductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_regrid_conductivity" << endl;
  }

  data_ops::set_value(m_conduct_cell, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>&   solver = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();
    
    const int idx = solver_it.get_solver();
    const int q   = species->get_charge();

    if(q != 0 && solver->is_mobile()){
      data_ops::set_value(m_particle_scratch1, 0.0);

      solver->deposit_particles(m_particle_scratch1, *m_conductivity_particles[idx]); // This deposit mu*mass

      // Copy to fluid realm and add to fluid stuff
      m_fluid_scratch1.copy(m_particle_scratch1);
      data_ops::incr(m_conduct_cell, m_fluid_scratch1, Abs(q));
    }
  }
  
  data_ops::scale(m_conduct_cell, units::s_Qe);

  m_amr->average_down(m_conduct_cell, m_fluid_realm, m_phase);
  m_amr->interp_ghost_pwl(m_conduct_cell, m_fluid_realm, m_phase);
  m_amr->interpolate_to_centroids(m_conduct_cell, m_fluid_realm, m_phase);

  // Now do the faces
  this->compute_face_conductivity();
}

void ito_plasma_godunov::compute_regrid_rho(){
  CH_TIME("ito_plasma_godunov::compute_regrid_rho");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_regrid_rho" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>&   solver = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();
    
    solver->deposit_particles(solver->get_state(), *m_rho_dagger_particles[idx]); 
  }
}

void ito_plasma_godunov::regrid_conductivity(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("ito_plasma_godunov::regrid_conductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid_conductivity" << endl;
  }
  
  const Interval interv(0,0);
    
  // Regrid the conductivity
  Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->get_eb_pwl_interp(m_fluid_realm, m_phase);

  // These levels have not changed
  for (int lvl = 0; lvl <= Max(0,a_lmin-1); lvl++){
    m_cache[lvl]->copyTo(*m_conduct_cell[lvl]); // Base level should never change, but ownership can.
  }

  // These levels have changed
  for (int lvl = a_lmin; lvl <= a_new_finest_level; lvl++){
    interpolator[lvl]->interpolate(*m_conduct_cell[lvl], *m_conduct_cell[lvl-1], interv);
    if(lvl <= Min(a_old_finest_level, a_new_finest_level)){
      m_cache[lvl]->copyTo(*m_conduct_cell[lvl]);
    }
  }
}

void ito_plasma_godunov::copy_conductivity_particles(){
  CH_TIME("ito_plasma_godunov::copy_conductivity_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::copy_conductivity_particles" << endl;
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
	List<godunov_particle>& gdnv_parts  = (*m_conductivity_particles[idx])[lvl][dit()].listItems();

	gdnv_parts.clear();

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

void ito_plasma_godunov::copy_rho_dagger_particles(){
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
	List<godunov_particle>& gdnv_parts  = (*m_rho_dagger_particles[idx])[lvl][dit()].listItems();

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

void ito_plasma_godunov::copy_particles_to_scratch(){
  CH_TIME("ito_plasma_godunov::copy_particles_to_scratch");
  if(m_verbosity > 5){
    pout() << m_name + "::copy_particles_to_scratch" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();

    particle_container<ito_particle>& scratch         = solver->get_scratch_particles();
    const particle_container<ito_particle>& particles = solver->get_particles();

    solver->clear(scratch);
    scratch.add_particles(particles);
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
  Real velo_time     = 0.0;
  Real diff_time     = 0.0;

  Real total_time    = -MPI_Wtime();
  
  // Particle algorithms
  particle_time = -MPI_Wtime();
  if(m_algorithm == which_algorithm::euler){
    this->advance_particles_euler(a_dt);
  }
  else if(m_algorithm == which_algorithm::semi_implicit){
    this->advance_particles_si(a_dt);
    m_prevDt = a_dt;
  }
  else if(m_algorithm == which_algorithm::split_semi_implicit){
    this->advance_particles_split_si(a_dt);
    m_prevDt = a_dt;
  }
  else if(m_algorithm == which_algorithm::split_pc){
    this->advance_particles_split_pc(a_dt);
  }
  particle_time += MPI_Wtime();

  // Compute current and relaxation time.
  relax_time = -MPI_Wtime();
  this->compute_J(m_J, a_dt);
  m_dt_relax = this->compute_relaxation_time(); // This is for the restricting the next step.
  relax_time += MPI_Wtime();

  // Move photons
  photon_time = -MPI_Wtime();
  this->advance_photons(a_dt);
  photon_time += MPI_Wtime();

  // If we are using the LEA, we must compute the Ohmic heating term. This must be done
  // BEFORE sorting the particles per cell. 
  if(m_physics->get_coupling() == ito_plasma_physics::coupling::LEA){
    this->compute_EdotJ_source();
  }
  
  // Sort the particles and photons per cell so we can call reaction algorithms
  sort_time = -MPI_Wtime();
  m_ito->sort_particles_by_cell();
  this->sort_bulk_photons_by_cell();
  this->sort_source_photons_by_cell();
  sort_time += MPI_Wtime();

  // Chemistry kernel.
  reaction_time = -MPI_Wtime();
  this->advance_reaction_network(a_dt);
  reaction_time += MPI_Wtime();

  // Make superparticles
  super_time = -MPI_Wtime();
  if((m_step+1) % m_merge_interval == 0 && m_merge_interval > 0){
    m_ito->make_superparticles(m_ppc);
  }
  super_time += MPI_Wtime();

  // Sort particles per patch.
  sort_time -= MPI_Wtime();
  m_ito->sort_particles_by_patch();
  this->sort_bulk_photons_by_patch();
  this->sort_source_photons_by_patch();
  sort_time += MPI_Wtime();

  // Clear other data holders for now. BC comes later...
  clear_time = -MPI_Wtime();
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->clear(solver_it()->get_eb_particles());
    solver_it()->clear(solver_it()->get_domain_particles());
  }
  clear_time += MPI_Wtime();

  // 
  m_ito->deposit_particles();

  // Prepare next step
  velo_time = -MPI_Wtime();
  this->compute_ito_velocities();
  velo_time += MPI_Wtime();
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
  velo_time     *= 100./total_time;
  diff_time     *= 100./total_time;

  if(m_profile){
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
	   << "velo time     = " << velo_time << "%" << endl
      	   << "diff time     = " << diff_time << "%" << endl
	   << "total time    = " << total_time << " (seconds)" << endl
	   << endl;
  }
  
  return a_dt;
}

void ito_plasma_godunov::advance_particles_euler(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advance_particles_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_particles_euler" << endl;
  }

  this->set_old_positions();
  this->advect_particles_euler(a_dt);
  this->diffuse_particles_euler(a_dt);
  
  this->intersect_particles(a_dt);
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->remove_eb_particles();
  }

  m_ito->remap();
  m_ito->deposit_particles();

  this->solve_poisson();
}

void ito_plasma_godunov::advance_particles_si(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advance_particles_si");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_particles_si" << endl;
  }

  Real time_old     = 0.0;
  Real time_setup   = 0.0;
  Real time_diffuse = 0.0;
  Real time_remap   = 0.0;
  Real time_deposit = 0.0;
  Real time_solve   = 0.0;
  Real time_swap    = 0.0;
  Real time_velo    = 0.0;
  Real time_advect  = 0.0;
  Real time_isect   = 0.0;
  Real time_total   = 0.0;

  time_total = -MPI_Wtime();

  // First, advect the current particles out of the domain. Then remove the EB particles and revert the particles
  time_old = -MPI_Wtime();
  this->set_old_positions();
  time_old += MPI_Wtime();

  // Need to copy the current particles because they will be used for computing the conductivity during regrids
  this->copy_conductivity_particles();

  // Compute conductivity and setup poisson
  time_setup = -MPI_Wtime();
  this->compute_conductivity();
  this->setup_semi_implicit_poisson(a_dt);
  time_setup += MPI_Wtime();

  // Diffuse the particles now
  time_diffuse = -MPI_Wtime();
  this->diffuse_particles_euler(a_dt);
  time_diffuse += MPI_Wtime();

  // Remap and deposit, only need to do this for diffusive solvers. 
  time_remap = -MPI_Wtime();
  this->remap_diffusive_particles(); 
  time_remap += MPI_Wtime();
  time_deposit -= MPI_Wtime();
  //  m_ito->deposit_particles();
  this->deposit_diffusive_particles();
  time_deposit += MPI_Wtime();

  // Need to copy the current particles because they will be used for the space charge during regrids. 
  this->copy_rho_dagger_particles();

  // Now compute the electric field
  time_solve -= MPI_Wtime();
  this->solve_poisson();
  time_solve += MPI_Wtime();
  // We have field at k+1 but particles have been diffused. The ones that are diffusive AND mobile are put back to X^k positions
  // and then we compute velocities with E^(k+1). 
  time_swap -= MPI_Wtime();
  this->swap_particle_positions();   // After this, oldPosition() holds X^\dagger, and position() holds X^k. 
  time_remap -= MPI_Wtime();
  this->remap_diffusive_particles(); // Only need to do this for the ones that were diffusive
  time_remap += MPI_Wtime();
  time_swap += MPI_Wtime();
  time_velo -= MPI_Wtime();
  this->compute_ito_velocities();
  time_velo += MPI_Wtime();
  time_advect -= MPI_Wtime();
  this->advect_particles_si(a_dt);  // This 
  time_advect += MPI_Wtime();
  
  // Remap, redeposit, store invalid particles, and intersect particles. Deposition is for relaxation time computation.
  time_remap -= MPI_Wtime();
  this->remap_mobile_or_diffusive_particles();
  time_remap += MPI_Wtime();

  // Do intersection test and remove EB particles. These particles are NOT allowed to react later.
  time_isect -= MPI_Wtime();
  this->intersect_particles(a_dt);
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->remove_eb_particles();
  }
  time_isect += MPI_Wtime();

  time_deposit -= MPI_Wtime();
  this->deposit_mobile_or_diffusive_particles();
  time_deposit += MPI_Wtime();

  time_total += MPI_Wtime();

  time_old     *= 100./time_total;
  time_setup   *= 100./time_total;
  time_diffuse *= 100./time_total;
  time_remap   *= 100./time_total;
  time_deposit *= 100./time_total;
  time_solve   *= 100./time_total;
  time_swap    *= 100./time_total;
  time_velo    *= 100./time_total;
  time_advect  *= 100./time_total;
  time_isect   *= 100./time_total;

  if(m_profile){
    pout() << endl
	   << "ito_plasma_godunov::advance_particles_si breakdown:" << endl
	   << "===================================================" << endl
	   << "set old pos   = " << time_old << "%" << endl
	   << "poisson setup = " << time_setup << "%" << endl
	   << "diffuse part  = " << time_diffuse << "%" << endl
	   << "remap time    = " << time_remap << "%" << endl
	   << "deposit time  = " << time_deposit << "%" << endl
	   << "poisson solve = " << time_solve << "%" << endl
	   << "swap time     = " << time_swap << "%" << endl
      	   << "velo time     = " << time_velo << "%" << endl
	   << "advect time   = " << time_advect << "%" << endl
      	   << "isect time    = " << time_isect << "%" << endl
	   << "total time    = " << time_total << " (seconds)" << endl
	   << endl;
  }
}

void ito_plasma_godunov::advance_particles_split_si(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advance_particles_split_si");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_particles_split_si" << endl;
  }

  // Set old positions and diffuse particles. 
  this->set_old_positions();
  this->diffuse_particles_euler(a_dt);
  
  // Do intersection test and remove EB particles. These particles are removed from the simulation. 
  this->intersect_particles(a_dt);
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->remove_eb_particles();
  }

  // Remap and deposit diffusive particles
  this->remap_diffusive_particles();
  this->deposit_diffusive_particles();

  // Explicit update of the electric field. Then recompute the advective velocities. 
  this->setup_standard_poisson();
  this->solve_poisson();

  // Recompute advective velocities. This is necessary because
  // we want the mobility at E^k
  this->compute_ito_velocities();

  // Compute conductivity and solve semi-implicit Poisson
  this->compute_conductivity();
  this->setup_semi_implicit_poisson(a_dt);
  this->solve_poisson();
  
  // We have field at k+1 and can now compute the velocities at E^(k+1) and move the particles
  // and then we compute velocities with E^(k+1).
  this->set_old_positions();
  this->compute_ito_velocities();
  this->advect_particles_euler(a_dt);  
  
  // Remap, redeposit, store invalid particles, and intersect particles. Deposition is for relaxation time computation.
  this->remap_mobile_particles();

  // Do intersection test and remove EB particles. These particles are NOT allowed to react later.
  this->intersect_particles(a_dt);
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->remove_eb_particles();
  }

  // Do final deposition. This is mostly for computing the current density. 
  this->deposit_mobile_particles();
}

void ito_plasma_godunov::advance_particles_split_pc(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advance_particles_split_pc");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_particles_split_pc" << endl;
  }

  // Set old positions and diffuse particles. 
  this->set_old_positions();
  this->diffuse_particles_euler(a_dt);
  
  // Do intersection test and remove EB particles. These particles are removed from the simulation. 
  this->intersect_particles(a_dt);
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->remove_eb_particles();
  }

  // Remap and deposit diffusive particles
  this->remap_diffusive_particles();
  this->deposit_diffusive_particles();

  // Explicit update of the electric field. Then recompute the advective velocities. 
  this->setup_standard_poisson();
  this->solve_poisson();

  // PREDICTOR-CORRECTOR BEGINS HERE
  this->set_old_positions();
  this->copy_particles_to_scratch();

  const int NPC = 10;
  for (int ipc =0; ipc < NPC; ipc++){
    // Recompute advective velocities. This is necessary because
    // we want the mobility at E^k
    this->compute_ito_velocities();
    this->advect_particles_euler(a_dt);

    // Remove particles that went through EB
    this->intersect_particles(a_dt);
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      solver_it()->remove_eb_particles();
    }

    m_ito->remap();
    m_ito->deposit_particles();

    // Solve Poisson. 
    this->solve_poisson();

    // Put particles back. 
    if(ipc < NPC - 1){
      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
	RefCountedPtr<ito_solver>& solver = solver_it();

	const particle_container<ito_particle>& scratch         = solver->get_scratch_particles();
	particle_container<ito_particle>& particles = solver->get_particles();

	solver->clear(particles);
	particles.add_particles(scratch);
      }
    }
  }
}

void ito_plasma_godunov::advect_particles_euler(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advect_particles_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advect_particles_euler" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    if(solver->is_mobile()){
      
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  // First step
	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];

	    // Update positions. 
	    p.position() += p.velocity()*a_dt;
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::advect_particles_si(const Real a_dt){
  CH_TIME("ito_plasma_godunov::advect_particles_si");
  if(m_verbosity > 5){
    pout() << m_name + "::advect_particles_si" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    if(solver->is_mobile()){
      
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  // First step
	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];

	    // Add diffusion hop again. The position after the diffusion hop is oldPosition() and X^k is in position()
	    const RealVect Xk = p.position();
	    p.position()      = p.oldPosition() + p.velocity()*a_dt;
	    p.oldPosition()   = Xk;
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::diffuse_particles_euler(const Real a_dt){
  CH_TIME("ito_plasma_godunov::diffuse_particles_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::diffuse_particles_euler" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    
    if(solver->is_diffusive()){
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  // Diffusion hop.
	  for (lit.rewind(); lit.ok(); ++lit){
	    ito_particle& p = particleList[lit];
	    const RealVect ran = solver->random_gaussian();
	    const RealVect hop = ran*sqrt(2.0*p.diffusion()*a_dt);

	    p.position() += hop;
	  }
	}
      }
    }
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

void ito_plasma_godunov::rewind_particles(){
  CH_TIME("ito_plasma_godunov::rewind_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::rewind_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    
    if(solver->is_diffusive()){
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  // Diffusion hop.
	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];
	    p.position() = p.oldPosition();
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::swap_particle_positions(){
  CH_TIME("ito_plasma_godunov::swap_particle_positions");
  if(m_verbosity > 5){
    pout() << m_name + "::swap_particle_positions" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();

    // No need to do this if solver is only mobile because the diffusion step didn't change the position.
    // Likewise, if the solver is only diffusive then advect_particles_si routine won't trigger so no need for that either. 
    if(mobile && diffusive){ 
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl          = m_amr->get_grids(m_particle_realm)[lvl];
	ParticleData<ito_particle>& particles = solver->get_particles()[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	  List<ito_particle>& particleList = particles[dit()].listItems();
	  ListIterator<ito_particle> lit(particleList);

	  for (lit.rewind(); lit; ++lit){
	    ito_particle& p = particleList[lit];

	    // We have made a diffusion hop, but we need p.position() to be X^k and p.oldPosition() to be the jumped position. 
	    const RealVect tmp = p.position();
	    
	    p.position()    = p.oldPosition();
	    p.oldPosition() = tmp;
	  }
	}
      }
    }
  }
}

void ito_plasma_godunov::intersect_particles(const Real a_dt){
  CH_TIME("ito_plasma_godunov::intersect_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::intersect_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();

    if(mobile || diffusive){
      solver->intersect_particles();
    }
  }
}

void ito_plasma_godunov::compute_conductivity(){
  CH_TIME("ito_plasma_godunov::compute_conductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_conductivity" << endl;
  }

  ito_plasma_stepper::compute_conductivity(m_conduct_cell);

  // Now do the faces
  this->compute_face_conductivity();
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
