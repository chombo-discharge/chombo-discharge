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

#define DEBUG 1

using namespace physics::ito_plasma;

ito_plasma_godunov::ito_plasma_godunov(){
  m_name = "ito_plasma_godunov";
  m_use_old_dt = false;

  m_dt_relax = 1.E99;

  ParmParse pp("ito_plasma_godunov");
  pp.get("particle_realm", m_particle_realm);
  pp.get("profile", m_profile);

  m_avg_cfl = 0.0;
}

ito_plasma_godunov::ito_plasma_godunov(RefCountedPtr<ito_plasma_physics>& a_physics){
  m_name    = "ito_plasma_godunov";
  m_physics = a_physics;
  m_use_old_dt = false;

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
  else{
    MayDay::Abort("ito_plasma_godunov::parse_options - unknown algorithm requested");
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

  m_amr->allocate(m_J,            m_fluid_realm, m_phase, SpaceDim);
  m_amr->allocate(m_scratch1,     m_fluid_realm, m_phase, 1);
  m_amr->allocate(m_scratch2,     m_fluid_realm, m_phase, 1);
  m_amr->allocate(m_conduct_cell, m_fluid_realm, m_phase, 1);
  m_amr->allocate(m_conduct_face, m_fluid_realm, m_phase, 1);
  m_amr->allocate(m_conduct_eb,   m_fluid_realm, m_phase, 1);
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

void ito_plasma_godunov::post_checkpoint_setup() {
  CH_TIME("ito_plasma_godunov::post_checkpoint_setup");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::post_checkpoint_setup" << endl;
  }

  if(m_algorithm == which_algorithm::euler){ // Default method is just fine. 

  }
  else if(m_algorithm == which_algorithm::semi_implicit) {
    // Strange but true thing. The semi_implicit algorithm needs particles that ended up inside the EB to deposit
    // their charge on the mesh. These were copied to a separate container earlier which was checkpointed. We move
    // those back to the solvers and redeposit on the mesh. 
    //    m_ito->move_invalid_particles_to_particles();
    m_ito->deposit_particles();
    
    // Recompute poisson
    this->solve_poisson();
    this->allocate_internals();

  
    this->compute_ito_velocities();
    this->compute_ito_diffusion();
  }
}

void ito_plasma_godunov::compute_dt(Real& a_dt, time_code& a_timecode){
  CH_TIME("ito_plasma_godunov::compute_dt");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::compute_dt" << endl;
  }
  
  a_dt = m_ito->compute_dt();
  a_dt = a_dt*m_max_cells_hop;
  a_timecode = time_code::cfl;

  if(m_algorithm == which_algorithm::euler){
    const Real dtRelax = m_relax_factor*m_dt_relax;
    if(dtRelax < a_dt){
      a_dt = dtRelax;
      a_timecode = time_code::relaxation_time;
    }
  }
  else if(m_algorithm == which_algorithm::semi_implicit && m_use_old_dt == true){
    // a_dt = m_dt;
    // m_use_old_dt = false;
  }


  if(a_dt < m_min_dt){
    a_dt = m_min_dt;
    a_timecode = time_code::hardcap;
  }

  if(a_dt > m_max_dt){
    a_dt = m_max_dt;
    a_timecode = time_code::hardcap;
  }

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

  // I don't know what this does. 
  if(m_algorithm == which_algorithm::semi_implicit){
    m_use_old_dt = true;

    // Copy conductivity to scratch storage
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();
    m_amr->allocate(m_cache,  m_fluid_realm, m_phase, ncomp);
    for (int lvl = 0; lvl <= a_old_finest_level; lvl++){
      m_conduct_cell[lvl]->localCopyTo(*m_cache[lvl]);
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

  // Allocate new memory
  this->allocate_internals();

  // Regrid solvers
  m_ito->regrid(a_lmin,     a_old_finest_level, a_new_finest_level);
  m_poisson->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  m_rte->regrid(a_lmin,     a_old_finest_level, a_new_finest_level);
  m_sigma->regrid(a_lmin,   a_old_finest_level, a_new_finest_level);

  // Recompute the conductivity
  this->regrid_conductivity(a_lmin, a_old_finest_level, a_new_finest_level);
  this->compute_face_conductivity();

  // Set up semi-implicit poisson again, with particles after diffusion jump (the rho^\dagger)
  this->setup_semi_implicit_poisson(m_prevDt);
  this->deposit_scratch_particles();
  
  // Compute the electric field
  const bool converged = this->solve_poisson();
  if(!converged){
    MayDay::Abort("ito_plasma_stepper::regrid - Poisson solve did not converge after regrid!!!");
  }

  // Recompute new velocities and diffusion coefficients
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
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

void ito_plasma_godunov::deposit_scratch_particles(){
  CH_TIME("ito_plasma_godunov::deposit_scratch_particles");
  if(m_verbosity > 5){
    pout() << m_name + "::deposit_scratch_particles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->deposit_particles(solver_it()->get_state(), solver_it()->get_scratch_particles().get_particles());
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

  // Clear other data holders for now. BC comes later
  clear_time = -MPI_Wtime();
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->clear(solver_it()->get_eb_particles());
    solver_it()->clear(solver_it()->get_domain_particles());
    solver_it()->remove_eb_particles();
  }
  clear_time += MPI_Wtime();

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
  
  m_ito->remap();
  m_ito->deposit_particles();
  this->intersect_particles(a_dt);
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
  Real time_copy    = 0.0;
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

  // Compute conductivity and setup poisson
  time_setup = -MPI_Wtime();
  this->compute_conductivity();
  this->setup_semi_implicit_poisson(a_dt);
  time_setup += MPI_Wtime();

  // Diffuse the particles now
  time_diffuse = -MPI_Wtime();
  this->diffuse_particles_euler(a_dt);
  time_diffuse += MPI_Wtime();

  // Store particles. These are the ones we need when we redo the semi-implicit Poisson solve during regrids.
  time_copy = -MPI_Wtime();
  this->copy_particles_to_scratch();
  time_copy += MPI_Wtime();

  // Remap and deposit
  time_remap = -MPI_Wtime();
  m_ito->remap();
  time_remap += MPI_Wtime();
  time_deposit -= MPI_Wtime();
  m_ito->deposit_particles();
  time_deposit += MPI_Wtime();


  // Now compute the electric field
  time_solve -= MPI_Wtime();
  this->solve_poisson();
  time_solve += MPI_Wtime();
  // We have field at k+1 but particles have been diffused. Put them back to X^k positions and compute
  // velocities with E^(k+1). Then move them.
  time_swap -= MPI_Wtime();
  this->swap_particle_positions();
  time_swap += MPI_Wtime();
  time_velo -= MPI_Wtime();
  this->compute_ito_velocities();
  time_velo += MPI_Wtime();
  time_advect -= MPI_Wtime();
  this->advect_particles_si(a_dt);
  time_advect += MPI_Wtime();
  
  // Remap, redeposit, store invalid particles, and intersect particles. Deposition is for relaxation time computation.
  time_remap -= MPI_Wtime();
  m_ito->remap();
  time_remap += MPI_Wtime();

  time_deposit -= MPI_Wtime();
  m_ito->deposit_particles();
  time_deposit += MPI_Wtime();

  time_isect -= MPI_Wtime();
  this->intersect_particles(a_dt);   // This moves particles that bumped into the EB into data holders for parsing BCs.
  time_isect += MPI_Wtime();

  time_total += MPI_Wtime();

  time_old     *= 100./time_total;
  time_setup   *= 100./time_total;
  time_diffuse *= 100./time_total;
  time_copy    *= 100./time_total;
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
	   << "copy part     = " << time_copy << "%" << endl
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
	    //	    p.oldPosition() = p.position();
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

	    // Add diffusion hop again. The position after the diffusion hop is oldPosition(). Go there,
	    // then advect. 
	    p.position() = p.oldPosition() + p.velocity()*a_dt;
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

	    //	    p.oldPosition() = p.position();
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
	    const RealVect tmp = p.position();
	    
	    p.position()    = p.oldPosition();
	    p.oldPosition() = p.position();
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
     solver->intersect_particles();
  }
}

void ito_plasma_godunov::compute_conductivity(){
  CH_TIME("ito_plasma_godunov::compute_conductivity");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_conductivity" << endl;
  }

  // Get handle to E gas-side E
  EBAMRCellData Egas;
  m_amr->allocate_ptr(Egas);
  m_amr->alias(Egas, m_phase, m_poisson->get_E());

  // Compute |E| and reset conductivity
  data_ops::vector_length(m_scratch1, Egas);
  data_ops::set_value(m_conduct_cell, 0.0);
  
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver   = solver_it();
    const RefCountedPtr<ito_species>& species = solver->get_species();

    const int q = species->get_charge();
    
    if(solver->is_mobile() &&  q != 0){

      // These things are defined on the particle realm. Copy them to the fluid realm. 
      const EBAMRCellData& velo  = solver->get_velo_cell();
      const EBAMRCellData& state = solver->get_state();

      m_fluid_scratch1.copy(state);
      m_fluid_scratchD.copy(velo);

      data_ops::vector_length(m_scratch2, m_fluid_scratchD);   // m_scratch2     = |v|
      data_ops::divide_scalar(m_scratch2, m_scratch1);         // m_scratch2     = |v|/|E| = mu
      data_ops::multiply(m_scratch2,      m_fluid_scratch1);   // m_scratch2     = mu*phi
      data_ops::incr(m_conduct_cell,      m_scratch2, Abs(q)); // m_conduct_cell = mu*phi*q
    }
  }

  // Scale by unit charge
  data_ops::scale(m_conduct_cell, units::s_Qe);

  m_amr->average_down(m_conduct_cell,     m_fluid_realm, m_phase);
  m_amr->interp_ghost_pwl(m_conduct_cell, m_fluid_realm, m_phase);

  // See if this helps...
  m_amr->interpolate_to_centroids(m_conduct_cell, m_fluid_realm, m_phase);


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
