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

  m_avg_cfl = 0.0;
}

ito_plasma_godunov::ito_plasma_godunov(RefCountedPtr<ito_plasma_physics>& a_physics){
  m_name    = "ito_plasma_godunov";
  m_physics = a_physics;
  m_use_old_dt = false;

  m_dt_relax = 1.E99;

  ParmParse pp("ito_plasma_godunov");
  pp.get("particle_realm", m_particle_realm);

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

  if(m_algorithm == which_algorithm::semi_implicit){
    m_use_old_dt = true;
  }
}

void ito_plasma_godunov::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("ito_plasma_godunov::regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_godunov::regrid" << endl;
  }

  if(false){//m_algorithm == which_algorithm::semi_implicit){ 

    this->allocate_internals();

    // Regrid solvers
    m_ito->regrid(a_lmin,     a_old_finest_level, a_new_finest_level);
    m_poisson->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
    m_rte->regrid(a_lmin,     a_old_finest_level, a_new_finest_level);
    m_sigma->regrid(a_lmin,   a_old_finest_level, a_new_finest_level);

    // Deposit particles
    m_ito->deposit_particles();
  }
  else{
    ito_plasma_stepper::regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  }
}

Real ito_plasma_godunov::advance(const Real a_dt) {
  CH_TIME("ito_plasma_godunov::advance");
  if(m_verbosity > 5){
    pout() << m_name + "::advance" << endl;
  }

  // Particle algorithms
  if(m_algorithm == which_algorithm::euler){
    this->advance_particles_euler(a_dt);
  }
  else if(m_algorithm == which_algorithm::semi_implicit){
    this->advance_particles_si(a_dt);
  }

  // Compute current and relaxation time.
  this->compute_J(m_J, a_dt);
  m_dt_relax = this->compute_relaxation_time(); // This is for the restricting the next step. 

  // Move photons
  this->advance_photons(a_dt);

  // Sort the particles and photons per cell so we can call reaction algorithms
  m_ito->sort_particles_by_cell();
  this->sort_bulk_photons_by_cell();
  this->sort_source_photons_by_cell();

  // Chemistry kernel.
  this->advance_reaction_network(a_dt);

  // Make superparticles
  if((m_step+1) % m_merge_interval == 0 && m_merge_interval > 0){
    m_ito->make_superparticles(m_ppc);
  }

  // Sort particles per patch. 
  m_ito->sort_particles_by_patch();
  this->sort_bulk_photons_by_patch();
  this->sort_source_photons_by_patch();

  // Clear other data holders for now. BC comes later
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->clear(solver_it()->get_eb_particles());
    solver_it()->clear(solver_it()->get_domain_particles());
    solver_it()->remove_eb_particles();
  }

  // Prepare next step
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
  
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

  this->set_old_positions();

  // Compute conductivity and setup poissonp
  this->compute_conductivity();
  this->setup_semi_implicit_poisson(a_dt);

  // Diffuse particles
  this->diffuse_particles_euler(a_dt);

  // Remap and deposit
  m_ito->remap();
  m_ito->deposit_particles();

  // Now compute the electric field
  this->solve_poisson();

  // We have field at k+1 but particles have been diffused. Put them back to X^k positions and compute
  // velocities with E^(k+1)
  this->swap_particle_positions();
  this->compute_ito_velocities();
  
  // Compute new ito velocities and advect the particles
  this->advect_particles_si(a_dt);

  // Remap, intersect, and redeposit
  m_ito->remap();
  this->intersect_particles(a_dt);
  m_ito->deposit_particles();
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
	    p.oldPosition() = p.position();
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

	    p.oldPosition() = p.position();
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
