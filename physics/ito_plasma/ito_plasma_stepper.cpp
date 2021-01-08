/*!
  @file   ito_plasmpa_stepper.cpp
  @brief  Implementation of ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
*/

#include "ito_plasma_stepper.H"
#include "data_ops.H"
#include "units.H"
#include "poisson_multifluid_gmg.H"

#include <EBArith.H>
#include <PolyGeom.H>

using namespace physics::ito_plasma;

ito_plasma_stepper::ito_plasma_stepper(){
  m_verbosity = -1;
  m_name      = "ito_plasma_stepper";
  m_phase     = phase::gas;

  m_dt   = 0.0;
  m_time = 0.0;
  
  m_halo_buffer = 1;
  m_pvr_buffer  = 0;
  m_load_ppc    = 1.0;

  m_nwo_reactions         = true;
  m_regrid_superparticles = true;

  m_fluid_realm    = realm::primal;
  m_particle_realm = realm::primal;
}

ito_plasma_stepper::ito_plasma_stepper(RefCountedPtr<ito_plasma_physics>& a_physics) : ito_plasma_stepper(){
  m_physics   = a_physics;
}

ito_plasma_stepper::~ito_plasma_stepper(){
}

void ito_plasma_stepper::set_verbosity(const int a_verbosity){
  CH_TIME("ito_plasma_stepper::set_verbosity");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_verbosity" << endl;
  }
  m_verbosity = a_verbosity;
}

void ito_plasma_stepper::setup_solvers(){
  CH_TIME("ito_plasma_stepper::setup_solver");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_solvers" << endl;
  }

  // Parse class options
  this->parse_options();

  // Set up solvers
  this->setup_ito();
  this->setup_poisson();
  this->setup_rte();
  this->setup_sigma();

  // Allocate internal stuff
  this->allocate_internals();

  // Set the particle buffers for the Ito solver
  this->set_particle_buffers();
}

void ito_plasma_stepper::setup_ito(){
  CH_TIME("ito_plasma_stepper::setup_ito");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_ito" << endl;
  }

  m_ito->set_verbosity(m_verbosity);
  m_ito->parse_options();
  m_ito->set_amr(m_amr);
  m_ito->set_phase(m_phase);
  m_ito->set_computational_geometry(m_compgeom);
  m_ito->set_realm(m_particle_realm);
}

void ito_plasma_stepper::setup_poisson(){
  CH_TIME("ito_plasma_stepper::setup_poisson");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_poisson" << endl;
  }

  m_poisson->set_verbosity(m_verbosity);
  m_poisson->parse_options();
  m_poisson->set_amr(m_amr);
  m_poisson->set_computational_geometry(m_compgeom);
  m_poisson->set_potential(m_potential); // Needs to happen AFTER set_poisson_wall_func
  m_poisson->set_realm(m_fluid_realm);
  m_poisson->sanity_check();
}

void ito_plasma_stepper::setup_rte(){
  CH_TIME("ito_plasma_stepper::setup_rte");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_rte" << endl;
  }

  m_rte->set_verbosity(m_verbosity);
  m_rte->parse_options();
  m_rte->set_phase(m_phase);
  m_rte->set_amr(m_amr);
  m_rte->set_computational_geometry(m_compgeom);
  m_rte->set_realm(m_particle_realm);
  m_rte->sanity_check();
}

void ito_plasma_stepper::setup_sigma(){
  CH_TIME("ito_plasma_stepper::setup_sigma");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_sigma" << endl;
  }

  m_sigma = RefCountedPtr<sigma_solver> (new sigma_solver());
  m_sigma->set_amr(m_amr);
  m_sigma->set_verbosity(m_verbosity);
  m_sigma->set_computational_geometry(m_compgeom);
  m_sigma->set_realm(m_fluid_realm);
}

void ito_plasma_stepper::set_particle_buffers(){
  CH_TIME("ito_plasma_stepper::set_particle_buffers");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_particle_buffers" << endl;
  }
  
  m_ito->set_halo_buffer(m_halo_buffer);
  m_ito->set_pvr_buffer(m_pvr_buffer);

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();

    solver->set_halo_buffer(m_halo_buffer);
    solver->set_pvr_buffer(m_pvr_buffer);
  }
}

void ito_plasma_stepper::allocate() {
  CH_TIME("ito_plasma_stepper::allocate");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::allocate" << endl;
  }

  m_ito->allocate_internals();
  m_rte->allocate_internals();
  m_poisson->allocate_internals();
  m_sigma->allocate_internals();
}

void ito_plasma_stepper::post_initialize(){

}

void ito_plasma_stepper::initial_data(){
  CH_TIME("ito_plasma_stepper::initial_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::initial_data" << endl;
  }

  m_ito->initial_data(); // This deposits, of course. 
  m_rte->initial_data();
  this->initial_sigma();

  m_ito->sort_particles_by_cell( ito_solver::which_container::bulk);
  m_ito->make_superparticles(    ito_solver::which_container::bulk, m_ppc);
  m_ito->sort_particles_by_patch(ito_solver::which_container::bulk);
  
  // Solve Poisson equation and compute the E-field
  this->solve_poisson();

  // Fill solvers with velocities and diffusion coefficients
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_stepper::initial_sigma(){
  CH_TIME("ito_plasma_stepper::initial_sigma");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::initial_sigma" << endl;
  }

  const RealVect origin  = m_amr->get_prob_lo();
  const int finest_level = m_amr->get_finest_level();

  EBAMRIVData& sigma = m_sigma->get_state();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->get_grids(m_fluid_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_fluid_realm, phase::gas)[lvl];
    const Real dx                = m_amr->get_dx()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      BaseIVFAB<Real>& state = (*sigma[lvl])[dit()];

      const EBISBox& ebisbox = ebisl[dit()];
      const IntVectSet& ivs  = state.getIVS();
      const EBGraph& ebgraph = state.getEBGraph();
      
      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();
	const RealVect pos  = origin + vof.gridIndex()*dx + 0.5*ebisbox.bndryCentroid(vof)*dx;
	
	for (int comp = 0; comp < state.nComp(); comp++){
	  state(vof, comp) = m_physics->initial_sigma(m_time, pos);
	}
      }
    }
  }

  m_amr->average_down(sigma, m_fluid_realm, phase::gas);
  m_sigma->reset_cells(sigma);
}

void ito_plasma_stepper::post_checkpoint_setup(){
  CH_TIME("ito_plasma_stepper::post_checkpoint_setup");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::post_checkpoint_setup" << endl;
  }

  //this->solve_poisson();
  this->allocate_internals();

  m_ito->remap();
  
  this->post_checkpoint_poisson();
  
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_stepper::post_checkpoint_poisson(){
  CH_TIME("ito_plasma_stepper::post_checkpoint_poisson");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::post_checkpoint_poisson" << endl;
  }

  // Do some post checkpointing stuff. 
  m_poisson->post_checkpoint();
  
  // Do ghost cells and then compute E
  MFAMRCellData& state = m_poisson->get_state();

  m_amr->average_down(state, m_fluid_realm);
  m_amr->interp_ghost(state, m_fluid_realm);

  m_poisson->compute_E();    // Solver checkpoints the potential. Now compute the field.

  // Interpolate the fields to centroids
  EBAMRCellData E;
  m_amr->allocate_ptr(E);
  m_amr->alias(E, m_phase, m_poisson->get_E());
  
  // Fluid realm
  m_fluid_E.copy(E);
  m_amr->average_down(m_fluid_E, m_fluid_realm, m_phase);
  m_amr->interp_ghost_pwl(m_fluid_E, m_fluid_realm, m_phase);
  m_amr->interpolate_to_centroids(m_fluid_E, m_fluid_realm, m_phase);

  // Particle realm
  m_particle_E.copy(E);
  m_amr->average_down(m_particle_E, m_particle_realm, m_phase);
  m_amr->interp_ghost_pwl(m_particle_E, m_particle_realm, m_phase);
  m_amr->interpolate_to_centroids(m_particle_E, m_fluid_realm, m_phase);

  // Compute maximum E
  // const Real Emax = this->compute_Emax(m_phase);
  // std::cout << Emax << std::endl;
}

void ito_plasma_stepper::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const {
  CH_TIME("ito_plasma_stepper::write_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_checkpoint_data" << endl;
  }

  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver = solver_it();
    solver->write_checkpoint_level(a_handle, a_lvl);
  }

  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<mc_photo>& solver = solver_it();
    solver->write_checkpoint_level(a_handle, a_lvl);
  }

  m_poisson->write_checkpoint_level(a_handle, a_lvl);
  m_sigma->write_checkpoint_level(a_handle, a_lvl);
}

void ito_plasma_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){
  CH_TIME("ito_plasma_stepper::read_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::read_checkpoint_data" << endl;
  }

  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    solver->read_checkpoint_level(a_handle, a_lvl);
  }

  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    solver->read_checkpoint_level(a_handle, a_lvl);
  }

  m_poisson->read_checkpoint_level(a_handle, a_lvl);
  m_sigma->read_checkpoint_level(a_handle, a_lvl);
}

void ito_plasma_stepper::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {
  CH_TIME("ito_plasma_stepper::write_plot_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_plot_data" << endl;
  }

  // Poisson solver copies over its output data
  a_plotvar_names.append(m_poisson->get_plotvar_names());
  m_poisson->write_plot_data(a_output, a_icomp);

  // Surface charge solver writes
  a_plotvar_names.append(m_sigma->get_plotvar_names());
  m_sigma->write_plot_data(a_output, a_icomp);

  // Ito solvers copy their output data
  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    a_plotvar_names.append(solver->get_plotvar_names());
    solver->write_plot_data(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    a_plotvar_names.append(solver->get_plotvar_names());
    solver->write_plot_data(a_output, a_icomp);
  }

  // Write the current to the output
  this->write_J(a_output, a_icomp);
  a_plotvar_names.push_back("x-J");
  a_plotvar_names.push_back("y-J");
  if(SpaceDim == 3){
    a_plotvar_names.push_back("z-J");
  }

  // Write the number of particles per patch
  this->write_num_particles_per_patch(a_output, a_icomp);
  a_plotvar_names.push_back("particles_per_patch");
}

void ito_plasma_stepper::write_J(EBAMRCellData& a_output, int& a_icomp) const{
  CH_TIME("ito_plasma_stepper::write_J");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_J" << endl;
  }

  const Interval src_interv(0, SpaceDim-1);
  const Interval dst_interv(a_icomp, a_icomp + SpaceDim -1);

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    if(m_J.get_realm() == a_output.get_realm()){
      m_J[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else {
      m_J[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }
  a_icomp += SpaceDim;
}

void ito_plasma_stepper::write_num_particles_per_patch(EBAMRCellData& a_output, int& a_icomp) const {
  CH_TIME("ito_plasma_stepper::write_num_particles_per_patch");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_num_particles_per_patch" << endl;
  }

  const Interval src_interv(0, 0);
  const Interval dst_interv(a_icomp, a_icomp);

  data_ops::set_value(m_particle_scratch1, 0.0);
  
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const particle_container<ito_particle>& particles = solver_it()->get_particles(ito_solver::which_container::bulk);

    for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[lvl];
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	(*m_particle_scratch1[lvl])[dit()] += particles[lvl][dit].numItems();
      }
    }
  }

  // Copy to output holder
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    if(m_particle_scratch1.get_realm() == a_output.get_realm()){
      m_particle_scratch1[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else {
      m_particle_scratch1[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }
  
  a_icomp += 1;
}

void ito_plasma_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("ito_plasma_stepper::synchronize_solver_times");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::synchronize_solver_times" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;

  m_ito->set_time(a_step,     a_time, a_dt);
  m_poisson->set_time(a_step, a_time, a_dt);
  m_rte->set_time(a_step,     a_time, a_dt);
  m_sigma->set_time(a_step,   a_time, a_dt);
}

void ito_plasma_stepper::print_step_report(){
  CH_TIME("ito_plasma_stepper::print_step_report");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::print_step_report" << endl;
  }

  const Real Emax = this->compute_Emax(m_phase);
  
  const size_t l_particles        = m_ito->get_num_particles(ito_solver::which_container::bulk, true);
  const size_t g_particles        = m_ito->get_num_particles(ito_solver::which_container::bulk, false);
  
  const size_t l_eb_particles     = m_ito->get_num_particles(ito_solver::which_container::eb, true);
  const size_t g_eb_particles     = m_ito->get_num_particles(ito_solver::which_container::eb, false);
  
  const size_t l_domain_particles = m_ito->get_num_particles(ito_solver::which_container::domain, true);
  const size_t g_domain_particles = m_ito->get_num_particles(ito_solver::which_container::domain, false);

  const size_t l_source_particles = m_ito->get_num_particles(ito_solver::which_container::source, true);
  const size_t g_source_particles = m_ito->get_num_particles(ito_solver::which_container::source, false);

  Real avg;
  Real sigma;
  
  int minRank;
  int maxRank;
  
  size_t minParticles;
  size_t maxParticles;
  
  // Compute some particle statistics
  this->get_particle_statistics(avg, sigma, minParticles, maxParticles, minRank, maxRank);

  // How was the time step restricted
  std::string str;
  switch(m_timecode){
  case time_code::physics:
    str = "dt restricted by 'physics'";
    break;
  case time_code::cfl:
    str = "dt restricted by 'cfl'";
    break;
  case time_code::relaxation_time:
    str = "dt restricted by 'relaxation time'";
    break;
  case time_code::hardcap:
    str = "dt restricted by 'hardcap'";
    break;
  default:
    str = "dt restricted by 'unspecified'";
    break;
  }
  pout() << "                                   " + str << endl;
  pout() << "                                   Emax      = " << Emax << endl
	 << "                                   #part     = " << l_particles << " (" << g_particles << ")" << endl
	 << "                                   #eb part  = " << l_eb_particles << " (" << g_eb_particles << ")" << endl
	 << "                                   #dom part = " << l_domain_particles << " (" << g_domain_particles << ")" << endl
    	 << "                                   #src part = " << l_source_particles << " (" << g_source_particles << ")" << endl
	 << "                                   #part min = " << minParticles << " (on rank = " << minRank << ")" << endl
    	 << "                                   #part max = " << maxParticles << " (on rank = " << maxRank << ")" << endl
    	 << "                                   #part avg = " << avg << endl
	 << "                                   #part dev = " << sigma << " (" << 100.*sigma/avg << "%)" << endl;
}

void ito_plasma_stepper::get_particle_statistics(Real& a_avg, Real& a_sigma, size_t& a_minPart, size_t& a_maxPart, int& a_minRank, int& a_maxRank){
  CH_TIME("ito_plasma_stepper::get_particle_statistics");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::get_particle_statistics" << endl;
  }

  const int srcProc = 0; 
  const int nProc   = numProc();

  const size_t numLocal = m_ito->get_num_particles(ito_solver::which_container::bulk, true);

  // Gather on source proc
  Vector<size_t> allCounts(nProc);
  gather(allCounts, numLocal, srcProc);


  // Compute and broadcast the average and the standard deviation
  if(procID() == srcProc){

    a_avg = 0.0;
    for (int i = 0; i < nProc; i++){
      a_avg += 1.0*allCounts[i];
    }
    a_avg *= 1./nProc;

    a_sigma = 0.0;
    for (int i = 0; i < nProc; i++){
      a_sigma += std::pow(1.0*allCounts[i]-a_avg,2);
    }
    a_sigma = sqrt(a_sigma/nProc);
  }
  
  broadcast(a_avg,   srcProc);
  broadcast(a_sigma, srcProc);

  
  // Get the minimum/maximum number of particles
  a_minRank = srcProc;
  a_maxRank = srcProc;
  
  a_minPart = std::numeric_limits<size_t>::max();
  a_maxPart = 0;

  if(procID() == srcProc){
    for (int i = 0; i < nProc; i++){
      if(allCounts[i] < a_minPart){
	a_minPart = allCounts[i];
	a_minRank = i;
      }

      if(allCounts[i] > a_maxPart){
	a_maxPart = allCounts[i];
	a_maxRank = i;
      }
    }
  }

  broadcast(a_minRank, srcProc);
  broadcast(a_maxRank, srcProc);
  
  broadcast(a_minPart, srcProc);
  broadcast(a_maxPart, srcProc);
}

void ito_plasma_stepper::print_timer_diagnostics(Real& a_timer, const std::string a_prefix){
  CH_TIME("ito_plasma_stepper::print_timer_diagnostics");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::print_timer_diagnostics" << endl;
  }

  const int srcProc = 0; 
  const int nProc   = numProc();

  const size_t numLocal = m_ito->get_num_particles(ito_solver::which_container::bulk, true);

  // Gather all timers on source proc
  Vector<Real> allTimers(nProc);
  gather(allTimers, a_timer, srcProc);

  Real averageTime;
  Real sigmaTime;

  Real minTime;
  Real maxTime;

  int minRank;
  int maxRank;

  // Compute and broadcast the average and the standard deviation
  if(procID() == srcProc){

    averageTime = 0.0;
    for (int i = 0; i < nProc; i++){
      averageTime += allTimers[i];
    }
    averageTime *= 1./nProc;

    sigmaTime = 0.0;
    for (int i = 0; i < nProc; i++){
      sigmaTime += std::pow(1.0*allTimers[i]-averageTime,2);
    }
    sigmaTime = sqrt(sigmaTime/nProc);
  }
  
  broadcast(averageTime, srcProc);
  broadcast(sigmaTime,   srcProc);
  
  // Get the minimum/maximum number of particles
  minRank = srcProc;
  maxRank = srcProc;
  
  minTime = std::numeric_limits<Real>::max();
  maxTime = 0;

  if(procID() == srcProc){
    for (int i = 0; i < nProc; i++){
      if(allTimers[i] < minTime){
	minTime = allTimers[i];
	minRank = i;
      }

      if(allTimers[i] > maxTime){
	maxTime = allTimers[i];
	maxRank = i;
      }
    }
  }

  broadcast(minRank, srcProc);
  broadcast(maxRank, srcProc);
  
  broadcast(minTime, srcProc);
  broadcast(maxTime, srcProc);

  // Fix formatting for the various fields
  std::stringstream ss_locTime;
  std::stringstream ss_minTime;
  std::stringstream ss_maxTime;
  std::stringstream ss_avgTime;
  std::stringstream ss_sigTime;

  ss_locTime << std::fixed << std::setprecision(2) << a_timer;
  ss_minTime << std::fixed << std::setprecision(2) << minTime;
  ss_maxTime << std::fixed << std::setprecision(2) << maxTime;
  ss_avgTime << std::fixed << std::setprecision(2) << averageTime;
  ss_sigTime << std::fixed << std::setprecision(2) << sigmaTime;


  pout() << std::left << std::setw(25) << a_prefix
	 << " | " << std::right << std::setw(8) << ss_locTime.str()
    	 << " | " << std::right << std::setw(8) << ss_minTime.str()
    	 << " | " << std::right << std::setw(8) << ss_maxTime.str()
    	 << " | " << std::right << std::setw(8) << ss_avgTime.str()
	 << " | " << std::right << std::setw(8) << ss_sigTime.str()
    	 << " | " << std::right << std::setw(8) << minRank
	 << " | " << std::right << std::setw(8) << maxRank
    	 << " | " << endl;
}

void ito_plasma_stepper::print_timer_head(){
  pout() << "--------------------------------------------------------------------------------------------------------"
	 << endl
	 << std::left << std::setw(25) << "Kernel"
    	 << " | " << std::right << std::setw(8) << "Loc."
    	 << " | " << std::right << std::setw(8) << "Min."
	 << " | " << std::right << std::setw(8) << "Max."
    	 << " | " << std::right << std::setw(8) << "Avg." 
	 << " | " << std::right << std::setw(8) << "Dev."
	 << " | " << std::right << std::setw(8)  << "Min rank"
	 << " | " << std::right << std::setw(8)  << "Max rank"
	 << " | " << endl
  	 << "-------------------------------------------------------------------------------------------------------|"
	 << endl;
}

void ito_plasma_stepper::print_timer_tail(){
  pout() << "--------------------------------------------------------------------------------------------------------\n";
}

void ito_plasma_stepper::parse_filters(){
  CH_TIME("ito_plasma_stepper::compute_dt");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::parse_filters" << endl;
  }

  ParmParse pp(m_name.c_str());

  // Build filters. Always uses a compensation step.
  for (int i = 0; i < 100; i++){

    const int ndigits = round(log10(1.0 + 1.0*i));
    char* cstr = new char[1+ndigits];
    sprintf(cstr, "%d", 1 + i);

    const std::string str = "filter_" + std::string(cstr);

    if(pp.contains(str.c_str())){
      Real alpha;
      int stride;
      int N;
      bool comp;
    
      pp.get(str.c_str(), alpha,    0);
      pp.get(str.c_str(), stride,   1);
      pp.get(str.c_str(), N,        2);
      pp.get(str.c_str(), comp,     3);

      m_filters.emplace_front(alpha, stride, N);
      if(comp){
	const Real alphaC = (N+1) - N*alpha;
	m_filters.emplace_front(alphaC, stride, 1);
      }
    }

    delete cstr;
  }
}

void ito_plasma_stepper::compute_dt(Real& a_dt, time_code& a_timecode){
  CH_TIME("ito_plasma_stepper::compute_dt");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_dt" << endl;
  }
  
  a_dt = m_ito->compute_dt();
  a_dt = a_dt*m_max_cells_hop;
  a_timecode = time_code::cfl;
  
}

void ito_plasma_stepper::register_realms(){
  CH_TIME("ito_plasma_stepper::register_realms");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::register_realms" << endl;
  }

  m_amr->register_realm(m_fluid_realm);
  m_amr->register_realm(m_particle_realm);
}

void ito_plasma_stepper::register_operators(){
  CH_TIME("ito_plasma_stepper::register_operators");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::register_operators" << endl;
  }

  m_ito->register_operators();
  m_poisson->register_operators();
  m_rte->register_operators();
  m_sigma->register_operators();
}
  
void ito_plasma_stepper::pre_regrid(const int a_lmin, const int a_old_finest_level){
  CH_TIME("ito_plasma_stepper::pre_regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::pre_regrid" << endl;
  }

  m_ito->pre_regrid(a_lmin,     a_old_finest_level);
  m_poisson->pre_regrid(a_lmin, a_old_finest_level);
  m_rte->pre_regrid(a_lmin,     a_old_finest_level);
  m_sigma->pre_regrid(a_lmin,   a_old_finest_level);
}

void ito_plasma_stepper::deallocate(){
  CH_TIME("ito_plasma_stepper::deallocate");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::deallocate" << endl;
  }

  // Don't deallocate anything. 
}

void ito_plasma_stepper::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("ito_plasma_stepper::regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::regrid" << endl;
  }

  // Allocate new memory
  this->allocate_internals();

  // Regrid solvers
  m_ito->regrid(a_lmin,     a_old_finest_level, a_new_finest_level);
  m_poisson->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  m_rte->regrid(a_lmin,     a_old_finest_level, a_new_finest_level);
  m_sigma->regrid(a_lmin,   a_old_finest_level, a_new_finest_level);

  if(m_regrid_superparticles){
    m_ito->sort_particles_by_cell( ito_solver::which_container::bulk);
    m_ito->make_superparticles(    ito_solver::which_container::bulk, m_ppc);
    m_ito->sort_particles_by_patch(ito_solver::which_container::bulk);
  }

  // Redeposit particles
  m_ito->deposit_particles();

  // Recompute the electric field
  const bool converged = this->solve_poisson();
  if(!converged){
    MayDay::Abort("ito_plasma_stepper::regrid - Poisson solve did not converge after regrid!!!");
  }

  // Recompute new velocities and diffusion coefficients
  this->compute_ito_velocities();
  this->compute_ito_diffusion();
}

void ito_plasma_stepper::post_regrid(){

}

int ito_plasma_stepper::get_num_plot_vars() const {
  CH_TIME("ito_plasma_stepper::get_num_plot_vars");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::get_num_plot_vars" << endl;
  }

  int ncomp = 0;
  
  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    ncomp += solver->get_num_plotvars();
  }
  
  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    ncomp += solver->get_num_plotvars();
  }

  ncomp += m_poisson->get_num_plotvars();
  ncomp += m_sigma->get_num_plotvars();
  ncomp += SpaceDim; // For plotting the current density
  ncomp += 1;        // For plotting the number of particles per cell

  return ncomp;
}

void ito_plasma_stepper::set_ito(RefCountedPtr<ito_layout<ito_solver> >& a_ito){
  CH_TIME("ito_plasma_stepper::set_ito");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_ito" << endl;
  }

  m_ito = a_ito;
}

void ito_plasma_stepper::set_poisson(RefCountedPtr<poisson_solver>& a_poisson){
  CH_TIME("ito_plasma_stepper::set_poisson");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_poisson" << endl;
  }

  m_poisson = a_poisson;
}

void ito_plasma_stepper::set_rte(RefCountedPtr<rte_layout<mc_photo> >& a_rte){
  CH_TIME("ito_plasma_stepper::set_rte");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_rte" << endl;
  }
  
  m_rte = a_rte;
}

void ito_plasma_stepper::set_potential(Real (*a_potential)(const Real a_time)){
  CH_TIME("ito_plasma_stepper::set_potential");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_potential" << endl;
  }

  m_potential = a_potential;
}

Real ito_plasma_stepper::compute_Emax(const phase::which_phase a_phase) {
  CH_TIME("ito_plasma_stepper::compute_Emax");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_Emax" << endl;
  }

  // Get a handle to the E-field
  EBAMRCellData Ephase;
  m_amr->allocate_ptr(Ephase);
  m_amr->alias(Ephase, m_phase, m_poisson->get_E());

  // Interpolate to centroids
  EBAMRCellData E;
  m_amr->allocate(E, m_fluid_realm, m_phase, SpaceDim);
  data_ops::copy(E, Ephase);
  m_amr->interpolate_to_centroids(E, m_fluid_realm, m_phase);

  Real max, min;
  data_ops::get_max_min_norm(max, min, E);

  return max;
}

Real ito_plasma_stepper::get_time() const{
  CH_TIME("ito_plasma_stepper::get_time");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::get_time" << endl;
  }

  return m_time;
}

void ito_plasma_stepper::compute_E(MFAMRCellData& a_E, const MFAMRCellData& a_potential){
  CH_TIME("ito_plasma_stepper::compute_E(mfamrcell,mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(mfamrcell, mfamrcell" << endl;
  }

  m_amr->compute_gradient(a_E, a_potential, m_fluid_realm);
  data_ops::scale(a_E, -1.0);

  m_amr->average_down(a_E, m_fluid_realm);
  m_amr->interp_ghost(a_E, m_fluid_realm);

}

void ito_plasma_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase){
  CH_TIME("ito_plasma_stepper::compute_E(ebamrcell, phase)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamrcell, phase" << endl;
  }

  this->compute_E(a_E, a_phase, m_poisson->get_state());
}

void ito_plasma_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential){
  CH_TIME("ito_plasma_stepper::compute_E(ebamrcell, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamrcell, phase mfamrcell" << endl;
  }
  
  EBAMRCellData pot_gas;
  m_amr->allocate_ptr(pot_gas);
  m_amr->alias(pot_gas, a_phase, a_potential);

  m_amr->compute_gradient(a_E, pot_gas, m_fluid_realm, a_phase);
  data_ops::scale(a_E, -1.0);

  m_amr->average_down(a_E, m_fluid_realm, a_phase);
  m_amr->interp_ghost(a_E, m_fluid_realm, a_phase);
}

void ito_plasma_stepper::compute_E(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("ito_plasma_stepper::compute_E(ebamrflux, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamrflux, phase mfamrcell" << endl;
  }

  CH_assert(a_E_face[0]->nComp() == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){

    const DisjointBoxLayout& dbl = m_amr->get_grids(m_fluid_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->get_ebisl(m_fluid_realm, a_phase)[lvl];
    const ProblemDomain& domain  = m_amr->get_domains()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& E_cell = (*a_E_cell[lvl])[dit()];
      const EBISBox& ebisbox  = ebisl[dit()];
      const EBGraph& ebgraph  = ebisbox.getEBGraph();
      const Box& box          = dbl.get(dit());
      
      for (int dir = 0; dir < SpaceDim; dir++){
      	EBFaceFAB& E_face = (*a_E_face[lvl])[dit()][dir];
	E_face.setVal(0.0);

      	EBLevelDataOps::averageCellToFace(E_face,
      					  E_cell,
      					  ebgraph,
      					  box,
      					  0,
      					  dir,
      					  domain,
      					  dir,
      					  dir);
      }

    }
    a_E_face[lvl]->exchange();
  }
}

void ito_plasma_stepper::compute_E(EBAMRIVData& a_E_eb,  const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("ito_plasma_stepper::compute_E(ebamriv, phase, ebamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamriv, phase ebamrcell)" << endl;
  }

  CH_assert(a_E_eb[0]->nComp()   == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const irreg_amr_stencil<eb_centroid_interp>& interp_stencil = m_amr->get_eb_centroid_interp_stencils(m_fluid_realm, a_phase);
  interp_stencil.apply(a_E_eb, a_E_cell);
}

void ito_plasma_stepper::compute_rho(){
  CH_TIME("ito_plasma_stepper::compute_rho()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_rho()" << endl;
  }
  
  this->compute_rho(m_poisson->get_source(), m_ito->get_densities());
}

void ito_plasma_stepper::compute_rho(MFAMRCellData& a_rho, const Vector<EBAMRCellData*>&  a_densities){
  CH_TIME("ito_plasma_stepper::compute_rho(rho, densities)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_rho(rho, densities)" << endl;
  }

  // TLDR: a_densities is from the ito solvers so it is defined over the particle realm. But a_rho is defined over
  //       the fluid realm so we need scratch storage we can copy into. We use m_fluid_scratch1 for that. 

  // Reset
  data_ops::set_value(a_rho, 0.0);

  // Make alias
  EBAMRCellData rhoPhase;
  m_amr->allocate_ptr(rhoPhase);
  m_amr->alias(rhoPhase, m_phase, a_rho);

  // Increment each solver
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver   = solver_it();
    const RefCountedPtr<ito_species>& species = solver->get_species();
    const int idx = solver_it.get_solver();
    const int q   = species->get_charge();

    if(species->get_charge() != 0){
      m_fluid_scratch1.copy(*a_densities[idx]);
      data_ops::incr(rhoPhase, m_fluid_scratch1, q);
    }
  }

  data_ops::scale(a_rho, units::s_Qe);

  m_amr->average_down(a_rho, m_fluid_realm);
  m_amr->interp_ghost(a_rho, m_fluid_realm);


  // Add potential filters.
  if(false){//m_filter_rho){
    for (const auto& f : m_filters){
      const Real alpha  = std::get<0>(f);
      const int stride  = std::get<1>(f);
      const int num_app = std::get<2>(f);

      for (int iapp = 0; iapp < num_app; iapp++){
	data_ops::set_value(m_fluid_scratch1, 0.0);
	m_fluid_scratch1.copy(rhoPhase);
	data_ops::set_covered_value(m_fluid_scratch1, 0.0, 0);
	data_ops::filter_smooth(rhoPhase, m_fluid_scratch1, stride, alpha);

	m_amr->average_down(rhoPhase, m_fluid_realm, m_phase);
	m_amr->interp_ghost(rhoPhase, m_fluid_realm, m_phase);
      }
    }
  }

  // Interpolate to centroids
  m_amr->interpolate_to_centroids(rhoPhase, m_fluid_realm, m_phase);
}

void ito_plasma_stepper::compute_conductivity(EBAMRCellData& a_conductivity){
  CH_TIME("ito_plasma_stepper::compute_conductivity(conductivity)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_conductivity(conductivity)" << endl;
  }

  this->compute_conductivity(a_conductivity, m_ito->get_particles(ito_solver::which_container::bulk));
  
}

void ito_plasma_stepper::compute_conductivity(EBAMRCellData& a_conductivity, const Vector<particle_container<ito_particle>* >& a_particles){
  CH_TIME("ito_plasma_stepper::compute_conductivity(conductivity)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_conductivity(conductivity)" << endl;
  }
  
  data_ops::set_value(a_conductivity, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>&   solver = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();
    
    const int idx = solver_it.get_solver();
    const int q   = species->get_charge();

    if(Abs(q) > 0 && solver->is_mobile()){
      data_ops::set_value(m_particle_scratch1, 0.0);

      solver->deposit_conductivity(m_particle_scratch1, *a_particles[idx]);

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

void ito_plasma_stepper::compute_J(EBAMRCellData& a_J, const Real a_dt){
  CH_TIME("ito_plasma_stepper::compute_J(J)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_J(J)" << endl;
  }

  // TLDR: a_J is defined over the fluid realm but the computation takes place on the particle realm.
  //       If the realms are different we compute on a scratch storage instead


  this->compute_conductivity(m_fluid_scratch1);
  data_ops::copy(a_J, m_fluid_E);

  data_ops::multiply_scalar(a_J, m_fluid_scratch1);
}

Real ito_plasma_stepper::compute_relaxation_time(){
  CH_TIME("ito_plasma_stepper::compute_relaxation_time()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_relaxation_time()" << endl;
  }

  Real dt = 1.E99;
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const Real thisDt = this->compute_relaxation_time(lvl);

    dt = Min(dt, thisDt);
  }

  return dt;
}

Real ito_plasma_stepper::compute_relaxation_time(const int a_level){
  CH_TIME("ito_plasma_stepper::compute_relaxation_time(level)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_relaxation_time(level)" << endl;
  }

  Real dt = 1.E99;

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_fluid_realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Real thisDt = this->compute_relaxation_time(a_level, dit());

    dt = Min(dt, thisDt);
  }

#ifdef CH_MPI
  Real tmp = dt;
  int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("cdr_solver::compute_cfl_dt() - communication error on norm");
  }
  dt = tmp;
#endif

  return dt;
}

Real ito_plasma_stepper::compute_relaxation_time(const int a_level, const DataIndex a_dit){
  CH_TIME("ito_plasma_stepper::compute_relaxation_time(level, dit)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_relaxation_time(level, dit)" << endl;
  }

  const int comp    = 0;
  const Real SAFETY = 1.E-10;

  const Box box = m_amr->get_grids(m_fluid_realm)[a_level].get(a_dit);
  const EBISBox& ebisbox = m_amr->get_ebisl(m_fluid_realm, m_phase)[a_level][a_dit];

  // Get a handle to the E-field
  EBAMRCellData amrE;
  m_amr->allocate_ptr(amrE);
  m_amr->alias(amrE, m_phase, m_poisson->get_E());
  
  const EBCellFAB& E = (*amrE[a_level])[a_dit];
  const EBCellFAB& J = (*m_J[a_level])[a_dit];

  EBCellFAB dt(ebisbox, box, 1);
  EBCellFAB e_magnitude(ebisbox, box, 1);
  EBCellFAB j_magnitude(ebisbox, box, 1);

  e_magnitude.setVal(0.0);
  j_magnitude.setVal(0.0);

  data_ops::vector_length(e_magnitude, E, box);
  data_ops::vector_length(j_magnitude, J, box);
  j_magnitude += SAFETY;

  dt.setVal(units::s_eps0);
  dt *= e_magnitude;
  dt /= j_magnitude;
  
  return dt.min(comp);
}

bool ito_plasma_stepper::solve_poisson(){
  CH_TIME("ito_plasma_stepper::solve_poisson()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::solve_poisson()" << endl;
  }

  // Computes rho
  this->compute_rho();

  // Solves the Poisson equation
  const bool converged = m_poisson->solve(m_poisson->get_state(),
					  m_poisson->get_source(),
					  m_sigma->get_state(),
					  false);

  // Computes cell-centered E onto storage in the field solver
  m_poisson->compute_E();

  // Code below here interpolates E to centroids on both realms
  EBAMRCellData E;
  m_amr->allocate_ptr(E);
  m_amr->alias(E, m_phase, m_poisson->get_E());
  
  // Fluid realm
  m_fluid_E.copy(E);
  m_amr->average_down(m_fluid_E, m_fluid_realm, m_phase);
  m_amr->interp_ghost_pwl(m_fluid_E, m_fluid_realm, m_phase);
  m_amr->interpolate_to_centroids(m_fluid_E, m_fluid_realm, m_phase);

  // Particle realm
  m_particle_E.copy(E);
  m_amr->average_down(m_particle_E, m_particle_realm, m_phase);
  m_amr->interp_ghost_pwl(m_particle_E, m_particle_realm, m_phase);
  m_amr->interpolate_to_centroids(m_particle_E, m_particle_realm, m_phase);
    
  return converged;
}

bool ito_plasma_stepper::solve_poisson(MFAMRCellData&                a_potential,
				       MFAMRCellData&                a_rho,
				       const Vector<EBAMRCellData*>& a_densities,
				       const EBAMRIVData&            a_sigma){
  CH_TIME("ito_plasma_stepper::solve_poisson(phi, rho, densities, sigma)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::solve_poisson(phi, rho, densities, sigma)" << endl;
  }

  this->compute_rho(a_rho, a_densities);

  const bool converged = m_poisson->solve(a_potential,
					  a_rho,
					  a_sigma,
					  false);
  m_poisson->compute_E();
  
  return converged;
}

void ito_plasma_stepper::intersect_particles(const which_particles a_which_particles, const EB_representation a_representation, const bool a_delete){
  CH_TIME("ito_plasma_stepper::intersect_particles(which_particles, EB_representation)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::intersect_particles(which_particles, EB_representation)" << endl;
  }

  this->intersect_particles(a_which_particles,
			    ito_solver::which_container::bulk,
			    ito_solver::which_container::eb,
			    ito_solver::which_container::domain,
			    a_representation,
			    a_delete);
}

void ito_plasma_stepper::intersect_particles(const which_particles             a_which_particles,
					     const ito_solver::which_container a_particles,
					     const ito_solver::which_container a_eb_particles,
					     const ito_solver::which_container a_domain_particles,
					     const EB_representation           a_representation,
					     const bool                        a_delete) {
  CH_TIME("ito_plasma_stepper::intersect_particles(which_particles, string, string, string, EB_representation, bool)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::intersect_particles(which_particles, string, string, string, EB_representation, bool)" << endl;
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
      solver->intersect_particles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::all_mobile:
      if(mobile) solver->intersect_particles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::all_diffusive:
      if(diffusive) solver->intersect_particles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) solver->intersect_particles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) solver->intersect_particles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) solver->intersect_particles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) solver->intersect_particles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) solver->intersect_particles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    default:
      MayDay::Abort("ito_plasma_stepper::intersect_particles_particles(which_particles, string, string, string, EB_representation, bool) - logic bust");
    }
  }  
}

void ito_plasma_stepper::remove_covered_particles(const which_particles a_which_particles, const EB_representation a_representation, const Real a_tolerance){
  CH_TIME("ito_plasma_stepper::remove_covered_particles(which_particles, representation, tolerance)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::remove_covered_particles(which_particles, representation, tolerance)" << endl;
  }

  this->remove_covered_particles(a_which_particles, ito_solver::which_container::bulk, a_representation, a_tolerance);
}

void ito_plasma_stepper::remove_covered_particles(const which_particles             a_which,
						  const ito_solver::which_container a_container,
						  const EB_representation           a_representation,
						  const Real                        a_tolerance){
  CH_TIME("ito_plasma_stepper::remove_covered_particles(which_particles, container, EB_representation, tolerance)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::remove_covered_particles(which_particles, container, EB_representation, tolerance)" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>&   solver = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();
    const bool charged   = species->get_charge() != 0;

    switch(a_which) {
    case which_particles::all:
      solver->remove_covered_particles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::all_mobile:
      if(mobile) solver->remove_covered_particles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::all_diffusive:
      if(diffusive) solver->remove_covered_particles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) solver->remove_covered_particles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) solver->remove_covered_particles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) solver->remove_covered_particles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) solver->remove_covered_particles(a_container, a_representation, a_tolerance);
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) solver->remove_covered_particles(a_container, a_representation, a_tolerance);
      break;
    default:
      MayDay::Abort("ito_plasma_stepper::remove_covered_particles_particles(which particles) - logic bust");
    }
  }  
}

void ito_plasma_stepper::remap_particles(const which_particles a_which_particles){
  CH_TIME("ito_plasma_stepper::remap_particles(which_particles)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::remap_particles(which_particles)" << endl;
  }

  this->remap_particles(a_which_particles, ito_solver::which_container::bulk);
}

void ito_plasma_stepper::remap_particles(const which_particles a_which_particles, const ito_solver::which_container a_container){
  CH_TIME("ito_plasma_stepper::remap_particles(which_particles)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::remap_particles(which_particles)" << endl;
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
      solver->remap(a_container);
      break;
    case which_particles::all_mobile:
      if(mobile) solver->remap(a_container);
      break;
    case which_particles::all_diffusive:
      if(diffusive) solver->remap(a_container);
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) solver->remap(a_container);
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) solver->remap(a_container);
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) solver->remap(a_container);
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) solver->remap(a_container);
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) solver->remap(a_container);
      break;
    default:
      MayDay::Abort("ito_plasma_stepper::remap_particles(which particles) - logic bust");
    }
  }
}

void ito_plasma_stepper::deposit_particles(const which_particles a_which_particles){
  CH_TIME("ito_plasma_stepper::deposit_particles(which_particles)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::deposit_particles(which_particles)" << endl;
  }

  this->deposit_particles(a_which_particles, ito_solver::which_container::bulk);

}

void ito_plasma_stepper::deposit_particles(const which_particles a_which_particles, const ito_solver::which_container a_container){
  CH_TIME("ito_plasma_stepper::deposit_particles(which_particles)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::deposit_particles(which_particles)" << endl;
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
      solver->deposit_particles(a_container);
      break;
    case which_particles::all_mobile:
      if(mobile) solver->deposit_particles(a_container);
      break;
    case which_particles::all_diffusive:
      if(diffusive) solver->deposit_particles(a_container);
      break;
    case which_particles::charged_mobile:
      if(charged && mobile) solver->deposit_particles(a_container);
      break;
    case which_particles::charged_diffusive:
      if(charged && diffusive) solver->deposit_particles(a_container);
      break;
    case which_particles::all_mobile_or_diffusive:
      if(mobile || diffusive) solver->deposit_particles(a_container);
      break;
    case which_particles::charged_and_mobile_or_diffusive:
      if(charged && (mobile || diffusive)) solver->deposit_particles(a_container);
      break;
    case which_particles::stationary:
      if(!mobile && !diffusive) solver->deposit_particles(a_container);
      break;
    default:
      MayDay::Abort("ito_plasma_stepper::deposit_particles(which_particles) - logic bust");
    }
  }
}

void ito_plasma_stepper::set_ito_velocity_funcs(){
  CH_TIME("ito_plasma_stepper::set_ito_velocity_funcs");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_ito_velocity_funcs" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver   = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    if(solver->is_mobile()){
      EBAMRCellData& velo_func = solver->get_velo_func();
      velo_func.copy(m_particle_E);
      
      const int q = species->get_charge();
      const int s = (q > 0) - (q < 0);
      
      data_ops::scale(velo_func, s);
    }
  }
}

void ito_plasma_stepper::compute_ito_velocities(){
  CH_TIME("ito_plasma_stepper::compute_ito_velocities()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_velocities()" << endl;
  }

  const ito_plasma_physics::coupling which_coupling = m_physics->get_coupling();

  // Set velocity functions
  this->set_ito_velocity_funcs();

  // Compute mobilities based on appropriate coupling
  switch(which_coupling){
  case ito_plasma_physics::coupling::LFA:
    this->compute_ito_mobilities_lfa();
    break;
  case ito_plasma_physics::coupling::LEA:
    this->compute_ito_mobilities_lea();
    break;
  }

  // Interpolate velocity function to particle position
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->interpolate_velocities(); // Interpolates v = +/- mu*E
  }
}

void ito_plasma_stepper::compute_ito_diffusion(){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion()" << endl;
  }

  const ito_plasma_physics::coupling which_coupling = m_physics->get_coupling();

  // Compute mobilities based on appropriate coupling
  switch(which_coupling){
  case ito_plasma_physics::coupling::LFA:
    this->compute_ito_diffusion_lfa();
    break;
  case ito_plasma_physics::coupling::LEA:
    this->compute_ito_diffusion_lea();
    break;
  }
}

void ito_plasma_stepper::compute_ito_mobilities_lfa(){
  CH_TIME("ito_plasma_stepper::compute_ito_mobilities_lfa()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_mobilities_lfa()" << endl;
  }

  Vector<EBAMRCellData*> meshMobilities = m_ito->get_mobility_funcs();
  this->compute_ito_mobilities_lfa(meshMobilities, m_fluid_E, m_time);
}

void ito_plasma_stepper::compute_ito_mobilities_lfa(Vector<EBAMRCellData*>& a_meshMobilities, const EBAMRCellData& a_E, const Real a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_mobilities_lfa(mobilities, E, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_mobilities_lfa(mobilities, E, time)" << endl;
  }


  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){

    // Computation is done on the fluid realm. 
    Vector<LevelData<EBCellFAB>* > meshMobilities;
    for (int i = 0; i < a_meshMobilities.size(); i++){
      meshMobilities.push_back(&(*(m_fscratch1[i])[lvl]));
    }
    
    this->compute_ito_mobilities_lfa(meshMobilities, *a_E[lvl], lvl, a_time);
  }

  // Average down and interpolate ghost cells. Then interpolate mobilities to particle positions. 
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<ito_solver>& solver = solver_it();

    if(solver->is_mobile()){

      
#if 0 // In principle, we should be able to average down and interpolate on the fluid realm and then copy directly to the particle realm.
      // But we need to make sure that EBAMRData::copy also gets ghost cells 
      m_amr->average_down(m_fscratch1[idx], m_fluid_realm, m_phase);
      m_amr->interp_ghost(m_fscratch1[idx], m_fluid_realm, m_phase);

      a_meshMobilities[idx]->copy(m_fscratch1[idx]);
#else
      // Copy to particle realm, build ghost cells and the interpolate the mobilities to particle positions. 
      a_meshMobilities[idx]->copy(m_fscratch1[idx]);
      
      m_amr->average_down(*a_meshMobilities[idx], m_particle_realm, m_phase);
      m_amr->interp_ghost(*a_meshMobilities[idx], m_particle_realm, m_phase);
#endif

      solver->interpolate_mobilities();
    }
  }
}

void ito_plasma_stepper::compute_ito_mobilities_lfa(Vector<LevelData<EBCellFAB>* >& a_meshMobilities,
						    const LevelData<EBCellFAB>&     a_E,
						    const int                       a_level,
						    const Real                      a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_mobilities_lfa(mobilities, E, level, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_mobilities_lfa(mobilities, E, level, time)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_fluid_realm)[a_level];
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const EBCellFAB& E = a_E[dit()];
    const Box bx       = dbl.get(dit());

    Vector<EBCellFAB*> meshMobilities;
    for (int i = 0; i < a_meshMobilities.size(); i++){
      meshMobilities.push_back(&((*a_meshMobilities[i])[dit()]));
    }

    this->compute_ito_mobilities_lfa(meshMobilities, E, a_level, dit(), bx, a_time);
  }
}

void ito_plasma_stepper::compute_ito_mobilities_lfa(Vector<EBCellFAB*>& a_meshMobilities,
						    const EBCellFAB&    a_E,
						    const int           a_level,
						    const DataIndex     a_dit,
						    const Box           a_box,
						    const Real          a_time) {
  CH_TIME("ito_plasma_stepper::compute_ito_mobilities_lfa(meshMobilities, E, level, dit, box, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_mobilities_lfa(meshMobilities, E, level, dit, box, time)" << endl;
  }

  const int comp         = 0;
  const Real dx          = m_amr->get_dx()[a_level];
  const RealVect prob_lo = m_amr->get_prob_lo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv   = bit();
    const RealVect pos = m_amr->get_prob_lo() + dx*(RealVect(iv) + 0.5*RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv,0), E(iv, 1), E(iv, 2)));

    // Call ito_physics and compute diffusion for each particle species
    const Vector<Real> mobilities = m_physics->compute_ito_mobilities_lfa(a_time, pos, e);
    
    // Put mobilities in data holder
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      (*a_meshMobilities[idx]).getSingleValuedFAB()(iv, comp) = mobilities[idx];
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->get_vofit(m_fluid_realm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, prob_lo);
    
    // Compute diffusion
    const Vector<Real> mobilities = m_physics->compute_ito_mobilities_lfa(a_time, pos, e);

    // Put diffusion in the appropriate place.
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      (*a_meshMobilities[idx])(vof, comp) = mobilities[idx];
    }
  }
  

  // Covered is bogus.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    a_meshMobilities[idx]->setCoveredCellVal(0.0, comp);
  }
}

void ito_plasma_stepper::compute_ito_mobilities_lea(){
  CH_TIME("ito_plasma_stepper");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_mobilities_lea()" << endl;
  }

  // This is really simple because the solvers do this directly...
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->update_mobilities();
  }
}

void ito_plasma_stepper::compute_ito_diffusion_lfa(){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion_lfa()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion_lfa()" << endl;
  }

  Vector<EBAMRCellData*> diffco_funcs = m_ito->get_diffco_funcs();
  Vector<EBAMRCellData*> densities    = m_ito->get_densities();

  this->compute_ito_diffusion_lfa(diffco_funcs, densities, m_particle_E, m_time);
}

void ito_plasma_stepper::compute_ito_diffusion_lfa(Vector<EBAMRCellData*>&       a_diffco_funcs,
						   const Vector<EBAMRCellData*>& a_densities,
						   const EBAMRCellData&          a_E,
						   const Real                    a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion_lfa(velo, E, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion_lfa(velo, E, time)" << endl;
  }

  // TLDR: In this routine we make m_fscratch1 hold the diffusion coefficients on the fluid realm and m_fscratch2 hold the particle densities on the fluid realm.
  //       This requires a couple of copies. 

  // 1. Copy particle realm densities to fluid realm scratch data. 
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();

    m_fscratch2[idx].copy(*a_densities[idx]);
  }

  // 2. Compute on each level. On the fluid realm. 
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){

    const int num_ito_species = m_physics->get_num_ito_species();

    Vector<LevelData<EBCellFAB>* > diffco_funcs(num_ito_species);
    Vector<LevelData<EBCellFAB>* > densities(num_ito_species);
    
    for (int idx = 0; idx < a_diffco_funcs.size(); idx++){
      diffco_funcs[idx] = &(*(m_fscratch1[idx])[lvl]);
      densities[idx]    = &(*(m_fscratch2[idx])[lvl]);
    }

    this->compute_ito_diffusion_lfa(diffco_funcs, densities, *m_fluid_E[lvl], lvl, a_time);
  }

  // Average down, interpolate ghost cells, and then interpolate to particle positions
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<ito_solver>& solver = solver_it();

    if(solver->is_diffusive()){

#if 0 // In principle, we should be able to average down and interpolate ghost cells on the fluid realm, and copy the entire result over to the particle realm.
      m_amr->average_down(m_fscratch1[idx], m_fluid_realm, m_phase);
      m_amr->interp_ghost(m_fscratch2[idx], m_fluid_realm, m_phase);
      a_diffco_funcs[idx]->copy(m_fluid_scratch1[idx]);
#else // Instead, we copy to the particle realm and average down there, then interpolate. 
      a_diffco_funcs[idx]->copy(m_fscratch1[idx]);
      
      m_amr->average_down(*a_diffco_funcs[idx], m_particle_realm, m_phase);
      m_amr->interp_ghost(*a_diffco_funcs[idx], m_particle_realm, m_phase);
#endif

      solver->interpolate_diffusion();
    }
  }
}

void ito_plasma_stepper::compute_ito_diffusion_lfa(Vector<LevelData<EBCellFAB>* >&       a_diffco,
						   const Vector<LevelData<EBCellFAB>* >& a_densities,
						   const LevelData<EBCellFAB>&           a_E,
						   const int                             a_level,
						   const Real                            a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion_lfa(velo, E, level, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion_lfa(velo, E, level, time)" << endl;
  }

  const int num_ito_species = m_physics->get_num_ito_species();

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_fluid_realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());

    Vector<EBCellFAB*> diffusion(num_ito_species);
    Vector<EBCellFAB*> densities(num_ito_species);;
    
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();

      if(solver_it()->is_diffusive()){
	diffusion[idx] = &(*a_diffco[idx])[dit()];
      }
      densities[idx] = &(*a_densities[idx])[dit()];
    }

    this->compute_ito_diffusion_lfa(diffusion, densities, a_E[dit()], a_level, dit(), box, a_time);
  }
}

void ito_plasma_stepper::compute_ito_diffusion_lfa(Vector<EBCellFAB*>&       a_diffco,
						   const Vector<EBCellFAB*>& a_densities,
						   const EBCellFAB&          a_E,
						   const int                 a_level,
						   const DataIndex           a_dit,
						   const Box                 a_box,
						   const Real                a_time){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion_lfa(velo, E, level, dit, time)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion_lfa(velo, E, level, dit, time)" << endl;
  }

  const int comp         = 0;
  const Real dx          = m_amr->get_dx()[a_level];
  const RealVect prob_lo = m_amr->get_prob_lo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv   = bit();
    const RealVect pos = m_amr->get_prob_lo() + dx*(RealVect(iv) + 0.5*RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv,0), E(iv, 1), E(iv, 2)));

    // Make grid densities
    Vector<Real> densities;
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      densities.push_back((*a_densities[idx]).getSingleValuedFAB()(iv, comp));
    }

    // Call ito_physics and compute diffusion for each particle species
    const Vector<Real> diffusion = m_physics->compute_ito_diffusion_lfa(a_time, pos, e, densities);
    
    // Put diffusion where they belong
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<ito_solver>& solver = solver_it();
      if(solver->is_diffusive()){
	const int idx = solver_it.get_solver();
	(*a_diffco[idx]).getSingleValuedFAB()(iv, comp) = diffusion[idx];
      }
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->get_vofit(m_fluid_realm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, prob_lo);
    
    // Get densities
    Vector<Real> densities;
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      densities.push_back((*a_densities[idx])(vof, comp));
    }
    
    // Compute diffusion
    const Vector<Real> diffusion = m_physics->compute_ito_diffusion_lfa(a_time, pos, e, densities);

    // Put diffusion in the appropriate place.
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      if(solver_it()->is_diffusive()){
	const int idx = solver_it.get_solver();
	(*a_diffco[idx])(vof, comp) = diffusion[idx];
      }
    }
  }

  // Covered is bogus.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    if(solver_it()->is_diffusive()){
      const int idx = solver_it.get_solver();
      a_diffco[idx]->setCoveredCellVal(0.0, comp);
    }
  }
}

void ito_plasma_stepper::compute_ito_diffusion_lea(){
  CH_TIME("ito_plasma_stepper::compute_ito_diffusion_lea()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_ito_diffusion_lea()" << endl;
  }
  
  // This is really simple because the solvers do this directly... No monkeying with interpolations or anything. 
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->update_diffusion();
  }
}

void ito_plasma_stepper::compute_reactive_particles_per_cell(EBAMRCellData& a_ppc){
  CH_TIME("ito_plasma_stepper::compute_reactive_particles_per_cell(ppc)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_particles_per_cell(ppc)" << endl;
  }

  data_ops::set_value(a_ppc, 0.0);
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    this->compute_reactive_particles_per_cell(*a_ppc[lvl], lvl);
  }
}

void ito_plasma_stepper::compute_reactive_particles_per_cell(LevelData<EBCellFAB>& a_ppc, const int a_level){
  CH_TIME("ito_plasma_stepper::compute_reactive_particles_per_cell(ppc, lvl)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_particles_per_cell(ppc, lvl)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[a_level];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_particle_realm, m_phase)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      
    const Box box = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];
    
    this->compute_reactive_particles_per_cell(a_ppc[dit()], a_level, dit(), box, ebisbox);
  }
}

void ito_plasma_stepper::compute_reactive_particles_per_cell(EBCellFAB& a_ppc, const int a_level, const DataIndex a_dit, const Box a_box, const EBISBox& a_ebisbox){
  CH_TIME("ito_plasma_stepper::compute_reactive_particles_per_cell(ppc, lvl, dit, box)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_particles_per_cell(ppc, lvl, dit, box)" << endl;
  }

  BaseFab<Real>& numFab = a_ppc.getSingleValuedFAB();

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    const int idx                     = solver_it.get_solver();

    const particle_container<ito_particle>& particles = solver->get_particles(ito_solver::which_container::bulk);
    const BinFab<ito_particle>& cellParticles         = particles.get_cell_particles(a_level, a_dit);


    // Regular cells. 
    for (BoxIterator bit(a_box); bit.ok(); ++bit){
      const IntVect iv = bit();

      Real num = 0.0;
      
      if(a_ebisbox.isRegular(iv)){

	
	const List<ito_particle>& listParticles = cellParticles(iv, 0);
	for (ListIterator<ito_particle> lit(listParticles); lit.ok(); ++lit){
	  num += lit().mass();
	}
      }
      
      numFab(iv, idx) = num;
    }

    // Irregular cells.
    VoFIterator& vofit = (*m_amr->get_vofit(m_particle_realm, m_phase)[a_level])[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof       = vofit();
      const IntVect iv          = vof.gridIndex();
      const RealVect normal     = a_ebisbox.normal(vof);
      const RealVect ebCentroid = a_ebisbox.bndryCentroid(vof);

      Real num = 0.0;

      const List<ito_particle>& listParticles = cellParticles(iv, 0);
      
      for (ListIterator<ito_particle> lit(listParticles); lit.ok(); ++lit){
	const RealVect& pos = lit().position();

	if(PolyGeom::dot((pos-ebCentroid), normal) >= 0.0){
	  num += lit().mass();
	}
      }
      
      numFab(iv, idx) = num;
    }
  }
}

void ito_plasma_stepper::compute_reactive_mean_energies_per_cell(EBAMRCellData& a_mean_energies){
  CH_TIME("ito_plasma_stepper::compute_mean-energies_per_cell(EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_particles_per_cell(EBAMRCellData)" << endl;
  }

  data_ops::set_value(a_mean_energies, 0.0);
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    this->compute_reactive_mean_energies_per_cell(*a_mean_energies[lvl], lvl);
  }
}

void ito_plasma_stepper::compute_reactive_mean_energies_per_cell(LevelData<EBCellFAB>& a_mean_energies, const int a_level){
  CH_TIME("ito_plasma_stepper::compute_reactive_mean_energies_per_cell(ppc, lvl)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_mean_energies_per_cell(ppc, lvl)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[a_level];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_particle_realm, m_phase)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      
    const Box box = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];
      
    this->compute_reactive_mean_energies_per_cell(a_mean_energies[dit()], a_level, dit(), box, ebisbox);
  }
}

void ito_plasma_stepper::compute_reactive_mean_energies_per_cell(EBCellFAB&      a_mean_energies,
								 const int       a_level,
								 const DataIndex a_dit,
								 const Box       a_box,
								 const EBISBox&  a_ebisbox){
  CH_TIME("ito_plasma_stepper::compute_reactive_mean_energies_per_cell(ppc, lvl, dit, box)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_reactive_mean_energies_per_cell(ppc, lvl, dit, box)" << endl;
  }

  BaseFab<Real>& numFab = a_mean_energies.getSingleValuedFAB();

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    const int idx                     = solver_it.get_solver();

    const particle_container<ito_particle>& particles = solver->get_particles(ito_solver::which_container::bulk);
    const BinFab<ito_particle>& cellParticles         = particles.get_cell_particles(a_level, a_dit);

    // Regular cells. 
    for (BoxIterator bit(a_box); bit.ok(); ++bit){
      const IntVect iv = bit();

      if (a_ebisbox.isRegular(iv)){
	Real m = 0.0;
	Real E = 0.0;
	
	const List<ito_particle>& listParticles = cellParticles(iv, 0);

	for (ListIterator<ito_particle> lit(listParticles); lit.ok(); ++lit){
	  m += lit().mass();
	  E += lit().mass()*lit().energy();
	}

	numFab(iv, idx) = E/m;
      }
    }

    // Irregular cells.
    VoFIterator& vofit = (*m_amr->get_vofit(m_particle_realm, m_phase)[a_level])[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof       = vofit();
      const IntVect iv          = vof.gridIndex();
      const RealVect normal     = a_ebisbox.normal(vof);
      const RealVect ebCentroid = a_ebisbox.bndryCentroid(vof);

      Real m = 0.0;
      Real E = 0.0;

      const List<ito_particle>& listParticles = cellParticles(iv, 0);
      
      for (ListIterator<ito_particle> lit(listParticles); lit.ok(); ++lit){
	const RealVect& pos = lit().position();

	if(PolyGeom::dot((pos-ebCentroid), normal) >= 0.0){
	  m += lit().mass();
	  E += lit().mass()*lit().energy();
	}
      }

      numFab(iv, idx) = E/m;
    }
  }
}

void ito_plasma_stepper::advance_reaction_network_nwo(const Real a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network_nwo(dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network_nwo(dt)" << endl;
  }

  this->advance_reaction_network_nwo(m_fluid_E, m_EdotJ, a_dt);
}

void ito_plasma_stepper::advance_reaction_network_nwo(const EBAMRCellData& a_E, const EBAMRCellData& a_EdotJ, const Real a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(ppc, ypc, E, sources, dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(ppc, ypc, E, sources, dt)" << endl;
  }

  // 1. Compute the number of particles per cell. Set the number of photons to be generated per cell to zero. 
  this->compute_reactive_particles_per_cell(m_particle_ppc);
  this->compute_reactive_mean_energies_per_cell(m_particle_eps);

  m_fluid_ppc.copy(m_particle_ppc);
  m_fluid_eps.copy(m_particle_eps);
  
  data_ops::set_value(m_fluid_ypc,    0.0);
  data_ops::set_value(m_particle_ypc, 0.0);
  data_ops::copy(m_particle_old, m_particle_ppc);

  // 2. Solve for the new number of particles per cell. This also obtains the number of photons to be generated in each cell. 
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    this->advance_reaction_network_nwo(*m_fluid_ppc[lvl], *m_fluid_ypc[lvl], *m_fluid_eps[lvl], *a_E[lvl], *a_EdotJ[lvl], lvl, a_dt);
  }

  // 3. Copy the results to the particle realm.
  m_particle_ppc.copy(m_fluid_ppc);
  m_particle_ypc.copy(m_fluid_ypc);
  m_particle_eps.copy(m_fluid_eps);

  // 4. Reconcile particles on the particle realm. Not implemented (yet).
  this->reconcile_particles(m_particle_ppc, m_particle_old, m_particle_eps, m_particle_ypc);
}


void ito_plasma_stepper::advance_reaction_network_nwo(LevelData<EBCellFAB>&       a_particlesPerCell,
						      LevelData<EBCellFAB>&       a_newPhotonsPerCell,
						      LevelData<EBCellFAB>&       a_meanParticleEnergies,
						      const LevelData<EBCellFAB>& a_E,
						      const LevelData<EBCellFAB>& a_EdotJ,
						      const int                   a_level,
						      const Real                  a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(ppc, ypc, energies, E, sources, level, dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(ppc, ypc, energies, E, sources, level, dt)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_fluid_realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    this->advance_reaction_network_nwo(a_particlesPerCell[dit()],
				       a_newPhotonsPerCell[dit()],
				       a_meanParticleEnergies[dit()],
				       a_E[dit()],
				       a_EdotJ[dit()],
				       a_level,
				       dit(),
				       dbl[dit()],
				       m_amr->get_dx()[a_level],
				       a_dt);
  }
}

void ito_plasma_stepper::advance_reaction_network_nwo(EBCellFAB&       a_particlesPerCell,
						      EBCellFAB&       a_newPhotonsPerCell,
						      EBCellFAB&       a_meanParticleEnergies,
						      const EBCellFAB& a_E,
						      const EBCellFAB& a_EdotJ,
						      const int        a_level,
						      const DataIndex  a_dit,
						      const Box        a_box,
						      const Real       a_dx,
						      const Real       a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network_nwo(ppc, ypc, E, sources, level, dit, box, dx, dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network_nwo(ppc, ypc, E, sources, level, dit, box, dx, dt)" << endl;
  }
  
  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();

  const RealVect prob_lo = m_amr->get_prob_lo();

  const EBISBox& ebisbox = m_amr->get_ebisl(m_fluid_realm, m_phase)[a_level][a_dit];
  const EBISBox& ebgraph = m_amr->get_ebisl(m_fluid_realm, m_phase)[a_level][a_dit];

  const BaseFab<Real>& Efab = a_E.getSingleValuedFAB();

  Vector<long long> particles(num_ito_species);
  Vector<long long> newPhotons(num_rte_species);
  Vector<Real>      meanEnergies(num_ito_species);
  Vector<Real>      energySources(num_ito_species);

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();

    if(ebisbox.isRegular(iv)){
      const Real kappa = 1.0;
      const Real dV    = kappa*pow(a_dx, SpaceDim);

      const RealVect pos = prob_lo + a_dx*(RealVect(iv) + 0.5*RealVect::Unit);
      const RealVect E   = RealVect(D_DECL(Efab(iv, 0), Efab(iv, 1), Efab(iv, 2)));

      // Initialize for this cell. 
      for (int i = 0; i < num_ito_species; i++){
	particles[i]     = llround(a_particlesPerCell.getSingleValuedFAB()(iv, i));
	meanEnergies[i]  = a_meanParticleEnergies.getSingleValuedFAB()(iv,i);
	energySources[i] = a_EdotJ.getSingleValuedFAB()(iv, i)*dV/units::s_Qe;
      }

      for (int i = 0; i < num_rte_species; i++){
	newPhotons[i]= 0LL;
      }
	   
      // Do the physics advance
      m_physics->advance_particles(particles, newPhotons, meanEnergies, energySources, a_dt, E, a_dx, kappa);

      // Set result
      for (int i = 0; i < num_ito_species; i++){
	a_particlesPerCell.getSingleValuedFAB()(iv, i)     = 1.0*particles[i];
	a_meanParticleEnergies.getSingleValuedFAB()(iv, i) = 1.0*meanEnergies[i];
      }

      for (int i = 0; i < num_rte_species; i++){
	a_newPhotonsPerCell.getSingleValuedFAB()(iv, i) = 1.0*newPhotons[i];
      }
    }
  }

  // Irregular cells
  VoFIterator& vofit = (*m_amr->get_vofit(m_fluid_realm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const Real kappa    = ebisbox.volFrac(vof);
    const Real dV       = kappa*pow(a_dx, SpaceDim);
    const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, prob_lo);
    const RealVect E    = RealVect(D_DECL(a_E(vof,0), a_E(vof,1), a_E(vof,2)));


    // Initialize for this cell. 
    for (int i = 0; i < num_ito_species; i++){
      particles[i]     = a_particlesPerCell(vof, i);
      meanEnergies[i]  = a_meanParticleEnergies(vof, i);
      energySources[i] = a_EdotJ(vof, i)*dV/units::s_Qe;
    }

    for (int i = 0; i < num_rte_species; i++){
      newPhotons[i]= 0LL;
    }

    m_physics->advance_particles(particles, newPhotons, meanEnergies, energySources, a_dt, E, a_dx, kappa);


    // Set result
    for (int i = 0; i < num_ito_species; i++){
      a_particlesPerCell(vof, i)     = 1.0*particles[i];
      a_meanParticleEnergies(vof, i) = 1.0*meanEnergies[i];
    }
    
    for (int i = 0; i < num_rte_species; i++){
      a_newPhotonsPerCell(vof, i) = 1.0*newPhotons[i];
    }
  }
}

void ito_plasma_stepper::reconcile_particles(const EBAMRCellData& a_newParticlesPerCell,
					     const EBAMRCellData& a_oldParticlesPerCell,
					     const EBAMRCellData& a_meanParticleEnergies,
					     const EBAMRCellData& a_newPhotonsPerCell){
  CH_TIME("ito_plasma_stepper::reconcile_particles(EBAMRCellData, EBAMRCellData, EBAMRCellData)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::reconcile_particles(EBAMRCellData, EBAMRCellData, EBAMRCellData)";
  }

  for(int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){

    this->reconcile_particles(*a_newParticlesPerCell[lvl], *a_oldParticlesPerCell[lvl], *a_meanParticleEnergies[lvl], *a_newPhotonsPerCell[lvl], lvl);
  }
}

void ito_plasma_stepper::reconcile_particles(const LevelData<EBCellFAB>& a_newParticlesPerCell,
					     const LevelData<EBCellFAB>& a_oldParticlesPerCell,
					     const LevelData<EBCellFAB>& a_meanParticleEnergies,
					     const LevelData<EBCellFAB>& a_newPhotonsPerCell,
					     const int                   a_level){
  CH_TIME("ito_plasma_stepper::reconcile_particles(LevelData<EBCellFAB>x3, int)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::reconcile_particles(LevelData<EBCellFAB>x3, int)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    this->reconcile_particles(a_newParticlesPerCell[dit()],
			      a_oldParticlesPerCell[dit()],
			      a_meanParticleEnergies[dit()],
			      a_newPhotonsPerCell[dit()],
			      a_level,
			      dit(),
			      dbl[dit()],
			      m_amr->get_dx()[a_level]);
  }
}

void ito_plasma_stepper::reconcile_particles(const EBCellFAB& a_newParticlesPerCell,
					     const EBCellFAB& a_oldParticlesPerCell,
					     const EBCellFAB& a_meanParticleEnergies,
					     const EBCellFAB& a_newPhotonsPerCell,
					     const int        a_level,
					     const DataIndex  a_dit,
					     const Box        a_box,
					     const Real       a_dx){
  CH_TIME("ito_plasma_stepper::reconcile_particles(EBCellFABx3, int, DataIndex, Box, Real)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::reconcile_particles(EBCellFABx3, int, DataIndex, Box, Real)" << endl;
  }

  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();

  const RealVect prob_lo = m_amr->get_prob_lo();

  const EBISBox& ebisbox = m_amr->get_ebisl(m_particle_realm, m_phase)[a_level][a_dit];
  const EBISBox& ebgraph = m_amr->get_ebisl(m_particle_realm, m_phase)[a_level][a_dit];

  Vector<BinFab<ito_particle>* > particlesFAB(num_ito_species);
  Vector<BinFab<photon>* >       sourcePhotonsFAB(num_rte_species);
  Vector<BinFab<photon>* >       bulkPhotonsFAB(num_rte_species);


  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    const int idx = solver_it.get_solver();
    
    particle_container<ito_particle>& solverParticles = solver->get_particles(ito_solver::which_container::bulk);
    
    particlesFAB[idx] = &(solverParticles.get_cell_particles(a_level, a_dit));
  }
  
  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    const int idx = solver_it.get_solver();
    
    particle_container<photon>& solverBulkPhotons = solver->get_bulk_photons();
    particle_container<photon>& solverSourPhotons = solver->get_source_photons();
    
    bulkPhotonsFAB[idx]   = &(solverBulkPhotons.get_cell_particles(a_level, a_dit));
    sourcePhotonsFAB[idx] = &(solverSourPhotons.get_cell_particles(a_level, a_dit));
  }

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    if(ebisbox.isRegular(iv)){
      const RealVect cellPos       = prob_lo + a_dx*(RealVect(iv) + 0.5*RealVect::Unit);
      const RealVect centroidPos   = cellPos;
      const RealVect lo            = -0.5*RealVect::Unit;
      const RealVect hi            = 0.5*RealVect::Unit;
      const RealVect bndryCentroid = RealVect::Zero;
      const RealVect bndryNormal   = RealVect::Zero;
      const Real     kappa         = 1.0;

      Vector<List<ito_particle>* > particles(num_ito_species);
      Vector<List<photon>* >       bulkPhotons(num_rte_species);
      Vector<List<photon>* >       sourcePhotons(num_rte_species);
      Vector<RefCountedPtr<rte_species> > photoSpecies(num_rte_species);

      Vector<Real>      particleMeanEnergies(num_ito_species);
      Vector<long long> numNewParticles(num_ito_species);
      Vector<long long> numOldParticles(num_ito_species);
      Vector<long long> numNewPhotons(num_rte_species);

      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.get_solver();

	particles[idx]            = &((*particlesFAB[idx])(iv, 0));
	particleMeanEnergies[idx] = a_meanParticleEnergies.getSingleValuedFAB()(iv, idx);
	numNewParticles[idx]      = llround(a_newParticlesPerCell.getSingleValuedFAB()(iv, idx));
	numOldParticles[idx]      = llround(a_oldParticlesPerCell.getSingleValuedFAB()(iv, idx));
      }

      for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.get_solver();

	bulkPhotons[idx]   = &((*bulkPhotonsFAB[idx])(iv, 0));
	sourcePhotons[idx] = &((*sourcePhotonsFAB[idx])(iv, 0));
	photoSpecies[idx]  = solver_it()->get_species();
	numNewPhotons[idx] = llround(a_newPhotonsPerCell.getSingleValuedFAB()(iv, idx));

	sourcePhotons[idx]->clear();
      }

      // Reconcile particles, photons, and photoionization
      m_physics->reconcile_particles(particles, numNewParticles, numOldParticles, cellPos, centroidPos, lo, hi, bndryCentroid, bndryNormal, a_dx, kappa);
      m_physics->reconcile_photons(sourcePhotons, numNewPhotons, cellPos, centroidPos, lo, hi, bndryCentroid, bndryNormal, a_dx, kappa);
      m_physics->reconcile_photoionization(particles, particleMeanEnergies, numNewParticles, bulkPhotons);
      m_physics->set_mean_particle_energy(particles, particleMeanEnergies);
      
      // Clear the bulk photons - they have now been absorbed on the mesh. 
      for (int i = 0; i < num_rte_species; i++){
	//	bulkPhotons[i]->clear();
      }
    }
  }

  // Irregular cells
  VoFIterator& vofit = (*m_amr->get_vofit(m_particle_realm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof          = vofit();
    const IntVect iv             = vof.gridIndex();
    const Real kappa             = ebisbox.volFrac(vof);
    const RealVect cellPos       = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, prob_lo);
    const RealVect centroidPos   = ebisbox.centroid(vof);
    const RealVect bndryCentroid = ebisbox.bndryCentroid(vof);
    const RealVect bndryNormal   = ebisbox.normal(vof);

    RealVect lo            = -0.5*RealVect::Unit;
    RealVect hi            = 0.5*RealVect::Unit;
    if(kappa < 1.0){
      data_ops::compute_min_valid_box(lo, hi, bndryNormal, bndryCentroid);
    }

    Vector<List<ito_particle>* > particles(num_ito_species);
    Vector<List<photon>* >       bulkPhotons(num_rte_species);
    Vector<List<photon>* >       sourcePhotons(num_rte_species);
    Vector<RefCountedPtr<rte_species> > photoSpecies(num_rte_species);

    Vector<Real>      particleMeanEnergies(num_ito_species);
    Vector<long long> numNewParticles(num_ito_species);
    Vector<long long> numOldParticles(num_ito_species);
    Vector<long long> numNewPhotons(num_rte_species);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();

      particles[idx]            = &((*particlesFAB[idx])(iv, 0));
      particleMeanEnergies[idx] = a_meanParticleEnergies(vof, idx);
      numNewParticles[idx]      = llround(a_newParticlesPerCell(vof, idx));
      numOldParticles[idx]      = llround(a_oldParticlesPerCell(vof, idx));
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();

      bulkPhotons[idx]   = &((*bulkPhotonsFAB[idx])(iv, 0));
      sourcePhotons[idx] = &((*sourcePhotonsFAB[idx])(iv, 0));
      photoSpecies[idx]  = solver_it()->get_species();
      numNewPhotons[idx] = llround(a_newPhotonsPerCell(vof, idx));

      sourcePhotons[idx]->clear();
    }

    // Reconcile particles, photons, and photoionization
    m_physics->reconcile_particles(particles, numNewParticles, numOldParticles, cellPos, centroidPos, lo, hi, bndryCentroid, bndryNormal, a_dx, kappa);
    m_physics->reconcile_photons(sourcePhotons, numNewPhotons, cellPos, centroidPos, lo, hi, bndryCentroid, bndryNormal, a_dx, kappa);
    m_physics->reconcile_photoionization(particles, particleMeanEnergies, numNewParticles, bulkPhotons);
    m_physics->set_mean_particle_energy(particles, particleMeanEnergies);
      
    // Clear the bulk photons - they have now been absorbed on the mesh. 
    for (int i = 0; i < num_rte_species; i++){
      //      bulkPhotons[i]->clear();
    }
  }
}

void ito_plasma_stepper::advance_reaction_network(const Real a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(a_dt)" << endl;
  }

  if(m_nwo_reactions){
    this->advance_reaction_network_nwo(a_dt);
  }
  else{
    const int num_ito_species = m_physics->get_num_ito_species();
    const int num_rte_species = m_physics->get_num_rte_species();
  
    Vector<particle_container<ito_particle>* > particles(num_ito_species);  // Current particles. 
    Vector<particle_container<photon>* > bulk_photons(num_rte_species);     // Photons absorbed on mesh
    Vector<particle_container<photon>* > new_photons(num_rte_species);      // Produced photons go here.

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      particles[solver_it.get_solver()] = &(solver_it()->get_particles(ito_solver::which_container::bulk));
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      bulk_photons[solver_it.get_solver()] = &(solver_it()->get_bulk_photons());
      new_photons[solver_it.get_solver()] = &(solver_it()->get_source_photons());
    }

    this->advance_reaction_network(particles, bulk_photons, new_photons, m_energy_sources, m_particle_E, a_dt);
  }
}

void ito_plasma_stepper::advance_reaction_network(Vector<particle_container<ito_particle>* >& a_particles,
						  Vector<particle_container<photon>* >&       a_photons,
						  Vector<particle_container<photon>* >&       a_newPhotons,
						  const Vector<EBAMRCellData>&                a_sources,
						  const EBAMRCellData&                        a_E,
						  const Real                                  a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(Vector<particle_container*> x3, Vector<EBAMRCellData*>, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(Vector<particle_container*> x3, Vector<EBAMRCellData*>, EBAMRCellData, Real)" << endl;
  }

  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();

  Vector<AMRCellParticles<ito_particle>* > particles(num_ito_species);
  Vector<AMRCellParticles<photon>* >       photons(num_ito_species);
  Vector<AMRCellParticles<photon>* >       newPhotons(num_ito_species);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    particles[idx] = &(a_particles[idx]->get_cell_particles());
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    photons[idx]    = &(a_photons[idx]->get_cell_particles());
    newPhotons[idx] = &(a_newPhotons[idx]->get_cell_particles());
  }
				
  //Advance reaction network
  this->advance_reaction_network(particles, photons, newPhotons, a_sources, m_particle_E, a_dt);

}

void ito_plasma_stepper::advance_reaction_network(Vector<AMRCellParticles<ito_particle>* >& a_particles,
						  Vector<AMRCellParticles<photon>* >&       a_photons,
						  Vector<AMRCellParticles<photon>* >&       a_newPhotons,
						  const Vector<EBAMRCellData>&              a_sources,
						  const EBAMRCellData&                      a_E,
						  const Real                                a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(Vector<AMRCellParticles*> x 3, Vector<EBAMRCellData*>, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(Vector<AMRCellParticles*> x 3, Vector<EBAMRCellData*>, EBAMRCellData, Real)" << endl;
  }

  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    Vector<LayoutData<BinFab<ito_particle> >* > particles(num_ito_species);
    Vector<LayoutData<BinFab<photon> >* >       photons(num_rte_species);
    Vector<LayoutData<BinFab<photon> >* >       newPhotons(num_rte_species);
    Vector<LevelData<EBCellFAB>* >              sources(num_ito_species);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      particles[idx] = &(*(*a_particles[idx])[lvl]);
      sources[idx]   = &(*(a_sources[idx])[lvl]);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      photons[idx]    = &(*(*a_photons[idx])[lvl]);
      newPhotons[idx] = &(*(*a_newPhotons[idx])[lvl]);
    }
      
    this->advance_reaction_network(particles, photons, newPhotons, sources, *a_E[lvl], lvl, a_dt);
  }
}

void ito_plasma_stepper::advance_reaction_network(Vector<LayoutData<BinFab<ito_particle> >* >& a_particles,
						  Vector<LayoutData<BinFab<photon> >* >&       a_photons,
						  Vector<LayoutData<BinFab<photon> >* >&       a_newPhotons,
						  const Vector<LevelData<EBCellFAB>* >&        a_sources,
						  const LevelData<EBCellFAB>&                  a_E,
						  const int                                    a_lvl,
						  const Real                                   a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(Vector<LD<BinFab>* > x 3, Vector<LD<EBCellFAB>*>, EBAMRCellData, level, dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(Vector<LD<BinFab>* > x 3, Vector<LD<EBCellFAB>*>, EBAMRCellData, level, dt)" << endl;
  }

  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[a_lvl];
  const Real dx = m_amr->get_dx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());

    Vector<BinFab<ito_particle>* > particles(num_ito_species);
    Vector<BinFab<photon>* >       photons(num_rte_species);;
    Vector<BinFab<photon>* >       newPhotons(num_rte_species);
    Vector<EBCellFAB*>             sources(num_ito_species);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      particles[idx] = &((*a_particles[idx])[dit()]);
      sources[idx]   = &((*a_sources[idx])[dit()]);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      photons[idx]    = &((*a_photons[idx])[dit()]);
      newPhotons[idx] = &((*a_newPhotons[idx])[dit()]);
    }

    this->advance_reaction_network(particles, photons, newPhotons, sources, a_E[dit()], a_lvl, dit(), box, dx, a_dt);
  }
}

void ito_plasma_stepper::advance_reaction_network(Vector<BinFab<ito_particle>* >& a_particles,
						  Vector<BinFab<photon>* >&       a_photons,
						  Vector<BinFab<photon>* >&       a_newPhotons,
						  const Vector<EBCellFAB*>&       a_sources,
						  const EBCellFAB&                a_E,
						  const int                       a_lvl,
						  const DataIndex                 a_dit,
						  const Box                       a_box,
						  const Real                      a_dx,
						  const Real                      a_dt){
  CH_TIME("ito_plasma_stepper::advance_reaction_network(Vector<BinFab*> x 3, Vector<EBCellFAB*>, EBCellFAB, level, dit, box, dx, dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_reaction_network(Vector<BinFab*> x 3, Vector<EBCellFAB*>, EBCellFAB, level, dit, box, dx, dt)" << endl;
  }

  const int comp = 0;

  const int num_ito_species = m_physics->get_num_ito_species();
  const int num_rte_species = m_physics->get_num_rte_species();

  const RealVect prob_lo = m_amr->get_prob_lo();
  const RealVect dx      = a_dx*RealVect::Unit;

  const EBISBox& ebisbox = m_amr->get_ebisl(m_particle_realm, m_phase)[a_lvl][a_dit];
  const EBISBox& ebgraph = m_amr->get_ebisl(m_particle_realm, m_phase)[a_lvl][a_dit];

  const BaseFab<Real>& Efab = a_E.getSingleValuedFAB();

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();
    
    if(ebisbox.isRegular(iv)){
      const Real kappa   = 1.0;
      const RealVect pos = prob_lo + a_dx*(RealVect(iv) + 0.5*RealVect::Unit);
      const RealVect e   = RealVect(D_DECL(Efab(iv, 0), Efab(iv, 1), Efab(iv, 2)));
      
      Vector<List<ito_particle>* > particles(num_ito_species);
      Vector<List<photon>* >       photons(num_rte_species);
      Vector<List<photon>* >       newPhotons(num_rte_species);
      Vector<Real>                 sources(num_ito_species);

      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.get_solver();
      
	List<ito_particle>& bp = (*a_particles[idx])(iv, comp);
	particles[idx] = &bp;

	const BaseFab<Real>& sourcesFAB = a_sources[idx]->getSingleValuedFAB();
	sources[idx] = sourcesFAB(iv, comp);
      }

      for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.get_solver();
      
	List<photon>& bp    = (*a_photons[idx])(iv, comp);
	List<photon>& bpNew = (*a_newPhotons[idx])(iv, comp);
      
	photons[idx]    = &bp;
	newPhotons[idx] = &bpNew;
      }

      // Dummy stuff for regular cells
      const RealVect lo = -0.5*RealVect::Unit;
      const RealVect hi =  0.5*RealVect::Unit;
      const RealVect n  = RealVect::Zero;
      const RealVect c  = RealVect::Zero;

      // Advance reactions
      m_physics->advance_reaction_network(particles, photons, newPhotons, sources, e, pos, c, c, n, lo, hi, a_dx, kappa, a_dt);
    }
  }

  // Now do the irregular cells
  VoFIterator& vofit = (*m_amr->get_vofit(m_particle_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex vof = vofit();
    const IntVect iv   = vof.gridIndex();
    const RealVect pos = prob_lo + a_dx*(RealVect(iv) + 0.5*RealVect::Unit);
    const RealVect cen = ebisbox.centroid(vof);
    const Real kappa   = ebisbox.volFrac(vof);
    const RealVect e   = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect n   = ebisbox.normal(vof);
    const RealVect ebc = ebisbox.bndryCentroid(vof);


    // Compute a small box that encloses the cut-cell volume
    RealVect lo = -0.5*RealVect::Unit;
    RealVect hi =  0.5*RealVect::Unit;
    if(kappa < 1.0){
      data_ops::compute_min_valid_box(lo, hi, n, ebc);
    }

    Vector<List<ito_particle>* > particles(num_ito_species);
    Vector<List<photon>* >       photons(num_rte_species);
    Vector<List<photon>* >       newPhotons(num_rte_species);
    Vector<Real>                 sources(num_ito_species);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      
      List<ito_particle>& bp = (*a_particles[idx])(iv, comp);
      particles[idx] = &bp;

      const BaseFab<Real>& sourcesFAB = a_sources[idx]->getSingleValuedFAB();
      sources[idx] = sourcesFAB(iv, comp);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      
      List<photon>& bp    = (*a_photons[idx])(iv, comp);
      List<photon>& bpNew = (*a_newPhotons[idx])(iv, comp);
      
      photons[idx]    = &bp;
      newPhotons[idx] = &bpNew;
    }

    // Advance reactions
    m_physics->advance_reaction_network(particles, photons, newPhotons, sources, e, pos, cen, ebc, n, lo, hi, a_dx, kappa, a_dt);
  }
}

Real ito_plasma_stepper::compute_physics_dt() const{
  CH_TIME("ito_plasma_stepper::compute_physics_dt()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_physics_dt()" << endl;
  }

  // TLDR: This is done on the particle realm because of the densities (which are defined on the particle realm). 

  const Real dt = this->compute_physics_dt(m_particle_E, m_ito->get_densities());

  return dt;
}

Real ito_plasma_stepper::compute_physics_dt(const EBAMRCellData& a_E, const Vector<EBAMRCellData*> a_densities) const {
  CH_TIME("ito_plasma_stepper::compute_physics_dt(EBAMRCellFAB, Vector<EBAMRCellFAB*>)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_physics_dt(EBAMRCellFAB, Vector<EBAMRCellFAB*>)" << endl;
  }



  const int num_ito_species = m_physics->get_num_ito_species();

  Real minDt = 1.E99;
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    
    Vector<LevelData<EBCellFAB>*> densities(num_ito_species);
    
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();

      densities[idx] = &(*(*a_densities[idx])[lvl]);
    }

    const Real levelDt = this->compute_physics_dt(*a_E[lvl], densities, lvl);

    minDt = Min(minDt, levelDt);
  }

  return minDt;
}

Real ito_plasma_stepper::compute_physics_dt(const LevelData<EBCellFAB>& a_E, const Vector<LevelData<EBCellFAB> *> a_densities, const int a_level) const {
  CH_TIME("ito_plasma_stepper::compute_physics_dt(LevelData<EBCellFAB>, Vector<LevelData<EBCellFAB> *>, int)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_physics_dt(LevelData<EBCellFAB>, Vector<LevelData<EBCellFAB> *>, int)" << endl;
  }

  const int num_ito_species = m_physics->get_num_ito_species();

  const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[a_level];

  Real minDt = 1.E99;
  
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

    Vector<EBCellFAB*> densities(num_ito_species);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      
      densities[idx] = &((*a_densities[idx])[dit()]);
    }

    const Real boxDt = this->compute_physics_dt(a_E[dit()], densities, a_level, dit(), dbl.get(dit()));

    minDt = Min(minDt, boxDt);
  }

  // MPI reduction....
#ifdef CH_MPI
  Real tmp = 1.;
  int result = MPI_Allreduce(&minDt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("ito_plasma_stepper::compute_physics_dt(lvl) - communication error on norm");
  }
  minDt = tmp;
#endif  

  return minDt;
}

Real ito_plasma_stepper::compute_physics_dt(const EBCellFAB& a_E, const Vector<EBCellFAB*> a_densities, const int a_level, const DataIndex a_dit, const Box a_box) const {
  CH_TIME("ito_plasma_stepper::compute_physics_dt(EBCellFAB, Vector<EBCellFAB*>, lvl, dit, box)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_physics_dt(EBCellFAB, Vector<EBCellFAB*>, lvl, dit, box)" << endl;
  }

  Real minDt = 1.E99;

  const int num_ito_species = m_physics->get_num_ito_species();
  
  const int comp         = 0;
  const Real dx          = m_amr->get_dx()[a_level];
  const RealVect prob_lo = m_amr->get_prob_lo();
  const BaseFab<Real>& E = a_E.getSingleValuedFAB();
  const EBISBox& ebisbox = m_amr->get_ebisl(m_particle_realm, m_phase)[a_level][a_dit];

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv   = bit();
    const RealVect pos = m_amr->get_prob_lo() + dx*(RealVect(iv) + 0.5*RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv,0), E(iv, 1), E(iv, 2)));

    if(ebisbox.isRegular(iv)){
      Vector<Real> densities(num_ito_species);
      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
	const int idx = solver_it.get_solver();
	
	const BaseFab<Real>& basefab = a_densities[idx]->getSingleValuedFAB();
	densities[idx] = basefab(iv, comp);
      }

      const Real cellDt = m_physics->compute_dt(e, pos, densities);

      minDt = Min(minDt, cellDt);
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->get_vofit(m_particle_realm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect e    = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect pos  = EBArith::getVofLocation(vof, dx*RealVect::Unit, prob_lo);

    Vector<Real> densities(num_ito_species);
    
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      densities[idx] = (*a_densities[idx])(vof, comp);
    }

    const Real cellDt = m_physics->compute_dt(e, pos, densities);

    minDt = Min(minDt, cellDt);
  }

  return minDt;
}

void ito_plasma_stepper::advance_photons(const Real a_dt){
  CH_TIME("ito_plasma_stepper::advance_photons(a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::advance_advance_photons(a_dt)" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    
    // Add source photons and move the photons
    particle_container<photon>& photons        = solver->get_photons();
    particle_container<photon>& bulkPhotons    = solver->get_bulk_photons();
    particle_container<photon>& ebPhotons      = solver->get_eb_photons();
    particle_container<photon>& domainPhotons  = solver->get_domain_photons();
    particle_container<photon>& sourcePhotons  = solver->get_source_photons();

    if(solver->is_instantaneous()){
      solver->clear(photons);

      // Add source photons
      photons.add_particles(sourcePhotons);
      solver->clear(sourcePhotons);

      // Instantaneous advance
      solver->advance_photons_stationary(bulkPhotons, ebPhotons, domainPhotons, photons);
    }
    else{
      // Add source photons
      photons.add_particles(sourcePhotons);
      solver->clear(sourcePhotons);

      // Stationary advance
      solver->advance_photons_transient(bulkPhotons, ebPhotons, domainPhotons, photons, a_dt);
    }
  }
}

void ito_plasma_stepper::sort_photons_by_cell(){
  CH_TIME("ito_plasma_stepper::sort_photons_by_cell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_photons_by_cell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_photons_by_cell();
  }
}

void ito_plasma_stepper::sort_photons_by_patch(){
  CH_TIME("ito_plasma_stepper::sort_photons_by_patch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_photons_by_patch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_photons_by_patch();
  }
}

void ito_plasma_stepper::sort_source_photons_by_cell(){
  CH_TIME("ito_plasma_stepper::sort_source_photons_by_cell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_source_photons_by_cell()" << endl;
  }
  
  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_source_photons_by_cell();
  }
}

void ito_plasma_stepper::sort_source_photons_by_patch(){
  CH_TIME("ito_plasma_stepper::sort_source_photons_by_patch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_source_photons_by_patch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_source_photons_by_patch();
  }
}

void ito_plasma_stepper::sort_bulk_photons_by_cell(){
  CH_TIME("ito_plasma_stepper::sort_bulk_photons_by_cell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_bulk_photons_by_cell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_bulk_photons_by_cell();
  }
}

void ito_plasma_stepper::sort_bulk_photons_by_patch(){
  CH_TIME("ito_plasma_stepper::sort_bulk_photons_by_patch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_bulk_photons_by_patch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_bulk_photons_by_patch();
  }
}

void ito_plasma_stepper::sort_eb_photons_by_cell(){
  CH_TIME("ito_plasma_stepper::sort_eb_photons_by_cell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_eb_photons_by_cell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_eb_photons_by_cell();
  }
}

void ito_plasma_stepper::sort_eb_photons_by_patch(){
  CH_TIME("ito_plasma_stepper::sort_eb_photons_by_patch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_eb_photons_by_patch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_eb_photons_by_patch();
  }
}

void ito_plasma_stepper::sort_domain_photons_by_cell(){
  CH_TIME("ito_plasma_stepper::sort_domain_photons_by_cell()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_domain_photons_by_cell()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_domain_photons_by_cell();
  }
}

void ito_plasma_stepper::sort_domain_photons_by_patch(){
  CH_TIME("ito_plasma_stepper::sort_domain_photons_by_patch()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::sort_domain_photons_by_patch()" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->sort_domain_photons_by_patch();
  }

}

bool ito_plasma_stepper::load_balance_realm(const std::string a_realm) const{
  CH_TIME("time_stepper::load_balance_realm");
  if(m_verbosity > 5){
    pout() << "time_stepper::load_balance_realm" << endl;
  }

  bool ret = false;

  if(a_realm == m_particle_realm && m_load_balance){
    ret = true;
  }
  
  return ret;
}

void ito_plasma_stepper::load_balance_boxes(Vector<Vector<int> >&            a_procs,
					    Vector<Vector<Box> >&            a_boxes,
					    const std::string                a_realm,
					    const Vector<DisjointBoxLayout>& a_grids,
					    const int                        a_lmin,
					    const int                        a_finest_level){
  CH_TIME("ito_plasma_stepper::load_balance_boxes");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper_stepper::load_balance_boxes" << endl;
  }

  if(m_load_balance && a_realm == m_particle_realm){
    this->load_balance_particle_realm(a_procs, a_boxes, a_realm, a_grids, a_lmin, a_finest_level);
  }
  else{
    MayDay::Abort("ito_plasma_stepper::load_balance_boxes - shouldn't happen, how did you get here..?");
  }
}

Vector<long int> ito_plasma_stepper::get_checkpoint_loads(const std::string a_realm, const int a_level) const {
  CH_TIME("ito_plasma_stepper::get_checkpoint_loads(...)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper_stepper::get_checkpoint_loads(...)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->get_grids(a_realm)[a_level];
  const int count = dbl.size();

  Vector<long int> loads(count, 0L);
  if(m_load_balance && a_realm == m_particle_realm){

    Vector<RefCountedPtr<ito_solver> > lb_solvers = this->get_lb_solvers();
    
    for (int i = 0; i < lb_solvers.size(); i++){
      Vector<long int> my_loads(count, 0L);
      lb_solvers[i]->compute_loads(my_loads, dbl, a_level);

      for (int i = 0; i < count; i++){
	loads[i] += my_loads[i];
      }
    }
  }
  else{
    loads = time_stepper::get_checkpoint_loads(a_realm, a_level);
  }

  return loads;
}

void ito_plasma_stepper::load_balance_particle_realm(Vector<Vector<int> >&            a_procs,
						     Vector<Vector<Box> >&            a_boxes,
						     const std::string                a_realm,
						     const Vector<DisjointBoxLayout>& a_grids,
						     const int                        a_lmin,
						     const int                        a_finest_level){
  CH_TIME("ito_plasma_stepper::load_balance_particle_realm(...)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper_stepper::load_balance_particle_realm(...)" << endl;
  }
  
  // Decompose the DisjointBoxLayout
  a_procs.resize(1 + a_finest_level);
  a_boxes.resize(1 + a_finest_level);
  
  for (int lvl = a_lmin; lvl <= a_finest_level; lvl++){
    a_procs[lvl] = a_grids[lvl].procIDs();
    a_boxes[lvl] = a_grids[lvl].boxArray();
  }

  // Get the particles that we will use for load balancing. 
  Vector<RefCountedPtr<ito_solver> > lb_solvers = this->get_lb_solvers();

  // Regrid particles onto the "dummy grids" a_grids
  for (int i = 0; i < lb_solvers.size(); i++){
    particle_container<ito_particle>& particles = lb_solvers[i]->get_particles(ito_solver::which_container::bulk);
    
    particles.regrid(a_grids, m_amr->get_domains(), m_amr->get_dx(), m_amr->get_ref_rat(), a_lmin, a_finest_level);

    // If we make superparticles during regrids, do it here so we can better estimate the computational loads for each patch. This way, if a grid is removed the realistic
    // load estimate of the underlying grid(s) is improved.
    if(m_regrid_superparticles){
      particles.sort_particles_by_cell();
      lb_solvers[i]->make_superparticles(ito_solver::which_container::bulk, m_ppc);
      particles.sort_particles_by_patch();
    }
  }

  // Do Load balance each level. 
  for (int lvl = a_lmin; lvl <= a_finest_level; lvl++){
    const int count = a_grids[lvl].size();

    Vector<long int> total_loads(count, 0L);

    for (int i = 0; i < lb_solvers.size(); i++){
      Vector<long int> loads(count, 0L);

      lb_solvers[i]->compute_loads(loads, a_grids[lvl], lvl);

      for (int i = 0; i < count; i++){
	total_loads[i] += loads[i];
      }
    }

    // Now add the "constant" loads.
    for (int ibox = 0; ibox < a_boxes[lvl].size(); ibox++){
      total_loads[ibox] += lround(m_load_ppc)*a_boxes[lvl][ibox].numPts();
    }

    // Do the friggin load balancing. 
    LoadBalance(a_procs[lvl], total_loads, a_boxes[lvl]);
  }

  // Go back to "pre-regrid" mode so we can get particles to the correct patches after load balancing. 
  for (int i = 0; i < lb_solvers.size(); i++){
    particle_container<ito_particle>& particles = lb_solvers[i]->get_particles(ito_solver::which_container::bulk);
    particles.pre_regrid(a_lmin);
  }
}

Vector<RefCountedPtr<ito_solver> > ito_plasma_stepper::get_lb_solvers() const {
  CH_TIME("ito_plasma_stepper::get_lb_solvers()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::get_lb_solvers()" << endl;
  }

  Vector<RefCountedPtr<ito_solver> > lb_solvers;

  if(m_load_balance_idx < 0){
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<ito_solver>& solver = solver_it();
      
      lb_solvers.push_back(solver);
    }
  }
  else {
    RefCountedPtr<ito_solver>& solver = m_ito->get_solvers()[m_load_balance_idx];
    lb_solvers.push_back(solver);
  }

  return lb_solvers;
}

void ito_plasma_stepper::compute_EdotJ_source(const Real a_dt){
  CH_TIME("ito_plasma_stepper::compute_EdotJ_source(a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_EdotJ_source(a_dt)" << endl;
  }

  // Swap between these two. 
  if(m_nwo_reactions){
    //    this->compute_EdotJ_source_nwo();
    this->compute_EdotJ_source_nwo2(a_dt);
  }
  else{
    this->compute_EdotJ_source();
  }
}

void ito_plasma_stepper::compute_EdotJ_source(){
  CH_TIME("ito_plasma_stepper::compute_EdotJ_source()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_EdotJ_source()" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver   = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();
    const int q   = species->get_charge();

    data_ops::set_value(m_energy_sources[idx], 0.0);

    // Do mobile contribution. 
    if(q != 0 && solver->is_mobile()){

      // Drift contribution
      solver->deposit_conductivity(m_particle_scratch1, solver->get_particles(ito_solver::which_container::bulk)); // Deposit mu*n
      data_ops::copy(m_particle_scratchD, m_particle_E); // Could use m_particle_E or solver's m_velo_func here, but m_velo_func = +/- E (depends on q)
      
      data_ops::multiply_scalar(m_particle_scratchD, m_particle_scratch1);        // m_particle_scratchD = mu*n*E
      data_ops::dot_prod(m_particle_scratch1, m_particle_E, m_particle_scratchD); // m_particle_scratch1 = mu*n*E*E
      data_ops::incr(m_energy_sources[idx], m_particle_scratch1, 1.0);            // a_source[idx] += mu*n*E*E
    }

    // Diffusive contribution
    if(q != 0 && solver->is_diffusive()){

      // Compute the negative gradient of the diffusion term
      solver->deposit_diffusivity(m_particle_scratch1, solver->get_particles(ito_solver::which_container::bulk));
      m_amr->compute_gradient(m_particle_scratchD, m_particle_scratch1, m_particle_realm, m_phase);
      data_ops::scale(m_particle_scratchD, -1.0); // scratchD = -grad(D*n)
      
      data_ops::dot_prod(m_particle_scratch1, m_particle_scratchD, m_particle_E); // m_particle_scratch1 = -E*grad(D*n)
      data_ops::incr(m_energy_sources[idx], m_particle_scratch1, 1.0);            // a_source[idx]
    }
    
    if (q != 0 && (solver->is_mobile() || solver->is_diffusive())){
      data_ops::scale(m_energy_sources[idx], Abs(q)*units::s_Qe);
    }
  }
}

void ito_plasma_stepper::compute_EdotJ_source_nwo(){
  CH_TIME("ito_plasma_stepper::compute_EdotJ_source_nwo()");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_EdotJ_source_nwo()" << endl;
  }

  data_ops::set_value(m_EdotJ, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver   = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();
    const int q   = species->get_charge();

    // Do mobile contribution. Computes Z*e*E*mu*n*E*E
    if(q != 0 && solver->is_mobile()){
      solver->deposit_conductivity(m_particle_scratch1, solver->get_particles(ito_solver::which_container::bulk)); // Deposit mu*n
      m_fluid_scratch1.copy(m_particle_scratch1);                                 // Copy mu*n to fluid realm
      data_ops::copy(m_fluid_scratchD, m_fluid_E);                                // m_fluid_scratchD = E
      data_ops::multiply_scalar(m_fluid_scratchD, m_fluid_scratch1);              // m_fluid_scratchD = E*mu*n
      data_ops::dot_prod(m_fluid_scratch1, m_fluid_E, m_fluid_scratchD);          // m_particle_scratch1 = E.dot.(E*mu*n)
      data_ops::scale(m_fluid_scratch1, Abs(q)*units::s_Qe);                      // m_particle_scratch1 = Z*e*mu*n*E*E

      m_amr->average_down(m_fluid_scratch1, m_fluid_realm, m_phase);
      m_amr->interp_ghost(m_fluid_scratch1, m_fluid_realm, m_phase);
      data_ops::plus(m_EdotJ, m_fluid_scratch1, 0, idx, 1);                       // a_source[idx] += Z*e*mu*n*E*E
    }

    // Diffusive contribution. Computes -Z*e*E*grad(D*n)
    if(q != 0 && solver->is_diffusive()){
      solver->deposit_diffusivity(m_particle_scratch1, solver->get_particles(ito_solver::which_container::bulk));            // Deposit D*n
      m_fluid_scratch1.copy(m_particle_scratch1);                                           // Copy D*n to fluid realm
      m_amr->compute_gradient(m_fluid_scratchD, m_fluid_scratch1, m_fluid_realm, m_phase);  // scratchD = grad(D*n)
      data_ops::scale(m_fluid_scratchD, -1.0);                                              // scratchD = -grad(D*n)
      data_ops::dot_prod(m_fluid_scratch1,  m_fluid_scratchD, m_fluid_E);                   // scratch1 = -E.dot.grad(D*n)
      data_ops::scale(m_fluid_scratch1, Abs(q)*units::s_Qe);                                // scratch1 = -Z*e*E*grad(D*n)

      m_amr->average_down(m_fluid_scratch1, m_fluid_realm, m_phase);
      m_amr->interp_ghost(m_fluid_scratch1, m_fluid_realm, m_phase);
      
      data_ops::plus(m_EdotJ, m_fluid_scratch1, 0, idx, 1);                                 // source  += -Z*e*E*grad(D*n)
    }
  }
}

void ito_plasma_stepper::compute_EdotJ_source_nwo2(const Real a_dt){
  CH_TIME("ito_plasma_stepper::compute_EdotJ_source_nwo2(a_dt)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_EdotJ_source_nwo2(a_dt)" << endl;
  }

  data_ops::set_value(m_EdotJ, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver   = solver_it();
    RefCountedPtr<ito_species>& species = solver->get_species();

    const int idx = solver_it.get_solver();
    const int q   = species->get_charge();

    const bool mobile    = solver->is_mobile();
    const bool diffusive = solver->is_diffusive();

    particle_container<ito_particle>& particles = solver->get_particles(ito_solver::which_container::bulk);

    const DepositionType::Which deposition = solver->get_deposition();

    if((mobile || diffusive) && q != 0){

      // We will interpolate m_particle_E onto particle velocity vectors.
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	  const EBCellFAB& E     = (*m_particle_E[lvl])[dit()];
	  const EBISBox& ebisbox = E.getEBISBox();
	  const FArrayBox& Efab  = E.getFArrayBox();
	  const RealVect dx      = m_amr->get_dx()[lvl]*RealVect::Unit;
	  const RealVect origin  = m_amr->get_prob_lo();
	  const Box box          = dbl[dit()];

	  List<ito_particle>& particleList = particles[lvl][dit()].listItems();

	  // This interpolates the velocity function on to the particle velocities
	  EBParticleInterp meshInterp(box, ebisbox, dx, origin, true);
	  meshInterp.interpolateVelocity(particleList, Efab, deposition);

	  // Go through the particles and set their mass to E.dot(X^new - X^old)
	  for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
	    ito_particle& p = lit();

	    const Real m         = p.mass();
	    const RealVect& v    = p.velocity(); // Actually = E(X^new)
	    const RealVect& Xnew = p.position();
	    const RealVect& Xold = p.oldPosition();

	    p.tmp()  = m;
	    p.mass() = m*PolyGeom::dot(v, Xnew-Xold);
	  } 
	}
      }

      // Deposit the result
      solver->deposit_particles(m_particle_scratch1, particles, deposition);
      m_fluid_scratch1.copy(m_particle_scratch1);

      // Scale by Qe/dt to make it Joule/dt. Then add to correct index
      data_ops::scale(m_fluid_scratch1, q*units::s_Qe/a_dt);
      data_ops::plus(m_EdotJ, m_fluid_scratch1, 0, idx, 1);

      // Set p.mass() back to the original value
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	const DisjointBoxLayout& dbl = m_amr->get_grids(m_particle_realm)[lvl];
	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	  List<ito_particle>& particleList = particles[lvl][dit()].listItems();

	  for (ListIterator<ito_particle> lit(particleList); lit.ok(); ++lit){
	    ito_particle& p = lit();
	    p.mass() = p.tmp();
	  }
	}
      }
    }
  }
}
