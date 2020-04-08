/*!
  @file   brownian_walker_stepper.cpp
  @brief  Implementation of brownian_walker_stepper.H
  @author Robert Marskar
  @data   March 2020
*/

#include "brownian_walker_stepper.H"
#include "brownian_walker_species.H"

#include <ParmParse.H>

using namespace physics::brownian_walker;

brownian_walker_stepper::brownian_walker_stepper(){
  ParmParse pp("brownian_walker");

  pp.get("diffco",   m_diffco);
  pp.get("omega",    m_omega);
  pp.get("cfl",      m_cfl);
}

brownian_walker_stepper::brownian_walker_stepper(RefCountedPtr<ito_solver>& a_solver) : brownian_walker_stepper() {
  m_solver = a_solver;
}

brownian_walker_stepper::~brownian_walker_stepper(){

}

void brownian_walker_stepper::initial_data(){
  CH_TIME("brownian_walker_stepper::initial_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::initial_data" << endl;
  }
  
  m_solver->initial_data();

  if(m_solver->is_diffusive()){
    m_solver->set_diffco(m_diffco);
  }
  if(m_solver->is_mobile()){
    this->set_velocity();
  }
}

void brownian_walker_stepper::set_velocity(){
  CH_TIME("brownian_walker_stepper::set_velocity");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::set_velocity" << endl;
  }
  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    this->set_velocity(lvl);
  }

  EBAMRCellData& vel = m_solver->get_velo_cell();
  m_amr->average_down(vel, phase::gas);
  m_amr->interp_ghost(vel, phase::gas);
}

void brownian_walker_stepper::set_velocity(const int a_level){
  CH_TIME("brownian_walker_stepper::set_velocity(level)");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::set_velocity(level)" << endl;
  }

  // TLDR: This code goes down to each cell on grid level a_level and sets the velocity to omega*r
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_level];
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box = dbl.get(dit());

    EBCellFAB& vel = (*(m_solver->get_velo_cell())[a_level])[dit()];
    BaseFab<Real>& vel_reg = vel.getSingleValuedFAB();

    vel.setVal(0.0);
      
    // Regular cells
    for (BoxIterator bit(box); bit.ok(); ++bit){
      const IntVect iv = bit();
      const RealVect pos = m_amr->get_prob_lo() + (RealVect(iv) + 0.5*RealVect::Unit)*m_amr->get_dx()[a_level];

      const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
      const Real theta = atan2(pos[1],pos[0]);

      vel_reg(iv,0) = -r*m_omega*sin(theta);
      vel_reg(iv,1) =  r*m_omega*cos(theta);
    }

    // Irregular and multicells
    const EBISBox& ebisbox = vel.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    for (VoFIterator vofit(ebisbox.getIrregIVS(box), ebgraph); vofit.ok(); ++vofit){

      const VolIndex vof = vofit();
      const IntVect iv   = vof.gridIndex();
      const RealVect pos = m_amr->get_prob_lo() + (RealVect(iv) + 0.5*RealVect::Unit)*m_amr->get_dx()[a_level];

      const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
      const Real theta = atan2(pos[1],pos[0]);

      vel(vof,0) = -r*m_omega*sin(theta);
      vel(vof,1) =  r*m_omega*cos(theta);
    }
  }
}

void brownian_walker_stepper::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const{
  CH_TIME("brownian_walker_stepper::write_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::write_checkpoint_data" << endl;
  }

  m_solver->write_checkpoint_level(a_handle, a_lvl);
}

void brownian_walker_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) {
  CH_TIME("brownian_walker_stepper::read_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::read_checkpoint_data" << endl;
  }
  
  m_solver->read_checkpoint_level(a_handle, a_lvl);
}

void brownian_walker_stepper::post_checkpoint_setup() {
  CH_TIME("brownian_walker_stepper::post_checkpoint_setup");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::post_checkpoint_setup" << endl;
  }

  MayDay::Abort("brownian_walker_stepper::post_checkpoint_setup - not implemented");
}

int brownian_walker_stepper::get_num_plot_vars() const {
  CH_TIME("brownian_walker_stepper::get_num_plot_vars");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::get_num_plot_vars" << endl;
  }

  return m_solver->get_num_plotvars();
}

void brownian_walker_stepper::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {
  CH_TIME("brownian_walker_stepper::write_plot_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::write_plot_data" << endl;
  }
  a_plotvar_names.append(m_solver->get_plotvar_names());
  m_solver->write_plot_data(a_output, a_icomp);
}

void brownian_walker_stepper::compute_dt(Real& a_dt, time_code::which_code& a_timecode) {
  CH_TIME("brownian_walker_stepper::compute_dt");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::compute_dt" << endl;
  }
  
  MayDay::Warning("brownian_walker_stepper::compute_dt - not implemented yet");

  a_dt = 0.1;
}

void brownian_walker_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("brownian_walker_stepper::synchronize_solver_times");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::synchronize_solver_times" << endl;
  }
  
  m_solver->set_time(a_step, a_time, a_dt);

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void brownian_walker_stepper::print_step_report() {
  CH_TIME("brownian_walker_stepper::print_step_report");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::print_step_report" << endl;
  }

  // Do nothing
}

bool brownian_walker_stepper::need_to_regrid() {
  return false;
}

void brownian_walker_stepper::cache() {
  CH_TIME("brownian_walker_stepper::deallocate");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::deallocate" << endl;
  }

  m_solver->cache_state();
}

void brownian_walker_stepper::deallocate() {
  CH_TIME("brownian_walker_stepper::deallocate");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::deallocate" << endl;
  }
}

void brownian_walker_stepper::setup_solvers() {
  CH_TIME("brownian_walker_stepper::setup_solvers");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::setup_solvers" << endl;
  }

  m_species = RefCountedPtr<ito_species> (new brownian_walker_species());

  m_solver->set_verbosity(-1);
  m_solver->parse_options();
  m_solver->set_amr(m_amr);
  m_solver->set_species(m_species);
  m_solver->set_amr(m_amr);
  m_solver->set_phase(phase::gas);
  m_solver->set_computational_geometry(m_compgeom);
  m_solver->allocate_internals(); // Allocate some internal storage
}

Real brownian_walker_stepper::advance(const Real a_dt) {
  CH_TIME("brownian_walker_stepper::advance");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::advance" << endl;
  }
  
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    const DisjointBoxLayout& dbl          = m_amr->get_grids()[lvl];
    ParticleData<ito_particle>& particles = *m_solver->get_particles()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      // Create a copy. 
      List<ito_particle>& particleList = particles[dit()].listItems();
      List<ito_particle>  particleCopy = List<ito_particle>(particleList);

      // The list iterator is NOT an indexing iterator but iterates over the list given
      // in the constructor. So, we need one for velocities and one for the copy
      ListIterator<ito_particle> lit(particleList);
      ListIterator<ito_particle> litC(particleCopy);

      
      // Compute velocities 
      m_solver->interpolate_velocities(lvl, dit()); 

      // Half Euler step and evaluate velocity at half step
      for (lit.rewind(); lit; ++lit){ 
	ito_particle& p = particleList[lit];
	p.position() += 0.5*p.velocity()*a_dt;
      }
      m_solver->interpolate_velocities(lvl, dit()); 

      // FInal stage
      for (lit.rewind(), litC.rewind(); lit, litC; ++lit, ++litC){
	ito_particle& p    = particleList[lit];
	ito_particle& oldP = particleCopy[litC];
	p.position() = oldP.position() + p.velocity()*a_dt;
      }
    }

    // Remap outcasts
    particles.gatherOutcast();
    particles.remapOutcast();
  }
  //  m_solver->move_particles_eulerf(a_dt);
  m_solver->deposit_particles();
}

void brownian_walker_stepper::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) {
  CH_TIME("brownian_walker_stepper::regrid");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::regrid" << endl;
  }

  m_solver->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  if(m_solver->is_diffusive()){
    m_solver->set_diffco(m_diffco);
  }
  if(m_solver->is_mobile()){
    this->set_velocity();
  }
}
