/*!
  @file   advection_diffusion_stepper.cpp
  @brief  Implementation of advection_diffusion_stepper
  @author Robert Marskar
  @data   March 2020
*/

#include <ParmParse.H>

#include "advection_diffusion_stepper.H"
#include "advection_diffusion_species.H"
#include "data_ops.H"

using namespace physics::advection_diffusion;

advection_diffusion_stepper::advection_diffusion_stepper(){
  ParmParse pp("advection_diffusion");

  m_phase = phase::gas;

  pp.get("realm",      m_realm);
  pp.get("verbosity",  m_verbosity); 
  pp.get("fhd",        m_fhd); 
  pp.get("diffco",     m_diffco);
  pp.get("omega",      m_omega);
  pp.get("cfl",        m_cfl);
  pp.get("integrator", m_integrator);
}

advection_diffusion_stepper::advection_diffusion_stepper(RefCountedPtr<cdr_solver>& a_solver) : advection_diffusion_stepper() {
  m_solver = a_solver;
}

advection_diffusion_stepper::~advection_diffusion_stepper(){
}

void advection_diffusion_stepper::parse_runtime_options() {
  ParmParse pp("advection_diffusion");

  pp.get("verbosity",  m_verbosity);
  pp.get("cfl",        m_cfl);
  pp.get("integrator", m_integrator);

  m_solver->parse_runtime_options();
}

void advection_diffusion_stepper::setup_solvers(){
  m_species = RefCountedPtr<cdr_species> (new advection_diffusion_species());

  // Solver setup
  m_solver->set_verbosity(m_verbosity);
  m_solver->set_species(m_species);
  m_solver->parse_options();
  m_solver->set_phase(m_phase);
  m_solver->set_amr(m_amr);
  m_solver->set_computational_geometry(m_compgeom);
  m_solver->sanity_check();
  m_solver->set_realm(m_realm);
  
  if(!m_solver->is_mobile() && !m_solver->is_diffusive()){
    MayDay::Abort("advection_diffusion_stepper::setup_solvers - can't turn off both advection AND diffusion");
  }
}

void advection_diffusion_stepper::post_initialize() {
}

void advection_diffusion_stepper::register_realms() {
  m_amr->register_realm(m_realm);
}

void advection_diffusion_stepper::register_operators(){
  m_solver->register_operators();
}

void advection_diffusion_stepper::allocate() {
  m_solver->allocate_internals();

  m_amr->allocate(m_tmp, m_realm, m_phase, 1);
  m_amr->allocate(m_k1,  m_realm, m_phase, 1);
  m_amr->allocate(m_k2,  m_realm, m_phase, 1);
}

void advection_diffusion_stepper::initial_data(){
  m_solver->initial_data();       // Fill initial through the cdr species

  m_solver->set_source(0.0);
  m_solver->set_ebflux(0.0);
  m_solver->set_domain_flux(0.0);
  if(m_solver->is_diffusive()){
    m_solver->set_diffco(m_diffco);
  }
  if(m_solver->is_mobile()){
    this->set_velocity();
  }

}

void advection_diffusion_stepper::set_velocity(){
  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    this->set_velocity(lvl);
  }

  EBAMRCellData& vel = m_solver->get_velo_cell();
  m_amr->average_down(vel, m_realm, m_phase);
  m_amr->interp_ghost(vel, m_realm, m_phase);
}

void advection_diffusion_stepper::set_velocity(const int a_level){
  // TLDR: This code goes down to each cell on grid level a_level and sets the velocity to omega*r
    
  const DisjointBoxLayout& dbl = m_amr->get_grids(m_realm)[a_level];
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
    VoFIterator& vofit = (*m_amr->get_vofit(m_realm, m_phase)[a_level])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){

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

void advection_diffusion_stepper::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const {
  m_solver->write_checkpoint_level(a_handle, a_lvl);
}

void advection_diffusion_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){
  m_solver->read_checkpoint_level(a_handle, a_lvl);
}

void advection_diffusion_stepper::post_checkpoint_setup(){
  m_solver->set_source(0.0);
  m_solver->set_ebflux(0.0);
  m_solver->set_domain_flux(0.0);
  if(m_solver->is_diffusive()){
    m_solver->set_diffco(m_diffco);
  }
  if(m_solver->is_mobile()){
    this->set_velocity();
  }
}

int advection_diffusion_stepper::get_num_plot_vars() const{
  return m_solver->get_num_plotvars();
}

void advection_diffusion_stepper::write_plot_data(EBAMRCellData&       a_output,
						  Vector<std::string>& a_plotvar_names,
						  int&                 a_icomp) const {
  a_plotvar_names.append(m_solver->get_plotvar_names());
  m_solver->write_plot_data(a_output, a_icomp);
}

void advection_diffusion_stepper::compute_dt(Real& a_dt, time_code& a_timecode){
  if(m_integrator == 0){      // Explicit diffusion. 
    const Real dt = m_solver->compute_advection_diffusion_dt();
    
    a_dt       = m_cfl*dt;
    a_timecode = time_code::advection_diffusion;
  }
  else if(m_integrator == 1){ // Implicit diffusion
    const Real dt = m_solver->compute_advection_dt();

    a_dt       = m_cfl*dt;
    a_timecode = time_code::advection;
  }
  else{
    MayDay::Abort("advection_diffusion_stepper::compute_dt - logic bust");
  }
}

Real advection_diffusion_stepper::advance(const Real a_dt){
  EBAMRCellData& state = m_solver->get_state();
  
  if(m_integrator == 0){ //   Use Heun's method
    m_solver->compute_divJ(m_k1, state, 0.0);
    data_ops::copy(m_tmp, state);
    data_ops::incr(m_tmp, m_k1, -a_dt); // m_tmp = phi - dt*div(J)

    m_solver->compute_divJ(m_k2, m_tmp, 0.0);
    data_ops::incr(state, m_k1, -0.5*a_dt);
    data_ops::incr(state, m_k2, -0.5*a_dt); // Done with deterministic update.

    //m_solver->make_non_negative(state);

    // Add random diffusion flux. This is equivalent to a 1st order. Godunov splitting
    if(m_fhd){
      m_solver->GWN_diffusion_source(m_k1, state); // k1 holds random diffusion
      data_ops::incr(state, m_k1, a_dt);
    }
  
    m_amr->average_down(state, m_realm, m_phase);
    m_amr->interp_ghost(state, m_realm, m_phase);
  }
  else if(m_integrator == 1){
    m_solver->compute_divF(m_k1, state, a_dt);
    data_ops::incr(state, m_k1, -a_dt);
    m_amr->average_down(state, m_realm, m_phase);

    if(m_solver->is_diffusive()){
      data_ops::copy(m_k2, state); // Now holds phiOld - dt*div(F)
      if(m_fhd){
	m_solver->GWN_diffusion_source(m_k1, state); // k1 holds random diffusion
	data_ops::incr(m_k2, m_k1, a_dt);
      }
      data_ops::set_value(m_k1, 0.0);
      m_solver->advance_euler(state, m_k2, m_k1, a_dt);
    }

    m_amr->average_down(state, m_realm, m_phase);
    m_amr->average_down(state, m_realm, m_phase);
  }
  else{
    MayDay::Abort("advection_diffusion_stepper - unknown integrator requested");
  }
  
  return a_dt;
}

void advection_diffusion_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt){
  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
  
  m_solver->set_time(a_step, a_time, a_dt);
}

void advection_diffusion_stepper::print_step_report(){

}

bool advection_diffusion_stepper::need_to_regrid(){
  return false;
}

void advection_diffusion_stepper::pre_regrid(const int a_lbase, const int a_old_finest_level){
  m_solver->pre_regrid(a_lbase, a_old_finest_level);
}

void advection_diffusion_stepper::deallocate(){
  m_solver->deallocate_internals();
}

void advection_diffusion_stepper::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){

  // Regrid CDR solver
  m_solver->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  m_solver->set_source(0.0);
  m_solver->set_ebflux(0.0);
  m_solver->set_domain_flux(0.0);
  if(m_solver->is_diffusive()){
    m_solver->set_diffco(m_diffco);
  }
  if(m_solver->is_mobile()){
    this->set_velocity();
  }


  // Allocate memory for RK steps
  m_amr->allocate(m_tmp, m_realm, m_phase, 1);
  m_amr->allocate(m_k1,  m_realm, m_phase, 1);
  m_amr->allocate(m_k2,  m_realm, m_phase, 1);
}

void advection_diffusion_stepper::post_regrid() {
}
