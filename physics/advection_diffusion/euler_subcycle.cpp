/*!
  @file   euler_subcycle.cpp
  @brief  Implementation of euler_subcycle
  @author Robert Marskar
  @data   March 2020
*/

#include <ParmParse.H>

#include "euler_subcycle.H"
#include "advection_diffusion_species.H"
#include "data_ops.H"

using namespace physics::advection_diffusion;

euler_subcycle::euler_subcycle(){
  ParmParse pp("advection_diffusion");

  pp.get("diffco",   m_diffco);
  pp.get("omega",    m_omega);
  pp.get("cfl",      m_cfl);
}

euler_subcycle::euler_subcycle(RefCountedPtr<cdr_solver>& a_solver) : euler_subcycle() {
  m_solver = a_solver;
}

euler_subcycle::~euler_subcycle(){
  
}

void euler_subcycle::setup_solvers(){
  advection_diffusion_stepper::setup_solvers();

  // Allocate memory for RK steps
  m_amr->allocate(m_tmp, phase::gas, 1);
  m_amr->allocate(m_k1,  phase::gas, 1);
  m_amr->allocate(m_k2,  phase::gas, 1);
}

Real euler_subcycle::advance(const Real a_dt){

  // Use Heun's method
  EBAMRCellData& state = m_solver->get_state();
  m_solver->compute_divJ(m_k1, state, 0.0);

  data_ops::copy(m_tmp, state);
  data_ops::incr(m_tmp, m_k1, -a_dt); // m_tmp = phi - dt*div(J)
  m_solver->make_non_negative(m_tmp);

  m_solver->compute_divJ(m_k2, m_tmp, 0.0);

  data_ops::incr(state, m_k1, -0.5*a_dt);
  data_ops::incr(state, m_k2, -0.5*a_dt);

  m_solver->make_non_negative(state);
  
  m_amr->average_down(state, phase::gas);
  m_amr->interp_ghost(state, phase::gas);

  
  return a_dt;
}

void euler_subcycle::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){

  // Regrid CDR solver
  m_solver->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  m_solver->set_source(0.0);
  m_solver->set_ebflux(0.0);
  if(m_solver->is_diffusive()){
    m_solver->set_diffco(m_diffco);
  }
  if(m_solver->is_mobile()){
    this->set_velocity();
  }

  // Allocate memory for RK steps
  m_amr->allocate(m_tmp, phase::gas, 1);
  m_amr->allocate(m_k1,  phase::gas, 1);
  m_amr->allocate(m_k2,  phase::gas, 1);
}
