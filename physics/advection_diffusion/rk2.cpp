/*!
  @file   rk2.cpp
  @brief  Implementation of rk2
  @author Robert Marskar
  @data   March 2020
*/

#include <ParmParse.H>

#include "rk2.H"
#include "advection_diffusion_species.H"
#include "data_ops.H"

using namespace physics::advection_diffusion;

rk2::rk2(){
  ParmParse pp("advection_diffusion");

  pp.get("diffco",   m_diffco);
  pp.get("omega",    m_omega);
  pp.get("cfl",      m_cfl);
}

rk2::rk2(RefCountedPtr<cdr_solver>& a_solver) : rk2() {
  m_solver = a_solver;
}

rk2::~rk2(){
  
}

void rk2::setup_solvers(){
  advection_diffusion_stepper::setup_solvers();

  m_amr->allocate(m_k1,  phase::gas, 1);
  m_amr->allocate(m_k2,  phase::gas, 1);
}

Real rk2::advance(const Real a_dt){
  EBAMRCellData& state = m_solver->get_state();
  
  m_solver->compute_divJ(m_k1, state, 0.0);  
  data_ops::incr(state, m_k1, -a_dt);        
  m_solver->compute_divJ(m_k2, state, 0.0);  
  data_ops::incr(state, m_k1,  0.5*a_dt);    
  data_ops::incr(state, m_k2, -0.5*a_dt);

  m_solver->make_non_negative(state);
  
  m_amr->average_down(state, phase::gas);
  m_amr->interp_ghost(state, phase::gas);
  
  return a_dt;
}

void rk2::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){

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

  m_amr->allocate(m_k1,  phase::gas, 1);
  m_amr->allocate(m_k2,  phase::gas, 1);
}
