/*!
  @file poisson_solver.cpp
  @brief Implementation of poisson_solver.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "poisson_solver.H"
#include "MFAliasFactory.H"


poisson_solver::poisson_solver(){
}

poisson_solver::~poisson_solver(){
}

void poisson_solver::setup_base(const computational_geometry* const a_compgeom){
  m_compgeom = a_compgeom;


  this->alias_internals();
}

void poisson_solver::solve() {
  this->solve(m_state, m_source);
}

void poisson_solver::solve(MFAMRCellData& a_state){
  this->solve(a_state, m_source);
}

MFAMRCellData& poisson_solver::get_state(){
  return m_state;
}

MFAMRCellData& poisson_solver::get_source(){
  return m_source;
}

EBAMRCellData& poisson_solver::get_state_phase(Phase::WhichPhase a_phase){
  this->alias_internals();
  
  if(a_phase == Phase::Gas){
    return m_state_gas;
  }
  else if(a_phase == Phase::Solid){
    return m_state_solid;
  }
  else{
    MayDay::Abort("poisson_solver::get_state_phase - unknown phase");
  }

}

EBAMRCellData& poisson_solver::get_source_phase(Phase::WhichPhase a_phase){
  this->alias_internals();

  if(a_phase == Phase::Gas){
    return m_state_gas;
  }
  else if(a_phase == Phase::Solid){
    return m_state_solid;
  }
  else{
    MayDay::Abort("poisson_solver::get_source_phase - unknown phase");
  }
}

EBAMRIVData& poisson_solver::get_jump(){
  return m_jump;
}

//
Real poisson_solver::get_time() const{
  return m_time;
}

//
void poisson_solver::set_time(const Real a_time) {
  m_time = a_time;
}

//
void poisson_solver::alias_internals(){
  aliasMF(m_state_gas,   Phase::Gas,   m_state);
  aliasMF(m_state_solid, Phase::Solid, m_state);

  aliasMF(m_source_gas,   Phase::Gas,   m_source);
  aliasMF(m_source_solid, Phase::Solid, m_source);
}
