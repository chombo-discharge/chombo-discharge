/*!
  @file poisson_solver.cpp
  @brief Implementation of poisson_solver.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "poisson_solver.H"
#include "MFAliasFactory.H"


poisson_solver::poisson_solver(){
  CH_TIME("poisson_solver::set_verbosity");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_verbosity" << endl;
  }
  
  this->set_verbosity(-1);
  this->allocate_wall_bc();
}

poisson_solver::~poisson_solver(){
  
}

void poisson_solver::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("poisson_solver::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_computational_geometry" << endl;
  }

  m_compgeom = a_compgeom;

  this->set_mfis(m_compgeom->get_mfis());
}

void poisson_solver::set_mfis(const RefCountedPtr<mfis>& a_mfis){
  CH_TIME("poisson_solver::set_mfis");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_mfis" << endl;
  }

  m_mfis = a_mfis;
}

void poisson_solver::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("poisson_solver::set_amr");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_amr" << endl;
  }

  m_amr = a_amr;
}

void poisson_solver::allocate_wall_bc(){
  CH_TIME("poisson_solver::poisson_solver(full)");
  if(m_verbosity > 5){
    pout() << "poisson_solver::poisson_solver(full)" << endl;
  }
  m_wallbc.resize(2*SpaceDim);
  for (int i = 0; i < 2*SpaceDim; i++){
    m_wallbc[i] = RefCountedPtr<wall_bc> (NULL);
  }
}

void poisson_solver::set_dirichlet_wall_bc(const int a_dir, Side::LoHiSide a_side, const Potential::GroundLive a_live){
  CH_TIME("poisson_solver::set_dirichlet_wall_bc");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_dirichlet_wall_bc" << endl;
  }

  const int idx = this->map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, WallBC::Dirichlet));
  m_wallbc[idx]->set_live(a_live);
}

void poisson_solver::set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value){
  CH_TIME("poisson_solver::set_neumann_wall_bc");
  if(m_verbosity > 5){
    pout() << "poisson_solver::set_neumann_wall_bc" << endl;
  }

  const int idx = this->map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, WallBC::Neumann));
  m_wallbc[idx]->set_value(a_value);
}

void poisson_solver::set_verbosity(const int a_verbosity){
  m_verbosity = a_verbosity;
}

void poisson_solver::solve() {
  this->solve(m_state, m_source);
}

void poisson_solver::solve(MFAMRCellData& a_state){
  this->solve(a_state, m_source);
}

void poisson_solver::set_time(const Real a_time) {
  m_time = a_time;
}

void poisson_solver::alias_internals(){
//   aliasMF(m_state_gas,   Phase::Gas,   m_state);
//   aliasMF(m_state_solid, Phase::Solid, m_state);

//   aliasMF(m_source_gas,   Phase::Gas,   m_source);
//   aliasMF(m_source_solid, Phase::Solid, m_source);
}

void poisson_solver::sanity_check(){
  CH_TIME("poisson_solver::sanity_check");
  if(m_verbosity > 4){
    pout() << "poisson_solver::sanity_check" << endl;
  }

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sideit; sideit.ok(); ++sideit){
      if(m_wallbc[map_bc(dir, sideit())].isNull()){
	pout() << "poisson_solver::sanity_check() - bc is null at coord = " << dir << ", side = " << sideit() << endl;
  	MayDay::Abort("poisson_solver::sanity_check() failed. Wall BC has not been set properly");
      }
    }
  }
}

int poisson_solver::map_bc(const int a_dir, const Side::LoHiSide a_side) const {
  CH_TIME("poisson_solver::get_wall_bc");
  if(m_verbosity > 999){
    pout() << "poisson_solver::get_wall_bc" << endl;
  }
  const int iside = (a_side == Side::Lo) ? 0 : 1;

  return 2*a_dir + iside;
}

Real poisson_solver::get_time() const{
  return m_time;
}

wall_bc& poisson_solver::get_wall_bc(const int a_dir, Side::LoHiSide a_side) const{
  CH_TIME("poisson_solver::get_wall_bc");
  if(m_verbosity > 5){
    pout() << "poisson_solver::get_wall_bc" << endl;
  }
  return *m_wallbc[this->map_bc(a_dir, a_side)];
}

MFAMRCellData& poisson_solver::get_state(){
  return m_state;
}

MFAMRCellData& poisson_solver::get_source(){
  return m_source;
}

EBAMRCellData& poisson_solver::get_state_phase(Phase::WhichPhase a_phase){
  //  this->alias_internals();
  
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
  //  this->alias_internals();

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
