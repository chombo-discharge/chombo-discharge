/*!
  @file plasma_engine.cpp
  @brief Implementation of plasma_engine.H
  @author Robert Marskar
  @date Nov. 2017
*/

#include "plasma_engine.H"

plasma_engine::plasma_engine(){
  CH_TIME("plasma_engine::plasma_engine(weak)");
  if(m_verbosity > 2){
    pout() << "plasma_engine::plasma_engine(weak)" << endl;
  }
}

plasma_engine::plasma_engine(const RefCountedPtr<computational_geometry>& a_compgeom,
			     const RefCountedPtr<plasma_kinetics>&        a_plaskin){
  CH_TIME("plasma_engine::plasma_engine(full)");
  if(m_verbosity > 2){
    pout() << "plasma_engine::plasma_engine(full)" << endl;
  }

  this->set_computational_geometry(a_compgeom);
  this->set_plasma_kinetics(a_plaskin);

  this->allocate_wall_bc();
}

plasma_engine::~plasma_engine(){
  CH_TIME("plasma_engine::~plasma_engine");
}

void plasma_engine::allocate_wall_bc(){
  m_wallbc.resize(2*SpaceDim);
  for (int i = 0; i < 2*SpaceDim; i++){
    m_wallbc[i] = RefCountedPtr<wall_bc> (NULL);
  }
}

void plasma_engine::set_verbosity(const int a_verbosity){
  CH_TIME("plasma_engine::set_verbosity");
  if(m_verbosity > 5){
    pout() << "plasma_engine::set_verbosity" << endl;
  }
  m_verbosity = a_verbosity;
}

void plasma_engine::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("plasma_engine::set_computational_geometry");
  if(m_verbosity > 3){
    pout() << "plasma_engine::set_computational_geometry" << endl;
  }
  m_compgeom = a_compgeom;
}

void plasma_engine::set_plasma_kinetics(const RefCountedPtr<plasma_kinetics>& a_plaskin){
  m_plaskin = a_plaskin;
}

void plasma_engine::setup_fresh(){
  CH_TIME("plasma_engine::setup_fresh");
  if(m_verbosity > 2){
    pout() << "plasma_engine::setup_fresh" << endl;
  }

  // Do a sanity check before doing anything actually expensive.
  this->sanity_check();

  // This stuff should come in through amr
  const int nCells = 512;
  ProblemDomain probdom(IntVect::Zero, (nCells - 1)*IntVect::Unit);
  const Real& finestdx = (m_physdom->get_prob_lo()[0] - m_physdom->get_prob_hi()[0])/nCells;
  
  m_compgeom->build_geometries(*m_physdom, probdom, finestdx, 8);
}

void plasma_engine::set_dirichlet_wall_bc(const int a_dir, Side::LoHiSide a_side, const Potential::GroundLive a_live){
  CH_TIME("plasma_engine::set_dirichlet_wall_bc");
  if(m_verbosity > 3){
    pout() << "plasma_engine::set_dirichlet_wall_bc" << endl;
  }

  const int idx = this->map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, WallBC::Dirichlet));
  m_wallbc[idx]->set_live(a_live);
}

void plasma_engine::set_neumann_wall_bc(const int a_dir, Side::LoHiSide a_side, const Real a_value){
  CH_TIME("plasma_engine::set_neumann_wall_bc");
  if(m_verbosity > 3){
    pout() << "plasma_engine::set_neumann_wall_bc" << endl;
  }

  const int idx = this->map_bc(a_dir, a_side);
  m_wallbc[idx] = RefCountedPtr<wall_bc> (new wall_bc(a_dir, a_side, WallBC::Neumann));
  m_wallbc[idx]->set_value(a_value);
}

void plasma_engine::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("plasma_engine::set_physical_domain");
  if(m_verbosity > 3){
    pout() << "plasma_engine::set_physical_domain" << endl;
  }
  m_physdom = a_physdom;
}

void plasma_engine::sanity_check(){
  CH_TIME("plasma_engine::sanity_check");
  if(m_verbosity > 2){
    pout() << "plasma_engine::sanity_check" << endl;
  }
  
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sideit; sideit.ok(); ++sideit){
      if(m_wallbc[map_bc(dir, sideit())].isNull()){
	MayDay::Abort("computational_geometry::sanity_check() failed. Wall BC has not been set properly");
      }
    }
  }
}

wall_bc& plasma_engine::get_wall_bc(const int a_dir, Side::LoHiSide a_side) const{
  CH_TIME("plasma_engine::get_wall_bc");
  if(m_verbosity > 2){
    pout() << "plasma_engine::get_wall_bc" << endl;
  }
  return *m_wallbc[this->map_bc(a_dir, a_side)];
}

int plasma_engine::map_bc(const int a_dir, const Side::LoHiSide a_side) const {
  const int iside = (a_side == Side::Lo) ? 0 : 1;

  return 2*a_dir + iside;
}
