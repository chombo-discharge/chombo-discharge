/*!
  @file   ito_plasma_physics.cpp
  @brief  Implementation of ito_plasma_physics.H
  @author Robert Marskar
  @date   June 2020
*/

#include "ito_plasma_physics.H"
#include "units.H"

#include <PolyGeom.H>

using namespace physics::ito_plasma;

ito_plasma_physics::ito_plasma_physics(){

  // Initialize RNGs
  m_udist11 = std::uniform_real_distribution<Real>(-1., 1.);
  m_udist01 = std::uniform_real_distribution<Real>(0., 1.);

  m_reactions.clear();
  m_photo_reactions.clear();

  m_Ncrit = 25;
  m_eps   = 0.03;
}

ito_plasma_physics::~ito_plasma_physics(){
}

const Vector<RefCountedPtr<ito_species> >& ito_plasma_physics::get_ito_species() const { 
  return m_ito_species; 
}

const Vector<RefCountedPtr<rte_species> >& ito_plasma_physics::get_rte_species() const {
  return m_rte_species;
}

int ito_plasma_physics::get_num_ito_species() const{
  return m_ito_species.size();
}

int ito_plasma_physics::get_num_rte_species() const {
  return m_rte_species.size();
}

Real ito_plasma_physics::initial_sigma(const Real a_time, const RealVect a_pos) const {
  return 0.0;
}



