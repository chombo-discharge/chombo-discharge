/*!
  @file   ito_plasma_physics.cpp
  @brief  Implementation of ito_plasma_physics.H
  @author Robert Marskar
  @date   June 2020
*/

#include "ito_plasma_physics.H"

#include <PolyGeom.H>

using namespace physics::ito_plasma;

ito_plasma_physics::ito_plasma_physics(){
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

RealVect ito_plasma_physics::random_position(const RealVect a_pos,
					     const RealVect a_lo,
					     const RealVect a_hi,
					     const RealVect a_bndryCentroid,
					     const RealVect a_normal,
					     const Real     a_dx,
					     const Real     a_kappa) const {

  RealVect pos;
  if(a_kappa < 1.0){ // Rejection sampling. 
    pos = this->random_position(a_lo, a_hi, a_bndryCentroid, a_normal);
  }
  else{ // Regular cell. Get a position. 
    pos = this->random_position(a_lo, a_hi);
  }

  pos = a_pos + pos*a_dx;

  return pos;
}

RealVect ito_plasma_physics::random_position(const RealVect a_lo,
					     const RealVect a_hi,
					     const RealVect a_bndryCentroid,
					     const RealVect a_normal) const {

  RealVect pos = this->random_position(a_lo, a_hi);
  bool valid   = PolyGeom::dot(pos-a_bndryCentroid, a_normal) >= 0.0;

  while(!valid){
    pos    = this->random_position(a_lo, a_hi);
    valid = PolyGeom::dot(pos-a_bndryCentroid, a_normal) > 0.0;
  }

  return pos;
}

RealVect ito_plasma_physics::random_position(const RealVect a_lo, const RealVect a_hi) const {

  RealVect pos = RealVect::Unit;

  for (int dir = 0; dir < SpaceDim; dir++){
    pos[dir] = a_lo[dir] + 0.5*(1.0 + m_udist11(m_rng))*(a_hi[dir] - a_lo[dir]);
  }

  return pos;
}

RealVect ito_plasma_physics::random_direction() const {
  const Real EPS = 1.E-8;
#if CH_SPACEDIM==2
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = m_udist11(m_rng);
    x2 = m_udist11(m_rng);
    r  = x1*x1 + x2*x2;
  }

  return RealVect(x1,x2)/sqrt(r);
#elif CH_SPACEDIM==3
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1*x1 + x2*x2;
  while(r >= 1.0 || r < EPS){
    x1 = m_udist11(m_rng);
    x2 = m_udist11(m_rng);
    r  = x1*x1 + x2*x2;
  }

  const Real x = 2*x1*sqrt(1-r);
  const Real y = 2*x2*sqrt(1-r);
  const Real z = 1 - 2*r;

  return RealVect(x,y,z);
#endif
}
