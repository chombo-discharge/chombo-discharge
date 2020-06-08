/*!
  @file   ito_plasma_air2.cpp
  @brief  Implementation of ito_plasma_air2.H
  @author Robert Marskar
  @date   June 2020
*/

#include "ito_plasma_air2.H"

using namespace physics::ito_plasma;

ito_plasma_air2::ito_plasma_air2(){
  m_num_ito_species = 2;
}

ito_plasma_air2::~ito_plasma_air2(){

}

Real ito_plasma_air2::compute_alpha(const RealVect a_E) const {
  return 0.0;
}

Vector<RealVect> ito_plasma_air2::compute_ito_velocities(const Real         a_time,
							 const RealVect     a_pos,
							 const RealVect     a_E,
							 const Vector<Real> a_cdr_densities) const {
  Vector<RealVect> velo(m_num_ito_species, RealVect::Unit);

  return velo;
}
