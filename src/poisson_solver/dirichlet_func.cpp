/*!
  @file    dirichlet_func.cpp
  @brief   Implementation of dirichlet_func.H
  @author  Robert Marskar
  @date    June 2018
*/

#include "dirichlet_func.H"

dirichlet_func::dirichlet_func(Real (*a_potential)(const Real a_time),
			       Real (*a_func)(const RealVect a_pos),
			       const RealVect a_origin){
  m_potential = a_potential;
  m_func      = a_func;
  m_origin    = a_origin;
}

dirichlet_func::~dirichlet_func(){

}

void dirichlet_func::set_time(const Real a_time){
  m_time = a_time;
}

Real dirichlet_func::value(const RealVect& a_point, const int& a_comp) const {
  // m_func wants physical coordinates but a_point is the computational coordinate
  return m_func(a_point + m_origin)*m_potential(m_time);
}

Real dirichlet_func::derivative(const RealVect& a_point, const int& a_comp, const int& a_dir) const {
  MayDay::Abort("poisson_multifluid_gmg::dirichlet_func::derivative - this should not be called. How did you get here?");

  return 0.0;
}
