/*!
  @file    potential_func.cpp
  @brief   Implementation of potential_func.H
  @author  Robert Marskar
  @date    June 2018
*/

#include "potential_func.H"

potential_func::potential_func(Real (*a_potential)(const Real a_time)){
  m_potential = a_potential;
  m_func = NULL;
}


potential_func::potential_func(Real (*a_potential)(const Real a_time), Real (*a_func)(const RealVect a_pos)){
  m_potential = a_potential;
}

potential_func::~potential_func(){

}

void potential_func::set_time(const Real a_time){
  m_time = a_time;
}

Real potential_func::value(const RealVect& a_point, const int& a_comp) const {
  if(m_func == NULL){
    return m_potential(m_time);
  }
  else {
    return m_func(a_point)*m_potential(m_time);
  }
}

Real potential_func::derivative(const RealVect& a_point, const int& a_comp, const int& a_dir) const {
  MayDay::Abort("poisson_multifluid_gmg::potential_func::derivative - this should not be called. How did you get here?");
}
