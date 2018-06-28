/*!
  @file    neumann_func.cpp
  @brief   Implementation of neumann_func.H
  @author  Robert Marskar
  @date    June 2018
*/

#include "neumann_func.H"

neumann_func::neumann_func(Real (*a_value)(const Real a_time),
			   Real (*a_func)(const RealVect a_pos),
			   const RealVect a_origin){
  m_value  = a_value;
  m_func   = a_func;
  m_origin = a_origin;
}

neumann_func::~neumann_func(){

}

void neumann_func::set_time(const Real a_time){
  m_time = a_time;
}

Real neumann_func::value(const RealVect& a_point, const int& a_comp) const {
  // m_func wants physical coordinates but a_point is the computational coordinate
  return m_func(a_point + m_origin)*m_value(m_time);
}

Real neumann_func::derivative(const RealVect& a_point, const int& a_comp, const int& a_dir) const {
  MayDay::Abort("neumann_func::derivative - this should not be called. How did you get here?");
}
