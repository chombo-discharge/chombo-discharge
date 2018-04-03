/*!
  @brief  bv_sphere.cpp
  @brief  Implementation of bv_sphere.H
  @author Robert Marskar
  @date   Apr. 2018
  @todo   Segregate implementation
*/

#include "bv_sphere.H"

bv_sphere::bv_sphere(){
  
}

bv_sphere::bv_sphere(const RealVect a_center, const Real a_radius){
  m_center = a_center;
  m_radius = a_radius;
}

bv_sphere::~bv_sphere(){
    
}

RealVect bv_sphere::get_center() const {
  return m_center;
}

Real bv_sphere::get_radius() const {
  return m_radius;
}
