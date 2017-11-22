/*!
  @file new_sphere_if.H
  @brief Implementationof new_sphere_if.H
  @date Nov. 2017
  @author Robert Marskar
*/

#include "new_sphere_if.H"

new_sphere_if::new_sphere_if(const RealVect& a_center, const Real& a_radius, const bool& a_inside){

  //
  m_center = a_center;
  m_radius = a_radius;
  m_inside = a_inside;
}

new_sphere_if::new_sphere_if(const new_sphere_if& a_inputIF){
  m_center = a_inputIF.m_center;
  m_radius = a_inputIF.m_radius;
  m_inside = a_inputIF.m_inside;
}

new_sphere_if::~new_sphere_if(){
}

Real new_sphere_if::value(const RealVect& a_point) const{

  const RealVect newPos = a_point - m_center;

  Real retval = m_radius - newPos.vectorLength();

  if(m_inside){
    retval = -retval;
  }
  
  return retval;
}

BaseIF* new_sphere_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new new_sphere_if(*this));
}

