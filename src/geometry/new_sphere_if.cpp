/*!
  @file   new_sphere_if.H
  @brief  Implementationof new_sphere_if.H
  @date   Nov. 2017
  @author Robert Marskar
*/

#include "new_sphere_if.H"

#include "CD_NamespaceHeader.H"

new_sphere_if::new_sphere_if(const RealVect& a_center, const Real& a_radius, const bool& a_fluidInside){
  m_center = a_center;
  m_radius = a_radius;
  m_fluidInside = a_fluidInside;
}

new_sphere_if::new_sphere_if(const new_sphere_if& a_inputIF){
  m_center      = a_inputIF.m_center;
  m_radius      = a_inputIF.m_radius;
  m_fluidInside = a_inputIF.m_fluidInside;
}

new_sphere_if::~new_sphere_if(){
}

Real new_sphere_if::value(const RealVect& a_point) const{

  const RealVect newPos = a_point - m_center;

  Real retval = newPos.vectorLength() - m_radius; // Currently negative inside

  if(!m_fluidInside){
    retval = -retval;
  }
  
  return retval;
}

BaseIF* new_sphere_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new new_sphere_if(*this));
}
#include "CD_NamespaceFooter.H"
