/*!
  @file   dcel_if.cpp
  @brief  Implementation of dcel_if.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_if.H"

dcel_if::dcel_if(){

}

dcel_if::dcel_if(const RefCountedPtr<dcel_mesh>& a_mesh, const bool a_inside){
  m_mesh   = a_mesh;
  m_inside = a_inside;
}

dcel_if::dcel_if(const dcel_if& a_object){
  m_mesh   = a_object.m_mesh;
  m_inside = a_object.m_inside;
}

dcel_if::~dcel_if(){

}

Real dcel_if::value(const RealVect& a_point) const{
  Real retval = m_mesh->signed_distance(a_point);
  if(m_inside){
    retval *= -1;
  }

  return retval;
}

BaseIF* dcel_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new dcel_if(*this));
}
