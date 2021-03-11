/*!
  @file   dcel_if.cpp
  @brief  Implementation of dcel_if.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_vec.H"
#include "dcel_if.H"

dcel_if::dcel_if(const std::shared_ptr<dcel::mesh>& a_mesh, const bool a_flipInside){
  m_mesh       = a_mesh;
  m_flipInside = a_flipInside;
}

dcel_if::dcel_if(const dcel_if& a_object){
  m_mesh       = a_object.m_mesh;
  m_flipInside = a_object.m_flipInside;
}

dcel_if::~dcel_if(){

}

Real dcel_if::value(const RealVect& a_point) const {
  dcel::Vec3<double> p(a_point[0], a_point[1], a_point[2]);
  
  Real retval = m_mesh->signedDistance(p); // dcel::mesh can return either positive or negative for outside. 
  
  if(m_flipInside){
    retval = -retval;
  }

  return retval;
}

BaseIF* dcel_if::newImplicitFunction() const {
  return static_cast<BaseIF*> (new dcel_if(*this));
}
