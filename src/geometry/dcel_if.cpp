/*!
  @file   dcel_if.cpp
  @brief  Implementation of dcel_if.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_if.H"

dcel_if::dcel_if(const std::shared_ptr<dcel::mesh>& a_mesh, const bool a_fluidInside){
  m_mesh        = a_mesh;
  m_fluidInside = a_fluidInside;
}

dcel_if::dcel_if(const dcel_if& a_object){
  m_mesh        = a_object.m_mesh;
  m_fluidInside = a_object.m_fluidInside;
}

dcel_if::~dcel_if(){

}

Real dcel_if::value(const RealVect& a_point) const {
  Real retval = m_mesh->signedDistance(a_point); // dcel::mesh returns positive for outside. 
  
  if(!m_fluidInside){
    retval = -retval;
  }

  return retval;
}

BaseIF* dcel_if::newImplicitFunction() const {
  return static_cast<BaseIF*> (new dcel_if(*this));
}
