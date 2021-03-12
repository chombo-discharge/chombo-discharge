/*!
  @file   dcel_ifI.H
  @brief  Implementation of dcel_if.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_IFI_H_
#define _DCEL_IFI_H_

#include "dcel_vec.H"
#include "dcel_if.H"

template <class T>
dcel_if<T>::dcel_if(const std::shared_ptr<mesh>& a_mesh, const bool a_flipInside){
  m_mesh       = a_mesh;
  m_flipInside = a_flipInside;
}

template <class T>
dcel_if<T>::dcel_if(const dcel_if& a_object){
  m_mesh       = a_object.m_mesh;
  m_flipInside = a_object.m_flipInside;
}

template <class T>
dcel_if<T>::~dcel_if(){

}

template <class T>
Real dcel_if<T>::value(const RealVect& a_point) const {
  
  dcel::Vec3T<T> p(a_point[0], a_point[1], a_point[2]);
  
  T retval = m_mesh->signedDistance(p); // dcel::mesh can return either positive or negative for outside. 
  
  if(m_flipInside){
    retval = -retval;
  }

  return Real(retval);
}

template <class T>
BaseIF* dcel_if<T>::newImplicitFunction() const {
  return static_cast<BaseIF*> (new dcel_if(*this));
}

#endif
