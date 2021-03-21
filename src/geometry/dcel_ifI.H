/*!
  @file   dcel_ifI.H
  @brief  Implementation of dcel_if.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_IFI_H_
#define _DCEL_IFI_H_

#include "Vec.H"

#include "dcel_if.H"

using namespace dcel;

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
  
  Vec3 p(a_point[0], a_point[1], a_point[2]);

  T retval = m_mesh->signedDistance(p); // Note that dcel::mesh can return either positive or negative for outside, depending on the orientation of the input normals. 
  
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

