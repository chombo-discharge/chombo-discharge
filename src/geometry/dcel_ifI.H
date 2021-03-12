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

using namespace dcel;

template <class T, class BV>
dcel_if<T, BV>::dcel_if(const std::shared_ptr<mesh>& a_mesh, const bool a_flipInside){
  m_mesh       = a_mesh;
  m_flipInside = a_flipInside;
  m_root       = nullptr;
}

template <class T, class BV>
dcel_if<T, BV>::dcel_if(const dcel_if& a_object){
  m_mesh       = a_object.m_mesh;
  m_flipInside = a_object.m_flipInside;
}

template <class T, class BV>
dcel_if<T, BV>::~dcel_if(){

}

template <class T, class BV>
Real dcel_if<T, BV>::value(const RealVect& a_point) const {
  
  Vec3 p(a_point[0], a_point[1], a_point[2]);
  
  T retval = m_mesh->signedDistance(p); // dcel::mesh can return either positive or negative for outside. 
  
  if(m_flipInside){
    retval = -retval;
  }

  return Real(retval);
}

template <class T, class BV>
BaseIF* dcel_if<T, BV>::newImplicitFunction() const {
  return static_cast<BaseIF*> (new dcel_if(*this));
}

template <class T, class BV>
inline
decltype(auto) dcel_if<T, BV>::buildBVH(SortingFunctor sortFunc, CostFunctor costFunc) noexcept{
  
  m_root = std::make_shared<Node>(m_mesh->getFaces());
  
  m_root->sortAndPartitionObjects(sortFunc, costFunc);
}

#endif

