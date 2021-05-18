/*!
  @file   bvh_ifI.H
  @brief  Implementation of bvh_if.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _BVH_IFI_H_
#define _BVH_IFI_H_

#include "Vec.H"

#include "bvh_if.H"

#include <chrono>

#include "CD_NamespaceHeader.H"

using namespace dcel;
using namespace BVH;
using namespace std::chrono;

template <class T, class BV>
bvh_if<T, BV>::bvh_if(){
  m_numCalled = 0L;
  m_timespan = duration<double>(0.0);
}

template <class T, class BV>
bvh_if<T, BV>::bvh_if(const std::shared_ptr<Node>& a_root, const bool a_flipInside) : bvh_if() {
  m_root       = a_root;
  m_flipInside = a_flipInside;
}

template <class T, class BV>
bvh_if<T, BV>::bvh_if(const bvh_if& a_primitive) : bvh_if() {
  m_root       = a_primitive.m_root;
  m_flipInside = a_primitive.m_flipInside;
}

template <class T, class BV>
bvh_if<T, BV>::~bvh_if(){
  if(m_numCalled > 0L){
    std::cout << "On destructor: Calls: " << m_numCalled << "\t Tot: " << m_timespan.count() << "\t Avg./Call = " << m_timespan.count()/(1.0*m_numCalled) << std::endl;
  }
}

template <class T, class BV>
Real bvh_if<T, BV>::value(const RealVect& a_point) const {
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  Vec3 p(a_point[0], a_point[1], a_point[2]); 

  auto d = m_root->pruneOrdered2(p);

  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

  if(m_flipInside){
    d = -d;
  }

  m_timespan += time_span;
  m_numCalled++;

  return Real(d);
}

template <class T, class BV>
BaseIF* bvh_if<T, BV>::newImplicitFunction() const {
  return static_cast<BaseIF*> (new bvh_if(*this));
}
#include "CD_NamespaceFooter.H"

#endif

