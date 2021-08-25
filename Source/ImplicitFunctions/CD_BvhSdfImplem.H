/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BvhSdfImplem.H
  @brief  Implementation of CD_BvhSdf.H
  @author Robert Marskar
*/

#ifndef CD_BvhSdfImplem_H
#define CD_BvhSdfImplem_H

// Std includes
#include <chrono>

// Our includes
#include <CD_Vec.H>
#include <CD_BvhSdf.H>
#include <CD_NamespaceHeader.H>

using namespace Dcel;
using namespace BVH;
using namespace std::chrono;

template <class T, class BV, int K>
BvhSdf<T, BV, K>::BvhSdf(const std::shared_ptr<Node>& a_root, const bool a_flipInside, const Real a_zCoord) {
  m_root       = a_root;
  m_zCoord     = a_zCoord;
  m_flipInside = a_flipInside;

  m_numCalled = 0L;
  m_timespan = std::chrono::duration<double>(0.0);    
}

template <class T, class BV, int K>
BvhSdf<T, BV, K>::BvhSdf(const BvhSdf& a_primitive) {
  m_root       = a_primitive.m_root;
  m_zCoord     = a_primitive.m_zCoord;  
  m_flipInside = a_primitive.m_flipInside;

  m_numCalled = 0L;
  m_timespan = std::chrono::duration<double>(0.0);    
}

template <class T, class BV, int K>
BvhSdf<T, BV, K>::~BvhSdf(){
  if(m_numCalled > 0L){
    pout() << "In file CD_BvhSdfImplem: BvhSdf::~BvhSdf() On destructor: Calls: " << m_numCalled << "\t Tot: " << m_timespan.count() << "\t Avg./Call = " << m_timespan.count()/(1.0*m_numCalled) << endl;
  }
}

template <class T, class BV, int K>
Real BvhSdf<T, BV, K>::value(const RealVect& a_point) const {

  // TLDR: In 2D we ignore z-variations, freezing that coordinate to some value. We then use the bounding volume hierarchy methodology for computing the signed distance. 
#if CH_SPACEDIM==2
  Vec3 p(a_point[0], a_point[1], m_zCoord);
#else
  Vec3 p(a_point[0], a_point[1], a_point[2]);
#endif

  high_resolution_clock::time_point t1 = high_resolution_clock::now();  
  auto d = m_root->pruneTree(p);
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

  if(m_flipInside){
    d = -d;
  }

  m_timespan += time_span;
  m_numCalled++;

  return Real(d);
}

template <class T, class BV, int K>
BaseIF* BvhSdf<T, BV, K>::newImplicitFunction() const {
  return static_cast<BaseIF*> (new BvhSdf(*this));
}

#include <CD_NamespaceFooter.H>

#endif

