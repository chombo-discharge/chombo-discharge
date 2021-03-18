/*!
  @file   bvh_ifI.H
  @brief  Implementation of bvh_if.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _BVH_IFI_H_
#define _BVH_IFI_H_

#include "dcel_vec.H"
#include "bvh_if.H"

#include <chrono>

using namespace dcel;
using namespace dcel::BVH;
using namespace std::chrono;

template <class T>
using PrimitiveList = std::vector<std::shared_ptr<const faceT<T> > >; 

template <class T, class BV>
StopFunctionT<T, faceT<T>, BV> bvh_if<T, BV>::defaultStopFunction = [](const NodeT<T, faceT<T>, BV>& a_node){
  const auto& primitives = a_node.getPrimitives();
  const int depth        = a_node.getDepth();

  return primitives.size() < 4;
};

template <class T, class BV>
BVConstructorT<faceT<T>, BV> bvh_if<T, BV>::defaultBVConstructor = [](const PrimitiveList<T>& a_primitives){
  std::vector<Vec3T<T> > coordinates;

  for (const auto& f : a_primitives){
    const auto faceCoordinates = f->getAllVertexCoordinates();
    coordinates.insert(coordinates.end(), faceCoordinates.begin(), faceCoordinates.end());
  }

  return BV(coordinates);
};

template <class T, class BV>
PartitionFunctionT<faceT<T> > bvh_if<T, BV>::defaultPartitionFunction = [](const PrimitiveList<T>& a_primitives){
  Vec3T<T> minCorner = a_primitives.front()->getCentroid();
  Vec3T<T> maxCorner = a_primitives.front()->getCentroid();

  for(const auto& p : a_primitives){
    minCorner = min(minCorner, p->getCentroid());
    maxCorner = max(minCorner, p->getCentroid());
  }

  const Vec3T<T> delta = maxCorner - minCorner;
  const int splitDir   = delta.maxDir(true);

  PrimitiveList<T> sortedPrimitives(a_primitives);
  auto comparator = [splitDir](const std::shared_ptr<const faceT<T> >& f1, const std::shared_ptr<const faceT<T> >& f2){
    return f1->getCentroid()[splitDir] < f2->getCentroid()[splitDir];
  };
  
  std::sort(sortedPrimitives.begin(), sortedPrimitives.end(), comparator);

  const int splitIndex = (sortedPrimitives.size()-1)/2;

  PrimitiveList<T> lPrims(sortedPrimitives.begin(), sortedPrimitives.begin() + splitIndex+1);
  PrimitiveList<T> rPrims(sortedPrimitives.begin() + splitIndex + 1, sortedPrimitives.end());
  
  return std::make_pair(lPrims, rPrims);
};

template <class T, class BV>
PartitionFunctionT<faceT<T> > bvh_if<T, BV>::partitionMinimumOverlap = [](const PrimitiveList<T>& a_primitives){

  const int splitIndex    = (a_primitives.size() - 1)/2;

  T minOverlap = std::numeric_limits<T>::infinity();

  int bestDir;

  PrimitiveList<T> sortedPrimitives(a_primitives);
  
  // Identify the best direction by the minimum overlap
  for (int dir = 0; dir < 3; dir++){

    std::sort(sortedPrimitives.begin(), sortedPrimitives.end(),
	      [dir](const std::shared_ptr<const faceT<T> >& f1, const std::shared_ptr<const faceT<T> >& f2){
		return f1->getCentroid()[dir] < f2->getCentroid()[dir];
	      });

    PrimitiveList<T> lPrims(sortedPrimitives.begin(), sortedPrimitives.begin() + splitIndex+1);
    PrimitiveList<T> rPrims(sortedPrimitives.begin() + splitIndex + 1, sortedPrimitives.end());

    const BV lBV = bvh_if<T, BV>::defaultBVConstructor(lPrims);
    const BV rBV = bvh_if<T, BV>::defaultBVConstructor(rPrims);

    const T curOverlap = getOverlappingVolume(lBV, rBV);

    if (curOverlap < minOverlap){
      minOverlap = curOverlap;
      bestDir = dir;
    }
  }

  // Build the two halves. 
  std::sort(sortedPrimitives.begin(), sortedPrimitives.end(),
	    [bestDir](const std::shared_ptr<const faceT<T> >& f1, const std::shared_ptr<const faceT<T> >& f2){
	      return f1->getCentroid()[bestDir] < f2->getCentroid()[bestDir];
	    });

  PrimitiveList<T> leftPrims(sortedPrimitives.begin(), sortedPrimitives.begin() + splitIndex+1);
  PrimitiveList<T> rightPrims(sortedPrimitives.begin() + splitIndex + 1, sortedPrimitives.end());

  return std::make_pair(leftPrims, rightPrims);
};

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


  // auto d = m_root->pruneOrdered(p);
  auto d = m_root->pruneOrdered2(p);
  // auto d = m_root->pruneUnordered(p);
  // auto d = m_root->pruneUnordered2(p);
  // auto d = m_root->prunePriorityQueue(p);
  // auto d = m_root->prunePriorityQueue2(p);

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

#endif

