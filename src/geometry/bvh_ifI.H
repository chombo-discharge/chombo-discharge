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

using namespace dcel;
using namespace dcel::BVH;

template <class T, class BV>
StopFunctionT<T, faceT<T>, BV> bvh_if<T, BV>::defaultStopFunction = [](const NodeT<T, faceT<T>, BV>& a_node){
  const auto& primitives = a_node.getPrimitives();
  const int depth        = a_node.getDepth();

  return primitives.size() < 2;
};

template <class T, class BV>
BVConstructorT<T, faceT<T>, BV> bvh_if<T, BV>::defaultBVConstructor = [](const NodeT<T, faceT<T>, BV>& a_node){
  const auto& a_primitives = a_node.getPrimitives();
  
  std::vector<Vec3T<T> > coordinates;

  for (const auto& f : a_primitives){
    const auto faceCoordinates = f->getAllVertexCoordinates();
    coordinates.insert(coordinates.end(), faceCoordinates.begin(), faceCoordinates.end());
  }

  return BV(coordinates);
};

template <class T, class BV>
PartitionFunctionT<T, faceT<T>, BV> bvh_if<T, BV>::defaultPartitionFunction = [](NodeT<T, faceT<T>, BV>& a_node){

  // TLDR: Very, very naive partitioning of primitives. We search for the spatial direction where the centroids are spaced furthest apart,
  //       and split the vector in half along that direction. 

  auto& primitives = a_node.getPrimitives();

  // Determine along which axis to split.
  Vec3T<T> minCorner = primitives.front()->getCentroid();
  Vec3T<T> maxCorner = primitives.front()->getCentroid();

  for(const auto& f : primitives){
    minCorner = std::min(minCorner, f->getCentroid());
    maxCorner = std::max(maxCorner, f->getCentroid());
  }

  const Vec3T<T> delta = maxCorner - minCorner;
  const int splitDir   = delta.maxDir(true);


  // Sort centroids along splitDir and split the vector in half.
  auto comparator = [splitDir](const std::shared_ptr<const faceT<T> >& f1, const std::shared_ptr<const faceT<T> >& f2){
    return f1->getCentroid()[splitDir] < f2->getCentroid()[splitDir];
  };
  
  std::sort(primitives.begin(), primitives.end(), comparator);
  
  return (primitives.size()-1)/2;
};

template <class T, class BV>
bvh_if<T, BV>::bvh_if(const std::shared_ptr<Node>& a_root, const bool a_flipInside){
  m_root       = a_root;
  m_flipInside = a_flipInside;

}

template <class T, class BV>
bvh_if<T, BV>::bvh_if(const bvh_if& a_primitive){
  m_root       = a_primitive.m_root;
  m_flipInside = a_primitive.m_flipInside;
}

template <class T, class BV>
bvh_if<T, BV>::~bvh_if(){

}

template <class T, class BV>
Real bvh_if<T, BV>::value(const RealVect& a_point) const {
  Vec3 p(a_point[0], a_point[1], a_point[2]); 


#if 1
  auto d = m_root->pruneBranches(p);
#else
  auto d = m_root->pruneQ(p);
#endif

  if(m_flipInside){
    d = -d;
  }

  return Real(d);
}

template <class T, class BV>
BaseIF* bvh_if<T, BV>::newImplicitFunction() const {
  return static_cast<BaseIF*> (new bvh_if(*this));
}

#endif

