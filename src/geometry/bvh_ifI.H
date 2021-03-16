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
  const auto& objects = a_node.getObjects();
  const int depth     = a_node.getDepth();

  return objects.size() < 2;
};

template <class T, class BV>
BVConstructorT<T, faceT<T>, BV> bvh_if<T, BV>::defaultBVConstructor = [](const NodeT<T, faceT<T>, BV>& a_node){
  const auto& a_objects = a_node.getObjects();
  
  std::vector<Vec3T<T> > coordinates;

  for (const auto& f : a_objects){
    const auto faceCoordinates = f->getAllVertexCoordinates();
    coordinates.insert(coordinates.end(), faceCoordinates.begin(), faceCoordinates.end());
  }

  return BV(coordinates);
};

template <class T, class BV>
PartitionFunctionT<T, faceT<T>, BV> bvh_if<T, BV>::defaultPartitionFunction = [](NodeT<T, faceT<T>, BV>& a_node){

  // TLDR: Very, very naive partitioning of objects. We search for the spatial direction where the centroids are spaced furthest apart,
  //       and split the vector in half along that direction. 

  auto& objects = a_node.getObjects();

  // Determine along which axis to split.
  Vec3T<T> minCorner = objects.front()->getCentroid();
  Vec3T<T> maxCorner = objects.front()->getCentroid();

  for(const auto& f : objects){
    minCorner = std::min(minCorner, f->getCentroid());
    maxCorner = std::max(maxCorner, f->getCentroid());
  }

  const Vec3T<T> delta = maxCorner - minCorner;
  const int splitDir   = delta.maxDir(true);


  // Sort centroids along splitDir and split the vector in half.
  auto comparator = [splitDir](const std::shared_ptr<faceT<T> >& f1, const std::shared_ptr<faceT<T> >& f2){
    return f1->getCentroid()[splitDir] < f2->getCentroid()[splitDir];
  };
    
  std::sort(objects.begin(), objects.end(), comparator);
  
  return (objects.size()-1)/2;
};

template <class T, class BV>
bvh_if<T, BV>::bvh_if(const std::shared_ptr<mesh>& a_mesh, const bool a_flipInside){
  m_mesh       = a_mesh;
  m_flipInside = a_flipInside;
  m_root       = nullptr;
}

template <class T, class BV>
bvh_if<T, BV>::bvh_if(const bvh_if& a_object){
  m_mesh       = a_object.m_mesh;
  m_root       = a_object.m_root;
  m_flipInside = a_object.m_flipInside;
}

template <class T, class BV>
bvh_if<T, BV>::~bvh_if(){

}

template <class T, class BV>
Real bvh_if<T, BV>::value(const RealVect& a_point) const {
  
  Vec3 p(a_point[0], a_point[1], a_point[2]);

  //  T retval = m_mesh->signedDistance(p); // dcel::mesh can return either positive or negative for outside.

  T pruneDistance = std::numeric_limits<T>::infinity();
  T retval = m_root->signedDistance(p, pruneDistance);
  
  if(m_flipInside){
    retval = -retval;
  }

  return Real(retval);
}

template <class T, class BV>
BaseIF* bvh_if<T, BV>::newImplicitFunction() const {
  return static_cast<BaseIF*> (new bvh_if(*this));
}

template <class T, class BV>
inline
decltype(auto) bvh_if<T, BV>::buildBVH(const StopFunctionT<T, face, BV>&      stopFunc,
					const PartitionFunctionT<T, face, BV>& partFunc,
					const BVConstructorT<T, face, BV>&     bvFunc) noexcept {
  
  m_root = std::make_shared<Node>(m_mesh->getFaces());
  m_root->setDepth(0);
  m_root->setNodeType(BVH::NodeType::Leaf);
  
  m_root->topDownSortAndPartitionObjects(stopFunc, partFunc, bvFunc);
}

#endif

