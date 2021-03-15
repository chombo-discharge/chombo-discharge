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
using namespace dcel::BVH;

template <class T, class BV>
StopFunctionT<faceT<T>, BV> dcel_if<T, BV>::defaultStopFunction = [](const NodeT<faceT<T>, BV>& a_node){
  const auto& objects = a_node.getObjects();
  const int depth     = a_node.getDepth();

  return objects.size() < 2;
};

template <class T, class BV>
BVConstructorT<faceT<T>, BV> dcel_if<T, BV>::defaultBVConstructor = [](const NodeT<faceT<T>, BV>& a_node){
  const auto& a_objects = a_node.getObjects();
  
  std::vector<Vec3T<T> > coordinates;

  for (const auto& f : a_objects){
    const auto faceCoordinates = f->getAllVertexCoordinates();
    coordinates.insert(coordinates.end(), faceCoordinates.begin(), faceCoordinates.end());
  }

  return BV(coordinates);
};

template <class T, class BV>
PartitionFunctionT<faceT<T>, BV> dcel_if<T, BV>::defaultPartitionFunction = [](NodeT<faceT<T>, BV>& a_node){

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
dcel_if<T, BV>::dcel_if(const std::shared_ptr<mesh>& a_mesh, const bool a_flipInside){
  m_mesh       = a_mesh;
  m_flipInside = a_flipInside;
  m_root       = nullptr;
}

template <class T, class BV>
dcel_if<T, BV>::dcel_if(const dcel_if& a_object){
  m_mesh       = a_object.m_mesh;
  m_root       = a_object.m_root;
  m_flipInside = a_object.m_flipInside;
}

template <class T, class BV>
dcel_if<T, BV>::~dcel_if(){

}

template <class T, class BV>
Real dcel_if<T, BV>::value(const RealVect& a_point) const {
  
  Vec3 p(a_point[0], a_point[1], a_point[2]);
  Vec3T<double> p2(a_point[0], a_point[1], a_point[2]);

  //  T retval = m_mesh->signedDistance(p); // dcel::mesh can return either positive or negative for outside.

  double distSoFar = std::numeric_limits<double>::infinity();
  double retval = m_root->signedDistance(p2, distSoFar);

  //  const auto candidates = m_root->getMinSignedDistanceToPointCandidates(p2);
  
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
decltype(auto) dcel_if<T, BV>::buildBVH(const StopFunctionT<face, BV>&      stopFunc,
					const PartitionFunctionT<face, BV>& partFunc,
					const BVConstructorT<face, BV>&     bvFunc) noexcept {
  
  m_root = std::make_shared<Node>(m_mesh->getFaces());
  m_root->setDepth(0);
  m_root->setNodeType(BVH::NodeType::Leaf);
  
  m_root->topDownSortAndPartitionObjects(stopFunc, partFunc, bvFunc);

}

#endif

