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
SortingFunction<faceT<T> > dcel_if<T, BV>::defaultSortingFunction = [](const faceT<T>& f1, const faceT<T>& f2, const int dir){
  return f1.getCentroid()[dir] < f1.getCentroid()[dir];
};

template <class T, class BV>
CostFunction<faceT<T> > dcel_if<T, BV>::defaultCostFunction = [](const std::vector<std::shared_ptr<faceT<T> > >& a_objects, const int a_dir){
  //  std::cout << "calling costfunctor" << std::endl;
  return std::make_pair(0,a_dir);
};

template <class T, class BV>
StopFunction<faceT<T>, BV> dcel_if<T, BV>::defaultStopFunction = [](const NodeT<faceT<T>, BV>& a_node){
  const auto& objects = a_node.getObjects();
  const auto& depth   = a_node.getDepth();
  
  if(depth > 0) std::cout << "stopping because depth > 0" << std::endl;
  return depth>0;
  //  return (a_objects.size() < 5) || a_depth > 0;
};

template <class T, class BV>
BVConstructor<faceT<T>, BV> dcel_if<T, BV>::defaultBVConstructor = [](const std::vector<std::shared_ptr<faceT<T> > >& a_objects){
  std::vector<Vec3T<T>> coordinates;

  for (const auto& f : a_objects){
    const auto faceCoordinates = f->getAllVertexCoordinates();
    coordinates.insert(coordinates.end(), faceCoordinates.begin(), faceCoordinates.end());
  }

  return BV(coordinates);
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
decltype(auto) dcel_if<T, BV>::buildBVH(SortingFunction<face>   sortFunc,
					CostFunction<face>      costFunc,
					StopFunction<face, BV>  stopFunc,
					BVConstructor<face, BV> bvFunc) noexcept {
  
  m_root = std::make_shared<Node>(m_mesh->getFaces());
  m_root->setDepth(0);
  m_root->setNodeType(BVH::NodeType::Leaf);
  m_root->computeBoundingVolume(bvFunc);
  
  m_root->topDownSortAndPartitionObjects(sortFunc, costFunc, stopFunc, bvFunc);
}

#endif

