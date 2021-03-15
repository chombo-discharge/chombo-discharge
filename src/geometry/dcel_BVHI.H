/*!
  @file   dcel_BVHI.H
  @brief  Implementation of dcel_BVH.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_BVHI_H_
#define _DCEL_BVHI_H_

#include "dcel_BVH.H"

using namespace dcel::BVH;

template <class T, class Obj, class BV>
inline
NodeT<T, Obj, BV>::NodeT() {
  m_parent = nullptr;
  m_left   = nullptr;
  m_right  = nullptr;

  m_depth    = 0;
  m_nodeType = NodeType::Regular;
}

template <class T, class Obj, class BV>
inline
NodeT<T, Obj, BV>::NodeT(NodePtr& a_parent) : NodeT<T, Obj, BV>() {
  m_parent   = a_parent;
  m_depth    = a_parent.m_depth + 1;
  m_nodeType = NodeType::Leaf;
}

template <class T, class Obj, class BV>
inline
NodeT<T, Obj, BV>::NodeT(const ObjectList& a_objects) : NodeT<T, Obj, BV>() {
  m_objects  = a_objects;
  m_nodeType = NodeType::Leaf;
}

template <class T, class Obj, class BV>
inline
NodeT<T, Obj, BV>::~NodeT() {
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::setParent(const NodePtr& a_parent) noexcept {
  m_parent = a_parent;
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::setLeft(const NodePtr& a_left) noexcept {
  m_left = a_left;
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::setRight(const NodePtr& a_right) noexcept{
  m_right = a_right;
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::setNodeType(const NodeType a_nodeType) noexcept {
  m_nodeType = a_nodeType;
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::setDepth(const int a_depth) noexcept {
  m_depth = a_depth;
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::setObjects(const ObjectList& a_objects) noexcept {
  m_objects = a_objects;
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::setToRegularNode() noexcept {
  m_nodeType = NodeType::Regular;
  m_objects.resize(0);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::computeBoundingVolume(const BVConstructor& bvFunction) noexcept {
  m_bv = bvFunction(m_objects);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getParent() noexcept {
  return (m_parent);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getParent() const noexcept {
  return (m_parent);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getLeft() noexcept {
  return (m_left);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getLeft() const noexcept {
  return (m_left);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getRight() noexcept {
  return (m_right);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getRight() const noexcept {
  return (m_right);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getNodeType() const noexcept {
  return m_nodeType;
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getDepth() const noexcept {
  return (m_depth);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getDistanceToBV(const Vec3& a_x0) const noexcept {
  return m_bv.getDistance(a_x0);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getDistanceToData(const Vec3& a_x0) const noexcept {
  T ret = std::numeric_limits<T>::infinity();

  for (const auto& o : m_objects){
    const auto cur = std::abs(o->unsignedDistance2(a_x0));

    ret = std::min(cur, ret);
  }

  return sqrt(ret);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getBoundingVolume() noexcept {
  return (m_bv);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getBoundingVolume() const noexcept {
  return (m_bv);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getObjects() noexcept {
  return (m_objects);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getObjects() const noexcept {
  return (m_objects);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::topDownSortAndPartitionObjects(const StopFunction&      a_stopFunc,
							      const PartitionFunction& a_partFunc,
							      const BVConstructor&     a_bvFunc) noexcept -> void {

  //  const bool stopSplitting = a_stopFunc(m_objects, m_depth);
  //const bool splitThisNode = !(a_stopFunc(*this));
  const bool stopSplitting = a_stopFunc(*this);

  this->computeBoundingVolume(a_bvFunc);
  
  if(!stopSplitting){

    const int splitIndex = a_partFunc(*this);

    // Create the left and right nodes and move further down the tree.
    ObjectList lData(m_objects.begin(), m_objects.begin() + splitIndex+1);
    ObjectList rData(m_objects.begin() + splitIndex + 1, m_objects.end());

    this->insertNode(m_left,  lData);
    this->insertNode(m_right, rData);

    m_left ->topDownSortAndPartitionObjects(a_stopFunc, a_partFunc, a_bvFunc);
    m_right->topDownSortAndPartitionObjects(a_stopFunc, a_partFunc, a_bvFunc);

    this->setToRegularNode();
  }
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::insertNode(NodePtr& a_node, const ObjectList& a_objects) noexcept {
  a_node = std::make_shared<NodeT<T, Obj, BV> >();

  a_node->setObjects(a_objects);
  a_node->setParent(std::make_shared<NodeT<T, Obj, BV> >(*this));
  a_node->setNodeType(NodeType::Leaf);
  a_node->setDepth(m_depth+1);
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::getLeafData() const noexcept -> ObjectList {

  ObjectList objects = ObjectList();

  if(m_nodeType == NodeType::Leaf){
    objects = m_objects;
  }
  else{
    const auto leftObjects  = m_left->getLeafData();
    const auto rightObjects = m_right->getLeafData();

    objects.insert(objects.end(), leftObjects.begin(),  leftObjects.end());
    objects.insert(objects.end(), rightObjects.begin(), rightObjects.end());
  }

  return objects;
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::signedDistance(const Vec3& a_point, T& a_pruneDistance) const noexcept {
  ObjectList candidates = this->pruneTree(a_point, a_pruneDistance);

  T minDist = std::numeric_limits<T>::infinity();

  for (const auto& c : candidates){
    const T curDist = (T) c->signedDistance(a_point);

    if(std::abs(curDist) < std::abs(minDist)){
      minDist = curDist;
    }
  }
    
  return minDist;
}

template <class T, class Obj, class BV>
inline
decltype(auto) NodeT<T, Obj, BV>::pruneTree(const Vec3& a_point, T& a_pruneDistance) const noexcept -> ObjectList {

  ObjectList cand = ObjectList();

  if(m_nodeType == NodeType::Leaf){
    const T distanceToData = this->getDistanceToData(a_point);

    if(distanceToData < a_pruneDistance){
      cand = m_objects;
      a_pruneDistance = distanceToData;
    }
  }
  else {
    const T minL = m_left ->getDistanceToBV(a_point);
    const T minR = m_right->getDistanceToBV(a_point);

    const bool leftFirst = minL < minR;
    
    std::shared_ptr<NodeT<T, Obj, BV> > first  = leftFirst ? m_left  : m_right;
    std::shared_ptr<NodeT<T, Obj, BV> > second = leftFirst ? m_right : m_left;

    const T minFirst  = std::min(minL, minR);
    const T minSecond = std::max(minL, minR);

    if(minFirst < a_pruneDistance){
      const ObjectList branchCand = first->pruneTree(a_point, a_pruneDistance);
      cand.insert(cand.end(), branchCand.begin(), branchCand.end());
    }

    if (minSecond < a_pruneDistance){
      const ObjectList branchCand = second->pruneTree(a_point, a_pruneDistance);
      cand.insert(cand.end(), branchCand.begin(), branchCand.end());
    }
  }
  
  return cand;
}

#endif
