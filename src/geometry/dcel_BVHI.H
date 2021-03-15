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

template <class Obj, class BV>
inline
NodeT<Obj, BV>::NodeT() {
  m_parent = nullptr;
  m_left   = nullptr;
  m_right  = nullptr;

  m_depth    = 0;
  m_nodeType = NodeType::Regular;
}

template <class Obj, class BV>
inline
NodeT<Obj, BV>::NodeT(NodePtr& a_parent) : NodeT<Obj, BV>() {
  m_parent = a_parent;
  m_depth  = a_parent.m_depth + 1;
}

template <class Obj, class BV>
inline
NodeT<Obj, BV>::NodeT(const ObjectList& a_objects) : NodeT<Obj, BV>() {
  m_objects  = a_objects;
  m_nodeType = NodeType::Leaf;
}

template <class Obj, class BV>
inline
NodeT<Obj, BV>::~NodeT() {
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::setParent(const NodePtr& a_parent) noexcept {
  m_parent = a_parent;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::setLeft(const NodePtr& a_left) noexcept {
  m_left = a_left;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::setRight(const NodePtr& a_right) noexcept{
  m_right = a_right;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::setNodeType(const NodeType a_nodeType) noexcept {
  m_nodeType = a_nodeType;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::setDepth(const int a_depth) noexcept {
  m_depth = a_depth;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::setObjects(const ObjectList& a_objects) noexcept {
  m_objects = a_objects;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::setToRegularNode() noexcept {
  m_nodeType = NodeType::Regular;
  m_objects.resize(0);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::computeBoundingVolume(const BVConstructor& bvFunction) noexcept {
  m_bv = bvFunction(m_objects);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getParent() noexcept {
  return (m_parent);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getParent() const noexcept {
  return (m_parent);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getLeft() noexcept {
  return (m_left);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getLeft() const noexcept {
  return (m_left);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getRight() noexcept {
  return (m_right);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getRight() const noexcept {
  return (m_right);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getNodeType() const noexcept {
  return m_nodeType;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getDepth() const noexcept {
  return (m_depth);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getSignedDistanceToPoint(const Vec3& a_x0) const noexcept {
  return m_bv.getSignedDistanceToPoint(a_x0);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getBoundingVolume() noexcept {
  return (m_bv);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getBoundingVolume() const noexcept {
  return (m_bv);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getObjects() noexcept {
  return (m_objects);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getObjects() const noexcept {
  return (m_objects);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::topDownSortAndPartitionObjects(const StopFunction&      a_stopFunc,
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

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::insertNode(NodePtr& a_node, const ObjectList& a_objects) noexcept {
  a_node = std::make_shared<NodeT<Obj, BV> >();

  a_node->setObjects(a_objects);
  a_node->setParent(std::make_shared<NodeT<Obj, BV> >(*this));
  a_node->setNodeType(NodeType::Leaf);
  a_node->setDepth(m_depth+1);
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getLeafData() const noexcept -> ObjectList {

  auto objects = m_objects;

  if(m_nodeType == NodeType::Regular){
    const auto leftObjects  = m_left->getLeafData();
    const auto rightObjects = m_right->getLeafData();

    objects.insert(objects.end(), leftObjects.begin(),  leftObjects.end());
    objects.insert(objects.end(), rightObjects.begin(), rightObjects.end());
  }

  return objects;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::signedDistance(const Vec3& a_point, double& a_pruneDistance) const noexcept {
  ObjectList candidates = this->getMinSignedDistanceToPointCandidates(a_point, a_pruneDistance);

  //  pout() << candidates.size() << endl;

  double minDist = std::numeric_limits<double>::infinity();

  for (const auto& c : candidates){
    const double curDist = (double) c->signedDistance(a_point);

    if(std::abs(curDist) < std::abs(minDist)){
      minDist = curDist;
    }
  }
    
  return minDist;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::getMinSignedDistanceToPointCandidates(const Vec3& a_point, double& a_pruneDistance) const noexcept -> ObjectList {

  ObjectList cand = ObjectList();

  if(m_nodeType == NodeType::Leaf){
    bool returnObjects = false;
    
    for (const auto& O : m_objects){
      const double cur = O->signedDistance(a_point);

      if(std::abs(cur) < a_pruneDistance){
	a_pruneDistance = std::abs(cur);
	returnObjects   = true;
      }
    }

    if(returnObjects){
      cand = m_objects;
    }
  }
  else {
    const BV& bvL = m_left ->getBoundingVolume();
    const BV& bvR = m_right->getBoundingVolume();
    
    const double minL = std::abs(bvL.getDistance(a_point));
    const double minR = std::abs(bvR.getDistance(a_point));

    const bool pruneLeft = (minL > a_pruneDistance);
    if(!pruneLeft){
      const ObjectList lcand = m_left->getMinSignedDistanceToPointCandidates(a_point, a_pruneDistance);
      cand.insert(cand.end(), lcand.begin(), lcand.end());
    }

    const bool pruneRight = (minR > a_pruneDistance);
    if (!pruneRight){
      const ObjectList rcand = m_right->getMinSignedDistanceToPointCandidates(a_point, a_pruneDistance);
      cand.insert(cand.end(), rcand.begin(), rcand.end());
    }
  }
  
  return cand;
}

#endif
