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

template <class T, class P, class BV>
inline
NodeT<T, P, BV>::NodeT() {
  m_parent = nullptr;
  m_left   = nullptr;
  m_right  = nullptr;

  m_depth    = 0;
  m_nodeType = NodeType::Regular;
}

template <class T, class P, class BV>
inline
NodeT<T, P, BV>::NodeT(NodePtr& a_parent) : NodeT<T, P, BV>() {
  m_parent   = a_parent;
  m_depth    = a_parent.m_depth + 1;
  m_nodeType = NodeType::Leaf;
}

template <class T, class P, class BV>
inline
NodeT<T, P, BV>::NodeT(const PrimitiveList& a_primitives) : NodeT<T, P, BV>() {
  m_primitives  = a_primitives;
  m_nodeType = NodeType::Leaf;
  m_depth    = 0;
}

template <class T, class P, class BV>
inline
NodeT<T, P, BV>::~NodeT() {
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::setParent(const NodePtr& a_parent) noexcept {
  m_parent = a_parent;
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::setLeft(const NodePtr& a_left) noexcept {
  m_left = a_left;
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::setRight(const NodePtr& a_right) noexcept{
  m_right = a_right;
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::setNodeType(const NodeType a_nodeType) noexcept {
  m_nodeType = a_nodeType;
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::setDepth(const int a_depth) noexcept {
  m_depth = a_depth;
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::setPrimitives(const PrimitiveList& a_primitives) noexcept {
  m_primitives = a_primitives;
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::setToRegularNode() noexcept {
  m_nodeType = NodeType::Regular;
  m_primitives.resize(0);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::computeBoundingVolume(const BVConstructor& bvFunction) noexcept {
  m_bv = bvFunction(m_primitives);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getParent() noexcept {
  return (m_parent);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getParent() const noexcept {
  return (m_parent);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getLeft() noexcept {
  return (m_left);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getLeft() const noexcept {
  return (m_left);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getRight() noexcept {
  return (m_right);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getRight() const noexcept {
  return (m_right);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getNodeType() const noexcept {
  return m_nodeType;
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getDepth() const noexcept {
  return (m_depth);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getDistanceToBV(const Vec3& a_x0) const noexcept {
  return m_bv.getDistance(a_x0);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getDistanceToData(const Vec3& a_x0) const noexcept {
  T ret = std::numeric_limits<T>::infinity();

  for (const auto& o : m_primitives){
    const auto cur = o->signedDistance(a_x0);

    if(std::abs(cur) < std::abs(ret)){
      ret = cur;
    }
    // const auto cur = std::abs(o->unsignedDistance2(a_x0));

    // ret = std::min(cur, ret);
  }

  return ret;
  //  return sqrt(ret);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getBoundingVolume() noexcept {
  return (m_bv);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getBoundingVolume() const noexcept {
  return (m_bv);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getPrimitives() noexcept {
  return (m_primitives);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getPrimitives() const noexcept {
  return (m_primitives);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::topDownSortAndPartitionPrimitives(const StopFunction&      a_stopFunc,
								 const PartitionFunction& a_partFunc,
								 const BVConstructor&     a_bvFunc) noexcept -> void {
  const bool stopSplitting = a_stopFunc(*this);

  this->computeBoundingVolume(a_bvFunc);
  
  if(!stopSplitting){

    const int splitIndex = a_partFunc(*this);

    // Create the left and right nodes and move further down the tree.
    PrimitiveList lData(m_primitives.begin(), m_primitives.begin() + splitIndex+1);
    PrimitiveList rData(m_primitives.begin() + splitIndex + 1, m_primitives.end());

    this->insertNode(m_left,  lData);
    this->insertNode(m_right, rData);

    m_left ->topDownSortAndPartitionPrimitives(a_stopFunc, a_partFunc, a_bvFunc);
    m_right->topDownSortAndPartitionPrimitives(a_stopFunc, a_partFunc, a_bvFunc);

    this->setToRegularNode();
  }
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::insertNode(NodePtr& a_node, const PrimitiveList& a_primitives) noexcept {
  a_node = std::make_shared<NodeT<T, P, BV> >();

  a_node->setPrimitives(a_primitives);
  a_node->setParent(std::make_shared<NodeT<T, P, BV> >(*this));
  a_node->setNodeType(NodeType::Leaf);
  a_node->setDepth(m_depth+1);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getLeafData() const noexcept -> PrimitiveList {

  PrimitiveList primitives = PrimitiveList();

  if(m_nodeType == NodeType::Leaf){
    primitives = m_primitives;
  }
  else{
    const auto leftPrimitives  = m_left->getLeafData();
    const auto rightPrimitives = m_right->getLeafData();

    primitives.insert(primitives.end(), leftPrimitives.begin(),  leftPrimitives.end());
    primitives.insert(primitives.end(), rightPrimitives.begin(), rightPrimitives.end());
  }

  return primitives;
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::pruneTree(const Vec3& a_point, T& a_pruneDistance) const noexcept -> PrimitiveList {

  PrimitiveList cand = PrimitiveList();

  if(m_nodeType == NodeType::Leaf){
    const T distanceToData = this->getDistanceToData(a_point);

    if(std::abs(distanceToData) < std::abs(a_pruneDistance)){
      cand = m_primitives;
      a_pruneDistance = distanceToData;
    }
  }
  else {
    const T minL = m_left ->getDistanceToBV(a_point);
    const T minR = m_right->getDistanceToBV(a_point);

    const bool leftFirst = minL < minR;
    
    std::shared_ptr<NodeT<T, P, BV> > first  = leftFirst ? m_left  : m_right;
    std::shared_ptr<NodeT<T, P, BV> > second = leftFirst ? m_right : m_left;

    const T minFirst  = std::min(minL, minR);
    const T minSecond = std::max(minL, minR);

    if(minFirst < std::abs(a_pruneDistance)){
      const PrimitiveList branchCand = first->pruneTree(a_point, a_pruneDistance);
      cand.insert(cand.end(), branchCand.begin(), branchCand.end());
    }

    if (minSecond < std::abs(a_pruneDistance)){
      const PrimitiveList branchCand = second->pruneTree(a_point, a_pruneDistance);
      cand.insert(cand.end(), branchCand.begin(), branchCand.end());
    }
  }
  
  return cand;
}

#endif
