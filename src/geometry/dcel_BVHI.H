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
  m_primitives.resize(0);

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
//NodeT<T, P, BV>::NodeT(const PrimitiveList& a_primitives) : NodeT<T, P, BV>() {
NodeT<T, P, BV>::NodeT(const std::vector<std::shared_ptr<P> >& a_primitives) : NodeT<T, P, BV>() {
  for (const auto& p : a_primitives){
    m_primitives.emplace_back(p);
  }
  
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

  m_bv = a_bvFunc(*this);
  
  if(!stopSplitting){

    const int splitIndex = a_partFunc(*this);

    // Create the left and right nodes and move further down the tree.
    PrimitiveList lPrimitives(m_primitives.begin(), m_primitives.begin() + splitIndex+1);
    PrimitiveList rPrimitives(m_primitives.begin() + splitIndex + 1, m_primitives.end());

    this->insertNode(m_left,  lPrimitives);
    this->insertNode(m_right, rPrimitives);

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
decltype(auto) NodeT<T, P, BV>::getDistanceToBoundingVolume(const Vec3& a_point) const noexcept{
  return m_bv.getDistance(a_point);
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::getDistanceToPrimitives(const Vec3& a_point) const noexcept {
  T minDist = std::numeric_limits<T>::max();

  for (const auto& p : m_primitives){
    const auto curDist = p->signedDistance(a_point);

    if(std::abs(curDist) < std::abs(minDist)){
      minDist = curDist;
    }
  }

  return minDist;
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::pruneBranches(const Vec3& a_point) const noexcept {

  T minDist = std::numeric_limits<T>::infinity();

  this->pruneBranches(minDist, a_point);

  return minDist;
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::pruneBranches(T& a_closest, const Vec3& a_point) const noexcept -> void {
					      
  if(m_nodeType == NodeType::Leaf){
    const T primDist = this->getDistanceToPrimitives(a_point);

    if(std::abs(primDist) < std::abs(a_closest)){
      a_closest = primDist;
    }
  }
  else {
    const T minL = m_left->getDistanceToBoundingVolume(a_point);
    const T minR = m_right->getDistanceToBoundingVolume(a_point);

    const auto leftFirst = (minL < minR);
    
    const auto first  = leftFirst ? m_left  : m_right;
    const auto second = leftFirst ? m_right : m_left;

    const auto minFirst  = std::min(minL, minR);
    const auto minSecond = std::max(minL, minR);

    if(minFirst  < std::abs(a_closest)) first ->pruneBranches(a_closest, a_point);
    if(minSecond < std::abs(a_closest)) second->pruneBranches(a_closest, a_point);
  }
}

template <class T, class P, class BV>
inline
decltype(auto) NodeT<T, P, BV>::pruneQ(const Vec3&                           a_point,
				       const PrimitivesDistanceFunction&     a_primFunc,
				       const BoundingVolumeDistanceFunction& a_bvFunc) const noexcept {

  using QueueElement = std::pair<T, std::shared_ptr<const Node> >;
  using QueueVector  = std::vector<QueueElement>;
  using QueueComp    = std::function<bool(const QueueElement&, const QueueElement&)>;
  using Queue        = std::priority_queue<QueueElement, QueueVector, QueueComp>;
      
  auto CompareElements = [](const QueueElement& a_q1, const QueueElement& a_q2){
    return a_q1.first > a_q2.first;
  };

  Queue Q(CompareElements);
  Q.push(QueueElement(a_bvFunc(m_bv, a_point), std::make_shared<const Node>(*this)));

  auto d = std::numeric_limits<T>::infinity();

  while(!Q.empty()) {
    const QueueElement cur = Q.top();
    Q.pop();

    const auto ds = cur.first;
    const auto N  = cur.second;
    
    if (std::abs(ds) < std::abs(d)){
      if(N->getNodeType() == NodeType::Leaf){
	const auto primDist = a_primFunc(N->getPrimitives(), a_point);
	
	if(std::abs(primDist) < std::abs(d)){
	  d = primDist;
	}
      }
      else{
	const auto& lN = N->getLeft();
	const auto& rN = N->getRight();
	
	const auto dl = a_bvFunc(lN->getBoundingVolume(), a_point);
	const auto dr = a_bvFunc(rN->getBoundingVolume(), a_point);
	
	if(dl < std::abs(d)) Q.emplace(QueueElement(dl, std::make_shared<const Node>(*lN)));
	if(dr < std::abs(d)) Q.emplace(QueueElement(dr, std::make_shared<const Node>(*rN)));
      }
    }
  }
  
  return d;
}


#endif
