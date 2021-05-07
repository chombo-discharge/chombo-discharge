/*!
  @file   BVHI.H
  @brief  Implementation of BVH.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _BVHI_H_
#define _BVHI_H_

#include "BVH.H"

namespace ChomboDischarge {

  namespace BVH {

    int reguCalls = 0;
    int leafCalls = 0;

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
    NodeT<T, P, BV>::NodeT(const std::vector<std::shared_ptr<P> >& a_primitives) : NodeT<T, P, BV>() {
      for (const auto& p : a_primitives){
	m_primitives.emplace_back(p);
      }
  
      m_nodeType = NodeType::Leaf;
      m_depth    = 0;
    }

    template <class T, class P, class BV>
    inline
    NodeT<T, P, BV>::NodeT(const PrimitiveList& a_primitives) : NodeT<T, P, BV>() {
      m_primitives = a_primitives;
  
      m_nodeType = NodeType::Leaf;
      m_depth    = 0;
    }

    template <class T, class P, class BV>
    inline
    NodeT<T, P, BV>::~NodeT() {
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::setParent(const NodePtr& a_parent) noexcept {
      m_parent = a_parent;
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::setLeft(const NodePtr& a_left) noexcept {
      m_left = a_left;
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::setRight(const NodePtr& a_right) noexcept{
      m_right = a_right;
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::setNodeType(const NodeType a_nodeType) noexcept {
      m_nodeType = a_nodeType;
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::setDepth(const int a_depth) noexcept {
      m_depth = a_depth;
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::setPrimitives(const PrimitiveList& a_primitives) noexcept {
      m_primitives = a_primitives;
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::setToRegularNode() noexcept {
      m_nodeType = NodeType::Regular;
      m_primitives.resize(0);
    }

    template <class T, class P, class BV>
    inline
    NodeT<T, P, BV>& NodeT<T, P, BV>::getParent() noexcept {
      return (m_parent);
    }

    template <class T, class P, class BV>
    inline
    const NodeT<T, P, BV>& NodeT<T, P, BV>::getParent() const noexcept {
      return (m_parent);
    }

    template <class T, class P, class BV>
    inline
    NodeT<T, P, BV>&  NodeT<T, P, BV>::getLeft() noexcept {
      return (m_left);
    }

    template <class T, class P, class BV>
    inline
    const NodeT<T, P, BV>& NodeT<T, P, BV>::getLeft() const noexcept {
      return (m_left);
    }

    template <class T, class P, class BV>
    inline
    NodeT<T, P, BV>& NodeT<T, P, BV>::getRight() noexcept {
      return (m_right);
    }

    template <class T, class P, class BV>
    inline
    const NodeT<T, P, BV>& NodeT<T, P, BV>::getRight() const noexcept {
      return (m_right);
    }

    template <class T, class P, class BV>
    inline
    NodeType NodeT<T, P, BV>::getNodeType() const noexcept {
      return m_nodeType;
    }

    template <class T, class P, class BV>
    inline
    int NodeT<T, P, BV>::getDepth() const noexcept {
      return m_depth;
    }

    template <class T, class P, class BV>
    inline
    BV& NodeT<T, P, BV>::getBoundingVolume() noexcept {
      return (m_bv);
    }

    template <class T, class P, class BV>
    inline
    const BV& NodeT<T, P, BV>::getBoundingVolume() const noexcept {
      return (m_bv);
    }

    template <class T, class P, class BV>
    inline
    PrimitiveListT<P>& NodeT<T, P, BV>::getPrimitives() noexcept {
      return (m_primitives);
    }

    template <class T, class P, class BV>
    inline
    const PrimitiveListT<P>& NodeT<T, P, BV>::getPrimitives() const noexcept {
      return (m_primitives);
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::topDownSortAndPartitionPrimitives(const StopFunction&      a_stopFunc,
							    const PartitionFunction& a_partFunc,
							    const BVConstructor&     a_bvFunc) noexcept {
      const bool stopSplitting = a_stopFunc(*this);

      m_bv = a_bvFunc(m_primitives);
  
      if(!stopSplitting){
    
	const auto partitions = a_partFunc(m_primitives);
    
	const auto leftPrims  = partitions.first;
	const auto rightPrims = partitions.second;

	this->insertNode(m_left,  leftPrims);
	this->insertNode(m_right, rightPrims);

	m_left ->topDownSortAndPartitionPrimitives(a_stopFunc, a_partFunc, a_bvFunc);
	m_right->topDownSortAndPartitionPrimitives(a_stopFunc, a_partFunc, a_bvFunc);

	this->setToRegularNode();
      }
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::insertNode(NodePtr& a_node, const PrimitiveList& a_primitives) noexcept {
      a_node = std::make_shared<NodeT<T, P, BV> >();

      a_node->setPrimitives(a_primitives);
      a_node->setParent(std::make_shared<NodeT<T, P, BV> >(*this));
      a_node->setNodeType(NodeType::Leaf);
      a_node->setDepth(m_depth+1);
    }

    template <class T, class P, class BV>
    inline
    T NodeT<T, P, BV>::getDistanceToBoundingVolume(const Vec3& a_point) const noexcept{
      return m_bv.getDistance(a_point);
    }

    template <class T, class P, class BV>
    inline
    T NodeT<T, P, BV>::getDistanceToBoundingVolume2(const Vec3& a_point) const noexcept{
      return m_bv.getDistance2(a_point);
    }

    template <class T, class P, class BV>
    inline
    T NodeT<T, P, BV>::getDistanceToPrimitives(const Vec3& a_point) const noexcept {
      T minDist = std::numeric_limits<T>::max();

      for (const auto& p : m_primitives){
	const auto curDist = p->signedDistance(a_point);

	if(curDist*curDist < minDist*minDist){
	  minDist = curDist;
	}
      }

      return minDist;
    }

    template <class T, class P, class BV>
    inline
    T NodeT<T, P, BV>::pruneOrdered(const Vec3& a_point) const noexcept {

      T minDist = std::numeric_limits<T>::infinity();

      this->pruneOrdered(minDist, a_point);

      return minDist;
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::pruneOrdered(T& a_closest, const Vec3& a_point) const noexcept  {
      if(m_nodeType == NodeType::Leaf){
	leafCalls++;
	const T primDist = this->getDistanceToPrimitives(a_point);

	if(primDist*primDist < a_closest*a_closest){
	  a_closest = primDist;
	}
      }
      else {
	reguCalls++;
	const T minL = m_left ->getDistanceToBoundingVolume(a_point);
	const T minR = m_right->getDistanceToBoundingVolume(a_point);

	const auto leftFirst = (minL < minR);
    
	const auto first  = leftFirst ? m_left  : m_right;
	const auto second = leftFirst ? m_right : m_left;

	const auto minFirst  = std::min(minL, minR);
	const auto minSecond = std::max(minL, minR);

	if(minFirst*minFirst   < a_closest*a_closest) first ->pruneOrdered(a_closest, a_point);
	if(minSecond*minSecond < a_closest*a_closest) second->pruneOrdered(a_closest, a_point);
      }
    }

    template <class T, class P, class BV>
    inline
    T NodeT<T, P, BV>::pruneOrdered2(const Vec3& a_point) const noexcept {

      T minDist2 = std::numeric_limits<T>::infinity();

      std::shared_ptr<const P> closest = nullptr;

      this->pruneOrdered2(minDist2, closest, a_point);

      const T minDist = closest->signedDistance(a_point);

      return minDist;
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::pruneOrdered2(T& a_minDist2, std::shared_ptr<const P>& a_closest, const Vec3& a_point) const noexcept  {
      if(m_nodeType == NodeType::Leaf){
	leafCalls++;
	for (const auto& p : m_primitives){
	  const auto curDist2 = p->unsignedDistance2(a_point);

	  if(curDist2 < a_minDist2){
	    a_minDist2 = curDist2;
	    a_closest  = p;
	  }
	}
      }
      else{
	reguCalls++;
	const auto minL2 = m_left ->getDistanceToBoundingVolume2(a_point);
	const auto minR2 = m_right->getDistanceToBoundingVolume2(a_point);

	const auto leftFirst = (minL2 < minR2);
    
	const auto first  = leftFirst ? m_left  : m_right;
	const auto second = leftFirst ? m_right : m_left;

	const auto minFirst2  = std::min(minL2, minR2);
	const auto minSecond2 = std::max(minL2, minR2);

	if(minFirst2  < a_minDist2) first ->pruneOrdered2(a_minDist2, a_closest, a_point);
	if(minSecond2 < a_minDist2) second->pruneOrdered2(a_minDist2, a_closest, a_point);
      }
    }

    template <class T, class P, class BV>
    inline
    T NodeT<T, P, BV>::pruneUnordered(const Vec3& a_point) const noexcept {

      T minDist = std::numeric_limits<T>::infinity();

      this->pruneUnordered(minDist, a_point);

      return minDist;
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::pruneUnordered(T& a_closest, const Vec3& a_point) const noexcept  {
					      
      if(m_nodeType == NodeType::Leaf){
	const T primDist = this->getDistanceToPrimitives(a_point);

	if(primDist*primDist < a_closest*a_closest){
	  a_closest = primDist;
	}
      }
      else {
	const T minL = m_left ->getDistanceToBoundingVolume(a_point);
	const T minR = m_right->getDistanceToBoundingVolume(a_point);

	if(minL*minL < a_closest*a_closest) m_left ->pruneUnordered(a_closest, a_point);
	if(minR*minR < a_closest*a_closest) m_right->pruneUnordered(a_closest, a_point);
      }
    }

    template <class T, class P, class BV>
    inline
    T NodeT<T, P, BV>::pruneUnordered2(const Vec3& a_point) const noexcept {

      T minDist2 = std::numeric_limits<T>::infinity();

      std::shared_ptr<const P> closest = nullptr;

      this->pruneUnordered2(minDist2, closest, a_point);

      const T minDist = closest->signedDistance(a_point);

      return minDist;
    }

    template <class T, class P, class BV>
    inline
    void NodeT<T, P, BV>::pruneUnordered2(T& a_minDist2, std::shared_ptr<const P>& a_closest, const Vec3& a_point) const noexcept  {

      if(m_nodeType == NodeType::Leaf){
	for (const auto& p : m_primitives){
	  const auto curDist2 = p->unsignedDistance2(a_point);

	  if(curDist2 < a_minDist2){
	    a_minDist2 = curDist2;
	    a_closest  = p;
	  }
	}
      }
      else{
	const auto minL2 = m_left ->getDistanceToBoundingVolume2(a_point);
	const auto minR2 = m_right->getDistanceToBoundingVolume2(a_point);

	if(minL2 < a_minDist2) m_left ->pruneOrdered2(a_minDist2, a_closest, a_point);
	if(minR2 < a_minDist2) m_right->pruneOrdered2(a_minDist2, a_closest, a_point);
      }
    }

    template <class T, class P, class BV>
    inline
    T NodeT<T, P, BV>::prunePriorityQueue(const Vec3& a_point) const noexcept {
  
      using QueueElement = std::pair<T, std::shared_ptr<const Node> >;
      
      auto CompareElements = [](const QueueElement& a_q1, const QueueElement& a_q2){
	return a_q1.first > a_q2.first;
      };

      std::vector<QueueElement> Q(0);

      const auto ds = this->getDistanceToBoundingVolume(a_point);
      const auto N  = std::make_shared<const Node>(*this);
      Q.emplace_back(QueueElement(ds, N));

      auto d = std::numeric_limits<T>::infinity();

      while(!Q.empty()) {

	std::sort(Q.begin(), Q.end(), CompareElements);
	const QueueElement cur = Q.back();
	Q.pop_back();

	const auto ds = cur.first;
	const auto N  = cur.second;
    
	if (ds*ds < d*d){
	  if(N->getNodeType() == NodeType::Leaf){
	    const T primDist = N->getDistanceToPrimitives(a_point);
	
	    if(std::abs(primDist) < std::abs(d)){
	      d = primDist;
	    }
	  }
	  else{
	    const auto& lN = N->getLeft();
	    const auto& rN = N->getRight();

	    const T dl = lN->getDistanceToBoundingVolume(a_point);
	    const T dr = rN->getDistanceToBoundingVolume(a_point);
	
	    if(dl*dl < d*d) Q.emplace_back(QueueElement(dl, std::make_shared<const Node>(*lN)));
	    if(dr*dr < d*d) Q.emplace_back(QueueElement(dr, std::make_shared<const Node>(*rN)));
	  }
	}
      }
  
      return d;
    }

    template <class T, class P, class BV>
    inline
    T NodeT<T, P, BV>::prunePriorityQueue2(const Vec3& a_point) const noexcept {
  
      using QueueElement = std::pair<T, std::shared_ptr<const Node> >;
      
      auto CompareElements = [](const QueueElement& a_q1, const QueueElement& a_q2){
	return a_q1.first > a_q2.first;
      };
  
      std::shared_ptr<const P> closestPrimitive = nullptr;
      auto minDist2 = std::numeric_limits<T>::infinity();

      std::vector<QueueElement> Q(0);

      // Push node onto queue.
      const auto d = this->getDistanceToBoundingVolume2(a_point);
      const auto N = std::make_shared<const Node>(*this);
      Q.emplace_back(QueueElement(d, N));


      while(!Q.empty()) {


	std::sort(Q.begin(), Q.end(), CompareElements);
	const QueueElement cur = Q.back();
	Q.pop_back();

	const auto dist2 = cur.first;
    
	if (dist2 < minDist2){
	  const auto curNode = cur.second;
      
	  if(curNode->getNodeType() == NodeType::Leaf){
	    for (const auto& p : curNode->getPrimitives()){
	      const auto curDist2 = p->unsignedDistance2(a_point);

	      if(curDist2 < minDist2){
		minDist2         = curDist2;
		closestPrimitive = p;
	      }
	    }
	  }
	  else{
	    const auto& leftNode  = curNode->getLeft();
	    const auto& rightNode = curNode->getRight();

	    const T distL2 = leftNode ->getDistanceToBoundingVolume2(a_point);
	    const T distR2 = rightNode->getDistanceToBoundingVolume2(a_point);
	
	    if(distL2 < minDist2) Q.emplace_back(QueueElement(distL2, std::make_shared<const Node>(*leftNode)));
	    if(distR2 < minDist2) Q.emplace_back(QueueElement(distR2, std::make_shared<const Node>(*rightNode)));
	  }
	}
      }

      return closestPrimitive->signedDistance(a_point);
    }
  }
}

#endif
