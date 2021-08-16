/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BVHImplem.H
  @brief  Implementation of CD_BVH.H
  @author Robert Marskar
*/

#ifndef CD_BVHImplem_H_
#define CD_BVHImplem_H_

#include <CD_BVH.H>
#include <CD_NamespaceHeader.H>

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
  NodeT<T, P, BV>::NodeT(const NodePtr& a_parent) : NodeT<T, P, BV>() {
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
  T NodeT<T, P, BV>::pruneTree(const Vec3& a_point, const Prune a_pruning) const noexcept {
    T ret;
    
    switch(a_pruning){
    case Prune::Ordered:
      ret = this->pruneOrdered(a_point);
      break;
    case Prune::Ordered2:
      ret = this->pruneOrdered2(a_point);
      break;
    case Prune::Unordered:
      ret = this->pruneUnordered(a_point);
      break;
    case Prune::Unordered2:
      ret = this->pruneUnordered2(a_point);
      break;
    default:
      std::cerr << "In file CD_BVHImplem function NodeT<T, P, BV>::pruneTree(Vec3, Prune) -- bad input enum for 'Prune'\n";
    };

    return ret;
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

    // TLDR: Beginning at some node, this routine descends the left/right branches in the tree. It always descends the branch with the shortest distance
    //       to the bounding volume first. The other branch is investigated only after the full sub-tree beneath the first branch has completed. Since the shortest
    //       distance to primitives is updated underway, there is a decent chance that the secondary subtree can be pruned. Hence why this routine is more efficient
    //       than prunedUnordered.
    if(m_nodeType == NodeType::Leaf){
      leafCalls++;
      const T primDist = this->getDistanceToPrimitives(a_point);

      if(primDist*primDist < a_closest*a_closest){
	a_closest = primDist;
      }
    }
    else {
      reguCalls++;
      const T minL = m_left ->getDistanceToBoundingVolume(a_point); // Distance from a_point to the left sub-trees bounding volume. 
      const T minR = m_right->getDistanceToBoundingVolume(a_point); // Distance from a_point to the right sub-trees bounding volume. 

      // Check the sub-tree with the shortest distance between a_point and the sub-trees bounding volume first. 
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

    // TLDR: Beginning at some node, this routine descends the left/right branches in the tree. It always descends the branch with the shortest distance
    //       to the bounding volume first. The other branch is investigated only after the full sub-tree beneath the first branch has completed. Since the shortest
    //       distance to primitives is updated underway, there is a decent chance that the secondary subtree can be pruned. Hence why this routine is more efficient
    //       than prunedUnordered. The only difference between this routine is that we prune based on the squared unsigned distance. This
    //       permits us to skip some computations involving a couple of square roots and only compute the signed distance at the end. 
    
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

      // Check the sub-tree with the shortest unsigned squared distance between a_point and the sub-trees bounding volume first.       
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
}

#include <CD_NamespaceFooter.H>

#endif
