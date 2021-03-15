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

template <class T, class BV>
inline
NodeT<T, BV>::NodeT() {
  m_parent = nullptr;
  m_left   = nullptr;
  m_right  = nullptr;

  m_depth    = 0;
  m_nodeType = NodeType::Regular;
}

template <class T, class BV>
inline
NodeT<T, BV>::NodeT(std::shared_ptr<Node>& a_parent) : NodeT<T, BV>() {
  m_parent = a_parent;
  m_depth  = a_parent.m_depth + 1;
}

template <class T, class BV>
inline
NodeT<T, BV>::NodeT(std::vector<std::shared_ptr<T> >& a_objects) : NodeT<T, BV>() {
  m_objects  = a_objects;
  m_nodeType = NodeType::Leaf;
}

template <class T, class BV>
inline
NodeT<T, BV>::~NodeT() {
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::setParent(const std::shared_ptr<Node>& a_parent) noexcept {
  m_parent = a_parent;
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::setLeft(const std::shared_ptr<Node>& a_left) noexcept {
  m_left = a_left;
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::setRight(const std::shared_ptr<Node>& a_right) noexcept{
  m_right = a_right;
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::setNodeType(const NodeType a_nodeType) noexcept {
  m_nodeType = a_nodeType;
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::setDepth(const int a_depth) noexcept {
  m_depth = a_depth;
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::setObjects(const std::vector<std::shared_ptr<T> >& a_objects) noexcept {
  m_objects = a_objects;
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getParent() noexcept {
  return (m_parent);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getParent() const noexcept {
  return (m_parent);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getLeft() noexcept {
  return (m_left);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getLeft() const noexcept {
  return (m_left);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getRight() noexcept {
  return (m_right);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getRight() const noexcept {
  return (m_right);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getNodeType() const noexcept {
  return m_nodeType;
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getDepth() const noexcept {
  return (m_depth);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getSignedDistanceToPoint(const Vec3& a_x0) const noexcept {
  return m_bv.getSignedDistanceToPoint(a_x0);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getBoundingVolume() noexcept {
  return (m_bv);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getBoundingVolume() const noexcept {
  return (m_bv);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getObjects() noexcept {
  return (m_objects);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getObjects() const noexcept {
  return (m_objects);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::topDownSortAndPartitionObjects(SortingFunction<T>   sortingFunc,
							    CostFunction<T>      costFunc,
							    StopFunction<T>      stopFunc,
							    BVConstructor<T, BV> bvFunc) noexcept -> void {

  const bool stopSplitting = stopFunc(m_objects, m_depth);

  if(procID() == 0){
    std::cout << "my data size = " << m_objects.size() << std::endl;
  }

  if(!stopSplitting){
    int splitDir     = 0;
    int splitIndex   = 0;
    double splitCost = std::numeric_limits<double>::max();

    // Sorting lambda, captures spatial direction in which objects are sorted. 
    int dir = 0;
    auto comp = [sortingFunc, dir](const std::shared_ptr<T>& t1, const std::shared_ptr<T>& t2) {
      return sortingFunc(*t1, *t2, dir);
    };


    // Determine in which spatial direction we partition the objects. 
    for (dir = 0; dir < 3; dir++){

      // Sort objects and compute a splitting direction and cost
      std::sort(m_objects.begin(), m_objects.end(), comp); 

      const std::pair<double, int> costAndIndex = costFunc(m_objects, dir);

      if(costAndIndex.first < splitCost){
	splitDir = dir;
      }
    }

    // Direction has been determined, partition objects along given direction. 
    dir = splitDir;
    std::sort(m_objects.begin(), m_objects.end(), comp);
    const std::pair<double, int> costAndIndex = costFunc(m_objects, dir);
    splitIndex = costAndIndex.second;

    // Create the left and right nodes. 
    std::vector<std::shared_ptr<T> > lData (m_objects.begin(), m_objects.begin() + splitIndex-1);
    std::vector<std::shared_ptr<T> > rData(m_objects.begin() + splitIndex, m_objects.end());

    m_nodeType = NodeType::Regular; 
    m_objects.resize(0);

    m_left  = std::make_shared<NodeT<T, BV> >();
    m_right = std::make_shared<NodeT<T, BV> >();

    m_left->setObjects(lData);
    m_right->setObjects(rData);

    m_left->computeBoundingVolume(bvFunc);
    m_right->computeBoundingVolume(bvFunc);

    m_left->setParent(std::make_shared<NodeT<T, BV> >(*this));
    m_right->setParent(std::make_shared<NodeT<T, BV> >(*this));

    m_left->setNodeType(NodeType::Leaf);
    m_right->setNodeType(NodeType::Leaf);

    m_left->setDepth(m_depth + 1);
    m_right->setDepth(m_depth + 1);

    m_left->topDownSortAndPartitionObjects(sortingFunc, costFunc, stopFunc, bvFunc);
    m_right->topDownSortAndPartitionObjects(sortingFunc, costFunc, stopFunc, bvFunc);
  }
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::computeBoundingVolume(BVConstructor<T, BV> bvFunction) noexcept {
  m_bv = bvFunction(m_objects);
}

#endif
