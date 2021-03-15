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
NodeT<Obj, BV>::NodeT(std::shared_ptr<Node>& a_parent) : NodeT<Obj, BV>() {
  m_parent = a_parent;
  m_depth  = a_parent.m_depth + 1;
}

template <class Obj, class BV>
inline
NodeT<Obj, BV>::NodeT(std::vector<std::shared_ptr<Obj> >& a_objects) : NodeT<Obj, BV>() {
  m_objects  = a_objects;
  m_nodeType = NodeType::Leaf;
}

template <class Obj, class BV>
inline
NodeT<Obj, BV>::~NodeT() {
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::setParent(const std::shared_ptr<Node>& a_parent) noexcept {
  m_parent = a_parent;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::setLeft(const std::shared_ptr<Node>& a_left) noexcept {
  m_left = a_left;
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::setRight(const std::shared_ptr<Node>& a_right) noexcept{
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
decltype(auto) NodeT<Obj, BV>::setObjects(const std::vector<std::shared_ptr<Obj> >& a_objects) noexcept {
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
decltype(auto) NodeT<Obj, BV>::computeBoundingVolume(const BVConstructor<Obj, BV> bvFunction) noexcept {
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
decltype(auto) NodeT<Obj, BV>::topDownSortAndPartitionObjects(SortingFunction<Obj>   sortingFunc,
							      CostFunction<Obj>      costFunc,
							      StopFunction<Obj>      stopFunc,
							      BVConstructor<Obj, BV> bvFunc) noexcept -> void {

  const bool stopSplitting = stopFunc(m_objects, m_depth);

  if(procID() == 0){
    std::cout << "my data size = " << m_objects.size() << std::endl;
  }

  this->computeBoundingVolume(bvFunc);

  if(!stopSplitting){
    int splitDir     = 0;
    int splitIndex   = 0;
    double splitCost = std::numeric_limits<double>::max();

    // Sorting lambda, captures spatial direction in which objects are sorted. 
    int dir = 0;
    auto comp = [sortingFunc, dir](const std::shared_ptr<Obj>& t1, const std::shared_ptr<Obj>& t2) {
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


    // Create the left and right nodes and move further down the tree. 
    std::vector<std::shared_ptr<Obj> > lData(m_objects.begin(), m_objects.begin() + splitIndex);
    std::vector<std::shared_ptr<Obj> > rData(m_objects.begin() + splitIndex+1, m_objects.end());

    this->insertNode(m_left,  lData);
    this->insertNode(m_right, rData);

    m_left->topDownSortAndPartitionObjects(sortingFunc, costFunc, stopFunc, bvFunc);
    m_right->topDownSortAndPartitionObjects(sortingFunc, costFunc, stopFunc, bvFunc);

    this->setToRegularNode();
  }
}

template <class Obj, class BV>
inline
decltype(auto) NodeT<Obj, BV>::insertNode(std::shared_ptr<Node>& a_node, const std::vector<std::shared_ptr<Obj> >& a_objects) noexcept {
  a_node = std::make_shared<NodeT<Obj, BV> >();

  a_node->setObjects(a_objects);
  a_node->setParent(std::make_shared<NodeT<Obj, BV> >(*this));
  a_node->setNodeType(NodeType::Leaf);
  a_node->setDepth(m_depth+1);
}





#endif
