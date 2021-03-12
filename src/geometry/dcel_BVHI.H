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
}

template <class T, class BV>
inline
NodeT<T, BV>::NodeT(std::shared_ptr<Node>& a_parent){
  m_parent = a_parent;
  m_left   = nullptr;
  m_right  = nullptr;
}

template <class T, class BV>
inline
NodeT<T, BV>::NodeT(std::vector<std::shared_ptr<T> >& a_objects) {
  m_objects = a_objects;
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
decltype(auto) NodeT<T, BV>::isRoot() const noexcept {
  return m_parent == nullptr;
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::isLeaf() const noexcept {
  return m_left == nullptr && m_right == nullptr;
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::isPointInside(const Vec3& a_x0) const noexcept {
  return m_bv.isPointInside(a_x0);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::isPointOutside(const Vec3& a_x0) const noexcept {
  return m_bv.isPointOutside(a_x0);
}

template <class T, class BV>
inline
decltype(auto) NodeT<T, BV>::getDistanceToPoint(const Vec3& a_x0) const noexcept {
  return m_bv.getDistanceToPoint(a_x0);
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
decltype(auto) NodeT<T, BV>::sortAndPartitionObjects(SortingFunctor<T> sortingFunc, CostFunctor<T> costFunc) noexcept {

  int splitDir   = 0;
  int splitIndex = 0;
  
  for (int dir = 0; dir < 3; dir++){
    auto comp = [sortingFunc, dir](const std::shared_ptr<T>& t1, const std::shared_ptr<T>& t2) {
      return sortingFunc(*t1, *t2, dir);
    };

    std::cout << "partitioning in dir = " << dir << std::endl;

    // Have now sorted the objects based on sortingFunction along dir
    std::sort(m_objects.begin(), m_objects.end(), comp); // How to pass into blargh?

  }
}

#endif
