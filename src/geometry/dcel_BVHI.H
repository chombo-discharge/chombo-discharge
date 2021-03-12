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
NodeT<T, BV>::NodeT(std::shared_ptr<Node >& a_parent){
  m_parent = a_parent;
  m_left   = nullptr;
  m_right  = nullptr;
}

template <class T, class BV>
inline
NodeT<T, BV>::~NodeT() {
}

template <class T, class BV>
inline
void NodeT<T, BV>::setParent(const std::shared_ptr<Node>& a_parent) noexcept {
  m_parent = a_parent;
}

template <class T, class BV>
inline
void NodeT<T, BV>::setLeft(const std::shared_ptr<Node>& a_left) noexcept {
  m_left = a_left;
}

template <class T, class BV>
inline
void NodeT<T, BV>::setRight(const std::shared_ptr<Node >& a_right) noexcept{
  m_right = a_right;
}

template <class T, class BV>
inline
bool NodeT<T, BV>::isRoot() const noexcept {
  return m_parent == nullptr;
}

template <class T, class BV>
inline
bool NodeT<T, BV>::isLeaf() const noexcept {
  return m_left == nullptr && m_right == nullptr;
}

template <class T, class BV>
inline
bool NodeT<T, BV>::isPointInside(const Vec3& a_x0) const noexcept {
  return m_bv.isPointInside(a_x0);
}

template <class T, class BV>
inline
bool NodeT<T, BV>::isPointOutside(const Vec3& a_x0) const noexcept {
  return m_bv.isPointOutside(a_x0);
}

template <class T, class BV>
inline
double NodeT<T, BV>::getDistanceToPoint(const Vec3& a_x0) const noexcept {
  return m_bv.getDistanceToPoint(a_x0);
}

template <class T, class BV>
inline
BV& NodeT<T, BV>::getBoundingVolume() noexcept {
  return m_bv;
}

template <class T, class BV>
inline
const BV& NodeT<T, BV>::getBoundingVolume() const noexcept {
  return m_bv;
}

template <class T, class BV>
inline
std::vector<std::shared_ptr<T> >& NodeT<T, BV>::getObjects() noexcept {
  return m_objects;
}

template <class T, class BV>
inline
const std::vector<std::shared_ptr<T> >& NodeT<T, BV>::getObjects() const noexcept {
  return m_objects;
}

template <class T, class BV>
inline
Tree<T, BV>::Tree() {
  m_root = std::make_shared<Node>();
}

template <class T, class BV>
Tree<T, BV>::Tree(std::vector<std::shared_ptr<T> >& a_objects){
  m_root = std::make_shared<Node>();
}

template <class T, class BV>
inline
Tree<T, BV>::~Tree() {

}

#endif
