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
Node<T, BV>::Node() {
  m_parent = nullptr;
  m_left   = nullptr;
  m_right  = nullptr;
}

template <class T, class BV>
Node<T, BV>::Node(std::shared_ptr<Node<T, BV> >& a_parent){
  m_parent = a_parent;
  m_left   = nullptr;
  m_right  = nullptr;
}

template <class T, class BV>
Node<T, BV>::~Node() {
}

template <class T, class BV>
void Node<T, BV>::setParent(const std::shared_ptr<Node<T, BV> >& a_parent) noexcept {
  m_parent = a_parent;
}

template <class T, class BV>
void Node<T, BV>::setLeft(const std::shared_ptr<Node<T, BV> >& a_left) noexcept {
  m_left = a_left;
}

template <class T, class BV>
void Node<T, BV>::setRight(const std::shared_ptr<Node<T, BV> >& a_right) noexcept{
  m_right = a_right;
}

template <class T, class BV>
bool Node<T, BV>::isRoot() const noexcept {
  return m_parent == nullptr;
}

template <class T, class BV>
bool Node<T, BV>::isLeaf() const noexcept {
  return m_left == nullptr && m_right == nullptr;
}

template <class T, class BV>
bool Node<T, BV>::isPointInside(const RealVect& a_x0) const noexcept {
  return m_bv.isPointInside(a_x0);
}

template <class T, class BV>
bool Node<T, BV>::isPointOutside(const RealVect& a_x0) const noexcept {
  return m_bv.isPointOutside(a_x0);
}

template <class T, class BV>
Real Node<T, BV>::getDistanceToPoint(const RealVect& a_x0) const noexcept {
  return m_bv.getDistanceToPoint(a_x0);
}

template <class T, class BV>
BV& Node<T, BV>::getBoundingVolume() noexcept {
  return m_bv;
}

template <class T, class BV>
const BV& Node<T, BV>::getBoundingVolume() const noexcept {
  return m_bv;
}

template <class T, class BV>
std::vector<std::shared_ptr<T> >& Node<T, BV>::getObjects() noexcept {
  return m_objects;
}

template <class T, class BV>
const std::vector<std::shared_ptr<T> >& Node<T, BV>::getObjects() const noexcept {
  return m_objects;
}

template <class T, class BV>
Tree<T, BV>::Tree() {

}

template <class T, class BV>
Tree<T, BV>::~Tree() {

}

#endif
