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

template <class BV>
Node<BV>::Node() {
  m_parent = nullptr;
  m_left   = nullptr;
  m_right  = nullptr;
}

template <class BV>
Node<BV>::Node(std::shared_ptr<Node<BV> >& a_parent){
  m_parent = a_parent;
  m_left   = nullptr;
  m_right  = nullptr;
}

template <class BV>
Node<BV>::~Node() {

}

template <class BV>
void Node<BV>::setParent(const std::shared_ptr<Node<BV> >& a_parent) noexcept {
  m_parent = a_parent;
}

template <class BV>
void Node<BV>::setLeft(const std::shared_ptr<Node<BV> >& a_left) noexcept {
  m_left = a_left;
}

template <class BV>
void Node<BV>::setRight(const std::shared_ptr<Node<BV> >& a_right) noexcept{
  m_right = a_right;
}

template <class BV>
bool Node<BV>::isRoot() const noexcept {
  return m_parent == nullptr;
}

template <class BV>
bool Node<BV>::isLeaf() const noexcept {
  return m_left == nullptr && m_right == nullptr;
}

template <class BV>
bool Node<BV>::isPointInside(const RealVect& a_x0) const noexcept {
  return m_bv.isPointInside(a_x0);
}

template <class BV>
bool Node<BV>::isPointOutside(const RealVect& a_x0) const noexcept {
  return m_bv.isPointOutside(a_x0);
}

template <class BV>
Real Node<BV>::getDistanceToPoint(const RealVect& a_x0) const noexcept {
  return m_bv.getDistanceToPoint(a_x0);
}

template <class BV>
BV& Node<BV>::getBoundingVolume() noexcept {
  return m_bv;
}

template <class BV>
const BV& Node<BV>::getBoundingVolume() const noexcept {
  return m_bv;
}

template <class BV>
Tree<BV>::Tree() {

}

template <class BV>
Tree<BV>::~Tree() {

}

#endif
