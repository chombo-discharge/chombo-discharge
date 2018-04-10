/*!
  @file   kd_node.cpp
  @brief  Implementation of kd_node.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "kd_node.H"

template<typename T> kd_node<T>::kd_node(){

  this->set_leaf(false);
}

template<typename T> kd_node<T>::kd_node(const RefCountedPtr<kd_node>& a_root){
  this->set_leaf(false);
  this->set_root(a_root);
}

template<typename T> kd_node<T>::~kd_node(){}

template<typename T> const RefCountedPtr<kd_node<T> >& kd_node<T>::get_root() const{
  return m_root;
}

template<typename T> const RefCountedPtr<kd_node<T> >& kd_node<T>::get_left() const{
  return m_left;
}

template<typename T> const RefCountedPtr<kd_node<T> >& kd_node<T>::get_right() const{
  return m_right;
}

template<typename T> bool kd_node<T>::is_leaf() const {
  return m_leaf;
}

template<typename T> bool kd_node<T>::inside(const RealVect a_x0) const {
  bool inside = true;
  for (int dir = 0; dir < SpaceDim; dir++){
    if(a_x0[dir] < m_lo[dir] || a_x0[dir] > m_hi[dir]){
      inside = false;
      break;
    }
  }
  
  return inside;
}

template<typename T> void kd_node<T>::set_bbox(const RealVect a_lo, const RealVect a_hi){
  m_lo = a_lo;
  m_hi = a_hi;
}

template<typename T> void kd_node<T>::set_leaf(const bool a_leaf){
  m_leaf = a_leaf;
}

template<typename T> void kd_node<T>::set_root(const RefCountedPtr<kd_node>& a_root){
  m_root = a_root;
}

template<typename T> void kd_node<T>::set_left(const RefCountedPtr<kd_node>& a_left){
  m_left = a_left;
}

template<typename T> void kd_node<T>::set_right(const RefCountedPtr<kd_node>& a_right){
  m_right = a_right;
}

template<typename T> void kd_node<T>::set_data(const Vector<RefCountedPtr<T> >& a_data){
  m_data = a_data;
}

template<typename T> void kd_node<T>::clear_data(){
  m_data.resize(0);
}
