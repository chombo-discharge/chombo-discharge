/*!
  @file   kd_node.cpp
  @brief  Implementation of kd_node.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "kd_node.H"
#include <ParmParse.H>

template<typename T> kd_node<T>::kd_node(){
  this->set_leaf(true);
}

template<typename T> kd_node<T>::kd_node(const RefCountedPtr<kd_node>& a_parent){
  this->set_leaf(true);
  this->set_parent(a_parent);
}

template<typename T> kd_node<T>::~kd_node(){

}

template<typename T> Vector<RefCountedPtr<T> >& kd_node<T>::get_data(){
  return m_data;
}

template<typename T> const RefCountedPtr<kd_node<T> >& kd_node<T>::get_parent() const{
  return m_parent;
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

template<typename T> int kd_node<T>::get_depth(){

  int depth = 0;

  RefCountedPtr<kd_node<T> > parent = m_parent;
  while(!parent.isNull()){
    parent = parent->get_parent();
    depth++;
  }

  return depth;
}

template<typename T> int kd_node<T>::get_num_elements(){
  return m_data.size();
}

template<typename T> Real kd_node<T>::dist(const RealVect a_x0) const{

  RealVect delta;
  for (int dir = 0; dir < SpaceDim; dir++){
    delta[dir] = Max(m_lo[dir] - a_x0[dir], Max(0.0, a_x0[dir] - m_hi[dir]));
  }

  return delta.vectorLength();
}

template<typename T> RealVect kd_node<T>::get_bbox_lo() const {
  return m_lo;
}

template<typename T> RealVect kd_node<T>::get_bbox_hi() const {
  return m_hi;
}

template<typename T> void kd_node<T>::compute_bbox(){

  RealVect lo =  1.234567E89*RealVect::Unit;
  RealVect hi = -1.234567E89*RealVect::Unit;

  for (int i = 0; i < m_data.size(); i++){
    for (int dir = 0; dir < SpaceDim; dir++){
      const RealVect cur_lo = m_data[i]->get_bbox_lo();
      const RealVect cur_hi = m_data[i]->get_bbox_hi();

      if(cur_lo[dir] < lo[dir]){
	lo[dir] = cur_lo[dir];
      }
      if(cur_hi[dir] > hi[dir]){
	hi[dir] = cur_hi[dir];
      }
    }
  }

  this->set_bbox(lo, hi);
}

template<typename T> void kd_node<T>::set_bbox(const RealVect a_lo, const RealVect a_hi){
  m_lo = a_lo;
  m_hi = a_hi;
}

template<typename T> void kd_node<T>::set_leaf(const bool a_leaf){
  m_leaf = a_leaf;
}

template<typename T> void kd_node<T>::set_parent(const RefCountedPtr<kd_node>& a_parent){
  m_parent = a_parent;
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

template<typename T> void kd_node<T>::split_data(const RealVect a_pos,
						 const int      a_dir){

  Vector<RefCountedPtr<T> > data_left;
  Vector<RefCountedPtr<T> > data_right;
  
  for (int i = 0; i < m_data.size(); i++){
    if(m_data[i]->get_coord()[a_dir] < a_pos[a_dir]){
      data_left.push_back(m_data[i]);
    }
    else{
      data_right.push_back(m_data[i]);
    }
  }



  if(!m_left.isNull()){
    m_left->set_data(data_left);
  }
  if(!m_right.isNull()){
    m_right->set_data(data_right);
  }

  pout() << "pos = " << a_pos << "\t dir =" << a_dir << endl;
  pout() << "data size = " << m_data.size() << endl;
  pout() << "left data size = " << data_left.size() << endl;
  pout() << "right data size = " << data_right.size() << endl;
}

template <typename T> void kd_node<T>::split_data(){
  Vector<RefCountedPtr<T> > data_right;
  Vector<RefCountedPtr<T> > data_left;

  Vector<RealVect> coords;
  for (int i = 0; i < m_data.size(); i++){
    coords.push_back(m_data[i]->get_coord());
  }

  RealVect min_coord =  1.234567E89*RealVect::Unit;
  RealVect max_coord = -1.234567E89*RealVect::Unit;
  for (int i = 0; i < m_data.size(); i++){
    for (int dir = 0; dir < SpaceDim; dir++){
      if(coords[i][dir] < min_coord[dir]){
	min_coord[dir] = coords[i][dir];
      }
      if(coords[i][dir] > max_coord[dir]){
	max_coord[dir] = coords[i][dir];
      }
    }
  }

  // Find the longest
  Real longest_len = 0.0;
  int longest_dir;
  for (int dir = 0; dir < SpaceDim; dir++){
    const Real cur_longest = max_coord[dir] - min_coord[dir];
    if(cur_longest > longest_len){
      longest_len = cur_longest;
      longest_dir = dir;
    }
  }

  // Split data on a left/right structure
  const RealVect split_pos = min_coord + 0.5*longest_len*RealVect(BASISV(longest_dir));
  for (int i = 0; i < m_data.size(); i++){
    if(m_data[i]->get_coord()[longest_dir] < split_pos[longest_dir]){
      data_left.push_back(m_data[i]);
    }
    else{
      data_right.push_back(m_data[i]);
    }
  }

  if(!m_left.isNull()){
    m_left->set_data(data_left);
    m_left->compute_bbox();
  }
  if(!m_right.isNull()){
    m_right->set_data(data_right);
    m_right->compute_bbox();
  }

#if 0
  pout() << "pos = " << split_pos << "\t dir =" << longest_dir << endl;
  pout() << "data size = " << m_data.size() << endl;
  pout() << "left data size = " << data_left.size() << endl;
  pout() << "right data size = " << data_right.size() << endl;

#endif
}
