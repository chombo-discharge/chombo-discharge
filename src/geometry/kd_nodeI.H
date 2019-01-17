/*!
  @file   kd_nodeI.H
  @brief  Implementation of kd_node.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "kd_node.H"
#include <ParmParse.H>

template<typename T> kd_node<T>::kd_node(){
  this->set_leaf(true);
}

template<typename T> kd_node<T>::kd_node(const RefCountedPtr<kd_node<T> >& a_parent){
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

template<typename T> Real kd_node<T>::dist(const RealVect a_x0) const {
  if(m_leaf){
    return data_dist(a_x0);
  }
  else{
    return aabb_dist(a_x0);
  }
}

template<typename T> Real kd_node<T>::aabb_dist(const RealVect a_x0) const{
  RealVect delta;
  const Real zero = 0.0;
  for (int dir = 0; dir < SpaceDim; dir++){
    delta[dir] = Max(m_lo[dir] - a_x0[dir], Max(zero, a_x0[dir] - m_hi[dir]));
  }

  return delta.vectorLength();
}

template<typename T> Real kd_node<T>::data_dist(const RealVect a_x0) const {
  Real dist = 1.234567E89;
  for (int i = 0; i < m_data.size(); i++){
    const Real cur_dist = Abs(m_data[i]->signed_distance(a_x0));

    dist = (cur_dist < dist) ? cur_dist : dist;
  }

  return dist;
}

template<typename T> void kd_node<T>::compute_bbox(){

  // Compute AABB (axis-aligned bounding box)
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

  m_lo = lo;
  m_hi = hi;

  
#if 1 // Safety test, not part of final code
  // Corner test
  for (int dir = 0; dir < SpaceDim; dir++){
    CH_assert(lo[dir] <= hi[dir]);
  }

  // Check that every data point is inside the bounding box
  for (int i = 0; i < m_data.size(); i++){
    Vector<RealVect> points = m_data[i]->get_points();

    for (int j = 0; j < points.size(); j++){
      for (int dir = 0; dir < SpaceDim; dir++){
	if(points[j][dir] < lo[dir] || points[j][dir] > hi[dir]){
	  MayDay::Abort("kd_node<T>::compute_bbox - points are not contained in bbox. Something is wrong.");
	}
      }
    }
  }
#endif

// #if 1 // Debug
//   std::cout << m_lo << "\t" << m_hi << std::endl;
// #endif
//   MayDay::Abort("stop");
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

template <typename T> void kd_node<T>::split_data(){
  CH_assert(m_leaf);
  CH_assert(!m_left.isNull());
  CH_assert(!m_right.isNull());
  
  Vector<RefCountedPtr<T> > data_right;
  Vector<RefCountedPtr<T> > data_left;

  // Find maximum and minimum coordinates (i.e. the bounding box of the data centroids)
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

  // Find the longest dimension
  Real longest_len = 0.0;
  int longest_dir = 0;
  for (int dir = 0; dir < SpaceDim; dir++){
    const Real cur_longest = max_coord[dir] - min_coord[dir];
    CH_assert(cur_longest >= 0.0);
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
  else{
    MayDay::Abort("kd_node<T>::split_data - left is NULL");
  }
  if(!m_right.isNull()){
    m_right->set_data(data_right);
    m_right->compute_bbox();
  }
  else{
    MayDay::Abort("kd_node<T>::split_data - right is NULL");
  }
}
