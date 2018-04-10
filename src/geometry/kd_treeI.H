/*!
  @file   kd_treeI.H
  @brief  Implementation of kd_tree.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "kd_tree.H"

template<typename T> kd_tree<T>::kd_tree(){

}

template<typename T> kd_tree<T>::kd_tree(const Vector<RefCountedPtr<T> >& a_data,
					 const int                       a_max_depth,
					 const int                       a_max_elements){
  this->build_tree(a_data, a_max_depth, a_max_elements);
}

template<typename T> kd_tree<T>::~kd_tree(){

}

template<typename T> void kd_tree<T>::build_tree(const Vector<RefCountedPtr<T> >& a_data,
						 const int                      a_max_depth,
						 const int                      a_max_elements){
  if(m_root.isNull()){ // Create root
    m_root = RefCountedPtr<kd_node<T> > (new kd_node<T>());


    const RealVect lo = this->compute_bbox_lo(a_data);
    const RealVect hi = this->compute_bbox_hi(a_data);

    m_root->set_data(a_data);
    m_root->set_leaf(true);
    m_root->set_bbox(lo, hi);

    m_nodes.push_back(m_root);
  }

  if(a_max_depth > 0){
    
  }
}

template<typename T> RealVect kd_tree<T>::compute_bbox_lo(const Vector<RefCountedPtr<T> >& a_data){
  RealVect ret = 1.234567E89*RealVect::Unit;

  for (int i = 0; i < a_data.size(); i++){
    for (int dir = 0; dir < SpaceDim; dir++){
      const RealVect cur = a_data[i]->get_bbox_lo();
      if(cur[dir] < ret[dir]){
	ret[dir] = cur[dir];
      }
    }
  }

  return ret;
}

template<typename T> RealVect kd_tree<T>::compute_bbox_hi(const Vector<RefCountedPtr<T> >& a_data){
  RealVect ret = -1.234567E89*RealVect::Unit;

  for (int i = 0; i < a_data.size(); i++){
    for (int dir = 0; dir < SpaceDim; dir++){
      const RealVect cur = a_data[i]->get_bbox_hi();
      if(cur[dir] > ret[dir]){
	ret[dir] = cur[dir];
      }
    }
  }

  return ret;
}

template<typename T> bool kd_tree<T>::left(const RefCountedPtr<T>& a_data, const RealVect a_pos, const int a_dir){
  bool is_left = true;
  if(a_data->get_coord()[a_dir] > a_pos[a_dir]){
    is_left = false;
  }

  return is_left;
}

template <typename T> Vector<RefCountedPtr<T> > kd_tree<T>::get_candidates(const RealVect a_x0){

}
