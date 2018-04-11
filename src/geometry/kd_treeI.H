/*!
  @file   kd_treeI.H
  @brief  Implementation of kd_tree.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "kd_tree.H"

#include <ParmParse.H>

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
    m_root  = RefCountedPtr<kd_node<T> > (new kd_node<T>());
    m_depth = 0;

    const RealVect lo = this->compute_bbox_lo(a_data);
    const RealVect hi = this->compute_bbox_hi(a_data);
    
    m_root->set_data(a_data);
    m_root->set_leaf(true);
    m_root->set_bbox(lo, hi);

    m_nodes.push_back(m_root);
  }


  
#if 1
  if(a_max_depth > 0){ 

    bool keep_splitting = true;
    while(keep_splitting){
      keep_splitting = false;

      // 1. Iterate through all leaves, split nodes if we have to
      //

      Vector<RefCountedPtr<kd_node<T> > > children(0);
      for (int i = 0; i < m_nodes.size(); i++){

	RefCountedPtr<kd_node<T> >& node = m_nodes[i];
	
	//	const bool split_node = (node->get_num_elements() > a_max_elements) && (node->get_depth()) < a_max_depth;
	const bool split_node =
	  (node->get_depth() < a_max_depth) &&
	  (node->is_leaf()) &&
	  (node->get_depth() < a_max_depth) &&
	  (node->get_num_elements() > a_max_elements);
	
	if(split_node){
	  keep_splitting = true; // Must do another pass to see if we should split more nodes

	  //	  pout() << "splitting node, depth = " << node->get_depth() << endl;
	  
	  const RealVect lo = node->get_bbox_lo();
	  const RealVect hi = node->get_bbox_hi();

	  // Split along dominant dimension
	  Real max = 0.;
	  int split_dir;
	  for (int dir = 0; dir < SpaceDim; dir++){
	    if(hi[dir] - lo[dir] > max){
	      max       = hi[dir] - lo[dir];
	      split_dir = dir;
	    }
	  }
	  const RealVect split_pos = lo + 0.5*max*RealVect(BASISV(split_dir));
	  
	  RefCountedPtr<kd_node<T> > left  = RefCountedPtr<kd_node<T> > (new kd_node<T>(node));
	  RefCountedPtr<kd_node<T> > right = RefCountedPtr<kd_node<T> > (new kd_node<T>(node));

	  // Create children and divide up data
	  node->set_left(left);
	  node->set_right(right);
	  node->split_data(); // This also computes bounding boxes
	  node->clear_data();

	  // Update who is the leaf
	  node->set_leaf(false);
	  left->set_leaf(true);
	  right->set_leaf(true);

	  // Update depth if we went deeper
	  if(left->get_depth() > m_depth){ 
	    m_depth = left->get_depth();
	  }

	  // Store new nodes somewhere
	  children.push_back(left);
	  children.push_back(right);

	  //	  pout() << endl;
	}

      }

      // After each pass, append children to list of nodes
      for (int i = 0; i < children.size(); i++){
	m_nodes.push_back(children[i]);
      }

      //      pout() << endl << endl;
    }
  }
#endif

  CH_assert(!m_root.isNull());
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

template<typename T> Vector<RefCountedPtr<T> > kd_tree<T>::get_candidates(const RealVect a_x0){

  Vector<RefCountedPtr<T> > data;

  if(m_root->is_leaf()){ // Only root
    return m_root->get_data();
  }

  // Get all leaf data
  for (int i = 0; i < m_nodes.size(); i++){
    if(m_nodes[i]->is_leaf()){
      data.append(m_nodes[i]->get_data()); // We will need to search this node
    }
  }

  // Search through nodes in which a_x0 is inside the bbox
#if 0
  data.resize(0);
  Real shortest = 1.E99;
  this->tree_search2(data, shortest, m_root, a_x0);

#else // This appears to be the fastest method by far
  data.resize(0);
  Real shortest = 1.E99;
  Vector<RefCountedPtr<kd_node<T> > > candidates;
  Vector<Real> distances;
  this->tree_search3(candidates, distances, shortest, m_root, a_x0);

  // Do another filtering
  Real min = 1.E99;

  for (int i = 0; i < distances.size(); i++){
    if(distances[i] < min){
      min = distances[i];
    }
  }

  Vector<RefCountedPtr<kd_node<T> > > new_candidates;
  for (int i = 0; i < candidates.size(); i++){
    if(distances[i] <= min){
      new_candidates.push_back(candidates[i]);

    }
  }

  for(int i = 0; i < new_candidates.size(); i++){
    data.append(new_candidates[i]->get_data());
  }
#endif

  return data;
}

template<typename T> void kd_tree<T>::tree_search(Vector<RefCountedPtr<T> >&        a_data,
						  const RefCountedPtr<kd_node<T> >& a_node,
						  const RealVect                    a_x0){

  if(a_node->inside(a_x0)){

    if(a_node->is_leaf()){
      a_data.append(a_node->get_data());
    }
    else{
      RefCountedPtr<kd_node<T> > left  = a_node->get_left();
      RefCountedPtr<kd_node<T> > right = a_node->get_right();

      this->tree_search(a_data, left,  a_x0);
      this->tree_search(a_data, right, a_x0);
    }
  }

  return;
}

template<typename T> void kd_tree<T>::tree_search2(Vector<RefCountedPtr<T> >&        a_data,
						   Real&                             a_shortest,
						   const RefCountedPtr<kd_node<T> >& a_node,
						   const RealVect                    a_x0){

  const Real len = a_node->dist(a_x0);

  if(len <= a_shortest){
    if(a_node->is_leaf()){
      a_shortest = len;
      a_data.append(a_node->get_data());
    }
    else{
      RefCountedPtr<kd_node<T> > left  = a_node->get_left();
      RefCountedPtr<kd_node<T> > right = a_node->get_right();

      this->tree_search2(a_data, a_shortest, left,  a_x0);
      this->tree_search2(a_data, a_shortest, right, a_x0);
    }
  }
  else if(len > a_shortest){

  }
  else{
    if(a_node->is_leaf()){
      a_data.append(a_node->get_data());
    }
  }

  return;
}

template<typename T> void kd_tree<T>::tree_search3(Vector<RefCountedPtr<kd_node<T> > >& a_candidates,
						   Vector<Real>&                        a_distances,
						   Real&                                a_shortest,
						   const RefCountedPtr<kd_node<T> >&    a_node,
						   const RealVect                       a_x0){

  const Real len = a_node->dist(a_x0);

  if(len <= a_shortest){
    if(a_node->is_leaf()){
      a_shortest = len;
      a_candidates.push_back(a_node);
      a_distances.push_back(len);
    }
    else{
      RefCountedPtr<kd_node<T> > left  = a_node->get_left();
      RefCountedPtr<kd_node<T> > right = a_node->get_right();


      const Real l_len = left->dist(a_x0);
      const Real r_len = right->dist(a_x0);
      const Real min   = Min(l_len, r_len);
      
      if(l_len <= min){
	this->tree_search3(a_candidates, a_distances, a_shortest, left,  a_x0);
      }
      if (r_len <= min){
	this->tree_search3(a_candidates, a_distances, a_shortest, right, a_x0);
      }
    }
  }

  return;
}
