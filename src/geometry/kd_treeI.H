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

  bool keep_splitting = a_max_depth > 0;
  while(keep_splitting){
    keep_splitting = false;

    Vector<RefCountedPtr<kd_node<T> > > children(0);
    for (int i = 0; i < m_nodes.size(); i++){

      RefCountedPtr<kd_node<T> >& node = m_nodes[i];
	
      const bool split_node =
	(node->get_depth() < a_max_depth) &&          // Can we go deeper?
	(node->is_leaf()) &&                          // Is this a leaf node?
	(node->get_num_elements() > a_max_elements);  // How many elements are in the leaf node?

      if(split_node){
	keep_splitting = true; // Must do another pass to see if we should split more nodes
	  
	const RealVect lo = node->get_bbox_lo();
	const RealVect hi = node->get_bbox_hi();
	  
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
      }
    }

    // After each pass, append children to list of nodes
    for (int i = 0; i < children.size(); i++){
      m_nodes.push_back(children[i]);
    }
  }


  // Make sure that we find a leaf
  for (int i = 0; i < m_nodes.size(); i++){
    if(m_nodes[i]->is_leaf()){
      m_leaf = m_nodes[i];
      break;
    }
  }

  for (int i = 0; i < m_nodes.size(); i++){
    if(m_nodes[i]->is_leaf()){
      m_leaves.push_back(m_nodes[i]);
    }
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

template<typename T> Vector<RefCountedPtr<T> > kd_tree<T>::get_candidates(const RealVect a_x0){
#define print_time 0
#define check_all 0
  Vector<RefCountedPtr<T> > data;

  Real shortest = m_leaf->dist(a_x0);
  shortest = 1.E99;
  Vector<RefCountedPtr<kd_node<T> > > candidates;
  Vector<Real> distances;

#if check_all
  // Look through all leaves
  shortest = 1.E99;
  data.resize(0);
  for(int i = 0; i < m_leaves.size(); i++){
    if(m_leaves[i]->dist(a_x0) <= shortest){
      shortest = m_leaves[i]->dist(a_x0);
    }
  }

  data.resize(0);
  for (int i = 0; i < m_leaves.size(); i++){
    if(m_leaves[i]->dist(a_x0) <= shortest){
      data.append(m_leaves[i]->get_data());
    }
  }
#else
#if print_time
  const Real t1 = MPI_Wtime();
#endif
  this->tree_search3(candidates, distances, shortest, m_root, a_x0); // Tree traversal
#if print_time
  const Real t2 = MPI_Wtime();
#endif

  // shortest might have changed, eliminate further candidates
  for (int i = 0; i < distances.size(); i++){
    if(candidates[i]->dist(a_x0) <= shortest){
      data.append(candidates[i]->get_data());
    }
  }
  
#if print_time
  const Real t3 = MPI_Wtime();
  pout() << "total time = << t4 - t1 " << "\tsearch time = " << t2-t1 << "\t filter time = " << t3-t2 << endl;
#endif

#endif

  return data;
}

template<typename T> void kd_tree<T>::tree_search(Vector<RefCountedPtr<T> >&        a_data,
						  const RefCountedPtr<kd_node<T> >& a_node,
						  const RealVect                    a_x0){

  MayDay::Abort("kd_tree<T>::tree_search - should not be called");
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
  MayDay::Abort("kd_tree<T>::tree_search2 - should not be called");
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

      this->tree_search3(a_candidates, a_distances, a_shortest, left,  a_x0); // Returns immediately if left->dist() >  a_shortest
      this->tree_search3(a_candidates, a_distances, a_shortest, right, a_x0); // Returns immediately if right->dist() > a_shortest
    }
  }

  return;
}
