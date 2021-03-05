/*!
  @file   kd_treeI.H
  @brief  Implementation of kd_tree.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "kd_tree.H"
#include <chrono>
#include <ctime> 
#include <iostream>
#include <random>

#include <ParmParse.H>
#include <PolyGeom.H>

#include <chrono>

#define TREE_EPSILON 1.E-6

template<typename T>
kd_tree<T>::kd_tree(){
  m_root = nullptr;
  m_isDefined = false;
}

template<typename T>
kd_tree<T>::kd_tree(const std::vector<std::shared_ptr<T> >& a_data,
		    const int                               a_max_depth,
		    const int                               a_max_elements){
  m_root = nullptr;
  
  this->build_tree(a_data, a_max_depth, a_max_elements);


}

template<typename T>
kd_tree<T>::~kd_tree(){

}

template<typename T>
void kd_tree<T>::build_tree(const std::vector<std::shared_ptr<T> >& a_data,
			    const int                             a_max_depth,
			    const int                             a_max_elements){

  int datasize = a_data.size();
  if(m_root == nullptr){ // Create root
    m_root  = std::shared_ptr<kd_node<T> > (new kd_node<T>());
    m_depth = 0;
    
    m_root->set_data(a_data);
    m_root->set_leaf(true);
    m_root->compute_bbox();

    m_nodes.push_back(m_root);
  }

  bool keep_splitting = a_max_depth > 0;
  while(keep_splitting){
    keep_splitting = false;
    
    const int num_nodes = m_nodes.size();
    std::vector<std::shared_ptr<kd_node<T> > > children(0);
    
    for (int i = 0; i < num_nodes; i++){

      std::shared_ptr<kd_node<T> >& node = m_nodes[i];
	
      const bool split_node =
	(node->get_depth() < a_max_depth) &&          // Can we go deeper?
	(node->is_leaf()) &&                          // Is this a leaf node?
	(node->get_num_elements() > a_max_elements);  // How many elements are in the leaf node?

      if(split_node){
	keep_splitting = true; // Must do another pass to see if we should split more nodes
	  
	std::shared_ptr<kd_node<T> > left  = std::shared_ptr<kd_node<T> > (new kd_node<T>(node));
	std::shared_ptr<kd_node<T> > right = std::shared_ptr<kd_node<T> > (new kd_node<T>(node));

	// Create children and divide up data
	node->set_left(left);
	node->set_right(right);
	node->split_data(); // This also computes bounding boxes for children
	node->clear_data(); // Clears data for this node

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
    m_nodes.insert(m_nodes.end(), children.begin(), children.end());
  }

  // Make sure that we find a leaf
  for (int i = 0; i < m_nodes.size(); i++){
    if(m_nodes[i]->is_leaf()){
      m_leaf = m_nodes[i];
      break;
    }
  }

  int dataleaf = 0;
  for (int i = 0; i < m_nodes.size(); i++){
    if(m_nodes[i]->is_leaf()){
      m_leaves.push_back(m_nodes[i]);
    }
  }

  for (int i = 0; i < m_leaves.size(); i++){
    dataleaf += m_leaves[i]->get_data().size();
  }

  m_isDefined = true;
}

template<typename T>
std::vector<std::shared_ptr<T> > kd_tree<T>::find_closest(const RealVect a_x0){

  std::vector<std::shared_ptr<kd_node<T> > > queue;
  std::vector<Real> dist;

  queue.push_back(m_root);
  dist.push_back(m_root->dist(a_x0));

#if 0  
  auto start = std::chrono::system_clock::now();
  const Real schme = 5*5 - PolyGeom::dot(a_x0, a_x0);
  auto end = std::chrono::system_clock::now();
  
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s"
	    << "\t #candiates = " << queue[0]->get_data().size() << std::endl;
#endif
  this->tree_search(queue, dist, a_x0);


  return queue[0]->get_data();
}

template<typename T>
void kd_tree<T>::tree_search(std::vector<std::shared_ptr<kd_node<T> > >& a_queue,
			     std::vector<Real>&                        a_distances,
			     const RealVect                            a_x0){

  // We assume that a_queue and a_distances is sorted on the way in
  const int Qmax = 10;

  if(!a_queue[0]->is_leaf()){ // Examine children
    std::shared_ptr<kd_node<T> > left  = a_queue[0]->get_left();
    std::shared_ptr<kd_node<T> > right = a_queue[0]->get_right();

    const Real left_dist  = left->dist(a_x0);
    const Real right_dist = right->dist(a_x0);

    a_queue.push_back(left);
    a_queue.push_back(right);

    a_distances.push_back(left_dist);
    a_distances.push_back(right_dist);

    // Remove the element that we opened. Shift the list and pop from the back
    const int N = a_queue.size(); 
    for (int i = 0; i < a_queue.size() - 1; i++){
      a_queue[i]     = a_queue[i+1];
      a_distances[i] = a_distances[i+1];
    }
    a_queue.pop_back();
    a_distances.pop_back();
    
    this->insertion_sort(a_queue, a_distances);
    this->trim(a_queue, a_distances, Qmax);

    // After sorting, if the first element is a node, we need to keep going
    // There's more to this because there can be multiple nodes/leaves with the same priority
    if(!a_queue[0]->is_leaf()){
      this->tree_search(a_queue, a_distances, a_x0);
    }
  }
}

template<typename T>
void kd_tree<T>::insertion_sort(std::vector<std::shared_ptr<kd_node<T> > >& a_queue,
				std::vector<Real>&                        a_distances){

  const int N = a_queue.size();
  
  for (int i = 1; i < N; i++){

    const std::shared_ptr<kd_node<T> > cur_elem  = a_queue[i];
    const Real cur_dist = a_distances[i];
    
    int j = i;
    while(j > 0 && a_distances[j-1] > a_distances[j]){

      const std::shared_ptr<kd_node<T> > cur_elem = a_queue[j];
      const Real cur_dist = a_distances[j];

      a_queue[j]   = a_queue[j-1];
      a_queue[j-1] = cur_elem;

      a_distances[j]   = a_distances[j-1];
      a_distances[j-1] = cur_dist;

      j--;
    }
  }
}

template<typename T>
void kd_tree<T>::trim(std::vector<std::shared_ptr<kd_node<T> > >& a_queue,
		      std::vector<Real>&                          a_distances,
		      const int                                   a_max_elem){
  while(a_queue.size() > a_max_elem){
    a_queue.pop_back();
    a_distances.pop_back();
  }
}
