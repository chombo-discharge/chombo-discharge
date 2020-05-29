/*!
  @file   bvhI.H
  @brief  Implementation of bvhI.H
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _BVHI_H_
#define _BVHI_H_

template <class T>
bvh_node<T>::bvh_node(){
}

template <class T>
bvh_node<T>::bvh_node(std::shared_ptr<bvh_node<T> > a_parent){

}

template <class T>
bool bvh_node<T>::is_leaf() const{
  return m_leaf;
}

template <class T>
void bvh_node<T>::set_leaf(const bool a_leaf){
  m_leaf = a_leaf;
}

template <class T>
void bvh_node<T>::set_parent(std::shared_ptr<bvh_node> a_parent){
  m_parent = a_parent;
}

template <class T>
void bvh_node<T>::set_left(std::shared_ptr<bvh_node> a_left){
  m_left = a_left;
}

template <class T>
void bvh_node<T>::set_right(std::shared_ptr<bvh_node> a_right){
  m_right = a_right;
}

template <class T>
void bvh_node<T>::set_data(const std::vector<T>& a_data){
  m_data = a_data;
}

template <class T>
void bvh_node<T>::sort(const int a_dir){
  MayDay::Abort("bvh_node::sort - not implemented");
}

template <class T>
void bvh_node<T>::split_data(){
  MayDay::Abort("bvh_node::split_data - not implemented");
}

template <class T>
bvh_tree<T>::bvh_tree(std::vector<T>& a_data){

  // Make the root node
  m_root = std::make_shared<bvh_node<T> >();
  m_root->set_leaf(true);
  m_root->set_data(a_data);
}

template <class T>
bvh_tree<T>::~bvh_tree(){
}

template <class T>
std::vector<std::shared_ptr<bvh_node<T> > >& bvh_tree<T>::get_leaves(){
  return m_leaves;
}

template <class T>
void bvh_tree<T>::build_tree(const int a_numLeaves){

  MayDay::Abort("bvh_tree<T>::build_tree - not implmemented");
}



#endif
