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
  this->set_parent(a_parent);
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
void bvh_node<T>::set_parent(std::shared_ptr<bvh_node<T> > a_parent){
  m_parent = a_parent;
}

template <class T>
void bvh_node<T>::set_left(std::shared_ptr<bvh_node<T> > a_left){
  m_left = a_left;
}

template <class T>
void bvh_node<T>::set_right(std::shared_ptr<bvh_node<T> > a_right){
  m_right = a_right;
}

template <class T>
void bvh_node<T>::set_data(const std::vector<T>& a_data){
  m_data = a_data;

  m_mass = 0.0;
  
  for (auto v : m_data){
    m_mass += v.mass();
  }
}

template <class T>
std::vector<T>& bvh_node<T>::get_data(){
  return m_data;
}

template <class T>
const std::vector<T>& bvh_node<T>::get_data() const{
  return m_data;
}

template <class T>
bool bvh_node<T>::can_split() const {
  bool ret = false;
  if(m_data.size() > 0){
    if(m_data.size() > 1 || m_data[0].can_split()){
      ret = true;
    }
  }

  return ret;
}

template <class T>
Real bvh_node<T>::mass() const {
  return m_mass;
}

template <class T>
void bvh_node<T>::split(const int a_dir){

  // 0. If we made it in here, we KNOW that we can split.
  
  // 1. Sort data. Use a lambda for comparing. God I love C++11. 
  std::sort(std::begin(m_data), std::end(m_data), [a_dir](const T& p1, const T& p2){return p1[a_dir] < p2[a_dir];});

  // 2. Find the median. Keep track of the total mass in the leaves. We want to balance the masses
  //    as best as we can. The left leaf will always have an equal or slightly bigger mass than the right leaf. 
  Real mass_left  = 0.0;
  Real mass_right = m_mass;
  int splitIndex = -1;
  for (int i = 0; i < m_data.size(); i++){
    mass_left  += m_data[i].mass();
    mass_right -= m_data[i].mass();

    // We found the point where we should split. 
    if(mass_left > mass_right){
      splitIndex = i;
      break;
    }
  }

  // 3. Everything to the left of splitIndex goes into left data; everything to the right of splitIndex goes into the right data
  std::vector<T> leftData;
  std::vector<T> rightData;
  for (int i = 0; i < m_data.size(); i++){
    if(i < splitIndex){
      leftData.push_back(m_data[i]);
    }
    else if (i > splitIndex){
      rightData.push_back(m_data[i]);
    }
  }

  // 4. If we can split up the median point into several ones, we do it. This would be the case if the median point
  //    is a superparticle; in this case we divide the weights as best we can and then assign them to the left/right halves
  if(m_data[splitIndex].can_split()){
    const Real m  = m_data[splitIndex].mass();
    const int  N  = (int) m*(1.0 + 1.E-6);
    const int  NR = N/2;
    const int  NL = N - NR;

    const RealVect& pos  = m_data[splitIndex].pos();
    const Real massLeft  = m*NL;
    const Real massRight = m*NR;

    leftData.push_back(T(pos, massLeft));
    rightData.push_back(T(pos, massRight));
  }

  // 5. We have the data, we can make the new leaves now.
  auto left = std::make_shared<bvh_node<T> >();
  auto righ = std::make_shared<bvh_node<T> >();

  left->set_data(leftData);
  righ->set_data(leftData);

  left->set_parent(std::shared_ptr<bvh_node<T> >(this));
  righ->set_parent(std::shared_ptr<bvh_node<T> >(this));

  // This is no longer a leaf node. Set the pointers and clear data. 
  this->set_left(left);
  this->set_right(righ);
  this->set_leaf(false);
  m_mass = 0.;
  m_data.clear();

  MayDay::Abort("bvh_node::split - not implemented");
}

template <class T>
bvh_tree<T>::bvh_tree(std::vector<T>& a_data){

  // Make the root node
  m_root = std::make_shared<bvh_node<T> >();
  m_root->set_leaf(true);
  m_root->set_data(a_data);

  m_nodes.push_back(m_root);
  m_leaves.push_back(m_root);
}

template <class T>
bvh_tree<T>::~bvh_tree(){
}

template <class T>
std::vector<std::shared_ptr<bvh_node<T> > >& bvh_tree<T>::get_leaves(){
  return m_leaves;
}

template <class T>
std::vector<std::shared_ptr<bvh_node<T> > >& bvh_tree<T>::get_nodes(){
  return m_nodes;
}

template <class T>
void bvh_tree<T>::build_tree(const int a_numLeaves){

  int splitDir = 0;

  bool keepSplitting = true;
  while(keepSplitting){

    // Sort leaves by their mass, and then split them
    std::sort(std::begin(m_leaves), std::end(m_leaves), [](const bvh_node<T>& n1, const bvh_node<T>& n2){
	return n1.mass() < n2.mass();});

    // Split all the leaves into new ones.
    for(auto v : m_leaves){
      if(v->can_split()){
	v->split(splitDir);
      }
    }
    keepSplitting = false;
  }

  
  // Build all the leaves. This is a stupid way of doing it but stupid is what I need right now. 
  for (int i = 0; i < m_nodes.size(); i++){
    if(m_nodes[i]->is_leaf()){
      m_leaves.push_back(m_nodes[i]);
    }
  }
  
  //  MayDay::Abort("bvh_tree<T>::build_tree - not implmemented");
}



#endif
