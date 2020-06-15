/*!
  @file   bvhI.H
  @brief  Implementation of bvhI.H
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _BVHI_H_
#define _BVHI_H_

#define DEBUG 0

template <class T>
bvh_node<T>::bvh_node(){
}

template <class T>
bvh_node<T>::~bvh_node(){
  m_data.resize(0);
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
void bvh_node<T>::set_data(const std::vector<T>& a_data, const Real a_mass){
  m_data = a_data;
  m_mass = a_mass;
}

template <class T>
std::shared_ptr<bvh_node<T> > bvh_node<T>::get_left(){
  return m_left;
}

template <class T>
std::shared_ptr<bvh_node<T> > bvh_node<T>::get_right(){
  return m_right;
}

template <class T>
std::shared_ptr<bvh_node<T> > bvh_node<T>::get_parent(){
  return m_parent;
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
#if DEBUG // Debug
  if(m_mass < 0.0){
    MayDay::Abort("should not get negative mass!");
  }
#endif
  
  // 1. Sort data. Use a lambda for comparing. God I love C++11.
  std::sort(std::begin(m_data), std::end(m_data), [a_dir](const T& p1, const T& p2){return p1[a_dir] < p2[a_dir];});

  // 2. Add left and right particles but don't add the particle on the median yet. 
  std::vector<T> leftData;
  std::vector<T> rightData;
  
  Real mass_left = 0.0;
  Real mass_right = m_mass;

  int splitIndex = -1;
  const int dataSize = m_data.size();

  // Left particles, keep track of total mass. 
  for (int i = 0; i < dataSize; i++){
    const Real mass = m_data[i].mass();
    if(mass_left < mass_right){
      leftData.push_back(m_data[i]);
      mass_left  += mass;
      mass_right -= mass;
    }
    else{
      splitIndex = i; // Don't add anything else
      break;
    }
  }

  // Add the right particles
  mass_right = 0.0;
  for (int i = splitIndex+1; i < dataSize; i++){
    rightData.push_back(m_data[i]);
    mass_right += m_data[i].mass();
  }

  // 4. If we can split up the median point into several ones, we do it. This would be the case if the median point
  //    is a superparticle; in this case we divide the weights as best we can and then assign them to the left/right halves
  const T splitMass = m_data[splitIndex];
  if(false){//splitMass.can_split()){
    const Real m  = splitMass.mass();
    const int  N  = round(m*(1.0 + 1.E-6));
    const int  NR = N/2;
    const int  NL = N - NR;

    const Real massLeft  = (m/N)*NL;
    const Real massRight = (m/N)*NR;

#if DEBUG
    if(m < 0.0){
      MayDay::Abort("got negative mass 2");
    }
    if(NL < 0.0){
      std::cout << m << "\t" << N << "\t" << NR << "\t" << NL << "\t" << m_data.size() << "\t idx = " << splitIndex << std::endl;
      MayDay::Abort("got negative NL");
    }
    if(massLeft < 0.0){
      std::cout << "massLeft = " << massLeft << std::endl;
      MayDay::Abort("stop left");
    }
    if(massRight < 0.0){
      std::cout << "massRight = " << massRight << std::endl;
      MayDay::Abort("stop right");
    }
#endif

    // Split the divided particle between left and right data
    const RealVect pos  = splitMass.pos();
    
    leftData.push_back(T(pos, massLeft));
    rightData.push_back(T(pos, massRight));

    mass_left  += massLeft;
    mass_right += massRight;
  }
  else{
    if(mass_left <= mass_right){
      leftData.push_back(splitMass);
      mass_left  += splitMass.mass();
    }
    else{
      rightData.push_back(splitMass);
      mass_right  += splitMass.mass();
    }
  }

#if DEBUG // Debug
  if(!(rightData.size() > 0)) {
    std::cout << "bad right data size = " << m_data.size() << "\t mass = " << m_mass << "\t idx = " << splitIndex << std::endl;
    MayDay::Abort("right data size must be > 0");
  }
  if(!(leftData.size()  > 0)) {
    std::cout << "bad left data size = " << m_data.size() << "\t mass = " << m_mass << "\t idx = " << splitIndex << std::endl;
    MayDay::Abort("left data size must be > 0");
  }
#endif

  // 5. We have the data, we can make the new leaves now.
  auto left = std::make_shared<bvh_node<T> >();
  auto righ = std::make_shared<bvh_node<T> >();

  left->set_data(leftData, mass_left);
  righ->set_data(rightData, mass_right);

#if 0 // This breaks. Don't know why. 
  left->set_parent(this->shared_from_this());
  righ->set_parent(this->shared_from_this());
#endif

  left->set_leaf(true);
  righ->set_leaf(true);

  // This is no longer a leaf node. Set the pointers and clear data. 
  this->set_left(left);
  this->set_right(righ);
  this->set_leaf(false);
  m_mass = 0.;
  m_data.clear();
}

template <class T>
bvh_tree<T>::bvh_tree(std::vector<T>& a_data, const Real a_mass){

  // Make the root node
  m_root = std::make_shared<bvh_node<T> >();
  m_root->set_leaf(true);
  
  m_root->set_data(a_data, a_mass);

  m_leaves.resize(0);
  m_leaves.push_back(m_root);
}

template <class T>
bvh_tree<T>::~bvh_tree(){
  m_leaves.resize(0);
}

template <class T>
std::vector<std::shared_ptr<bvh_node<T> > >& bvh_tree<T>::get_leaves(){
  return m_leaves;
}

template <class T>
void bvh_tree<T>::build_tree(const int a_numLeaves){
  int splitDir = 0;

  bool keepSplitting = m_leaves.size() < a_numLeaves;
  
  while(keepSplitting){

    std::vector<std::shared_ptr<bvh_node<T> > > new_leaves;
    bool can_split = false;
    
    // Sort the leaves by their mass. We will split the ones with the largest masses.
    std::sort(std::begin(m_leaves), std::end(m_leaves),
    	      [](const std::shared_ptr<bvh_node<T> >& n1,
    		 const std::shared_ptr<bvh_node<T> >& n2){
    		return n1->mass() < n2->mass();
    	      });

    int leavesNeeded = a_numLeaves - m_leaves.size();
    for(int i = 0; i < m_leaves.size(); i++){
      if(m_leaves[i]->can_split() && leavesNeeded > 0){ // Only split if we actually need more leaves. 
	can_split = true;
	  
	m_leaves[i]->split(splitDir);
	new_leaves.push_back(m_leaves[i]->get_left());
	new_leaves.push_back(m_leaves[i]->get_right());

	leavesNeeded -= 1; // We only get one extra leaf when we split. 
      }
      else{
	new_leaves.push_back(m_leaves[i]);
      }
    }

    splitDir = (splitDir + 1) % SpaceDim; // New splitting coordinate
    m_leaves = new_leaves;                // New leaves.

    // Did we end up in a convergent situation? OK, then exit. 
    if(m_leaves.size() >= a_numLeaves){
      keepSplitting = false;
    }
    if(!can_split){
      keepSplitting = false;
    }
  }
}

#endif
