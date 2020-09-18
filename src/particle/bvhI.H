/*!
  @file   bvhI.H
  @brief  Implementation of bvhI.H
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _BVHI_H_
#define _BVHI_H_

#include <chrono>
#include <ctime>
#include <ratio>

using namespace std::chrono;

#define BVH_RESERVE_SIZE 2048
#define BVH_DEBUG 0

template <class T>
bvh_node<T>::bvh_node(){
  m_data.reserve(BVH_RESERVE_SIZE);
}

template <class T>
bvh_node<T>::~bvh_node(){
  m_data.resize(0);
}

template <class T>
bvh_node<T>::bvh_node(std::shared_ptr<bvh_node<T> >& a_parent){
  this->set_parent(a_parent);
}

template <class T>
inline bool bvh_node<T>::is_leaf() const{
  return m_leaf;
}

template <class T>
inline void bvh_node<T>::set_leaf(const bool a_leaf){
  m_leaf = a_leaf;
}

template <class T>
inline void bvh_node<T>::set_parent(std::shared_ptr<bvh_node<T> >& a_parent){
  m_parent = a_parent;
}

template <class T>
inline void bvh_node<T>::set_left(std::shared_ptr<bvh_node<T> >& a_left){
  m_left = a_left;
}

template <class T>
inline void bvh_node<T>::set_right(std::shared_ptr<bvh_node<T> >& a_right){
  m_right = a_right;
}

template <class T>
inline void bvh_node<T>::set_data(const std::vector<T>& a_data, const Real a_mass){
  m_data = a_data;
  m_mass = a_mass;
}

template <class T>
inline void bvh_node<T>::set_mass(const Real a_mass){
  m_mass = a_mass;
}

template <class T>
inline std::shared_ptr<bvh_node<T> >& bvh_node<T>::get_left(){
  return m_left;
}

template <class T>
inline std::shared_ptr<bvh_node<T> >& bvh_node<T>::get_right(){
  return m_right;
}


template <class T>
inline std::shared_ptr<bvh_node<T> >& bvh_node<T>::get_parent(){
  return m_parent;
}


template <class T>
inline std::vector<T>& bvh_node<T>::get_data(){
  return m_data;
}

template <class T>
inline const std::vector<T>& bvh_node<T>::get_data() const{
  return m_data;
}

template <class T>
inline bool bvh_node<T>::can_split() const {
  bool ret = false;
  if(m_data.size() > 0){
    if(m_data.size() > 1 || m_data[0].can_split()){
      ret = true;
    }
  }

  return ret;
}

template <class T>
inline Real bvh_node<T>::mass() const {
  return m_mass;
}

template <class T>
inline void bvh_node<T>::split(const int a_dir){

  // 0. If we made it in here, we KNOW that we can split.
#if BVH_DEBUG // Debug
  if(m_mass < 0.0){
    MayDay::Abort("should not get negative mass!");
  }
#endif

  const int data_size = m_data.size();
  
  // 1. Sort data. Use a lambda for comparing. God I love C++11.
  std::sort(std::begin(m_data), std::end(m_data), [a_dir](const T& p1, const T& p2){return p1[a_dir] < p2[a_dir];});

  // Init median on i = 0
  int splitIndex  = 0;
  Real mass_left  = 0.0;
  Real mass_right = m_mass - m_data[0].mass();

  for (int i = 1; i < data_size; i++){
    const Real& m1 = m_data[splitIndex].mass(); // This is the previous median. 
    const Real& m2 = m_data[splitIndex+1].mass();   // This is, maybe, our new median.

    if(mass_left + m1 < mass_right){ // This is equivalent to shifting the median. 
      mass_left  += m1;              
      mass_right = m_mass - mass_left - m2; // Total mass - left mass - median mass
      
      splitIndex++;
    }
    else{
      break;
    }
  }

  // 2. Make leaves and get hooks for data
  if(!m_left)  m_left  = std::make_shared<bvh_node<T> >();
  if(!m_right) m_right = std::make_shared<bvh_node<T> >();

  std::vector<T>& leftData  = m_left->get_data();
  std::vector<T>& rightData = m_right->get_data();

  leftData.assign(m_data.begin(), m_data.begin() + splitIndex);
  rightData.assign(m_data.begin() + splitIndex + 1, m_data.end());
  
  // If we can split up the median point into several ones, we do it. This would be the case if the median point
  // is a superparticle; in this case we divide the weights as best we can and then assign them to the left/right halves
  const T& splitMass = m_data[splitIndex];
  if(splitMass.can_split()){
    const RealVect& p   = splitMass.pos();
    const Real& mass    = splitMass.mass();       // Always positive
    const Real& energy  = splitMass.energy();     // Energy of split mass particle
    const Real massDiff = mass_right - mass_left; // Not necessarily positive


    // This hook is for when we don't have enough mass to distribute to the two parts. 
    if(mass <= Abs(massDiff)){ 
      if(massDiff > 0){ // Goes into left
	mass_left += mass; 
	leftData.emplace_back(T(p, mass, energy));
      }
      else{
	mass_right += mass; // Goes into right
	rightData.emplace_back(T(p, mass, energy));
      }
    }
    else {
      // Distribute mass as best we can. This might seem weird because we can have mass_left > mass_right in which
      // case massDiff < 0. This still works since we are simply subtracting the difference from the left leaf first and
      // then divide a larger total mass between the two leaves. Essentially, we are subtracting the part of the median
      // particle that gives mass_left > mass_right from the left particle first, and then we reassign a larger total mass.
      Real massLeft  = massDiff;
      Real massRight = 0.0;
      Real massRem   = mass - massDiff;

      const int N  = round(massRem*(1.0 + 1.E-6));
      if(N > 0){ 
	const int NR = N/2;
	const int NL = N - NR;

	massLeft  += (massRem/N)*NL;
	massRight += (massRem/N)*NR;
      }

      const Real energyLeft  = energy*massLeft/mass;
      const Real energyRight = energy*massRight/mass;

      leftData.emplace_back(T(p, massLeft, energyLeft));
      mass_left  += massLeft;

      rightData.emplace_back(T(p, massRight, energyRight));
      mass_right += massRight;
    }
  }
  else{
    //    We should be able to simplify this since it always goes into the left mass(?)
    //    if(mass_left <= mass_right){
      leftData.emplace_back(splitMass);
      mass_left += splitMass.mass();
    // }
    // else{
    //   MayDay::Abort("wtf");
    //   rightData.emplace_back(splitMass);
    //   mass_right += splitMass.mass();
    // }
  }

#if BVH_DEBUG // Make sure we don't break mass
  if(Abs(m_mass - (mass_left + mass_right)) > 1.E-6) MayDay::Abort("bvh_node::split - broke mass");
#endif

  m_left->set_mass(mass_left);
  m_right->set_mass(mass_right);

#if 0 // This breaks. Don't know why. 
  left->set_parent(this->shared_from_this());
  righ->set_parent(this->shared_from_this());
#endif

  // Update who is leaf
  m_left->set_leaf(true);
  m_right->set_leaf(true);
  this->set_leaf(false);
  
}

template <class T>
bvh_tree<T>::bvh_tree(){
  m_leaves.reserve(BVH_RESERVE_SIZE);
  m_new_leaves.reserve(BVH_RESERVE_SIZE);
}

template <class T>
bvh_tree<T>::bvh_tree(std::vector<T>& a_data, const Real a_mass) : bvh_tree() {
  this->define(a_data, a_mass);
}

template <class T>
bvh_tree<T>::~bvh_tree(){
  m_leaves.resize(0);
}

template <class T>
inline void bvh_tree<T>::define(std::vector<T>& a_data, const Real a_mass){

  // Make the root node
  if(!m_root) m_root = std::make_shared<bvh_node<T> >();

  m_root->set_leaf(true);
  m_root->set_data(a_data, a_mass);

  m_leaves.resize(1);
  m_leaves[0] = m_root;
}

template <class T>
inline std::vector<std::shared_ptr<bvh_node<T> > >& bvh_tree<T>::get_leaves(){
  return m_leaves;
}

template <class T>
inline void bvh_tree<T>::build_tree(const int a_numLeaves){

  bool keepSplitting = a_numLeaves > 1;

  // Reset leaves and root
  m_leaves.resize(a_numLeaves, nullptr);
  m_new_leaves.resize(a_numLeaves, nullptr);
  m_leaves[0] = m_root;

  bool can_split;
  int leavesNeeded;
  int cur_size = 1;
  int splitDir = 0;
  int N;
  
  while(keepSplitting){
    can_split = false;

    // Sort the leaves by their mass(?). We will split the ones with the largest masses.
    // std::sort(m_leaves.begin(), m_leaves.begin() + cur_size,
    // 	      [](const std::shared_ptr<bvh_node<T> >& n1,
    // 		 const std::shared_ptr<bvh_node<T> >& n2){
    // 		return n1->mass() < n2->mass();
    // 	      });

    N = 0;
    leavesNeeded = a_numLeaves - cur_size;
    for (int i = 0; i < cur_size; i++){
      
      if(m_leaves[i]->can_split() && leavesNeeded > 0){ // Only split if we actually need more leaves. 
	can_split = true;
#if 1 // Original code
	m_leaves[i]->split(splitDir);
#else
	for (int k = 0; k < 10; k++){
	  m_leaves[i]->split(splitDir);

	  splitDir = (splitDir+1)%SpaceDim;
	}
#endif
	
	m_new_leaves[N]   = m_leaves[i]->get_left(); 
	m_new_leaves[N+1] = m_leaves[i]->get_right();

	N += 2;

	leavesNeeded -= 1; // We only get one extra leaf when we split. 
      }
      else{
	m_new_leaves[N] = m_leaves[i]; N++;
      }
    }

    splitDir = (splitDir + 1) % SpaceDim; // New splitting coordinate

    m_leaves.assign(m_new_leaves.begin(), m_new_leaves.begin() + N);

    cur_size = N;

    // Did we end up in a convergent situation? OK -- then exit.
    keepSplitting = N < a_numLeaves && can_split;
  }
}

#endif
