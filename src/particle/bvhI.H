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

#define BVH_DEBUG 0

template <class T>
bvh_node<T>::bvh_node(){
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
bool bvh_node<T>::is_leaf() const{
  return m_leaf;
}

template <class T>
void bvh_node<T>::set_leaf(const bool a_leaf){
  m_leaf = a_leaf;
}

template <class T>
void bvh_node<T>::set_parent(std::shared_ptr<bvh_node<T> >& a_parent){
  m_parent = a_parent;
}

template <class T>
void bvh_node<T>::set_left(std::shared_ptr<bvh_node<T> >& a_left){
  m_left = a_left;
}

template <class T>
void bvh_node<T>::set_right(std::shared_ptr<bvh_node<T> >& a_right){
  m_right = a_right;
}

template <class T>
void bvh_node<T>::set_data(const std::vector<T>& a_data, const Real a_mass){
  m_data = a_data;
  m_mass = a_mass;
}

template <class T>
void bvh_node<T>::set_mass(const Real a_mass){
  m_mass = a_mass;
}

template <class T>
std::shared_ptr<bvh_node<T> >& bvh_node<T>::get_left(){
  return m_left;
}

template <class T>
std::shared_ptr<bvh_node<T> >& bvh_node<T>::get_right(){
  return m_right;
}

template <class T>
std::shared_ptr<bvh_node<T> >& bvh_node<T>::get_parent(){
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
#if BVH_DEBUG // Debug
  if(m_mass < 0.0){
    MayDay::Abort("should not get negative mass!");
  }
#endif
  
  // 1. Sort data. Use a lambda for comparing. God I love C++11.
  std::sort(std::begin(m_data), std::end(m_data), [a_dir](const T& p1, const T& p2){return p1[a_dir] < p2[a_dir];});

  // 2. Make leaves and get hooks for data
  m_left  = std::make_shared<bvh_node<T> >();
  m_right = std::make_shared<bvh_node<T> >();

  std::vector<T>& leftData  = m_left->get_data();
  std::vector<T>& rightData = m_right->get_data();

  leftData.resize(0);
  rightData.resize(0);

  // Init median on i = 0
  int splitIndex  = 0;
  Real mass_left  = 0.0;
  Real mass_right = m_mass - m_data[0].mass();

  // Left particles, keep track of total mass. 
  for (int i = 1; i < m_data.size(); i++){
    const Real& m1 = m_data[splitIndex].mass(); // This is the previous median. 
    const Real& m2 = m_data[splitIndex+1].mass();   // This is, maybe, our new median.

    if(mass_left + m1 < mass_right){ // This is equivalent to shifting the median. The previous median goes into the left
      leftData.push_back(m_data[i-1]); // list and the median is shifted to the 
      mass_left  += m1;
      mass_right = m_mass - mass_left - m2; // Total mass - left mass - median mass
      
      splitIndex++;
    }
    else{
      break;
    }
  }
 
  // Add the right particles
  for (int i = splitIndex+1; i < m_data.size(); i++){
    rightData.push_back(m_data[i]);
  }

  
  // If we can split up the median point into several ones, we do it. This would be the case if the median point
  // is a superparticle; in this case we divide the weights as best we can and then assign them to the left/right halves
  const T& splitMass = m_data[splitIndex];
  if(splitMass.can_split()){
    const RealVect& p   = splitMass.pos();
    const Real& mass    = splitMass.mass();
    const Real massDiff = mass_right - mass_left;

    // This hook is for when we don't have enough mass to distribute to the two parts. 
    if(mass <= Abs(massDiff)){ 
      if(massDiff > 0){ // Goes into left
	mass_left += mass; 
	leftData.push_back(T(p, mass));
      }
      else{
	mass_right += mass; // Goes into right
	rightData.push_back(T(p, mass));
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

      if(massLeft > 0.0){
	leftData.push_back(T(p, massLeft));
	mass_left  += massLeft;
      }

      if(massRight > 0.0){
	rightData.push_back(T(p, massRight));
	mass_right += massRight;
      }

#if ITO_DEBUG
      if(massLeft <= 0.0) MayDay::Abort("bvh_node<T>::split - should not get negative left mass");
      if(massRight <= 0.0) MayDay::Abort("bvh_node<T>::split - should not get negative right mass");
#endif
    }
  }
  else{
    //    We should be able to simplify this since it always goes into the left mass(?)
    if(mass_left <= mass_right){
      leftData.push_back(splitMass);
      mass_left += splitMass.mass();
    }
    else{
      MayDay::Abort("how?????");
      rightData.push_back(splitMass);
      mass_right += splitMass.mass();
    }
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
  
  // m_mass = 0.;
  // m_data.clear();
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

#if BVH_DEBUG
  const int numParts = m_leaves[0]->get_data().size();
  duration<double> tSort(0.0);//   = duration_cast<duration<double> >(0.0);
  duration<double> tSplit(0.0);//   = duration_cast<duration<double> >(0.0);
  duration<double> tAssign(0.0);//   = duration_cast<duration<double> >(0.0);
  duration<double> tTotal(0.0);//   = duration_cast<duration<double> >(0.0);
  auto t0 = std::chrono::high_resolution_clock::now();
#endif
  
  while(keepSplitting){

    std::vector<std::shared_ptr<bvh_node<T> > > new_leaves;
    bool can_split = false;
    
    // Sort the leaves by their mass. We will split the ones with the largest masses.
#if BVH_DEBUG
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    std::sort(std::begin(m_leaves), std::end(m_leaves),
    	      [](const std::shared_ptr<bvh_node<T> >& n1,
    		 const std::shared_ptr<bvh_node<T> >& n2){
    		return n1->mass() < n2->mass();
    	      });
#if BVH_DEBUG
    auto t2 = std::chrono::high_resolution_clock::now();
#endif

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
#if BVH_DEBUG
    auto t3 = std::chrono::high_resolution_clock::now();
#endif

    splitDir = (splitDir + 1) % SpaceDim; // New splitting coordinate

    // New leaves
    m_leaves.resize(new_leaves.size());
    for (int i = 0; i < new_leaves.size(); i++){
      m_leaves[i] = new_leaves[i];
    }

    // Did we end up in a convergent situation? OK, then exit. 
    if(m_leaves.size() >= a_numLeaves){
      keepSplitting = false;
    }
    if(!can_split){
      keepSplitting = false;
    }

#if BVH_DEBUG
    auto t4 = std::chrono::high_resolution_clock::now();
    
    tSort   += duration_cast<duration<double> > (t2-t1);
    tSplit  += duration_cast<duration<double> > (t3-t2);
    tAssign += duration_cast<duration<double> > (t4-t3);
#endif

  }

#if BVH_DEBUG
  auto tf = std::chrono::high_resolution_clock::now();
  tTotal = duration_cast<duration<double> >(tf-t0);

  pout() << "sort time = "   << 100.*tSort.count()/tTotal.count() << "%" << endl
	 << "split time = "  << 100.*tSplit.count()/tTotal.count() << "%" << endl
	 << "assign time = " << 100.*tAssign.count()/tTotal.count() << "%" << endl
    	 << "Total time = " << tTotal.count() << endl
	 << "Particles = " << numParts << endl
	 << "Time per part = " << tTotal.count()/numParts << endl
	 << endl;
#endif
}

#endif
