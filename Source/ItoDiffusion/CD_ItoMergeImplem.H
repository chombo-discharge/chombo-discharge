/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ItoMergeImplem.H
  @brief  Implementation of CD_ItoMerge.H
  @author Robert Marskar
  @date   May 2020
*/

#ifndef CD_ItoMergeImplem_H
#define CD_ItoMergeImplem_H

// Std includes
#include <chrono>
#include <ctime>
#include <ratio>

using namespace std::chrono;

#define BVH_RESERVE_SIZE 2048
#define BVH_DEBUG 0

#include <CD_NamespaceHeader.H>

namespace ItoMerge {

  template <class T>
  Node<T>::Node(){
    m_data.reserve(BVH_RESERVE_SIZE);
  }

  template <class T>
  Node<T>::~Node(){
    m_data.resize(0);
  }

  template <class T>
  Node<T>::Node(std::shared_ptr<Node<T> >& a_parent){
    this->setParent(a_parent);
  }

  template <class T>
  inline bool Node<T>::isLeaf() const{
    return m_leaf;
  }

  template <class T>
  inline void Node<T>::setLeaf(const bool a_leaf){
    m_leaf = a_leaf;
  }

  template <class T>
  inline void Node<T>::setParent(std::shared_ptr<Node<T> >& a_parent){
    m_parent = a_parent;
  }

  template <class T>
  inline void Node<T>::setLeft(std::shared_ptr<Node<T> >& a_left){
    m_left = a_left;
  }

  template <class T>
  inline void Node<T>::setRight(std::shared_ptr<Node<T> >& a_right){
    m_right = a_right;
  }

  template <class T>
  inline void Node<T>::setData(const std::vector<T>& a_data, const Real a_mass){
    m_data = a_data;
    m_mass = a_mass;
  }

  template <class T>
  inline void Node<T>::setMass(const Real a_mass){
    m_mass = a_mass;
  }

  template <class T>
  inline std::shared_ptr<Node<T> >& Node<T>::getLeft(){
    return m_left;
  }

  template <class T>
  inline std::shared_ptr<Node<T> >& Node<T>::getRight(){
    return m_right;
  }


  template <class T>
  inline std::shared_ptr<Node<T> >& Node<T>::getParent(){
    return m_parent;
  }


  template <class T>
  inline std::vector<T>& Node<T>::getData(){
    return m_data;
  }

  template <class T>
  inline const std::vector<T>& Node<T>::getData() const{
    return m_data;
  }

  template <class T>
  inline bool Node<T>::canSplit() const {
    bool ret = false;
    if(m_data.size() > 0){
      if(m_data.size() > 1 || m_data[0].canSplit()){
	ret = true;
      }
    }

    //  return m_mass >= 2.0;
    return m_mass >= 2.0;
    return ret;
  }

  template <class T>
  inline Real Node<T>::mass() const {
    return m_mass;
  }

  template <class T>
  inline void Node<T>::split(const int a_dir){
#if BVH_DEBUG
    for (const auto& d : m_data){
      if(d.mass() < 1.0) {
	std::cout << d.mass() << std::endl;
	MayDay::Warning("Node<T>::split - how did mass be come <= 1??");
      }
    }
#endif

    // 0. If we made it in here, we KNOW that we can split.
#if BVH_DEBUG // Debugging flag to see if we're doing something bad. 
    if(m_mass <= 1.1){
      std::cout << m_mass << std::endl;
      MayDay::Warning("Node::split - m_mass <= 1.0, splitting should not be called");
    }
#endif

    const int data_size = m_data.size();
  
    // 1. Sort data. Use a lambda for comparing. God I love C++11.
    std::sort(std::begin(m_data), std::end(m_data), [a_dir](const T& p1, const T& p2){return p1[a_dir] < p2[a_dir];});

    // Init median on i = 0. mass_left and mass_right are the masses on each side of the median mass. 
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
    if(!m_left)  m_left  = std::make_shared<Node<T> >();
    if(!m_right) m_right = std::make_shared<Node<T> >();

    std::vector<T>& leftData  = m_left->getData();
    std::vector<T>& rightData = m_right->getData();

    leftData.assign(m_data.begin(), m_data.begin() + splitIndex);
    rightData.assign(m_data.begin() + splitIndex + 1, m_data.end());
  
    // If we can split up the median point into several ones, we do it. This would be the case if the median point
    // is a superparticle; in this case we divide the weights as best we can and then assign them to the left/right halves
    const T& splitMass = m_data[splitIndex];

#if 0 // Print direction and position
    std::cout << splitMass.pos()[0] << "\t" << splitMass.pos()[1] << "\t" << a_dir << std::endl;
#endif
    if(splitMass.canSplit()){
      const RealVect& p   = splitMass.pos();
      const Real& mass    = splitMass.mass();         // Always positive
      const Real& energy  = splitMass.energy();       // Total energy of split mass particle
      const Real massDiff = mass_right - mass_left;   // Not necessarily positive


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

	const long long N  = llround(massRem);//*(1.0 + 1.E-6));
	if(N > 0LL){ 
	  const long long NR = N/2;
	  const long long NL = N - NR;

	  massLeft  += (massRem/N)*NL;
	  massRight += (massRem/N)*NR;
	}

	if(massLeft > 0.0){
	  leftData.emplace_back(T(p, massLeft, energy));
	  mass_left  += massLeft;
	}

	if(massRight > 0.0){
	  rightData.emplace_back(T(p, massRight, energy));
	  mass_right += massRight;
	}
      }
    }
    else{
      //    We should be able to simplify this since it always goes into the left mass(?)
      if(mass_left <= mass_right){
	leftData.emplace_back(splitMass);
	mass_left += splitMass.mass();
      }
      else{
	rightData.emplace_back(splitMass);
	mass_right += splitMass.mass();
      }
    }

#if BVH_DEBUG // Make sure we don't break mass
    if(Abs(m_mass - (mass_left + mass_right)) > 1.E-6) MayDay::Abort("Node::split - broke mass");
    if(mass_left <= 0.1) {
      std::cout << m_mass << "\t" << mass_left << "\t" << mass_right << std::endl;
      MayDay::Abort("Node::split - got zero left mass");
    }
    if(mass_right <= 0.1) {
      std::cout << m_mass << "\t" << mass_left << "\t" << mass_right << std::endl;
      MayDay::Abort("Node::split - got zero right mass");
    }
#endif

#if BVH_DEBUG
    for (const auto& l : leftData){
      if(l.mass() <= 0.0){
	MayDay::Abort("Node::split - left data got a zero mass particle... :(");
      }
    }
    for (const auto& l : rightData){
      if(l.mass() <= 0.0){
	MayDay::Abort("Node::split - right data got a zero mass particle... :(");
      }
    }
#endif

    m_left->setMass(mass_left);
    m_right->setMass(mass_right);

#if 0 // This breaks. Don't know why. 
    left->setParent(this->shared_from_this());
    righ->setParent(this->shared_from_this());
#endif

    // Update who is leaf
    m_left->setLeaf(true);
    m_right->setLeaf(true);
    this->setLeaf(false);
  
  }

  template <class T>
  Tree<T>::Tree(){
    m_leaves.reserve(BVH_RESERVE_SIZE);
    m_new_leaves.reserve(BVH_RESERVE_SIZE);
  }

  template <class T>
  Tree<T>::Tree(std::vector<T>& a_data, const Real a_mass) : Tree() {
    this->define(a_data, a_mass);
  }

  template <class T>
  Tree<T>::~Tree(){
    m_leaves.resize(0);
  }

  template <class T>
  inline void Tree<T>::define(std::vector<T>& a_data, const Real a_mass){

    // Make the root node
    if(!m_root) m_root = std::make_shared<Node<T> >();

    m_root->setLeaf(true);
    m_root->setData(a_data, a_mass);

    m_leaves.resize(1);
    m_leaves[0] = m_root;
  }

  template <class T>
  inline std::vector<std::shared_ptr<Node<T> > >& Tree<T>::getLeaves(){
    return m_leaves;
  }

  template <class T>
  inline void Tree<T>::buildTree(const int a_firstDir, const int a_numLeaves){
  
#if BVH_DEBUG
    for (const auto& d : m_root->getData()){
      if(d.mass() < 1.0) MayDay::Abort("Node<T>::buildTree - how did mass become <= 0??");
    }
#endif

    bool keepSplitting = a_numLeaves > 1;

    // Reset leaves and root
    m_leaves.resize(a_numLeaves, nullptr);
    m_new_leaves.resize(a_numLeaves, nullptr);
    m_leaves[0] = m_root;

    bool canSplit;
    int leavesNeeded;
    int cur_size = 1;
    int splitDir = a_firstDir;
    int N;
  
    while(keepSplitting){
      canSplit = false;

      // Sort the leaves by their mass(?). We will split the ones with the largest masses.
      // std::sort(m_leaves.begin(), m_leaves.begin() + cur_size,
      // 	      [](const std::shared_ptr<Node<T> >& n1,
      // 		 const std::shared_ptr<Node<T> >& n2){
      // 		return n1->mass() < n2->mass();
      // 	      });

      N = 0;
      leavesNeeded = a_numLeaves - cur_size;
      for (int i = 0; i < cur_size; i++){
      
	if(m_leaves[i]->canSplit() && leavesNeeded > 0){ // Only split if we actually need more leaves. 
	  canSplit = true;

	  m_leaves[i]->split(splitDir);
	
	  m_new_leaves[N]   = m_leaves[i]->getLeft(); 
	  m_new_leaves[N+1] = m_leaves[i]->getRight();

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
      keepSplitting = N < a_numLeaves && canSplit;
    }
  }

}

#include <CD_NamespaceFooter.H>

#endif
