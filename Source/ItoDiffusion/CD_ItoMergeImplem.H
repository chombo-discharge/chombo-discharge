/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_NamespaceHeader.H>

#define ITOMERGE_VECTOR_CAPACITY 2048

namespace ItoMerge {

  template <class T>
  Node<T>::Node(){
    m_parent = nullptr;
    m_left   = nullptr;
    m_right  = nullptr;
    
    m_data.reserve(ITOMERGE_VECTOR_CAPACITY);
  }

  template <class T>
  Node<T>::~Node(){
    m_data.resize(0);
  }

  template <class T>
  Node<T>::Node(std::shared_ptr<Node<T> >& a_parent) : Node<T>(){
    this->setParent(a_parent);
  }

  template <class T>
  inline
  bool Node<T>::isLeafNode() const {
    return m_left == nullptr && m_right == nullptr;
  }

  template <class T>
  inline
  bool Node<T>::isRegularNode() const {
    return !(this->isLeafNode());
  }

  template <class T>
  inline
  bool Node<T>::isRootNode() const {
    return m_parent == nullptr;
  }    

  template <class T>
  inline
  void Node<T>::setParent(std::shared_ptr<Node<T> >& a_parent){
    m_parent = a_parent;
  }

  template <class T>
  inline
  void Node<T>::setLeft(std::shared_ptr<Node<T> >& a_left){
    m_left = a_left;
  }

  template <class T>
  inline
  void Node<T>::setRight(std::shared_ptr<Node<T> >& a_right){
    m_right = a_right;
  }

  template <class T>
  inline
  void Node<T>::setData(const std::vector<T>& a_data, const Real a_mass){
    m_data = a_data;
    m_mass = a_mass;
  }

  template <class T>
  inline
  void Node<T>::setMass(const Real a_mass){
    m_mass = a_mass;
  }

  template <class T>
  inline
  std::shared_ptr<Node<T> >& Node<T>::getLeft(){
    return m_left;
  }

  template <class T>
  inline
  std::shared_ptr<Node<T> >& Node<T>::getRight(){
    return m_right;
  }

  template <class T>
  inline
  std::shared_ptr<Node<T> >& Node<T>::getParent(){
    return m_parent;
  }

  template <class T>
  inline
  std::vector<T>& Node<T>::getData(){
    return m_data;
  }

  template <class T>
  inline
  const std::vector<T>& Node<T>::getData() const{
    return m_data;
  }

  template <class T>
  inline
  bool Node<T>::canSplit() const {
    bool ret = false;
    if(m_data.size() > 0){
      if(m_data.size() > 1 || m_data[0].canSplit()){
	ret = true;
      }
    }

    return m_mass >= 2.0;
    return ret;
  }

  template <class T>
  inline
  Real Node<T>::mass() const {
    return m_mass;
  }

  template <class T>
  inline
  void Node<T>::split(const int a_dir){
    CH_TIME("ItoMerge::Node<T>::split");

    // Some basic assertions -- the loop should evaluate to an empty loop when DEBUG=TRUE.
    CH_assert(m_data.size() > 0  );
    CH_assert(m_mass        > 1.0);
    
    for (const auto& n : m_data){
      CH_assert(n.mass() >= 1.0);
    }    

    // TLDR: This function splits the primitives in the node into two subvolumes. A splitting plane is made along the input coordinate and it is positioned
    //       such that the total mass on the left/right side of the splitting plane differ by at most one computational particle. The guts of the code mainly
    //       consists of tests for how to divide the "median" particle -- there are various cases that we need to cover here (such as the median particle being
    //       a physical particle or a superparticle). 
  
    // 1. Sort data along the input coordinate. We need this because we need to have something sensible to iterate through
    //    when we place the splitting plane. 
    std::sort(std::begin(m_data), std::end(m_data), [a_dir](const T& p1, const T& p2){return p1[a_dir] < p2[a_dir];});

    // 2. Figure out where to place the splitting plane -- this is not the physical coordinate but index in m_data (the list of primitives). We are
    //    actually looking for the particle on some median where the masses between the two half-planes can be made approximately the same. We begin
    //    by setting the median on i = 0 in the primitives vector (m_data). We then go through the primitives and check if adding the median particle
    //    ends up placing MORE mass in the left node than the right node. Once that happens we have the inflection point where we need to split. The median
    //    particle will be assigned later because it may be a superparticle and we can equilibriate the masses even further. 
    //
    //
    // In all of the below. massLeft and massRight are the masses on each side of the median particle. As remarked above, the median particle could be split
    // later. 
    int  medianParticleIndex = 0;
    Real massLeft            = 0.0;                                         // Zero mass, as it should be. 
    Real massRight           = m_mass - m_data[medianParticleIndex].mass(); // So, entire mass minus the median particle mass. 

    for (int i = 1; i < m_data.size(); i++){
      const Real& medianParticleMass = m_data[medianParticleIndex].mass(); // This is the current median particle mass. 

      // If adding the median particle to the mass in the left node still yields a lower mass than the right now, we just shift the median particle index. So, the previous
      // "median particle" goes to the left node but we then have a new median particle. 
      if(massLeft + medianParticleMass < massRight){
	massLeft += medianParticleMass;
	
	medianParticleIndex = i;

	massRight = m_mass - massLeft - m_data[medianParticleIndex].mass(); // Right mass = Total mass - left mass - median particle mass
      }
      else{
	break;
      }
    }

    // 2. Make the two leaves and push data into them. Currently, we know that all particles to the left of the "median" particle go into the left half-plane
    //    and all particles to the right of the "median" particle go into the right half-plane. We do not yet know how to place the median particle -- it could
    //    go into either the left node or the right node, or it could be split into and placed into either node. 
    m_left  = std::make_shared<Node<T> >();
    m_right = std::make_shared<Node<T> >();

    std::vector<T>& leftData  = m_left ->getData();
    std::vector<T>& rightData = m_right->getData();

    leftData. assign(m_data.begin()                          , m_data.begin() + medianParticleIndex);
    rightData.assign(m_data.begin() + medianParticleIndex + 1, m_data.end()                        );
  
    // 3. We next need to assign the median particle. If we can split up the median point into several ones, we do it. This would be the case if the median point
    //    is a superparticle; in this case we divide the weights as best we can and then assign them to the left/right halves. If we can't then we just assign the
    //    particle to whichever left/right child node has the lowest mass. 
    const T& medianParticle = m_data[medianParticleIndex];

    if(medianParticle.canSplit()){
      const RealVect& medianPos      = medianParticle.pos   ();
      const Real&     medianMass     = medianParticle.mass  (); // Always positive
      const Real&     medianEnergy   = medianParticle.energy(); // Note: This is average energy so it remain the same for subparticles. 
      const Real      massDiff = massRight - massLeft;          // Not necessarily positive.

      // This hook is for when all mass goes into one part. 
      if(medianMass < std::abs(massDiff)){
	if(massDiff > 0){ // Right mass is bigger than left so put entire particle in left mass. 
	  massLeft += medianMass; 
	  leftData.emplace_back(T(medianPos, medianMass, medianEnergy));
	}
	else{ // Left mass is bigger than right so put everything into the right node. 
	  massRight += medianMass; // Goes into right
	  rightData.emplace_back(T(medianPos, medianMass, medianEnergy));
	}
      }
      else {
	// Distribute mass as best we can. This might seem weird because we can have massLeft > massRight in which
	// case massDiff < 0. This still works since we are simply subtracting the difference from the left leaf first and
	// then divide a larger total mass between the two leaves. Essentially, we are subtracting the part of the median
	// particle that gives massLeft > massRight from the left particle first, and then we reassign a larger total mass.
	Real tmpMassLeft  = massDiff;
	Real tmpMassRight = 0.0;
	Real tmpMassRem   = medianMass - massDiff;

	const long long N  = llround(tmpMassRem);
	if(N > 0LL){ 
	  const long long NR = N/2;
	  const long long NL = N - NR;

	  tmpMassLeft  += (tmpMassRem/N)*NL;
	  tmpMassRight += (tmpMassRem/N)*NR;
	}

	if(tmpMassLeft > 0.0){
	  leftData.emplace_back(T(medianPos, tmpMassLeft, medianEnergy));
	  massLeft  += tmpMassLeft;
	}

	if(tmpMassRight > 0.0){
	  rightData.emplace_back(T(medianPos, tmpMassRight, medianEnergy));
	  massRight += tmpMassRight;
	}
      }
    }
    else{ 
      // Not a computational particle -- put the particle in the child node with the lowest mass. 
      if(massLeft <= massRight){
	leftData.emplace_back(medianParticle);
	massLeft += medianParticle.mass();
      }
      else{
	rightData.emplace_back(medianParticle);
	massRight += medianParticle.mass();
      }
    }

    // Update the tree structure. 
    auto thisPtr = this->shared_from_this();
    
    m_left ->setMass(massLeft );
    m_right->setMass(massRight);
    
    m_left ->setParent(thisPtr);
    m_right->setParent(thisPtr);    

    m_data.resize(0);        

    // Some assertions in case this code every breaks.
    CH_assert((std::abs(m_mass - (massLeft + massRight)) <= 1.E-12                         ));
    CH_assert(massLeft                                   >= std::numeric_limits<Real>::min());
    CH_assert(massRight                                  >= std::numeric_limits<Real>::min());

    // These should be evaluated to empty loopps by the compiler if DEBUG=FALSE.
    for (const auto& l : leftData){
      CH_assert(l.mass() >= std::numeric_limits<Real>::min());
    }

    for (const auto& r : rightData){
      CH_assert(r.mass() >= std::numeric_limits<Real>::min());
    }    
  }

  template <class T>
  Tree<T>::Tree(){
    m_isDefined = false;
  }

  template <class T>
  Tree<T>::Tree(const std::vector<T>& a_data, const Real a_mass) : Tree() {
    this->define(a_data, a_mass);
  }

  template <class T>
  Tree<T>::~Tree(){
    m_leaves.resize(0);
  }

  template <class T>
  inline
  void Tree<T>::define(const std::vector<T>& a_data, const Real a_mass){

    // Make the root node if we don't have it. 
    if(!m_root) {
      m_root = std::make_shared<Node<T> >();
    }

    m_root->setData(a_data, a_mass);

    m_leaves.resize(1);
    m_leaves[0] = m_root;

    m_isDefined = true;
  }

  template <class T>
  inline
  std::vector<std::shared_ptr<Node<T> > >& Tree<T>::getLeaves(){
    CH_assert(m_isDefined);
    
    return m_leaves;
  }

  template <class T>
  inline  
  const std::vector<std::shared_ptr<Node<T> > >& Tree<T>::getLeaves() const {
    CH_assert(m_isDefined);
    
    return m_leaves;
  }  

  template <class T>
  inline
  void Tree<T>::buildTree(const int a_firstDir, const int a_numTargetLeaves){
    CH_TIME("ItoMerge::Tree<T>::buildTree");

    CH_assert(m_isDefined                    );
    CH_assert(a_numTargetLeaves >= 1         );
    CH_assert(a_firstDir        >= 0         );
    CH_assert(a_firstDir        <= SpaceDim-1);    
  
    // Should evaluate to an empty loop and be purged by the compiler if DEBUG=FALSE
    for (const auto& d : m_root->getData()){
      CH_assert(d.mass() >= 1.0);
    }

    // TLDR: This function splits the primitives recursively until we either have a_numTargetLeaves leaf nodes OR none of the leaf nodes can be
    //       partitioned further. 

    std::vector<std::shared_ptr<Node<T> > > newLeaves;

    // Reset leaves and root. I'm just resizing here so I get the capacity
    // and avoid vector resizing during the splitting process.
    m_leaves. resize(a_numTargetLeaves, nullptr);
    newLeaves.resize(a_numTargetLeaves, nullptr);
    m_leaves[0] = m_root;

    // We always begin with one leaf node. 
    int numCurrentLeaves = 1;

    int  splitDir      = a_firstDir;    
    bool keepSplitting = (a_numTargetLeaves > 1);
    
    while(keepSplitting){

      // Flag which determines if we can split any of the nodes further. This will remain false if none of the
      // leaf nodes could be split OR if we didn't need more leaf nodes (i.e. the tree has finished building). 
      bool canSplit = false; 

#if 0 // I'm not sure this matters so I'm commenting it out. R.M. Nov. 2021. 
      // Sort the leaves by their mass(?). We will split the ones with the largest masses.      
      std::sort(m_leaves.begin(), m_leaves.begin() + numCurrentLeaves,
		[](const std::shared_ptr<Node<T> >& n1,
		   const std::shared_ptr<Node<T> >& n2){
		  return n1->mass() < n2->mass();
		});
#endif

      // These are the number of leaves that we need in order to reach a_numTargetLeaves. 
      int leavesNeeded = a_numTargetLeaves - numCurrentLeaves;

      int N = 0;      
      for (int i = 0; i < numCurrentLeaves; i++){

	// If we needed a new leaf AND we could split this node, split it and go further. This will give us two leaf nodes
	// in the "new" list of leaf nodes.
	if(m_leaves[i]->canSplit() && leavesNeeded > 0){ 
	  m_leaves[i]->split(splitDir);
	
	  newLeaves[N]   = m_leaves[i]->getLeft (); 
	  newLeaves[N+1] = m_leaves[i]->getRight();

	  N += 2;

	  leavesNeeded--;

	  canSplit = true;
	}
	else{ // We did not split the leaf node so just append the current leaf to the list of "new" leaves.
	  newLeaves[N] = m_leaves[i];
	  
	  N++;
	}
      }

      // Assign new leaf nodes. 
      m_leaves.assign(newLeaves.begin(), newLeaves.begin() + N);

      // Update the current number of leaf nodes. 
      numCurrentLeaves = N;

      // Did the tree finish building or are all the nodes unsplittable leaf nodes? If so, exit. 
      keepSplitting = (N < a_numTargetLeaves) && canSplit;

      // Update the splitting coordinate -- the next leaves are split along a different coordinate.
      splitDir = (splitDir + 1) % SpaceDim;       
    }
  }

  /*!
    @brief Default partitioning functino for splitting a node. 
  */
  template <class T>
  NodePartitioner<T> NodePartitionEqualMass = [](std::shared_ptr<Node<T> >& a_node, const int a_dir) -> void {
    return a_node->split(a_dir);
  };
  
}

#include <CD_NamespaceFooter.H>

#endif
