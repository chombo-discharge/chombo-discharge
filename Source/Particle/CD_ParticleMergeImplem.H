/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParticleMergeImplem.H
  @brief  Implementation of CD_ParticleMerge.H
  @author Robert Marskar
  @date   May 2020
*/

#ifndef CD_ParticleMergeImplem_H
#define CD_ParticleMergeImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_ParticleMerge.H>
#include <CD_ParticlePartitioners.H>
#include <CD_NamespaceHeader.H>

#define ParticleMergeVectorCapacity 2048

namespace ParticleMerge {

  template <class T>
  Node<T>::Node()
  {
    m_parent = nullptr;
    m_left   = nullptr;
    m_right  = nullptr;

    m_data.reserve(ParticleMergeVectorCapacity);
  }

  template <class T>
  Node<T>::~Node()
  {
    m_data.resize(0);
  }

  template <class T>
  Node<T>::Node(std::shared_ptr<Node<T>>& a_parent) : Node<T>()
  {
    this->setParent(a_parent);
  }

  template <class T>
  inline bool
  Node<T>::isLeafNode() const
  {
    return m_left == nullptr && m_right == nullptr;
  }

  template <class T>
  inline bool
  Node<T>::isRegularNode() const
  {
    return !(this->isLeafNode());
  }

  template <class T>
  inline bool
  Node<T>::isRootNode() const
  {
    return m_parent == nullptr;
  }

  template <class T>
  inline void
  Node<T>::setParent(std::shared_ptr<Node<T>>& a_parent)
  {
    m_parent = a_parent;
  }

  template <class T>
  inline void
  Node<T>::setLeft(std::shared_ptr<Node<T>>& a_left)
  {
    m_left = a_left;
  }

  template <class T>
  inline void
  Node<T>::setRight(std::shared_ptr<Node<T>>& a_right)
  {
    m_right = a_right;
  }

  template <class T>
  inline void
  Node<T>::setData(const std::vector<T>& a_data, const Real a_mass)
  {
    m_data = a_data;
    m_mass = a_mass;
  }

  template <class T>
  inline void
  Node<T>::setMass(const Real a_mass)
  {
    m_mass = a_mass;
  }

  template <class T>
  inline std::shared_ptr<Node<T>>&
  Node<T>::getLeft()
  {
    return m_left;
  }

  template <class T>
  inline std::shared_ptr<Node<T>>&
  Node<T>::getRight()
  {
    return m_right;
  }

  template <class T>
  inline std::shared_ptr<Node<T>>&
  Node<T>::getParent()
  {
    return m_parent;
  }

  template <class T>
  inline std::vector<T>&
  Node<T>::getData()
  {
    return m_data;
  }

  template <class T>
  inline const std::vector<T>&
  Node<T>::getData() const
  {
    return m_data;
  }

  template <class T>
  inline bool
  Node<T>::canSplit() const
  {
    bool ret = false;
    if (m_data.size() > 0) {
      if (m_data.size() > 1 || m_data[0].canSplit()) {
        ret = true;
      }
    }

    return m_mass >= 2.0;
    return ret;
  }

  template <class T>
  inline Real
  Node<T>::mass() const
  {
    return m_mass;
  }

  template <class T>
  Tree<T>::Tree()
  {
    m_isDefined = false;
  }

  template <class T>
  Tree<T>::Tree(const std::vector<T>& a_data) : Tree()
  {
    this->define(a_data);
  }

  template <class T>
  Tree<T>::~Tree()
  {
    m_leaves.resize(0);
  }

  template <class T>
  inline void
  Tree<T>::define(const std::vector<T>& a_data)
  {

    // Make the root node if we don't have it.
    if (!m_root) {
      m_root = std::make_shared<Node<T>>();
    }

    // Compute total mass.
    Real mass = 0.0;
    for (const auto& d : a_data) {
      mass += d.mass();
    }

    CH_assert(mass >= 1.0);

    m_root->setData(a_data, mass);

    m_leaves.resize(1);
    m_leaves[0] = m_root;

    m_isDefined = true;
  }

  template <class T>
  inline std::vector<std::shared_ptr<Node<T>>>&
  Tree<T>::getLeaves()
  {
    CH_assert(m_isDefined);

    return m_leaves;
  }

  template <class T>
  inline const std::vector<std::shared_ptr<Node<T>>>&
  Tree<T>::getLeaves() const
  {
    CH_assert(m_isDefined);

    return m_leaves;
  }

  template <class T>
  inline void
  Tree<T>::buildTree(const int a_firstDir, const int a_numTargetLeaves, const NodePartitioner<T>& a_partitioner)
  {
    CH_TIME("ParticleMerge::Tree<T>::buildTree");

    CH_assert(m_isDefined);
    CH_assert(a_numTargetLeaves >= 1);
    CH_assert(a_firstDir >= 0);
    CH_assert(a_firstDir <= SpaceDim - 1);

    // Should evaluate to an empty loop and be purged by the compiler if DEBUG=FALSE
    for (const auto& d : m_root->getData()) {
      CH_assert(d.mass() >= 1.0);
    }

    // TLDR: This function splits the primitives recursively until we either have a_numTargetLeaves leaf nodes OR none of the leaf nodes can be
    //       partitioned further.

    std::vector<std::shared_ptr<Node<T>>> newLeaves;

    // Reset leaves and root. I'm just resizing here so I get the capacity
    // and avoid vector resizing during the splitting process.
    m_leaves.resize(a_numTargetLeaves, nullptr);
    newLeaves.resize(a_numTargetLeaves, nullptr);
    m_leaves[0] = m_root;

    // We always begin with one leaf node.
    int numCurrentLeaves = 1;

    int  splitDir      = a_firstDir;
    bool keepSplitting = (a_numTargetLeaves > 1);

    while (keepSplitting) {

      // Flag which determines if we can split any of the nodes further. This will remain false if none of the
      // leaf nodes could be split OR if we didn't need more leaf nodes (i.e. the tree has finished building).
      bool canSplit = false;

#if 0 // I'm not sure this matters so I'm commenting it out. R.M. Nov. 2021. \
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
      for (int i = 0; i < numCurrentLeaves; i++) {

        // If we needed a new leaf AND we could split this node, split it and go further. This will give us two leaf nodes
        // in the "new" list of leaf nodes.
        if (m_leaves[i]->canSplit() && leavesNeeded > 0) {
          a_partitioner(m_leaves[i], splitDir);

          newLeaves[N]     = m_leaves[i]->getLeft();
          newLeaves[N + 1] = m_leaves[i]->getRight();

          N += 2;

          leavesNeeded--;

          canSplit = true;
        }
        else { // We did not split the leaf node so just append the current leaf to the list of "new" leaves.
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
} // namespace ParticleMerge

#include <CD_NamespaceFooter.H>

#endif
