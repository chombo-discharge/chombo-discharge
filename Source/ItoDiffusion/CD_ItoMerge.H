/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ItoMerge.H
  @brief  Basic merging utility for ItoSolver
  @author Robert Marskar
*/

#ifndef CD_ItoMerge_H
#define CD_ItoMerge_H

// Chombo includes
#include <RealVect.H>
#include <memory>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Namespace for merging utilities for ItoSolver
*/
namespace ItoMerge {
  /*!
    @brief Class the describes a single kd-tree node.
  */
  template <class T>
  class Node {// : public std::enable_shared_from_this<Node<T> > {
  public:

    Node();
    Node(std::shared_ptr<Node<T> >& a_parent);
    ~Node();

    inline bool canSplit() const;
    inline bool isLeaf() const;

    inline void setLeaf(const bool a_leaf);
    inline void setParent(std::shared_ptr<Node<T> >& a_parent);
    inline void setLeft(std::shared_ptr<Node<T> >& a_left);
    inline void setRight(std::shared_ptr<Node<T> >& a_right);
    inline void setData(const std::vector<T>& a_data, const Real a_mass);
    inline void setMass(const Real a_mass);
    inline void sort(const int a_dir);
    inline void split(const int a_dir);

    inline std::vector<T>& getData();
    inline const std::vector<T>& getData() const;

    inline Real mass() const;

    inline std::shared_ptr<Node<T> >& getParent();
    inline std::shared_ptr<Node<T> >& getLeft();
    inline std::shared_ptr<Node<T> >& getRight();

    //protected:

    std::shared_ptr<Node<T> > m_parent;
    std::shared_ptr<Node<T> > m_left;
    std::shared_ptr<Node<T> > m_right;

    bool m_leaf;

    RealVect m_lo;
    RealVect m_hi;

    Real m_mass;
    std::vector<T> m_data;
  };

  /*!
    @brief Class for making a kd-tree with approximately the same weight in each leaf
  */
  template <class T>
  class Tree {
  public:

    Tree();
  
    Tree(std::vector<T>& a_points, const Real a_mass);
    ~Tree();

    inline void buildTree(const int a_firstDir, const int a_numLeaves); // Desired number of leaves
    inline std::vector<std::shared_ptr<Node<T> > >& getLeaves();
    inline void define(std::vector<T>& a_points, const Real a_mass);
  
  protected:

    std::shared_ptr<Node<T> > m_root; // Root node

    std::vector<std::shared_ptr<Node<T> > > m_leaves;
    std::vector<std::shared_ptr<Node<T> > > m_new_leaves;
  };
}

#include <CD_NamespaceFooter.H>

#include <CD_ItoMergeImplem.H>

#endif
