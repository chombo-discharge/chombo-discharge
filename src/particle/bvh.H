/*!
  @file   bvh.H
  @brief  Declaration of a class for making a kd-tree hiearchy when generating superparticles
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _BVH_H_
#define _BVH_H_

#include <RealVect.H>
#include <memory>

/*!
  @brief Class the describes a single kd-tree node.
*/
template <class T>
class bvh_node {
public:

  bvh_node();
  bvh_node(std::shared_ptr<bvh_node<T> > a_parent);

  inline bool is_leaf() const;

  inline void set_leaf(const bool a_leaf);
  inline void set_parent(std::shared_ptr<bvh_node> a_parent);
  inline void set_left(std::shared_ptr<bvh_node> a_left);
  inline void set_right(std::shared_ptr<bvh_node> a_right);
  inline void set_data(const std::vector<T>& a_data);
  inline void sort(const int a_dir);
  inline void split_data();

  std::vector<T>& get_data();
  const std::vector<T>& get_data() const;

  const std::shared_ptr<T> get_parent() const;
  const std::shared_ptr<T> get_left() const;
  const std::shared_ptr<T> get_right() const;

protected:

  std::shared_ptr<T> m_parent;
  std::shared_ptr<T> m_left;
  std::shared_ptr<T> m_right;

  bool m_leaf;

  RealVect m_lo;
  RealVect m_hi;

  std::vector<T> m_data;
};

/*!
  @brief Class for making a kd-tree with approximately the same weight in each leaf
*/
template <class T>
class bvh_tree {
public:

  bvh_tree() = delete;
  
  bvh_tree(std::vector<T>& a_points);
  ~bvh_tree();

  void build_tree(const int a_numLeaves); // Desired number of leaves

  std::vector<std::shared_ptr<bvh_node<T> > >& get_nodes(); 
  std::vector<std::shared_ptr<bvh_node<T> > >& get_leaves(); 
  
protected:

  std::shared_ptr<bvh_node<T> > m_root; // Root node

  std::vector<std::shared_ptr<bvh_node<T> > > m_nodes;
  std::vector<std::shared_ptr<bvh_node<T> > > m_leaves;

  int m_depth;
};

#include "bvhI.H"

#endif
