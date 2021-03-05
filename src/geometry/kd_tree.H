/*!
  @file   kd_tree.H
  @brief  Declaration of a kd-tree space partitioning class for fast polygon search
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _KD_TREE_
#define _KD_TREE_

#include "kd_node.H"

#include <memory>
#include <vector>

/*!
  @brief Class for generating k-D binary search trees for fast polygon search
*/
template <typename T>
class kd_tree {
public:

  kd_tree();

  /*!
    @brief Full constructor
  */
  kd_tree(const std::vector<std::shared_ptr<T> >& a_data, const int a_max_depth, const int a_max_elements);

  /*!
    @brief Destructor
  */
  ~kd_tree();

  /*!
    @brief Build kd-tree
  */
  void build_tree(const std::vector<std::shared_ptr<T> >& a_data, const int a_max_depth, const int a_max_elements);

  /*!
    @brief Find closest data 
  */
  std::vector<std::shared_ptr<T> > find_closest(const RealVect a_x0);

  bool isDefined(){return m_isDefined;}
  int get_depth(){return m_depth;}

protected:

  bool m_isDefined;

  /*!
    @brief Root node
  */
  std::shared_ptr<kd_node<T> > m_root;

  /*!
    @brief One of the leaf nodes
  */
  std::shared_ptr<kd_node<T> > m_leaf;
  
  /*!
    @brief Nodes
  */
  std::vector<std::shared_ptr<kd_node<T> > > m_nodes;

  /*!
    @brief Leaves
  */
  std::vector<std::shared_ptr<kd_node<T> > > m_leaves;

  /*!
    @brief Depth
  */
  int m_depth;

  /*!
    @brief Tree search
  */
  
  void tree_search(std::vector<std::shared_ptr<kd_node<T> > >& a_candidates,
		   std::vector<Real>&                          a_distances,
		   const RealVect                              a_x0);

  /*!
    @brief Insertion sort
  */
  
  void insertion_sort(std::vector<std::shared_ptr<kd_node<T> > >& a_queue,
		      std::vector<Real>&                          a_distances);

  /*!
    @brief Trim std::vector
  */
  
  void trim(std::vector<std::shared_ptr<kd_node<T> > >& a_queue, std::vector<Real>& a_distances, const int max_elem);
		    
};

#include "kd_treeI.H"

#endif

