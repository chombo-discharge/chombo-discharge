/*!
  @file   kd_tree.H
  @brief  Declaration of a kd-tree space partitioning class for fast polygon search
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _KD_TREE_
#define _KD_TREE_

#include "kd_node.H"

/*!
  @brief Class for generating k-D binary search trees for fast polygon search
*/
template <typename T> class kd_tree {
public:

  kd_tree();

  /*!
    @brief Full constructor
  */
  kd_tree(const Vector<RefCountedPtr<T> >& a_data, const int a_max_depth, const int a_max_elements);

  /*!
    @brief Destructor
  */
  virtual ~kd_tree();

    /*!
    @brief Build kd-tree
  */
  virtual void build_tree(const Vector<RefCountedPtr<T> >& a_data, const int a_max_depth, const int a_max_elements);

  /*!
    @brief Get candidates from tree search
  */
  virtual Vector<RefCountedPtr<T> > get_candidates(const RealVect a_x0);

protected:

  /*!
    @brief Root node
  */
  RefCountedPtr<kd_node<T> > m_root;
  
  /*!
    @brief Nodes
  */
  Vector<RefCountedPtr<kd_node<T> > > m_nodes;

  /*!
    @brief Depth
  */
  int m_depth;

  /*!
    @brief Compute low corner for bounding box
  */
  virtual RealVect compute_bbox_lo(const Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Compute high corner for bounding box
  */
  virtual RealVect compute_bbox_hi(const Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Recursive search down the tree starting from a_parent
  */
  virtual void tree_search(Vector<RefCountedPtr<T> >&        a_data,
			   const RefCountedPtr<kd_node<T> >& a_parent,
			   const RealVect                    a_x0);

  /*!
    @brief Recursive search down the tree starting from a_parent
  */
  virtual void tree_search2(Vector<RefCountedPtr<T> >&        a_data,
			    Real&                             a_shortest,
			    const RefCountedPtr<kd_node<T> >& a_node,
			    const RealVect                    a_x0);
};

#include "kd_treeI.H"

#endif

