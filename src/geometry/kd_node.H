/*!
  @file   kd_node.H
  @brief  Declaration of a kd-node for space partioning
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _KD_NODE_
#define _KD_NODE_

#include "dcel_poly.H"

#include <RefCountedPtr.H>
#include <Vector.H>

/*!
  @brief Node class for kd_tree. This should probably be a template class. 
*/
template<typename T> class kd_node {
public:

  /*!
    @brief Weak constructor
  */
  kd_node();

  /*!
    @brief Full constructor
  */
  kd_node(const RefCountedPtr<kd_node<T> >& a_root);

  /*!
    @brief Destructor
  */
  virtual ~kd_node();

  /*!
    @brief Check if node is leaf
  */
  virtual bool is_leaf() const;

  /*!
    @brief Check if point intersects with bounding box
  */
  virtual bool inside(const RealVect a_x0) const;

  /*!
    @brief Set leaf or not
  */
  virtual void set_leaf(const bool a_leaf);

  /*!
    @brief Set root
  */
  virtual void set_root(const RefCountedPtr<kd_node<T> >& a_root);

  /*!
    @brief Set left
  */
  virtual void set_left(const RefCountedPtr<kd_node<T> >& a_left);

  /*!
    @brief Set right
  */
  virtual void set_right(const RefCountedPtr<kd_node<T> >& a_right);

  /*!
    @brief Set bounding box
  */
  virtual void set_bbox(const RealVect a_lo, const RealVect a_hi);

  /*!
    @brief Set data
  */
  virtual void set_data(const Vector<RefCountedPtr<T> >& a_data);

  /*!
    @brief Clear data
  */
  virtual void clear_data();

  /*!
    @brief Get root
  */
  virtual const RefCountedPtr<kd_node<T> >& get_root() const;

  /*!
    @brief Get left branch
  */
  virtual const RefCountedPtr<kd_node<T> >& get_left() const;

    /*!
    @brief Get left branch
  */
  virtual const RefCountedPtr<kd_node<T> >& get_right() const;

protected:

  /*!
    @brief Root node
  */
  RefCountedPtr<kd_node<T> > m_root;

  /*!
    @brief Left branch
  */
  RefCountedPtr<kd_node<T> > m_left;

  /*!
    @brief Right branch
  */
  RefCountedPtr<kd_node<T> > m_right;

  /*!
    @brief Set leaf or not
  */
  bool m_leaf;

  /*!
    @brief Low corner of bounding box
  */
  RealVect m_lo;

  /*!
    @brief High corner of bounding box
  */
  RealVect m_hi;
  
  /*!
    @brief Leaf data
  */
  Vector<RefCountedPtr<T> > m_data;
};

#include "kd_nodeI.H"

#endif

