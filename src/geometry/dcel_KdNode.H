/*!
  @file   dcel_KdNode.H
  @brief  Declaration of a kd-node for space partitioning
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _KD_NODE_
#define _KD_NODE_

#include <vector>
#include <memory>

namespace dcel {
  
  /*!
    @brief Template class for kd-trees
  */
  template<typename T>
  class kd_node {
  public:

    /*!
      @brief Weak constructor
    */
    kd_node();

    /*!
      @brief Full constructor
    */
    kd_node(const std::shared_ptr<kd_node<T> >& a_parent);

    /*!
      @brief Destructor
    */
    ~kd_node();

    /*!
      @brief Check if node is leaf
    */
    bool is_leaf() const;

    /*!
      @brief Check if point intersects with bounding box
    */
    bool inside(const RealVect a_x0) const;

    /*!
      @brief Get node depth
    */
    int get_depth();

    /*!
      @brief Get number of data elements
    */
    int get_num_elements();

    /*!
      @brief Compute distance. This is the bounding box distance if node, and data distance if leaf
    */
    Real dist(const RealVect a_x0) const;

    /*!
      @brief Compute distance to bounding box
    */
    Real aabb_dist(const RealVect a_x0) const;

    /*!
      @brief Compute distance to data
    */
    Real data_dist(const RealVect a_x0) const;

    /*!
      @brief Compute bounding box
    */
    void compute_bbox();

    /*!
      @brief Set leaf or not
    */
    void set_leaf(const bool a_leaf);

    /*!
      @brief Set parent
    */
    void set_parent(const std::shared_ptr<kd_node<T> >& a_parent);

    /*!
      @brief Set left
    */
    void set_left(const std::shared_ptr<kd_node<T> >& a_left);

    /*!
      @brief Set right
    */
    void set_right(const std::shared_ptr<kd_node<T> >& a_right);

    /*!
      @brief Set data
    */
    void set_data(const std::vector<std::shared_ptr<T> >& a_data);

    /*!
      @brief Clear data
    */
    void clear_data();

    /*!
      @brief Split data
    */
    void split_data();

    /*!
      @brief Return data
    */
    std::vector<std::shared_ptr<T> >& get_data();

    /*!
      @brief Get root
    */
    const std::shared_ptr<kd_node<T> >& get_parent() const;

    /*!
      @brief Get left branch
    */
    const std::shared_ptr<kd_node<T> >& get_left() const;

    /*!
      @brief Get left branch
    */
    const std::shared_ptr<kd_node<T> >& get_right() const;

  protected:

    /*!
      @brief Root node
    */
    std::shared_ptr<kd_node<T> > m_parent;

    /*!
      @brief Left branch
    */
    std::shared_ptr<kd_node<T> > m_left;

    /*!
      @brief Right branch
    */
    std::shared_ptr<kd_node<T> > m_right;

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
    std::vector<std::shared_ptr<T> > m_data;
  };
}

#include "dcel_KdNodeI.H"

#endif

