/*!
  @file   dcel_iterator.H
  @brief  Declaration of an iterators for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_ITERATOR_
#define _DCEL_ITERATOR_

#include "RealVect.H"

#include <memory>

class dcel_vert;
class dcel_edge;
class dcel_poly;

/*!
  @brief Class that iterates over edges
*/
class edge_iterator {
public:

  /*!
    @brief Weak constructor
  */
  edge_iterator();
  
  /*!
    @brief Iterator for all half edges of a polygon
  */
  edge_iterator(dcel_poly& a_poly);

  /*!
    @brief Iterator for outgoing half edges from the vertex
  */
  edge_iterator(dcel_vert& a_vert);

  /*!
    @brief Get edge
  */
  
  std::shared_ptr<dcel_edge>& operator() ();

  /*!
    @brief Reset iterator
  */
  
  void reset();

  /*!
    @brief Increment edge
  */
  
  void operator ++();

  /*!
    @brief Check if ok or not
  */
  
  bool ok();
  
protected:

  /*!
    @brief Flag for encapsulating hwo we iterator
  */
  bool m_polymode;

 /*!
    @brief Check if this is the first iteration
  */
  bool m_full_loop;
  
  /*!
    @brief Starting edge
  */
  std::shared_ptr<dcel_edge> m_begin;

  /*!
    @brief Current edge
  */
  std::shared_ptr<dcel_edge> m_current;
};

#endif
