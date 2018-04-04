/*!
  @file   dcel_iterator.H
  @brief  Declaration of an iterators for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_ITERATOR_
#define _DCEL_ITERATOR_

#include "RealVect.H"

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
  edge_iterator(const dcel_poly* const a_poly);

  /*!
    @brief Iterator for half edges that end on the input vertex
  */
  edge_iterator(const dcel_vert* const a_vert);

  /*!
    @brief Get edge
  */
  virtual const dcel_edge* operator() ();

  /*!
    @brief Reset iterator
  */
  virtual void reset();

  /*!
    @brief Increment edge
  */
  virtual void operator ++();

  /*!
    @brief Check if ok or not
  */
  virtual bool ok();
  
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
  const dcel_edge* m_begin;

  /*!
    @brief Current edge
  */
  const dcel_edge* m_current;
};

#endif
