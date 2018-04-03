/*!
  @file   dcel_edge.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_EDGE_
#define _DCEL_EDGE_

#include "RealVect.H"

class dcel_vert;
class dcel_poly;

/*!
  @brief Edge class for handling surface tesselations with DCEL meshes
*/
class dcel_edge {
public:

  /*!
    @brief Constructor
  */
  dcel_edge();

  /*!
    @brief Destructor
  */
  ~dcel_edge();

  /*!
    @brief Get vertex
  */
  virtual const dcel_vert* get_vert() const;

  /*!
    @brief Get other vertex
  */
  virtual const dcel_vert* get_other_vert() const;

  /*!
    @brief Get the other half edge
  */
  virtual const dcel_edge* get_pair() const;

  /*!
    @brief Get the previous half edge
  */
  virtual const dcel_edge* get_prev() const;

  /*!
    @brief Get the next half edge
  */
  virtual const dcel_edge* get_next() const;

  /*!
    @brief Get pseudo-normal
  */
  virtual RealVect get_normal() const;

  /*!
    @brief Get the signed distance between this edge and a point
  */
  virtual Real signed_distance(const RealVect a_x0) const;

protected:

  /*!
    @brief Edge pseudo-normal
  */
  RealVect m_normal;

  /*!
    @brief Vertex at end of edge
  */
  dcel_vert* m_vert; 

  /*!
    @brief Edge pair
  */
  dcel_edge* m_pair;

  /*!
    @brief Previous edge
  */
  dcel_edge* m_prev; 

  /*!
    @brief Next edge
  */
  dcel_edge* m_next;
};

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
    @brief Starting edge
  */
  const dcel_edge* m_begin;

  /*!
    @brief Current edge
  */
  const dcel_edge* m_current;
};

#endif
