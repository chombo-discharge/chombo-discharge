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
    @brief Define function
  */
  virtual void define(const dcel_vert* const a_vert,
		      const dcel_edge* const a_pair,
		      const dcel_edge* const a_next,
		      const dcel_edge* const a_prev,
		      const RealVect         a_normal = RealVect::Zero);

  /*!
    @brief Set polygon
  */
  virtual void set_poly(const dcel_poly* const a_poly);

  /*!
    @brief Set normal vector
  */
  virtual void set_normal(const RealVect a_normal);
  
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
    @brief Get polygon
  */
  virtual const dcel_poly* get_poly() const;
  
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
  const dcel_vert* m_vert; 

  /*!
    @brief Edge pair
  */
  const dcel_edge* m_pair;

  /*!
    @brief Previous edge
  */
  const dcel_edge* m_prev; 

  /*!
    @brief Next edge
  */
  const dcel_edge* m_next;

  /*!
    @brief Reference to polygon
  */
  const dcel_poly* m_poly;
};

#endif
