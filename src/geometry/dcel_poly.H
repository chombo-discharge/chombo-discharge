/*!
  @file   dcel_poly.cpp
  @brief  Declaration of a polygon class for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_POLY_
#define _DCEL_POLY_

#include "RealVect.H" 

class dcel_vert;
class dcel_edge;

/*!
  @brief Abstract class for declaring polygons with DCEL surface tesselations
*/
class dcel_poly {
public:
  
  /*!
    @brief Constructor
  */
  dcel_poly();

  /*!
    @brief Destructor
  */
  ~dcel_poly();

  /*!
    @brief Define edge
  */
  virtual void define(const RealVect a_normal, const dcel_edge* const a_edge);

  /*!
    @brief Set edge
  */
  virtual void set_edge(const dcel_edge* const a_edge);

    /*!
    @brief Set normal vector
  */
  virtual void set_normal(const RealVect a_normal);

  /*!
    @brief Compute polygon area
  */
  virtual void compute_area();

  /*!
    @brief Normalize normal vector
  */
  virtual void normalize();

  /*!
    @brief Compute centroid. 
    @note This currently computes the point centroid and not the mass centroid. This probably doesn't matter for performance. 
  */
  virtual void compute_centroid();

  /*!
    @brief Compute normal vector from vertices
  */
  virtual void compute_normal();

  /*!
    @brief Get edge
  */
  virtual const dcel_edge* get_edge() const;

  /*!
    @brief Get all vertices
  */
  virtual Vector<const dcel_vert*> get_vertices() const;

  /*!
    @brief Get all edges
  */
  virtual Vector<const dcel_edge*> get_edges() const;

  /*!
    @brief Get centroid
  */
  virtual RealVect get_centroid() const;
  
  /*!
    @brief Get face normal vector
  */
  virtual RealVect get_normal() const;

  /*!
    @brief Get polygon area
  */
  virtual Real get_area() const;

  /*!
    @brief Signed distance. 
  */
  virtual Real signed_distance(const RealVect a_x0) const;
  
protected:

  /*!
    @brief One of the edges of the polygon
  */
  const dcel_edge* m_edge;

  /*!
    @brief Face normal
  */
  RealVect m_normal;

  /*!
    @brief Polygon centroid
  */
  RealVect m_centroid;

  /*!
    @brief Polygon area
  */
  Real m_area;

  
};

#endif
