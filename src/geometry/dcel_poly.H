/*!
  @file   dcel_poly.H
  @brief  Declaration of a polygon class for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_POLY_
#define _DCEL_POLY_

#include <ritter_sphere.H>

#include "RealVect.H"
#include "RefCountedPtr.H"

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
  virtual void define(const RealVect a_normal, const RefCountedPtr<dcel_edge>& a_edge);

  /*!
    @brief Set edge
  */
  virtual void set_edge(const RefCountedPtr<dcel_edge>& a_edge);

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
  virtual void compute_normal(const bool a_outward_normal);

  /*!
    @brief Compute bounding box
  */
  virtual void compute_bbox();

  /*!
    @brief Get edge
  */
  virtual const RefCountedPtr<dcel_edge>& get_edge() const;

  /*!
    @brief Get edge
  */
  virtual RefCountedPtr<dcel_edge>& get_edge();

  /*!
    Get all positions
  */
  virtual Vector<RealVect> get_points();

  /*!
    @brief Get all vertices
  */
  virtual Vector<RefCountedPtr<dcel_vert> > get_vertices();

  /*!
    @brief Get all edges
  */
  virtual Vector<RefCountedPtr<dcel_edge> > get_edges();

  /*!
    @brief Get centroid
  */
  virtual RealVect get_centroid() const;

  /*!
    @brief Get the centroid
  */
  virtual RealVect get_coord() const;

  /*!
    @brief Get lower corner of bounding box
  */
  virtual RealVect get_bbox_lo() const;

    /*!
    @brief Get high corner of bounding box
  */
  virtual RealVect get_bbox_hi() const;
  
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
#if 0 // This is the actual signature
  virtual Real signed_distance(const RealVect a_x0) const;
#else // Just so that I can make things compile
  virtual Real signed_distance(const RealVect a_x0);
#endif
  
  
protected:
  
  /*!
    @brief One of the edges of the polygon
  */
  RefCountedPtr<dcel_edge> m_edge;

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

  /*!
    @brief Bounding box low corner
  */
  RealVect m_lo;

  /*!
    @brief Bounding box high corner
  */
  RealVect m_hi;

  /*!
    @brief Ritter sphere
  */
  ritter_sphere m_ritter;
};

#endif
