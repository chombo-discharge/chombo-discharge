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

#include <memory>
#include <vector>

namespace dcel {

  class vertex;
  class edge;

  /*!
    @brief Abstract class for declaring polygons with DCEL surface tesselations
  */
  class polygon {
  public:
  
    /*!
      @brief Constructor
    */
    polygon();

    /*!
      @brief Destructor
    */
    ~polygon();

    /*!
      @brief Define edge
    */
  
    void define(const RealVect& a_normal, const std::shared_ptr<edge>& a_edge);

    /*!
      @brief Set edge
    */
  
    void set_edge(const std::shared_ptr<edge>& a_edge);

    /*!
      @brief Set normal vector
    */
  
    void set_normal(const RealVect& a_normal);

    /*!
      @brief Compute polygon area
    */
  
    void compute_area();

    /*!
      @brief Normalize normal vector
    */
  
    void normalize();

    /*!
      @brief Compute centroid. 
      @note This currently computes the point centroid and not the mass centroid. This probably doesn't matter for performance. 
    */
  
    void compute_centroid();

    /*!
      @brief Compute normal vector from vertices
    */
  
    void compute_normal(const bool a_outward_normal);

    /*!
      @brief Compute bounding box
    */
  
    void compute_bbox();

    /*!
      @brief Get edge
    */
  
    const std::shared_ptr<edge>& get_edge() const;

    /*!
      @brief Get edge
    */
  
    std::shared_ptr<edge>& get_edge();

    /*!
      Get all positions
    */
  
    std::vector<RealVect> get_points();

    /*!
      @brief Get all vertices
    */
  
    std::vector<std::shared_ptr<vertex> > get_vertices();

    /*!
      @brief Get all edges
    */
  
    std::vector<std::shared_ptr<edge> > get_edges();

    /*!
      @brief Get centroid
    */
  
    RealVect get_centroid() const;

    /*!
      @brief Get the centroid
    */
  
    RealVect get_coord() const;

    /*!
      @brief Get lower corner of bounding box
    */
  
    RealVect get_bbox_lo() const;

    /*!
      @brief Get high corner of bounding box
    */
  
    RealVect get_bbox_hi() const;
  
    /*!
      @brief Get face normal vector
    */
  
    RealVect normal() const;

    /*!
      @brief Get polygon area
    */
  
    Real get_area() const;

    /*!
      @brief Signed distance. 
    */
#if 0 // This is the actual signature
  
    Real signed_distance(const RealVect a_x0) const;
#else // Just so that I can make things compile
  
    Real signed_distance(const RealVect a_x0);
#endif
  
  protected:
  
    /*!
      @brief One of the edges of the polygon
    */
    std::shared_ptr<edge> m_edge;

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

}

#endif
