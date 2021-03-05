/*!
  @file   dcel_poly.H
  @brief  Declaration of a polygon class for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_POLY_
#define _DCEL_POLY_

#include <dcel_BoundingVolumes.H>

#include "RealVect.H"

#include <memory>
#include <vector>

namespace dcel {

  class vertex;
  class edge;

  class polygon {
  public:
  
    polygon();
    polygon(const polygon& a_otherPolygon);
    ~polygon();

    void define(const RealVect& a_normal, const std::shared_ptr<edge>& a_edge) noexcept;

    void setEdge(const std::shared_ptr<edge>& a_edge) noexcept;
    void setNormal(const RealVect& a_normal) noexcept;

    void normalizeNormalVector() noexcept;
    void computeArea() noexcept;
    void computeCentroid() noexcept;
    void computeNormal(const bool a_outwardNormal) noexcept;
    void computeBoundingBox() noexcept;
    void computeBoundingSphere() noexcept;

    std::shared_ptr<edge>& getEdge() noexcept;
    const std::shared_ptr<edge>& getEdge() const noexcept;

    const std::vector<RealVect> getAllVertexCoordinates() const noexcept;
    const std::vector<std::shared_ptr<vertex> > getVertices() const noexcept;
    const std::vector<std::shared_ptr<edge> > getEdges() const noexcept;

    RealVect& getCentroid() noexcept;
    const RealVect& getCentroid() const noexcept;

    RealVect& getNormal() noexcept;
    const RealVect& getNormal() const noexcept;

    Real& getArea() noexcept;
    const Real& getArea() const noexcept;

    RealVect& getBoundingBoxLo() noexcept;
    const RealVect& getBoundingBoxLo() const noexcept;

    RealVect& getBoundingBoxHi() noexcept;
    const RealVect& getBoundingBoxHi() const noexcept;
    
    // Real signed_distance(const RealVect a_x0) const; Should be this signature. 
    Real signedDistance(const RealVect& a_x0) const noexcept;

#if 1 // Just so kd will shut up for now. 
    RealVect get_coord() const {
      return m_centroid;
    }
#endif
  
  protected:
  
    std::shared_ptr<edge> m_edge;

    Real m_area;
    
    RealVect m_normal;
    RealVect m_centroid;
    RealVect m_lo;
    RealVect m_hi;

    // Bounding volumes
    BoundingSphere m_boundingSphere;
    AABB m_boundingBox;
  };

}

#endif
