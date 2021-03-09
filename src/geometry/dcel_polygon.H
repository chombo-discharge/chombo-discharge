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

  struct Point2D {
    Point2D() = default;
    Point2D(const Real a_x, const Real a_y);
    ~Point2D() = default;
    Real x;
    Real y;
  };

  class polygon {
  public:
  
    polygon();
    polygon(const std::shared_ptr<edge>& a_edge);
    polygon(const polygon& a_otherPolygon);
    ~polygon();

    void define(const RealVect& a_normal, const std::shared_ptr<edge>& a_edge) noexcept;

    void setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept;
    void setNormal(const RealVect& a_normal) noexcept;

    void normalizeNormalVector() noexcept;
    void computeArea() noexcept;
    void computeCentroid() noexcept;
    void computeNormal(const bool a_outwardNormal) noexcept;
    void computeBoundingBox() noexcept;
    void computeBoundingSphere() noexcept;
    void computeVerticesAndEdges() noexcept;


    std::shared_ptr<edge>& getHalfEdge() noexcept;
    const std::shared_ptr<edge>& getHalfEdge() const noexcept;

    const std::vector<RealVect> getAllVertexCoordinates() const noexcept;

    std::vector<std::shared_ptr<vertex> >& getVertices() noexcept;
    const std::vector<std::shared_ptr<vertex> >& getVertices() const noexcept;

    std::vector<std::shared_ptr<edge> >& getEdges() noexcept;
    const std::vector<std::shared_ptr<edge> >& getEdges() const noexcept;

    const std::vector<std::shared_ptr<vertex> > gatherVertices() const noexcept;
    const std::vector<std::shared_ptr<edge> > gatherEdges() const noexcept;

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
    
    Real signedDistance(const RealVect& a_x0) const noexcept;
    Real unsignedDistance2(const RealVect& a_x0) const noexcept;

#if 1 // Just so kd will shut up for now. 
    RealVect get_coord() const {
      return m_centroid;
    }
#endif
  
  protected:


  
    std::shared_ptr<edge> m_halfEdge; // One of the edges. This must always be != nullptr

    std::vector<std::shared_ptr<vertex> > m_vertices;
    std::vector<std::shared_ptr<edge> > m_edges;

    Real     m_area;
    RealVect m_normal;
    RealVect m_centroid;
    RealVect m_lo;
    RealVect m_hi;

    // Potential bounding volumes
    BoundingSphere m_boundingSphere;
    AABB m_boundingBox;

    // Projects three-dimensional point a_p into the plane of the polygon
    RealVect projectPointIntoPolygonPlane(const RealVect& a_p) const noexcept;

    // Function that computes the subtended angles and checks for inside/outside
    bool isPointInsidePolygonAngleSum(const RealVect& a_p) const noexcept;
    bool isPointInsidePolygonWindingNumber(const RealVect& a_p) const noexcept;

    // These functions are for computing inside/outside tests using the winding number algorithm
    int m_ignoreDir;
    int m_xDir;
    int m_yDir;
    std::vector<Point2D> m_polygon2D;
    void computePolygon2D() noexcept; // Computes m_ignoreDir, m_xDir, m_yDir and m_polygon2D
    Point2D projectPointTo2D(const RealVect& a_x) const noexcept; // Must call computePolygon2D first
    int isLeft(const Point2D& a_P0, const Point2D& a_P1, const Point2D& a_P2) const noexcept; // Returns > 0 if left, < 0 if right, = 0 if on line.
    int wn_PnPoly(const Point2D& P, const std::vector<Point2D>& a_vertices) const noexcept;

  };

}

#endif
