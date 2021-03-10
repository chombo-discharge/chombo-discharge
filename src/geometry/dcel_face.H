/*!
  @file   dcel_face.H
  @brief  Declaration of a polygon class for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_FACE__
#define _DCEL_FACE__

#include "dcel_BoundingVolumes.H"
#include "dcel_algorithms.H"

#include "RealVect.H"

#include <memory>
#include <vector>

namespace dcel {

  class vertex;
  class edge;
  class face;
  class Polygon2D;

  // Class that describes one face of a watertight surface tesselation. This class
  // can also embed a 2D version of itself as a two-dimensional polygon, projected onto the
  // cardinal plane which maximizes the projected face's area. We use this embedding
  // for speeding up tests when checking if a point projected onto the plane of the face
  // lies inside the face or not. This takes storage, but is required for doing those
  // intersection tests. 
  class face {
  public:
    
    face();
    face(const std::shared_ptr<edge>& a_edge);
    face(const face& a_otherFace);
    ~face();

    void define(const RealVect& a_normal, const std::shared_ptr<edge>& a_edge) noexcept;

    void setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept;
    void setNormal(const RealVect& a_normal) noexcept;
    void setInsideOutsideAlgorithm(const InsideOutsideAlgorithm& a_algorithm) noexcept;
    void normalizeNormalVector() noexcept;
    void computeArea() noexcept;
    void computeCentroid() noexcept;
    void computeNormal() noexcept;
    void computeBoundingBox() noexcept;
    void computeBoundingSphere() noexcept;
    void computeVerticesAndEdges() noexcept;
    void computePolygon2D() noexcept;

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

#if 1 // Just so kd will shut up for now. Will redo the BVH interface when all bugs are squashed in the direct brute-force method. 
    RealVect get_coord() const {
      return m_centroid;
    }
#endif
  
  protected:
    
    std::shared_ptr<edge> m_halfEdge; // One of the edges. This must always be != nullptr for a valid face. 

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

    // 2D polygon
    std::shared_ptr<Polygon2D> m_poly2;
    InsideOutsideAlgorithm m_poly2Algorithm;

    // Projects three-dimensional point a_p into the plane of the face
    RealVect projectPointIntoFacePlane(const RealVect& a_p) const noexcept;
    bool isPointInsideFace(const RealVect& a_p) const noexcept;
  };

}
#endif
