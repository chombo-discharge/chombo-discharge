/*!
  @file   dcel_face.H
  @brief  Implementation of dcel_face.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_FACE__
#define _DCEL_FACE__

#include "dcel_BoundingVolumes.H"
#include "dcel_algorithms.H"
#include "dcel_poly.H"

#include "RealVect.H"

#include <memory>
#include <vector>

namespace dcel {

  class vertex;
  class edge;
  class face;

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

    inline
    void define(const RealVect& a_normal, const std::shared_ptr<edge>& a_edge) noexcept;

    inline
    void setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept;
    
    inline
    void setNormal(const RealVect& a_normal) noexcept;

    inline
    void normalizeNormalVector() noexcept;

    inline
    void setInsideOutsideAlgorithm(const InsideOutsideAlgorithm& a_algorithm) noexcept;

    inline
    void computeArea() noexcept;

    inline
    void computeCentroid() noexcept;

    inline
    void computeNormal() noexcept;

    inline
    void computeBoundingBox() noexcept;

    inline
    void computeBoundingSphere() noexcept;

    inline
    void computeVerticesAndEdges() noexcept;

    inline
    void computePolygon2D() noexcept;

    inline
    std::shared_ptr<edge>& getHalfEdge() noexcept;

    inline
    const std::shared_ptr<edge>& getHalfEdge() const noexcept;

    inline
    const std::vector<RealVect> getAllVertexCoordinates() const noexcept;

    inline
    std::vector<std::shared_ptr<vertex> >& getVertices() noexcept;

    inline
    const std::vector<std::shared_ptr<vertex> >& getVertices() const noexcept;

    inline
    std::vector<std::shared_ptr<edge> >& getEdges() noexcept;
    
    inline
    const std::vector<std::shared_ptr<edge> >& getEdges() const noexcept;

    inline
    const std::vector<std::shared_ptr<vertex> > gatherVertices() const noexcept;

    inline
    const std::vector<std::shared_ptr<edge> > gatherEdges() const noexcept;

    inline
    RealVect& getCentroid() noexcept;

    inline
    const RealVect& getCentroid() const noexcept;

    inline
    RealVect& getNormal() noexcept;

    inline
    const RealVect& getNormal() const noexcept;

    inline
    Real& getArea() noexcept;

    inline
    const Real& getArea() const noexcept;

    inline
    RealVect& getBoundingBoxLo() noexcept;

    inline
    const RealVect& getBoundingBoxLo() const noexcept;

    inline
    RealVect& getBoundingBoxHi() noexcept;

    inline
    const RealVect& getBoundingBoxHi() const noexcept;

    inline
    Real signedDistance(const RealVect& a_x0) const noexcept;

    inline
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
    std::shared_ptr<Polygon2D<double> > m_poly2;
    InsideOutsideAlgorithm m_poly2Algorithm;

    // Projects three-dimensional point a_p into the plane of the face
    RealVect projectPointIntoFacePlane(const RealVect& a_p) const noexcept;
    bool isPointInsideFace(const RealVect& a_p) const noexcept;
  };

}

#include "dcel_faceI.H"

#endif
