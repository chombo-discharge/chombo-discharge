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

#include "dcel_vec.H"
#include "dcel_vertex.H"
#include "dcel_edge.H"

#include <memory>
#include <vector>

namespace dcel {

  template <class T> class vertexT;
  template <class T> class edgeT;
  template <class T> class faceT;
  template <class T> class edge_iteratorT;
  
  template <class T>
  class faceT {
  public:

    using Vec3   = Vec3T<T>;
    using vertex = vertexT<T>;
    using edge   = edgeT<T>;
    using face   = faceT<T>;

    using edge_iterator = edge_iteratorT<T>;
    
    faceT();
    
    faceT(const std::shared_ptr<edge>& a_edge);
    
    faceT(const face& a_otherFace);
    
    ~faceT();

    inline
    auto define(const Vec3& a_normal, const std::shared_ptr<edge>& a_edge) noexcept -> void;

    inline
    auto reconcile() noexcept -> void;

    inline
    auto setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept -> void;

    inline
    auto setInsideOutsideAlgorithm(const InsideOutsideAlgorithm& a_algorithm) noexcept -> void;

    inline
    auto getHalfEdge() noexcept -> std::shared_ptr<edge>&;

    inline
    auto getHalfEdge() const noexcept -> const std::shared_ptr<edge>&;

    inline
    auto getCentroid() noexcept -> Vec3&;

    inline
    auto getCentroid() const noexcept -> const Vec3&;

    inline
    auto getNormal() noexcept -> Vec3&;

    inline
    auto getNormal() const noexcept -> const Vec3&;

    inline
    auto getBoundingBoxLo() noexcept -> Vec3&;

    inline
    auto getBoundingBoxLo() const noexcept -> const Vec3&;

    inline
    auto getBoundingBoxHi() noexcept -> Vec3&;

    inline
    auto getBoundingBoxHi() const noexcept -> const Vec3&;

    inline
    auto signedDistance(const Vec3& a_x0) const noexcept -> T;

    inline
    auto unsignedDistance2(const Vec3& a_x0) const noexcept -> T;
  
  protected:
    
    std::shared_ptr<edge> m_halfEdge; // One of the edges. This must always be != nullptr for a valid face. 

    std::vector<std::shared_ptr<vertex> > m_vertices; 
    std::vector<std::shared_ptr<edge> > m_edges;

    T m_area;
    
    Vec3 m_normal;
    Vec3 m_centroid;
    Vec3 m_lo;
    Vec3 m_hi;

    // Potential bounding volumes
    BoundingSphere<T> m_boundingSphere;
    AABB<T> m_boundingBox;

    // 2D polygon
    std::shared_ptr<Polygon2D<T> > m_poly2;
    InsideOutsideAlgorithm m_poly2Algorithm;

    inline
    auto computeArea() noexcept -> void;

    inline
    auto computeCentroid() noexcept -> void;

    inline
    auto computeNormal() noexcept -> void;

    inline
    auto computeBoundingBox() noexcept -> void;

    inline
    auto computeBoundingSphere() noexcept -> void;

    inline
    auto computePolygon2D() noexcept -> void;
    
    inline
    auto normalizeNormalVector() noexcept -> void;

    inline
    auto getArea() noexcept -> T&;

    inline
    auto getArea() const noexcept -> const T& ;

    inline
    auto getAllVertexCoordinates() const noexcept -> std::vector<Vec3>;

    inline
    auto gatherVertices() const noexcept -> std::vector<std::shared_ptr<vertex> >;

    inline
    auto projectPointIntoFacePlane(const Vec3& a_p) const noexcept -> Vec3;

    inline
    auto isPointInsideFace(const Vec3& a_p) const noexcept -> bool;
  };

}

#include "dcel_faceI.H"

#endif
