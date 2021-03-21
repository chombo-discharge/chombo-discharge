/*!
  @file   dcel_face.H
  @brief  Implementation of dcel_face.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_FACE__
#define _DCEL_FACE__

#include "Vec.H"

#include "dcel_algorithms.H"
#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_poly.H"

#include <memory>
#include <vector>

namespace dcel {

  template <class T> class vertexT;
  template <class T> class edgeT;
  template <class T> class faceT;
  template <class T> class edgeIteratorT;
  
  template <class T>
  class faceT {
  public:

    using Vec3   = Vec3T<T>;
    
    using vertex = vertexT<T>;
    using edge   = edgeT<T>;
    using face   = faceT<T>;

    using vertexPtr = std::shared_ptr<vertex>;
    using edgePtr   = std::shared_ptr<edge>;
    using facePtr   = std::shared_ptr<face>;

    using edgeIterator = edgeIteratorT<T>;
    
    faceT();
    
    faceT(const edgePtr& a_edge);
    
    faceT(const face& a_otherFace);
    
    ~faceT();

    inline
    decltype(auto) define(const Vec3& a_normal, const edgePtr& a_edge) noexcept;

    inline
    decltype(auto) reconcile() noexcept;

    inline
    decltype(auto) setHalfEdge(const edgePtr& a_halfEdge) noexcept;

    inline
    decltype(auto) setInsideOutsideAlgorithm(const InsideOutsideAlgorithm& a_algorithm) noexcept;

    inline
    decltype(auto) getHalfEdge() noexcept;

    inline
    decltype(auto) getHalfEdge() const noexcept;

    inline
    decltype(auto) getCentroid() noexcept;

    inline
    decltype(auto) getCentroid() const noexcept;

    inline
    decltype(auto) getNormal() noexcept;

    inline
    decltype(auto) getNormal() const noexcept;

    inline
    decltype(auto) signedDistance(const Vec3& a_x0) const noexcept;

    inline
    decltype(auto) unsignedDistance2(const Vec3& a_x0) const noexcept;

    inline
    decltype(auto) getAllVertexCoordinates() const noexcept;
  
  protected:
    
    edgePtr m_halfEdge; // One of the edges. This must always be != nullptr for a valid face.
    
    std::vector<edgePtr > m_edges; // Exists because of performance reasons. 

    T m_area;
    
    Vec3 m_normal;
    Vec3 m_centroid;

    // 2D polygon
    std::shared_ptr<Polygon2D<T> > m_poly2;
    InsideOutsideAlgorithm m_poly2Algorithm;

    inline
    decltype(auto) computeArea() noexcept;

    inline
    decltype(auto) computeCentroid() noexcept;

    inline
    decltype(auto) computeNormal() noexcept;

    inline
    decltype(auto) computePolygon2D() noexcept;
    
    inline
    decltype(auto) normalizeNormalVector() noexcept;

    inline
    decltype(auto) getArea() noexcept;

    inline
    decltype(auto) getArea() const noexcept;

    inline
    decltype(auto) gatherVertices() const noexcept;

    inline
    decltype(auto) computeAndStoreEdges() noexcept;

    inline
    decltype(auto) projectPointIntoFacePlane(const Vec3& a_p) const noexcept;

    inline
    decltype(auto) isPointInsideFace(const Vec3& a_p) const noexcept;
  };

}

#include "dcel_faceI.H"

#endif
