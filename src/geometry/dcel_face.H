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

#include <memory>
#include <vector>

namespace dcel {

  class vertex;
  class edge;
  class face;

  using Vec3 = Vec3T<double>;

  class face {
  public:
    
    face();
    
    face(const std::shared_ptr<edge>& a_edge);
    
    face(const face& a_otherFace);
    
    ~face();

    inline
    void define(const Vec3T<double>& a_normal, const std::shared_ptr<edge>& a_edge) noexcept;

    inline
    void reconcile() noexcept;

    inline
    void setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept;

    inline
    void setInsideOutsideAlgorithm(const InsideOutsideAlgorithm& a_algorithm) noexcept;

    inline
    std::shared_ptr<edge>& getHalfEdge() noexcept;

    inline
    const std::shared_ptr<edge>& getHalfEdge() const noexcept;

    inline
    Vec3T<double>& getCentroid() noexcept;

    inline
    const Vec3T<double>& getCentroid() const noexcept;

    inline
    Vec3T<double>& getNormal() noexcept;

    inline
    const Vec3T<double>& getNormal() const noexcept;

    inline
    Vec3T<double>& getBoundingBoxLo() noexcept;

    inline
    const Vec3T<double>& getBoundingBoxLo() const noexcept;

    inline
    Vec3T<double>& getBoundingBoxHi() noexcept;

    inline
    const Vec3T<double>& getBoundingBoxHi() const noexcept;

    inline
    double signedDistance(const Vec3T<double>& a_x0) const noexcept;

    inline
    double unsignedDistance2(const Vec3T<double>& a_x0) const noexcept;
  
  protected:
    
    std::shared_ptr<edge> m_halfEdge; // One of the edges. This must always be != nullptr for a valid face. 

    std::vector<std::shared_ptr<vertex> > m_vertices; 
    std::vector<std::shared_ptr<edge> > m_edges;

    double       m_area;
    Vec3T<double> m_normal;
    Vec3T<double> m_centroid;
    Vec3T<double> m_lo;
    Vec3T<double> m_hi;

    // Potential bounding volumes
    BoundingSphere m_boundingSphere;
    AABB m_boundingBox;

    // 2D polygon
    std::shared_ptr<Polygon2D<double> > m_poly2;
    InsideOutsideAlgorithm m_poly2Algorithm;

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
    void computePolygon2D() noexcept;
    
    inline
    void normalizeNormalVector() noexcept;

    inline
    double& getArea() noexcept;

    inline
    const double& getArea() const noexcept;

    inline
    const std::vector<Vec3T<double> > getAllVertexCoordinates() const noexcept;

    inline
    const std::vector<std::shared_ptr<vertex> > gatherVertices() const noexcept;

    // Projects three-dimensional point a_p into the plane of the face
    inline
    Vec3T<double> projectPointIntoFacePlane(const Vec3T<double>& a_p) const noexcept;

    inline
    bool isPointInsideFace(const Vec3T<double>& a_p) const noexcept;
  };

}

#include "dcel_faceI.H"

#endif
