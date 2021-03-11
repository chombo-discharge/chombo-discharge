/*!
  @file   dcel_mesh.H
  @brief  Declaration of a mesh class for handling surface tesselations
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_MESH_
#define _DCEL_MESH_

#include "dcel_BoundingVolumes.H"
#include "dcel_BVH.H"
#include "dcel_algorithms.H"

#include <vector>
#include <memory>

namespace dcel {

  class vertex;
  class edge;
  class face;

  class mesh {
  public:

    mesh();
    mesh(const mesh& a_otherMesh) = delete;
    mesh(std::vector<std::shared_ptr<face> >&   a_faces,
	 std::vector<std::shared_ptr<edge> >&   a_edges,
	 std::vector<std::shared_ptr<vertex> >& a_vertices);

    ~mesh();



    inline
    void define(std::vector<std::shared_ptr<face> >&   a_faces,
		std::vector<std::shared_ptr<edge> >&   a_edges,
		std::vector<std::shared_ptr<vertex> >& a_vertices) noexcept;
    
    inline
    void sanityCheck() const;

    inline
    void setSearchAlgorithm(const SearchAlgorithm a_algorithm) noexcept;

    inline
    void setInsideOutsideAlgorithm(const dcel::InsideOutsideAlgorithm a_algorithm) noexcept;

    inline
    void reconcile(VertexNormalWeight a_weight = VertexNormalWeight::Angle) noexcept;

    inline
    void computeBoundingSphere() noexcept;

    inline
    void computeBoundingBox() noexcept;

    inline
    std::vector<std::shared_ptr<vertex> >& getVertices() noexcept;

    inline
    const std::vector<std::shared_ptr<vertex> >& getVertices() const noexcept;

    inline
    std::vector<std::shared_ptr<edge> >& getEdges() noexcept ;

    inline
    const std::vector<std::shared_ptr<edge> >& getEdges() const noexcept;

    inline
    std::vector<std::shared_ptr<face> >& getFaces() noexcept;

    inline
    const std::vector<std::shared_ptr<face> >& getFaces() const noexcept;

    inline
    double signedDistance(const Vec3<double>& a_x0) const noexcept ;

    inline
    double signedDistance(const Vec3<double>& a_x0, SearchAlgorithm a_algorithm) const noexcept;

  protected:

    SearchAlgorithm m_algorithm;
    bool m_use_tree;

    BoundingSphere m_boundingSphere;
    AABB m_boundingBox;

    std::vector<std::shared_ptr<vertex> >  m_vertices;
    std::vector<std::shared_ptr<edge> >    m_edges;
    std::vector<std::shared_ptr<face> > m_faces;

    inline
    std::vector<Vec3<double> > getAllVertexCoordinates() const noexcept;

    inline
    void reconcileFaces() noexcept;

    inline
    void reconcileEdges() noexcept;

    inline
    void reconcileVertices(VertexNormalWeight a_weight) noexcept;

    inline
    double DirectSignedDistance(const Vec3<double>& a_point) const noexcept;

    inline
    double DirectSignedDistance2(const Vec3<double>& a_point) const noexcept;
  };
}

#include "dcel_meshI.H"
#endif
