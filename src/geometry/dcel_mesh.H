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

  template <class T> class vertexT;
  template <class T> class edgeT;
  template <class T> class faceT;

  template <class T>
  class mesh {
  public:

    using Vec3   = Vec3T<T>;
    using vertex = vertexT<T>;
    using edge   = edgeT<T>;
    using face   = faceT<T>;

    mesh();
    mesh(const mesh& a_otherMesh) = delete;
    mesh(std::vector<std::shared_ptr<face> >&   a_faces,
	 std::vector<std::shared_ptr<edge> >&   a_edges,
	 std::vector<std::shared_ptr<vertex> >& a_vertices);

    ~mesh();



    inline
    auto define(std::vector<std::shared_ptr<face> >&   a_faces,
		std::vector<std::shared_ptr<edge> >&   a_edges,
		std::vector<std::shared_ptr<vertex> >& a_vertices) noexcept -> void;
    
    inline
    auto sanityCheck() const noexcept -> void;

    inline
    auto setSearchAlgorithm(const SearchAlgorithm a_algorithm) noexcept -> void;

    inline
    auto setInsideOutsideAlgorithm(const dcel::InsideOutsideAlgorithm a_algorithm) noexcept -> void;

    inline
    auto reconcile(VertexNormalWeight a_weight = VertexNormalWeight::Angle) noexcept -> void;

    inline
    auto computeBoundingSphere() noexcept -> void;

    inline
    auto computeBoundingBox() noexcept -> void;

    inline
    auto getVertices() noexcept -> std::vector<std::shared_ptr<vertex> >&;

    inline
    auto getVertices() const noexcept -> const std::vector<std::shared_ptr<vertex> >&;

    inline
    auto getEdges() noexcept -> std::vector<std::shared_ptr<edge> >&;

    inline
    auto getEdges() const noexcept -> const std::vector<std::shared_ptr<edge> >&;

    inline
    auto getFaces() noexcept -> std::vector<std::shared_ptr<face> >&;

    inline
    auto getFaces() const noexcept -> const std::vector<std::shared_ptr<face> >&;

    inline
    auto signedDistance(const Vec3& a_x0) const noexcept -> T;

    inline
    auto signedDistance(const Vec3& a_x0, SearchAlgorithm a_algorithm) const noexcept -> T;

  protected:

    SearchAlgorithm m_algorithm;
    bool m_use_tree;

    BoundingSphere<T> m_boundingSphere;
    AABB<T> m_boundingBox;

    std::vector<std::shared_ptr<vertex> >  m_vertices;
    std::vector<std::shared_ptr<edge> >    m_edges;
    std::vector<std::shared_ptr<face> > m_faces;

    inline
    auto getAllVertexCoordinates() const noexcept -> std::vector<Vec3 >;

    inline
    auto reconcileFaces() noexcept -> void;

    inline
    auto reconcileEdges() noexcept -> void;

    inline
    auto reconcileVertices(VertexNormalWeight a_weight) noexcept -> void;

    inline
    auto DirectSignedDistance(const Vec3& a_point) const noexcept -> T;

    inline
    auto DirectSignedDistance2(const Vec3& a_point) const noexcept -> T;
  };
}

#include "dcel_meshI.H"
#endif
