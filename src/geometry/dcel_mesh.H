/*!
  @file   dcel_mesh.H
  @brief  Declaration of a mesh class for handling surface tesselations
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_MESH_
#define _DCEL_MESH_

#include "dcel_BoundingVolumes.H"
#include "dcel_KdTree.H"
#include "dcel_BVH.H"

#include <vector>
#include <memory>

namespace dcel {


  class vertex;
  class edge;
  class polygon;

  class mesh {
  public:

    enum class VertexNormalWeight {
      None,
      Angle,
    };

    enum class SearchAlgorithm{
      Direct,
      KdTree
    };

    mesh();
    mesh(const mesh& a_otherMesh) = delete;
    mesh(std::vector<std::shared_ptr<polygon> >& a_polygons,
	 std::vector<std::shared_ptr<edge> >&    a_edges,
	 std::vector<std::shared_ptr<vertex> >&  a_vertices);

    ~mesh();

    void sanityCheck() const;
    void define(std::vector<std::shared_ptr<polygon> >& a_polygons,
		std::vector<std::shared_ptr<edge> >&    a_edges,
		std::vector<std::shared_ptr<vertex> >&  a_vertices) noexcept;
    void setAlgorithm(SearchAlgorithm a_algorithm) noexcept;

    std::vector<std::shared_ptr<vertex> >& getVertices() noexcept;
    const std::vector<std::shared_ptr<vertex> >& getVertices() const noexcept;

    std::vector<std::shared_ptr<edge> >& getEdges() noexcept ;
    const std::vector<std::shared_ptr<edge> >& getEdges() const noexcept;

    std::vector<std::shared_ptr<polygon> >& getPolygons() noexcept;
    const std::vector<std::shared_ptr<polygon> >& getPolygons() const noexcept;

    std::vector<RealVect> getAllVertexCoordinates() const noexcept;

    void computeBoundingSphere() noexcept;
    void computeBoundingBox() noexcept;
    void reconcilePolygons(const bool a_outwardNormal) noexcept;
    void computeVertexNormals(VertexNormalWeight a_weight) noexcept;
    void computeEdgeNormals() noexcept;

    void buildKdTree(const int a_max_depth, const int a_max_elements) noexcept;

    Real signedDistance(const RealVect& a_x0) const noexcept ;
    Real signedDistance(const RealVect& a_x0, SearchAlgorithm a_algorithm) const noexcept;

  protected:

    SearchAlgorithm m_algorithm;
    bool m_reconciled;
    bool m_use_tree;

    BoundingSphere m_boundingSphere;
    AABB m_boundingBox;

    std::shared_ptr<kd_tree<polygon> > m_tree;

    std::vector<std::shared_ptr<vertex> >  m_vertices;
    std::vector<std::shared_ptr<edge> >    m_edges;
    std::vector<std::shared_ptr<polygon> > m_polygons;

    Real DirectSignedDistance(const RealVect& a_point) const noexcept;
    Real KdTreeSignedDistance(const RealVect& a_point) const noexcept;

    void computeVertexNormalAngleWeighted(std::shared_ptr<vertex>& a_vert) noexcept;
    void computeVertexNormalAverage(std::shared_ptr<vertex>& a_vert) noexcept;
  };
}

#endif
