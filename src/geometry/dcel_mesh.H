/*!
  @file   dcel_mesh.H
  @brief  Declaration of a mesh class for handling surface tesselations
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_MESH_
#define _DCEL_MESH_

#include "dcel_BoundingVolumes.H"
#include "kd_tree.H"

#include <vector>
#include <memory>

namespace dcel {

  class polygon;
  class edge;
  class vertex;

  class mesh {
  public:

    enum class SearchAlgorithm{
      BruteForce,
      KdTree
    };

    static bool s_angle_weighted;

    mesh();
    mesh(const mesh& a_otherMesh) = delete;
    mesh(std::vector<std::shared_ptr<polygon> >& a_polygons,
	 std::vector<std::shared_ptr<edge> >&    a_edges,
	 std::vector<std::shared_ptr<vertex> >&  a_vertices);

    ~mesh();

    bool sanityCheck() const;

    void define(std::vector<std::shared_ptr<polygon> >& a_polygons,
		std::vector<std::shared_ptr<edge> >&    a_edges,
		std::vector<std::shared_ptr<vertex> >&  a_vertices);

    void computeBoundingSphere();
    void computeBoundingBox();
    void computeVertexNormals();
    void computeEdgeNormals();
    void reconcilePolygons(const bool a_outward_normal, const bool recompute_vnormal);
    void buildKdTree(const int a_max_depth, const int a_max_elements);

    Real signedDistance(const RealVect a_x0);
    Real signedDistance(const RealVect& a_x0, SearchAlgorithm a_algorithm) const;

    std::vector<std::shared_ptr<vertex> >& getVertices();
    const std::vector<std::shared_ptr<vertex> >& getVertices() const;

    std::vector<std::shared_ptr<edge> >& getEdges();
    const std::vector<std::shared_ptr<edge> >& getEdges() const;

    std::vector<std::shared_ptr<polygon> >& getPolygons();
    const std::vector<std::shared_ptr<polygon> >& getPolygons() const;

    std::vector<RealVect> getAllVertexCoordinates() const;

  protected:


    bool m_reconciled;

    bool m_use_tree;

    BoundingSphere m_boundingSphere;
    AABB m_boundingBox;

    std::shared_ptr<kd_tree<polygon> > m_tree;

    std::vector<std::shared_ptr<vertex> >  m_vertices;
    std::vector<std::shared_ptr<edge> >    m_edges;
    std::vector<std::shared_ptr<polygon> > m_polygons;


    Real BruteForceSignedDistance(const RealVect& a_point) const;
    Real KdTreeSignedDistance(const RealVect& a_point) const;
  };
}

#endif
