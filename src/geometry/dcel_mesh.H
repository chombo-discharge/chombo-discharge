/*!
  @file   dcel_mesh.H
  @brief  Declaration of a mesh class for handling surface tesselations
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_MESH_
#define _DCEL_MESH_

#include "ritter_sphere.H"
#include "kd_tree.H"

#include <vector>
#include <memory>

namespace dcel {

  class polygon;
  class edge;
  class vertex;

  class mesh {
  public:

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


    void reconcilePolygons(const bool a_outward_normal, const bool recompute_vnormal);
    void buildKdTree(const int a_max_depth, const int a_max_elements);

    Real signedDistance(const RealVect a_x0);

    std::vector<std::shared_ptr<vertex> >& getVertices();
    const std::vector<std::shared_ptr<vertex> >& getVertices() const;

    std::vector<std::shared_ptr<edge> >& getEdges();
    const std::vector<std::shared_ptr<edge> >& getEdges() const;

    std::vector<std::shared_ptr<polygon> >& getPolygons();
    const std::vector<std::shared_ptr<polygon> >& getPolygons() const;

  protected:


    bool m_reconciled;

    bool m_use_tree;

    ritter_sphere m_sphere;

    std::shared_ptr<kd_tree<polygon> > m_tree;

    std::vector<std::shared_ptr<vertex> >  m_vertices;
    std::vector<std::shared_ptr<edge> >    m_edges;
    std::vector<std::shared_ptr<polygon> > m_polygons;

    void computeVertexNormals();
    void computeEdgeNormals();
    void computeBoundingSphere();
  };
}

#endif
