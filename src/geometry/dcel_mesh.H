/*!
  @file   dcel_mesh.H
  @brief  Declaration of a dcel_mesh class for handling surface tesselations
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

  /*!
    @brief Class for handling surface tesselations
  */
  class dcel_mesh {
  public:

    static bool s_angle_weighted;

    /*!
      @brief Constructor
    */
    dcel_mesh();

    /*!
      @brief Full constructor
    */
    dcel_mesh(std::vector<std::shared_ptr<polygon> >& a_polygons,
	      std::vector<std::shared_ptr<edge> >& a_edges,
	      std::vector<std::shared_ptr<vertex> >& a_vertices);

    /*!
      @brief Destructor
    */
    ~dcel_mesh();

    /*!
      @brief Check that mesh is OK
    */
  
    bool sanity_check() const;

    /*!
      @brief Define a mesh from existing polygons and edges
    */
  
    void define(std::vector<std::shared_ptr<polygon> >& a_polygons,
		std::vector<std::shared_ptr<edge> >& a_edges,
		std::vector<std::shared_ptr<vertex> >& a_vertices);

    /*!
      @brief Reconcile polygon edges. This gives each edge a reference to the polygon they circulate, and also computes the 
      polygon area
    */
  
    void reconcile_polygons(const bool a_outward_normal, const bool recompute_vnormal);

    /*!
      @brief Build kd-tree
    */
  
    void build_tree(const int a_max_depth, const int a_max_elements);

    /*!
      @brief Compute signed distance
    */
  
    Real signed_distance(const RealVect a_x0);

    /*!
      @brief Get vertices
    */
  
    std::vector<std::shared_ptr<vertex> >& get_vertices();

    /*!
      @brief Get edges
    */
  
    std::vector<std::shared_ptr<edge> >& get_edges();

    /*!
      @brief Get polygons
    */
  
    std::vector<std::shared_ptr<polygon> >& get_polygons();

  protected:

    /*!
      @brief Reconciled polygons or not
    */
    bool m_reconciled;

    /*!
      @brief Use kd-tree for searching
    */
    bool m_use_tree;

    /*!
      @brief Bounding sphere
    */
    ritter_sphere m_sphere;

    /*!
      @brief kd-tree for fast polysearch
    */
    std::shared_ptr<kd_tree<polygon> > m_tree;

    /*!
      @brief Mesh polygons
    */
    std::vector<std::shared_ptr<polygon> > m_polygons;

    /*!
      @brief Mesh edges
    */
    std::vector<std::shared_ptr<edge> > m_edges;

    /*!
      @brief Mesh vertices
    */
    std::vector<std::shared_ptr<vertex> > m_vertices;

    /*!
      @brief Compute vertex normals
    */
  
    void compute_vertex_normals();

    /*!
      @brief Compute edge normals
    */
  
    void compute_edge_normals();

    /*!
      @brief Compute bounding sphere for mesh
    */
    void compute_bounding_sphere();
  };
}

#endif
