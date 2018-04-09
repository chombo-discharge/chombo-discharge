/*!
  @file dcel_mesh.H
  @brief Declaration of a dcel_mesh class for handling surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_MESH_
#define _DCEL_MESH_

#include "ritter_sphere.H"
#include "dcel_poly.H"
#include "dcel_edge.H" 
#include "dcel_vert.H"

#include <Vector.H>
#include <RefCountedPtr.H>

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
  dcel_mesh(Vector<RefCountedPtr<dcel_poly> >& a_polygons,
	    Vector<RefCountedPtr<dcel_edge> >& a_edges,
	    Vector<RefCountedPtr<dcel_vert> >& a_vertices);

  /*!
    @brief Destructor
  */
  ~dcel_mesh();

  /*!
    @brief Check that mesh is OK
  */
  virtual bool sanity_check() const;

  /*!
    @brief Define a mesh from existing polygons and edges
  */
  virtual void define(Vector<RefCountedPtr<dcel_poly> >& a_polygons,
		      Vector<RefCountedPtr<dcel_edge> >& a_edges,
		      Vector<RefCountedPtr<dcel_vert> >& a_vertices);

  /*!
    @brief Reconcile polygon edges. This gives each edge a reference to the polygon they circulate, and also computes the 
    polygon area
    @param[in] a_area_weighted Use area weighted vertex normals
  */
  virtual void reconcile_polygons(const bool a_outward_normal);

  /*!
    @brief Compute signed distance
  */
  virtual Real signed_distance(const RealVect a_x0);

  /*!
    @brief Get vertices
  */
  virtual Vector<RefCountedPtr<dcel_vert> >& get_vertices();

  /*!
    @brief Get edges
  */
  virtual Vector<RefCountedPtr<dcel_edge> >& get_edges();

  /*!
    @brief Get polygons
  */
  virtual Vector<RefCountedPtr<dcel_poly> >& get_polygons();


protected:

  bool m_reconciled;

  /*!
    @brief Bounding sphere
  */
  ritter_sphere m_sphere;

  /*!
    @brief Mesh polygons
  */
  Vector<RefCountedPtr<dcel_poly> > m_polygons;

  /*!
    @brief Mesh edges
  */
  Vector<RefCountedPtr<dcel_edge> > m_edges;

  /*!
    @brief Mesh vertices
  */
  Vector<RefCountedPtr<dcel_vert> > m_vertices;

  /*!
    @brief Compute vertex normals
  */
  virtual void compute_vertex_normals();

  /*!
    @brief Compute edge normals
  */
  virtual void compute_edge_normals();

  /*!
    @brief Compute bounding sphere for mesh
  */
  virtual void compute_bounding_sphere();
};
#endif
