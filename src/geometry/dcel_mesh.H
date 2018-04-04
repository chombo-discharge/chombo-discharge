/*!
  @file dcel_mesh.H
  @brief Declaration of a dcel_mesh class for handling surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_MESH_
#define _DCEL_MESH_

#include "Vector.H"

#include "dcel_poly.H"
#include "dcel_edge.H" 
#include "dcel_vert.H"

/*!
  @brief Class for handling surface tesselations
*/
class dcel_mesh {
public:

  /*!
    @brief Constructor
  */
  dcel_mesh();

  /*!
    @brief Full constructor
  */
  dcel_mesh(Vector<dcel_poly*> a_polygons, Vector<dcel_edge*> a_edges, Vector<dcel_vert*> a_vertices);

  /*!
    @brief Destructor
  */
  ~dcel_mesh();

  /*!
    @brief Define a mesh from existing polygons and edges
  */
  virtual void define(Vector<dcel_poly*> a_polygons, Vector<dcel_edge*> a_edges, Vector<dcel_vert*> a_vertices);

  /*!
    @brief Compute signed distance
  */
  virtual Real signed_distance(const RealVect a_x0);

protected:

  /*!
    @brief Mesh polygons
  */
  Vector<dcel_poly*> m_polygons;

  /*!
    @brief Mesh edges
  */
  Vector<dcel_edge*> m_edges;

  /*!
    @brief Mesh vertices
  */
  Vector<dcel_vert*> m_vertices;
};
#endif
