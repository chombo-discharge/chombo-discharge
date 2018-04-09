/*!
  @file   ply_reader.H
  @brief  Static class for reading .ply-files and putting them in the half-edge mesh data structure
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _PLY_READER_
#define _PLY_READER_

#include "dcel_mesh.H"

class ply_reader {
public:

  /*!
    @brief Read an ASCII .ply file and put it in a mesh
  */
  static void read_ascii(dcel_mesh& a_mesh, const std::string a_filename);

  /*!
    @brief Read a binary .ply file and put it in a mesh. 
  */
  static void read_binary(dcel_mesh& a_mesh, const std::string a_filename);

protected:

  /*!
    @brief Read an ASCII header
    @param[out] a_first_vertex  The line number containing the first vertex
    @param[out] a_first_polygon The line number containing the first polygon
    @param[out] a_num_vertices  Number of vertices
    @param[out] a_num_polygons  Number of polygons
  */
  static void read_ascii_header(int& a_first_vertex,
				int& a_first_polygon,
				int& a_num_vertices,
				int& a_num_polygons,
				std::ifstream& a_inputstream);

  /*!
    @brief Read ASCII vertices. The next getline() from inputstream must read the first vertex (i.e. don't rewind the stream)
  */
  static void read_ascii_vertices(Vector<RefCountedPtr<dcel_vert> >& vertices,
				  const int a_num_vertices,
				  std::ifstream& a_inputstream);

  /*!
    @brief Read ASCII polygons and create mesh connectivity. The next getline() from inputstream must read the first polygon,
    i.e. we assume that read_ascii_vertices was called IMMEDIATELY before this function. That function will center the fstream 
    on the correct line
  */
  static void read_ascii_polygons(Vector<RefCountedPtr<dcel_poly> >& a_polygons,
				  Vector<RefCountedPtr<dcel_edge> >& a_edges,
				  Vector<RefCountedPtr<dcel_vert> >& a_vertices,
				  const int a_num_polygons,
				  std::ifstream& a_inputstream);
};
#endif
