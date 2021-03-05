/*!
  @file   dcel_parser.H
  @brief  Various classes for reading mesh files and putting them in a mesh
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_PARSER_
#define _DCEL_PARSER_

#include "dcel_mesh.H"

#include <vector>
#include <memory>

namespace dcel {

  namespace parser {
  
    /*!
      @class PLY
      @brief Class for generation a dcel::mesh from the Stanford PLY file format
    */
    class PLY {
    public:

      /*!
	@brief Read an ASCII .ply file and put it in a mesh
      */
      static void read_ascii(dcel::mesh& a_mesh, const std::string a_filename);

    protected:

      /*!
	@brief Read an ASCII header
	@param[out]   a_num_vertices  Number of vertices
	@param[out]   a_num_polygons  Number of polygons
	@param[inout] a_inputstream   File
      */
      static void read_ascii_header(int&           a_num_vertices,
				    int&           a_num_polygons,
				    std::ifstream& a_inputstream);

      /*!
	@brief Read ASCII vertices. The next getline() from inputstream must read the first vertex (i.e. don't rewind the stream
	before entering this routine)
      */
      static void read_ascii_vertices(std::vector<std::shared_ptr<dcel::vertex> >& a_vertices,
				      const int                                 a_num_vertices,
				      std::ifstream&                            a_inputstream);

      /*!
	@brief Read ASCII polygons and create mesh connectivity. The next getline() from inputstream must read the first polygon,
	i.e. we assume that read_ascii_vertices was called IMMEDIATELY before this function. That function will center the fstream 
	on the correct line
      */
      static void read_ascii_polygons(std::vector<std::shared_ptr<dcel::polygon> >& a_polygons,
				      std::vector<std::shared_ptr<dcel::edge> >& a_edges,
				      std::vector<std::shared_ptr<dcel::vertex> >& a_vertices,
				      const int                                 a_num_polygons,
				      std::ifstream&                            a_inputstream);
    };
  }
}

#endif
