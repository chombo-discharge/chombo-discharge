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
      static void readASCII(dcel::mesh& a_mesh, const std::string a_filename);

    protected:

      /*!
	@brief Read an ASCII header
	@param[out]   a_numVertices  Number of vertices
	@param[out]   a_numFaces  Number of faces
	@param[inout] a_inputstream   File
      */
      static void readHeaderASCII(int&           a_numVertices,
				  int&           a_numFaces,
				  std::ifstream& a_inputstream);

      /*!
	@brief Read ASCII vertices. The next getline() from inputstream must read the first vertex (i.e. don't rewind the stream
	before entering this routine)
      */
      static void readVerticesASCII(std::vector<std::shared_ptr<dcel::vertex> >& a_vertices,
				    const int                                    a_numVertices,
				    std::ifstream&                               a_inputstream);

      /*!
	@brief Read ASCII faces and create mesh connectivity. The next getline() from inputstream must read the first face,
	i.e. we assume that read_ascii_vertices was called IMMEDIATELY before this function. That function will center the fstream 
	on the correct line
      */
      static void readFacesASCII(std::vector<std::shared_ptr<dcel::face> >& a_faces,
				    std::vector<std::shared_ptr<dcel::edge> >&    a_edges,
				    std::vector<std::shared_ptr<dcel::vertex> >&  a_vertices,
				    const int                                     a_numFaces,
				    std::ifstream&                                a_inputstream);

      /*!
	@brief This sets all pair edges
      */
      static void reconcilePairEdges(std::vector<std::shared_ptr<dcel::edge> >& a_edges);

      /*!
	@brief Clean up memory
      */
      static void clearFaceCache(std::vector<std::shared_ptr<dcel::vertex> >& a_vertices);

      /*!
	@brief Give each edge a length and a vector for parametrization
      */
      static void computeEdgeLengths(std::vector<std::shared_ptr<dcel::edge> >& a_edges);

      /*!
	@brief Give each face access to edges and vertices so they don't have to be looked up every time. 
      */
      static void computeVerticesAndEdges(std::vector<std::shared_ptr<dcel::face> >& a_faces);
    };
  }
}

#endif
