/*!
  @file   dcel_parser.H
  @brief  Various classes for reading mesh files and putting them in a mesh
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_PARSER_
#define _DCEL_PARSER_

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"
#include "dcel_mesh.H"

#include <vector>
#include <memory>
#include <map>

namespace ChomboDischarge {
  namespace dcel {
    namespace parser {
  
      /*!
	@class PLY
	@brief Class for generation a dcel::mesh from the Stanford PLY file format
      */
      template <class T>
      class PLY {
      public:

	using vertex = vertexT<T>;
	using edge   = edgeT<T>;
	using face   = faceT<T>;
	using mesh   = meshT<T>;

	using edgeIterator = edgeIteratorT<T>;

	/*!
	  @brief Read an ASCII .ply file and put it in a mesh
	*/
	inline
	static void readASCII(mesh& a_mesh, const std::string a_filename);

      protected:

	/*!
	  @brief Read an ASCII header
	  @param[out]   a_numVertices  Number of vertices
	  @param[out]   a_numFaces  Number of faces
	  @param[inout] a_inputstream   File
	*/
	inline
	static void readHeaderASCII(int&           a_numVertices,
				    int&           a_numFaces,
				    std::ifstream& a_inputstream);

	/*!
	  @brief Read ASCII vertices. The next getline() from inputstream must read the first vertex (i.e. don't rewind the stream
	  before entering this routine)
	*/
	inline
	static void readVerticesASCII(std::vector<std::shared_ptr<vertex> >& a_vertices,
				      const int                              a_numVertices,
				      std::ifstream&                         a_inputstream);

	/*!
	  @brief Read ASCII faces and create mesh connectivity. The next getline() from inputstream must read the first face,
	  i.e. we assume that read_ascii_vertices was called IMMEDIATELY before this function. That function will center the fstream 
	  on the correct line
	*/
	inline
	static void readFacesASCII(std::vector<std::shared_ptr<face> >&   a_faces,
				   std::vector<std::shared_ptr<edge> >&   a_edges,
				   std::vector<std::shared_ptr<vertex> >& a_vertices,
				   const int                              a_numFaces,
				   std::ifstream&                         a_inputstream);

	/*!
	  @brief This sets all pair edges
	*/
	inline
	static void reconcilePairEdges(std::vector<std::shared_ptr<edge> >& a_edges);

	/*!
	  @brief Clean up memory
	*/
	inline
	static void clearFaceCache(std::vector<std::shared_ptr<vertex> >& a_vertices);
      };
    }
  }
}

#include "dcel_parserI.H"

#endif
