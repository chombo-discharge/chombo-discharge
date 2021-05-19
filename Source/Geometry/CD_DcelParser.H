/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DcelParser.H
  @brief  Utilities for passing data into dcel structures. 
  @author Robert Marskar
*/

#ifndef CD_DcelParser_H
#define CD_DcelParser_H

// Std includes
#include <vector>
#include <memory>
#include <map>

// Our includes
#include <CD_DcelVertex.H>
#include <CD_DcelEdge.H>
#include <CD_DcelFace.H>
#include <CD_DcelMesh.H>
#include "CD_NamespaceHeader.H"

namespace Dcel {
  namespace parser {
  
    /*!
      @class PLY
      @brief Class for generation a Dcel::mesh from the Stanford PLY file format
    */
    template <class T>
    class PLY {
    public:

      using vertex = vertexT<T>;
      using edge   = edgeT<T>;
      using Face   = FaceT<T>;
      using Mesh   = MeshT<T>;

      using edgeIterator = edgeIteratorT<T>;

      /*!
	@brief Read an ASCII .ply file and put it in a mesh
      */
      inline
      static void readASCII(Mesh& a_mesh, const std::string a_filename);

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
      static void readFacesASCII(std::vector<std::shared_ptr<Face> >&   a_faces,
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

#include <CD_NamespaceFooter.H>

#include <CD_DcelParserImplem.H>

#endif
