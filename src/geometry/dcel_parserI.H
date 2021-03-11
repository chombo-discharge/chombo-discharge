/*!
  @file   dcel_parserI.H
  @brief  Implementation of dcel_parser.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_PARSERI_H_
#define _DCEL_PARSERI_H_

#include "dcel_parser.H"
#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"
#include "dcel_mesh.H"
#include "dcel_iterator.H"

#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

namespace dcel {

  template <class T>
  inline
  auto parser::PLY<T>::readASCII(mesh<T>& a_mesh, const std::string a_filename) -> void {
    std::ifstream filestream(a_filename);

    if(filestream.is_open()){
      std::vector<std::shared_ptr<vertex> >& vertices  = a_mesh.getVertices();
      std::vector<std::shared_ptr<edge> >& edges       = a_mesh.getEdges();
      std::vector<std::shared_ptr<face> >& faces       = a_mesh.getFaces();

      vertices.resize(0);
      edges.resize(0);
      faces.resize(0);

      int numVertices;  // Number of vertices
      int numFaces;  // Number of faces

      dcel::parser::PLY<T>::readHeaderASCII(numVertices, numFaces, filestream); 
      dcel::parser::PLY<T>::readVerticesASCII(vertices, numVertices, filestream);
      dcel::parser::PLY<T>::readFacesASCII(faces, edges, vertices, numFaces, filestream);
      dcel::parser::PLY<T>::reconcilePairEdges(edges);
      dcel::parser::PLY<T>::clearFaceCache(vertices);
    
      a_mesh.sanityCheck();
    
      filestream.close();
    }
    else{
      const std::string error = "dcel::parser::PLY::readASCII - ERROR! Could not open file " + a_filename;
      std::cerr << error + "\n";
    }
  }

  template <class T>
  inline
  auto parser::PLY<T>::readHeaderASCII(int&           a_numVertices,
				       int&           a_numFaces,
				       std::ifstream& a_inputstream) -> void{

    std::string str1;
    std::string str2;
    std::string line;

    // Get number of vertices
    a_inputstream.clear();
    a_inputstream.seekg(0);
    while (getline(a_inputstream, line)){
      std::stringstream sstream(line);
      sstream >> str1 >> str2 >> a_numVertices;
      if(str1 == "element" && str2 == "vertex"){
	break;
      }
    }

    // Get number of faces
    a_inputstream.clear();
    a_inputstream.seekg(0);
    while (getline(a_inputstream, line)){
      std::stringstream sstream(line);
      sstream >> str1 >> str2 >> a_numFaces;
      if(str1 == "element" && str2 == "face"){
	break;
      }
    }

    // Find the line # containing "end_header" halt the input stream there
    a_inputstream.clear();
    a_inputstream.seekg(0);
    while (getline(a_inputstream, line)){
      std::stringstream sstream(line);
      sstream >> str1;
      if(str1 == "end_header"){
	break;
      }
    }
  }

  template <class T>
  inline
  auto parser::PLY<T>::readVerticesASCII(std::vector<std::shared_ptr<vertex> >& a_vertices,
					 const int                              a_numVertices,
					 std::ifstream&                         a_inputstream) -> void {

    Vec3T<double> pos;
    double& x = pos[0];
    double& y = pos[1];
    double& z = pos[2];

    Vec3T<double> norm;
    double& nx = norm[0];
    double& ny = norm[1];
    double& nz = norm[2];
  
    int num = 0;

    std::string line;
    while(std::getline(a_inputstream, line)){
      std::stringstream sstream(line);
      sstream >> x >> y >> z >> nx >> ny >> nz;

      a_vertices.emplace_back(std::make_shared<vertex>(pos, norm));

      // We have read all the vertices we should read. Exit now.
      num++;
      if(num == a_numVertices) break;
    }
  }

  template <class T>
  inline
  auto dcel::parser::PLY<T>::readFacesASCII(std::vector<std::shared_ptr<face> >&   a_faces,
					    std::vector<std::shared_ptr<edge> >&   a_edges,
					    std::vector<std::shared_ptr<vertex> >& a_vertices,
					    const int                              a_numFaces,
					    std::ifstream&                         a_inputstream) -> void {
    int numVertices;
    std::vector<int> vertexIndices;

    std::string line;
    int counter = 0;
    while(std::getline(a_inputstream, line)){
      counter++;
    
      std::stringstream sstream(line);

      sstream >> numVertices;
      vertexIndices.resize(numVertices);
      for (int i = 0; i < numVertices; i++){
	sstream >> vertexIndices[i];
      }

      if(numVertices < 3) std::cerr << "dcel::parser::PLY::readFacesASCII - a face must have at least three vertices!\n";
    
      // Get the vertices that make up this face. 
      std::vector<std::shared_ptr<vertex> > curVertices;
      for (int i = 0; i < numVertices; i++){
	const int vertexIndex = vertexIndices[i];
	curVertices.emplace_back(a_vertices[vertexIndex]);
      }

      // Build inside half edges and give each vertex an outgoing half edge. This may get overwritten later,
      // but the outgoing edge is not unique so it doesn't matter. 
      std::vector<std::shared_ptr<edge> > halfEdges;
      for (const auto& v : curVertices){
	halfEdges.emplace_back(std::make_shared<edge>(v));
	v->setEdge(halfEdges.back());
      }

      a_edges.insert(a_edges.end(), halfEdges.begin(), halfEdges.end());

      // Associate next/previous for the half edges inside the current face. Wish we had a circular iterator
      // but this will have to do. 
      for (int i = 0; i < halfEdges.size(); i++){
	auto& curEdge  = halfEdges[i];
	auto& nextEdge = halfEdges[(i+1)%halfEdges.size()];

	curEdge->setNextEdge(nextEdge);
	nextEdge->setPreviousEdge(curEdge);
      }

      // Construct a new face
      a_faces.emplace_back(std::make_shared<face>(halfEdges.front()));
      auto& curFace = a_faces.back();

      // Half edges get a reference to the currently created face
      for (auto& e : halfEdges){
	e->setFace(curFace);
      }

      // Must give vertices access to all faces associated with them since PLY files do not give any edge association. 
      for (auto& v : curVertices){
	v->addFaceToCache(curFace);
      }

    
      if(counter == a_numFaces) break;
    }
  }

  template <class T>
  inline
  auto parser::PLY<T>::reconcilePairEdges(std::vector<std::shared_ptr<edge> >& a_edges) -> void {
    for (auto& curEdge : a_edges){
      const auto& nextEdge = curEdge->getNextEdge();
    
      const auto& vertexStart = curEdge->getVertex();
      const auto& vertexEnd   = nextEdge->getVertex();

      for (const auto& p : vertexStart->getFaceCache()){
	for (edge_iterator edgeIt(*p); edgeIt.ok(); ++edgeIt){
	  const auto& polyCurEdge  = edgeIt();
	  const auto& polyNextEdge = polyCurEdge->getNextEdge();

	  const auto& polyVertexStart = polyCurEdge->getVertex();
	  const auto& polyVertexEnd   = polyNextEdge->getVertex();

	  if(vertexStart == polyVertexEnd && polyVertexStart == vertexEnd){ // Found the pair edge
	    curEdge->setPairEdge(polyCurEdge);
	    polyCurEdge->setPairEdge(curEdge);
	  }
	}
      }
    }
  }

  template <class T>
  inline
  auto dcel::parser::PLY<T>::clearFaceCache(std::vector<std::shared_ptr<vertex> >& a_vertices) -> void {
    for (auto& v : a_vertices){
      v->clearFaceCache();
    }
  }

}
#endif
