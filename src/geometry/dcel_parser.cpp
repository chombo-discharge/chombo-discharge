/*!
  @file   dcel_parser.cpp
  @brief  Implementation of dcel_parser.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_parser.H"
#include "dcel_mesh.H"
#include "dcel_polygon.H"
#include "dcel_edge.H"
#include "dcel_vertex.H"
#include "dcel_iterator.H"

#include <iostream>
#include <fstream>
#include <iterator>

void dcel::parser::PLY::readASCII(dcel::mesh& a_mesh, const std::string a_filename){
  std::ifstream filestream(a_filename);

  if(filestream.is_open()){
    std::vector<std::shared_ptr<dcel::vertex> >& vertices  = a_mesh.getVertices();
    std::vector<std::shared_ptr<dcel::edge> >& edges       = a_mesh.getEdges();
    std::vector<std::shared_ptr<dcel::polygon> >& polygons = a_mesh.getPolygons();

    vertices.resize(0);
    edges.resize(0);
    polygons.resize(0);

    int num_vertices;  // Number of vertices
    int num_polygons;  // Number of polygons

    dcel::parser::PLY::readHeaderASCII(num_vertices, num_polygons, filestream); 
    dcel::parser::PLY::readVerticesASCII(vertices, num_vertices, filestream);
    dcel::parser::PLY::readPolygonsASCII(polygons, edges, vertices, num_polygons, filestream);

    a_mesh.sanityCheck();
  
    filestream.close();
  }
  else{
    const std::string error = "dcel::parser::PLY::readASCII - ERROR! Could not open file " + a_filename;
    std::cerr << error + "\n";
  }
}

void dcel::parser::PLY::readHeaderASCII(int&           a_num_vertices,
					int&           a_num_polygons,
					std::ifstream& a_inputstream){

  std::string str1;
  std::string str2;
  std::string line;

  // Get number of vertices
  a_inputstream.clear();
  a_inputstream.seekg(0);
  while (getline(a_inputstream, line)){
    std::stringstream sstream(line);
    sstream >> str1 >> str2 >> a_num_vertices;
    if(str1 == "element" && str2 == "vertex"){
      break;
    }
  }

  // Get number of polygons
  a_inputstream.clear();
  a_inputstream.seekg(0);
  while (getline(a_inputstream, line)){
    std::stringstream sstream(line);
    sstream >> str1 >> str2 >> a_num_polygons;
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

void dcel::parser::PLY::readVerticesASCII(std::vector<std::shared_ptr<dcel::vertex> >& a_vertices,
					  const int                                    a_num_vertices,
					  std::ifstream&                               a_inputstream){

  RealVect pos;
  Real& x = pos[0];
  Real& y = pos[1];
  Real& z = pos[2];

  RealVect norm;
  Real& nx = norm[0];
  Real& ny = norm[1];
  Real& nz = norm[2];
  
  int num = 0;

  std::string line;
  while(std::getline(a_inputstream, line)){
    std::stringstream sstream(line);
    sstream >> x >> y >> z >> nx >> ny >> nz;

    a_vertices.emplace_back(std::make_shared<dcel::vertex>(pos, norm));

    // We have read all the vertices we should read. Exit now.
    num++;
    if(num == a_num_vertices) break;
  }
}

void dcel::parser::PLY::readPolygonsASCII(std::vector<std::shared_ptr<dcel::polygon> >& a_polygons,
					  std::vector<std::shared_ptr<dcel::edge> >&    a_edges,
					  std::vector<std::shared_ptr<dcel::vertex> >&  a_vertices,
					  const int                                     a_num_polygons,
					  std::ifstream&                                a_inputstream){
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
    
    // Get the vertices that make up this polygon. 
    std::vector<std::shared_ptr<dcel::vertex> > curVertices;
    for (int i = 0; i < numVertices; i++){
      const int vertexIndex = vertexIndices[i];
      curVertices.emplace_back(a_vertices[vertexIndex]);
    }

    // Build inside half edges and give each vertex an outgoing half edge. This may get overwritten later,
    // but the outgoing edge is not unique so it doesn't matter. 
    std::vector<std::shared_ptr<dcel::edge> > halfEdges;
    for (const auto& v : curVertices){
      halfEdges.emplace_back(std::make_shared<dcel::edge>(v));
      v->setEdge(halfEdges.back());
    }

    a_edges.insert(a_edges.end(), halfEdges.begin(), halfEdges.end());

    // Associate next/previous for the half edges inside the current polygon. Wish we had a circular iterator
    // but this will have to do. 
    for (int i = 0; i < halfEdges.size(); i++){
      auto& curEdge  = halfEdges[i];
      auto& nextEdge = halfEdges[(i+1)%halfEdges.size()];

      curEdge->setNextEdge(nextEdge);
      nextEdge->setPreviousEdge(curEdge);
    }

    // Construct a new polygon
    a_polygons.emplace_back(std::make_shared<dcel::polygon>(halfEdges.front()));
    auto& curPolygon = a_polygons.back();

    // Half edges get a reference to the currently created polygon
    for (auto& e : halfEdges){
      e->setPolygon(curPolygon);
    }

    // Check for pairs
    for (int i = 0; i < halfEdges.size(); i++){

      std::shared_ptr<dcel::edge>& edge   = halfEdges[i];
      std::shared_ptr<dcel::vertex>& vert = edge->getVertex();

      // Get all polygons connected to the current vertex and look for edge pairs
      std::vector<std::shared_ptr<dcel::polygon> >& polygons = vert->getPolycache();

      for (int j = 0; j < polygons.size(); j++){
	std::shared_ptr<dcel::edge>& other_polygon_edge = polygons[j]->getEdge();

	for (dcel::edge_iterator iter(*polygons[j]); iter.ok(); ++iter){
	  std::shared_ptr<dcel::edge>& other_polygon_edge = iter();

	  if(other_polygon_edge->getVertex() == edge->getPreviousEdge()->getVertex()){
	    edge->setPairEdge(other_polygon_edge);
	    other_polygon_edge->setPairEdge(edge);
	  }
	}
      }
    }

    // Must give vertices access to all polygons associated with them since PLY files do not give any edge association. 
    for (auto& v : curVertices){
      v->addPolygonToCache(curPolygon);
    }


    // Add edges and polygons
    for (int i = 0; i < halfEdges.size(); i++){
      //      halfEdges[i]->setPolygon(polygon);
      //      a_edges.push_back(halfEdges[i]);
    }
    //    a_polygons.push_back(polygon);
    
    if(counter == a_num_polygons) break;
  }
}
