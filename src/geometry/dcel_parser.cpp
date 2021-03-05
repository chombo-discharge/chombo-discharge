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

void dcel::parser::PLY::read_ascii(dcel::dcel_mesh& a_mesh, const std::string a_filename){
  std::ifstream filestream(a_filename);

  if(filestream.is_open()){
    std::vector<std::shared_ptr<dcel::vertex> >& vertices = a_mesh.get_vertices();
    std::vector<std::shared_ptr<dcel::edge> >& edges    = a_mesh.get_edges();
    std::vector<std::shared_ptr<dcel::polygon> >& polygons = a_mesh.get_polygons();

    vertices.resize(0);
    edges.resize(0);
    polygons.resize(0);

    int num_vertices;  // Number of vertices
    int num_polygons;  // Number of polygons

    dcel::parser::PLY::read_ascii_header(num_vertices, num_polygons, filestream); 
    dcel::parser::PLY::read_ascii_vertices(vertices, num_vertices, filestream);
    dcel::parser::PLY::read_ascii_polygons(polygons, edges, vertices, num_polygons, filestream);

    a_mesh.sanity_check();
  
    filestream.close();
  }
  else{
    const std::string error = "dcel::parser::PLY::read_ascii - ERROR! Could not open file " + a_filename;
    MayDay::Abort(error.c_str());
  }
}

void dcel::parser::PLY::read_ascii_header(int&           a_num_vertices,
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

void dcel::parser::PLY::read_ascii_vertices(std::vector<std::shared_ptr<dcel::vertex> >& a_vertices,
					    const int                                 a_num_vertices,
					    std::ifstream&                            a_inputstream){

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
    num++;
    std::stringstream sstream(line);
    sstream >> x >> y >> z >> nx >> ny >> nz;

    std::shared_ptr<dcel::vertex> vert = std::shared_ptr<dcel::vertex> (new dcel::vertex(pos));
    vert->set_normal(norm);
    a_vertices.push_back(vert);
    if(num == a_num_vertices){
      break;
    }
  } 
}

void dcel::parser::PLY::read_ascii_polygons(std::vector<std::shared_ptr<dcel::polygon> >& a_polygons,
					    std::vector<std::shared_ptr<dcel::edge> >& a_edges,
					    std::vector<std::shared_ptr<dcel::vertex> >& a_vertices,
					    const int a_num_polygons,
					    std::ifstream& a_inputstream){
  int num_vert;
  std::vector<int> which_vertices;

  std::string line;
  int counter = 0;
  while(std::getline(a_inputstream, line)){
    counter++;
    
    std::stringstream sstream(line);

    sstream >> num_vert;
    which_vertices.resize(num_vert);
    for (int i = 0; i < num_vert; i++){
      sstream >> which_vertices[i];
    }

    // Build polygon and inside edges
    std::shared_ptr<dcel::polygon> polygon = std::shared_ptr<dcel::polygon> (new dcel::polygon());

    // Get vertices. Add a reference to the newly created polygon
    std::vector<std::shared_ptr<dcel::vertex> > poly_vertices(num_vert);
    for (int i = 0; i < num_vert; i++){
      poly_vertices[i] = a_vertices[which_vertices[i]];
    }

    // Build inside edges. Polygon gets a reference to the edge
    std::vector<std::shared_ptr<dcel::edge> > poly_edges(num_vert);
    for (int i = 0; i < num_vert; i++){
      poly_edges[i] = std::shared_ptr<dcel::edge> (new dcel::edge());
      poly_edges[i]->set_vert(poly_vertices[(i+1)%num_vert]);
    }
    polygon->set_edge(poly_edges[0]);

    // Associate prev/next
    for (int i = 0; i < num_vert; i++){
      poly_edges[i]->set_next(poly_edges[(i+1)%num_vert]);
      poly_edges[(i+1)%num_vert]->set_prev(poly_edges[i]);
    }

    // Set edges emanating from vertices if that hasn't been done already
    for (int i = 0; i < poly_vertices.size(); i++){
      if(poly_vertices[i]->get_edge() == nullptr){
	poly_vertices[i]->set_edge(poly_edges[i]);
      }
    }

    // Check for pairs
    for (int i = 0; i < poly_edges.size(); i++){

      std::shared_ptr<dcel::edge>& edge = poly_edges[i];
      std::shared_ptr<dcel::vertex>& vert = edge->get_vert();

      // Get all polygons connected to the current vertex and look for edge pairs
      std::vector<std::shared_ptr<dcel::polygon> >& polygons = vert->get_polycache();

      for (int j = 0; j < polygons.size(); j++){
	std::shared_ptr<dcel::edge>& other_polygon_edge = polygons[j]->get_edge();

	for (dcel::edge_iterator iter(*polygons[j]); iter.ok(); ++iter){
	  std::shared_ptr<dcel::edge>& other_polygon_edge = iter();

	  if(other_polygon_edge->get_vert() == edge->get_prev()->get_vert()){
	    edge->set_pair(other_polygon_edge);
	    other_polygon_edge->set_pair(edge);
	  }
	}
      }
    }

    // Add reference to newly created polygon
    for (int i = 0; i < poly_vertices.size(); i++){
      poly_vertices[i]->add_polygon(polygon);
      CH_assert(!poly_vertices[i]->get_edge().isNull());
    }

    // Add edges and polygons
    for (int i = 0; i < poly_edges.size(); i++){
      a_edges.push_back(poly_edges[i]);
    }
    a_polygons.push_back(polygon);
    
    if(counter == a_num_polygons){
      break;
    }
  }
}
