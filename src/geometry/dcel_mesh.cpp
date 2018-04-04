/*!
  @file   dcel_mesh.cpp
  @brief  Implementation of dcel_mesh.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_mesh.H"

dcel_mesh::dcel_mesh(){

}

dcel_mesh::~dcel_mesh(){

}

dcel_mesh::dcel_mesh(Vector<dcel_poly*> a_polygons, Vector<dcel_edge*> a_edges, Vector<dcel_vert*> a_vertices){
  this->define(a_polygons, a_edges, a_vertices);
}

void dcel_mesh::define(Vector<dcel_poly*> a_polygons, Vector<dcel_edge*> a_edges, Vector<dcel_vert*> a_vertices){
  m_polygons = a_polygons;
  m_edges    = a_edges;
  m_vertices = a_vertices;
}

Real dcel_mesh::signed_distance(const RealVect a_x0){

  Real min_dist = 1.E99;

  // This is a very slow version of doing this; this should be accelerated by using a BVH-tree or kD-tree, but that'll have to wait
  // until we get the basics right. 
  for (int i = 0; i < m_polygons.size(); i++){
    const Real cur_dist = m_polygons[i]->signed_distance(a_x0);
    if(Abs(cur_dist) < Abs(min_dist)){
      min_dist = cur_dist;
    }
  }

  return min_dist;
}
