/*!
  @file   dcel_mesh.cpp
  @brief  Implementation of dcel_mesh.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_mesh.H"
#include "dcel_iterator.H"

#if 1
#include <ParmParse.H>
#endif

dcel_mesh::dcel_mesh(){
  m_reconciled = false;
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

void dcel_mesh::reconcile_polygons(){

  // Reconcile polygons; compute polygon area and provide edges explicit access
  // to their polygons
  for (int i = 0; i < m_polygons.size(); i++){
    dcel_poly* poly = m_polygons[i];

    // Every edge gets a reference to this polygon
    for (edge_iterator iter(const_cast<dcel_poly*> (poly)); iter.ok(); ++iter){
      dcel_edge* edge = const_cast<dcel_edge*> (iter());
      edge->set_poly(poly);
    }
    poly->compute_area();
    poly->normalize();
  }

  // Compute pseudonormals for vertices and edges. 
  this->compute_vertex_normals();
  this->compute_edge_normals();

  m_reconciled = true;
}

void dcel_mesh::compute_vertex_normals(){
  for (int i = 0; i < m_vertices.size(); i++){
    Vector<const dcel_poly*> polygons = m_vertices[i]->get_polygons();

    CH_assert(polygons.size() == 3);

    // Compute area-weighted normal vector
    RealVect normal = RealVect::Zero;
    for (int i = 0; i < polygons.size(); i++){
#if DCEL_ANGLE_WEIGHTED
      normal += polygons[i]->get_area()*polygons[i]->get_normal();
#else
      normal += polygons[i]->get_normal();
#endif
    }
    normal *= 1./normal.vectorLength();

    //    pout() << normal << endl;

    m_vertices[i]->set_normal(normal);
  }
}

void dcel_mesh::compute_edge_normals(){
  for (int i = 0; i < m_edges.size(); i++){

    const dcel_edge* edge      = m_edges[i];
    const dcel_edge* pair_edge = edge->get_pair();

    const dcel_poly* poly      = edge->get_poly();
    const dcel_poly* pair_poly = pair_edge->get_poly();

    const RealVect n1 = poly->get_normal();
    const RealVect n2 = pair_poly->get_normal();
    const RealVect n  = (n1+n2)/(n1+n2).vectorLength();

    dcel_edge* e = const_cast<dcel_edge*> (edge);
    e->set_normal(n);
  }
}


Real dcel_mesh::signed_distance(const RealVect a_x0){
  CH_assert(m_reconciled);
  
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
