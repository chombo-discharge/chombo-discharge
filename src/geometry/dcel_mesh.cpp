/*!
  @file   dcel_mesh.cpp
  @brief  Implementation of dcel_mesh.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_mesh.H"
#include "dcel_iterator.H"

#include <PolyGeom.H>

#if 1
#include <ParmParse.H>
#endif

bool dcel_mesh::s_angle_weighted = false; // Use angle-weighted vertex normal vectors

dcel_mesh::dcel_mesh(){
  m_reconciled = false;
  m_use_tree   = false;
}

dcel_mesh::~dcel_mesh(){

}

dcel_mesh::dcel_mesh(Vector<RefCountedPtr<dcel_poly> >& a_polygons,
		     Vector<RefCountedPtr<dcel_edge> >& a_edges,
		     Vector<RefCountedPtr<dcel_vert> >& a_vertices){
  m_reconciled = false;
  m_use_tree   = false;
  
  this->define(a_polygons, a_edges, a_vertices);
}

Vector<RefCountedPtr<dcel_vert> >& dcel_mesh::get_vertices(){
  return m_vertices;
}

Vector<RefCountedPtr<dcel_edge> >& dcel_mesh::get_edges(){
  return m_edges;
}

Vector<RefCountedPtr<dcel_poly> >& dcel_mesh::get_polygons(){
  return m_polygons;
}

bool dcel_mesh::sanity_check() const {
  for (int i = 0; i < m_edges.size(); i++){
    if(m_edges[i].isNull()){
      MayDay::Abort("dcel_mesh::sanity_check - edge is NULL");
    }
    else{
      const RefCountedPtr<dcel_edge>& edge = m_edges[i];
    
      if(edge->get_pair().isNull()){
	MayDay::Abort("dcel_mesh::sanity_check - pair edge is NULL");
      }
      else if(edge->get_next().isNull()){
	MayDay::Abort("dcel_mesh::sanity_check - next edge is NULL");
      }
      else if(edge->get_prev().isNull()){
	MayDay::Abort("dcel_mesh::sanity_check - prev edge is NULL");
      }
      else if(edge->get_vert().isNull()){
	MayDay::Abort("dcel_mesh::sanity_check - vertex is NULL");
      }
    }
  }

  for (int i = 0; i < m_vertices.size(); i++){
    if(m_vertices[i].isNull()){
      MayDay::Abort("dcel_mesh::sanity_check - m_vertices[i] is NULL");
    }
    else{
      if(m_vertices[i]->get_edge().isNull()){
	MayDay::Abort("dcel_mesh::sanity_check - vertex edge is NULL");
      }
    }
  }
}

void dcel_mesh::define(Vector<RefCountedPtr<dcel_poly> >& a_polygons,
		       Vector<RefCountedPtr<dcel_edge> >& a_edges,
		       Vector<RefCountedPtr<dcel_vert> >& a_vertices){
  m_polygons = a_polygons;
  m_edges    = a_edges;
  m_vertices = a_vertices;
}

void dcel_mesh::compute_bounding_sphere(){
  Vector<RealVect> pos;
  for (int i = 0; i < m_vertices.size(); i++){
    pos.push_back(m_vertices[i]->get_pos());
  }
  
  m_sphere.define(pos);
}

void dcel_mesh::reconcile_polygons(const bool a_outward_normal){

  // Reconcile polygons; compute polygon area and provide edges explicit access
  // to their polygons
  for (int i = 0; i < m_polygons.size(); i++){
    RefCountedPtr<dcel_poly>& poly = m_polygons[i];

    // Every edge gets a reference to this polygon
    for (edge_iterator iter(*poly); iter.ok(); ++iter){
      RefCountedPtr<dcel_edge>& edge = iter();
      edge->set_poly(poly);
    }
    poly->compute_normal(a_outward_normal);
    poly->compute_centroid();
    poly->compute_area();
    poly->normalize();
    poly->compute_bbox();
  }


  // Compute pseudonormals for vertices and edges. 
  this->compute_vertex_normals();
  this->compute_edge_normals();
  this->compute_bounding_sphere();

  m_reconciled = true;
}

void dcel_mesh::compute_vertex_normals(){
  for (int i = 0; i < m_vertices.size(); i++){
#if 0 // This doesn't work, why?!?
    const Vector<RefCountedPtr<dcel_poly> > polygons = m_vertices[i]->get_polygons();
#else
    const Vector<RefCountedPtr<dcel_poly> > polygons = m_vertices[i]->get_polycache();
#endif

    // Mean or area weighted
    if(!s_angle_weighted){
      RealVect normal = RealVect::Zero;
      for (int i = 0; i < polygons.size(); i++){
	//normal += polygons[i]->get_area()*polygons[i]->get_normal(); // Area weighted
	normal += polygons[i]->get_normal(); // Mean
      }
      normal *= 1./normal.vectorLength();
      m_vertices[i]->set_normal(normal);
    }
    else { // Angle-weighted normal vector
      RealVect normal = RealVect::Zero;
      for (edge_iterator iter(*m_vertices[i]); iter.ok(); ++iter){
	const RefCountedPtr<dcel_edge>& outgoing = iter();
	const RefCountedPtr<dcel_edge>& incoming = outgoing->get_prev();

	const RealVect origin = incoming->get_vert()->get_pos();
	const RealVect x2     = outgoing->get_vert()->get_pos();
	const RealVect x1     = incoming->get_other_vert()->get_pos();
	const Real len1       = (x1-origin).vectorLength();
	const Real len2       = (x2-origin).vectorLength();

	const RealVect norm = PolyGeom::cross(x2-origin, x1-origin)/(len1*len2);
	const Real alpha = asin(norm.vectorLength());

	normal += alpha*norm;
      
      }
      normal *= 1./normal.vectorLength();

      m_vertices[i]->set_normal(normal);
    }
  }
}

void dcel_mesh::compute_edge_normals(){
  for (int i = 0; i < m_edges.size(); i++){

    RefCountedPtr<dcel_edge>& cur_edge         = m_edges[i];
    const RefCountedPtr<dcel_edge>& pair_edge = cur_edge->get_pair();

    const RefCountedPtr<dcel_poly>& poly      = cur_edge->get_poly();
    const RefCountedPtr<dcel_poly>& pair_poly = pair_edge->get_poly();

    const RealVect n1 = poly->get_normal();
    const RealVect n2 = pair_poly->get_normal();
    const RealVect n  = (n1+n2)/(n1+n2).vectorLength();

    cur_edge->set_normal(n);
  }
}

void dcel_mesh::build_tree(){
  m_tree     = RefCountedPtr<kd_tree<dcel_poly> > (new kd_tree<dcel_poly>(m_polygons, 20, 1));
  m_use_tree = true;
}

Real dcel_mesh::signed_distance(const RealVect a_x0){
  CH_assert(m_reconciled);
  
  Real min_dist = 1.E99;


  if(m_sphere.inside(a_x0)){ // Bounding sphere contains point
    if(m_use_tree){ // Fast kd-tree search
      Vector<RefCountedPtr<dcel_poly> > candidates = m_tree->get_candidates(a_x0);

      if(candidates.size() > 0){
	for (int i = 0; i < candidates.size(); i++){
	  const Real cur_dist = candidates[i]->signed_distance(a_x0);
	  if(Abs(cur_dist) < Abs(min_dist)){
	    min_dist = cur_dist;
	  }
	}
      }
      else { // I think we can use any triangle for this
	//	return m_polygons[0]->signed_distance(a_x0);
      }
    }
    else{ // Brute force search
      for (int i = 0; i < m_polygons.size(); i++){
	const Real cur_dist = m_polygons[i]->signed_distance(a_x0);
	if(Abs(cur_dist) < Abs(min_dist)){
	  min_dist = cur_dist;
	}
      }
    }
  }
  else{ // We are outside every bounding box, simply return the distance to the bounding sphere
    min_dist = (a_x0 - m_sphere.get_center()).vectorLength() - m_sphere.get_radius();
  }

  return min_dist;
}
