/*!
  @file   dcel_poly.cpp
  @brief  Implementation of dcel_poly.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_vert.H"
#include "dcel_edge.H"
#include "dcel_iterator.H"
#include "dcel_poly.H"

#include <PolyGeom.H>

#define EPSILON 1.E-8
#define TWOPI 6.283185307179586476925287

dcel_poly::dcel_poly(){
  m_normal = RealVect::Zero;
}

dcel_poly::~dcel_poly(){

}

const RefCountedPtr<dcel_edge>& dcel_poly::get_edge() const{
  return m_edge;
}

RefCountedPtr<dcel_edge>& dcel_poly::get_edge(){
  return m_edge;
}

void dcel_poly::define(const RealVect a_normal, const RefCountedPtr<dcel_edge>& a_edge){
  this->set_normal(a_normal);
  this->set_edge(a_edge);
}

void dcel_poly::set_edge(const RefCountedPtr<dcel_edge>& a_edge){
  m_edge = a_edge;
}

void dcel_poly::set_normal(const RealVect a_normal){
  m_normal = a_normal;
}

void dcel_poly::normalize(){
  m_normal *= 1./m_normal.vectorLength();
}

void dcel_poly::compute_area() {
  const Vector<RefCountedPtr<dcel_vert> > vertices = this->get_vertices();

  Real area = 0.0;

  for (int i = 0; i < vertices.size() - 1; i++){
    const RealVect v1 = vertices[i]->get_pos();
    const RealVect v2 = vertices[i+1]->get_pos();
    area += PolyGeom::dot(PolyGeom::cross(v2,v1), m_normal);
  }

  m_area = Abs(0.5*area);
}

void dcel_poly::compute_centroid() {

  m_centroid = RealVect::Zero;
  const Vector<RefCountedPtr<dcel_vert> > vertices = this->get_vertices();

  for (int i = 0; i < vertices.size(); i++){
    m_centroid += vertices[i]->get_pos();
  }
  m_centroid = m_centroid/vertices.size();
}

void dcel_poly::compute_normal(const bool a_outward_normal){
  
  // We assume that the normal is defined by right-hand rule where the rotation direction is along the half edges
  const RefCountedPtr<dcel_vert>& v0 = m_edge->get_prev()->get_vert();
  const RefCountedPtr<dcel_vert>& v1 = m_edge->get_vert();
  const RefCountedPtr<dcel_vert>& v2 = m_edge->get_next()->get_vert();

  const RealVect x0 = v0->get_pos();
  const RealVect x1 = v1->get_pos();
  const RealVect x2 = v2->get_pos();
  
  m_normal = PolyGeom::cross(x2-x1,x1-x0);
  if(m_normal.vectorLength() < EPSILON){
    MayDay::Abort("dcel_poly::compute_normal - vertices lie on a line. Cannot compute normal vector");
  }
  else{
    m_normal = m_normal/m_normal.vectorLength();
  }

  if(!a_outward_normal){
    m_normal = -m_normal;
  }

}

Real dcel_poly::get_area() const{
  return m_area;
}

Real dcel_poly::signed_distance(const RealVect a_x0) {
  Real retval = 1.234567E89;

  Vector<RefCountedPtr<dcel_vert> > vertices = this->get_vertices();

  // Compute projection of x0 on the plane
  const RealVect x1 = vertices[0]->get_pos();
  const RealVect xp = a_x0 - PolyGeom::dot(m_normal, a_x0 - x1)*m_normal;

  // Use angle rule to check if projected point lies inside the polygon
  Real anglesum = 0.0;
  const int n = vertices.size();
  for(int i = 0; i < n; i++){
    const RealVect p1 = vertices[i]->get_pos() - xp;
    const RealVect p2 = vertices[(i+1)%n]->get_pos() - xp;

    const Real m1 = p1.vectorLength();
    const Real m2 = p2.vectorLength();

    if(m1*m2 <= EPSILON) {// Projected point hits a vertex or edge, return early. 
      anglesum = 0.;
      break;
    }
    else {
      const Real cos_theta = PolyGeom::dot(p1, p2)/(m1*m2);
      anglesum += acos(cos_theta);
    }
  }

  // Projected point is inside if angles sum to 2*pi
  bool inside = false;
  if(Abs(anglesum - TWOPI) < EPSILON){
    inside = true;
  }
#if 0
  inside = true;
#endif

  // If projection is inside, shortest distance is the normal component of the point
  if(inside){
    retval = PolyGeom::dot(a_x0-x1, m_normal);
  }
  else{ // The projected point lies outside the triangle. Check distance to edges/vertices
    const Vector<RefCountedPtr<dcel_edge> > edges = this->get_edges();
    for (int i = 0; i < edges.size(); i++){

      const Real cur_dist = edges[i]->signed_distance(a_x0);
      if(Abs(cur_dist) < Abs(retval)){
	retval = cur_dist;
      }
    }
  }

  return retval;
}

RealVect dcel_poly::get_normal() const {
  return m_normal;
}

RealVect dcel_poly::get_centroid() const {
  return m_centroid;
}

Vector<RefCountedPtr<dcel_vert> > dcel_poly::get_vertices(){
  Vector<RefCountedPtr<dcel_vert> > vertices;

  for (edge_iterator iter(*this); iter.ok(); ++iter){
    RefCountedPtr<dcel_edge>& edge = iter();
    vertices.push_back(edge->get_vert());
  }

  return vertices;
}

Vector<RefCountedPtr<dcel_edge> > dcel_poly::get_edges(){
  Vector<RefCountedPtr<dcel_edge> > edges;

  for (edge_iterator iter(*this); iter.ok(); ++iter){
    edges.push_back(iter());
  }

  return edges;
}
