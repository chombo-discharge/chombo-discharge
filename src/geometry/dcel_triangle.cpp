/*!
  @file   dcel_triangle.cpp
  @brief  Declaration of a polygon class for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_triangle.H"
#include "dcel_edge.H"
#include "dcel_vert.H"

#include <PolyGeom.H>

#if 1
#include <ParmParse.H>
#endif

dcel_triangle::dcel_triangle(){

}

dcel_triangle::~dcel_triangle(){

}

void dcel_triangle::define(const RealVect a_normal, const dcel_edge* const a_edge){
  dcel_poly::define(a_normal, a_edge);

  Vector<const dcel_vert*> vertices = this->get_vertices();
  m_x1 = vertices[0]->get_pos();
  m_x2 = vertices[1]->get_pos();
  m_x3 = vertices[2]->get_pos();
}

void dcel_triangle::compute_area(){
  const RealVect v = PolyGeom::cross(m_x3-m_x1, m_x3-m_x2);
  m_area = 0.5*v.vectorLength();
}

Real dcel_triangle::signed_distance(const RealVect a_x0) const {
  Real retval = 1.234567E89;

  // Projection of xp onto the triangle plane
  const RealVect xp = a_x0 - PolyGeom::dot(m_normal, a_x0-m_x1)*m_normal;

  // Outside test. This is taken from W. Heidrich 2005, Computing the Barycentric Coordinates of a Projected Point
  bool inside = false;
  const RealVect u = m_x2 - m_x1;
  const RealVect v = m_x3 - m_x1;
  const RealVect w = a_x0 - m_x1;
  const RealVect n = PolyGeom::cross(u,v);
  const Real alpha = PolyGeom::dot(PolyGeom::cross(u,w), n)/PolyGeom::dot(n, n);
  const Real beta  = PolyGeom::dot(PolyGeom::cross(w,v), n)/PolyGeom::dot(n, n);
  const Real gamma = 1.0 - alpha - beta;
  if(alpha >= 0.0 && alpha <= 1.0){
    if(beta >= 0.0 && beta <= 1.0){
      if(gamma >= 0.0 && gamma <= 1.0){
	inside = true;
      }
    }
  }

  if(inside){ // Distance is just the normal component of a_x0
    retval = PolyGeom::dot(a_x0-m_x1, m_normal);
  }
  else{ // The projected point lies outside the triangle. Check distance to edges/vertices
    Vector<const dcel_edge*> edges = this->get_edges();
    for (int i = 0; i < edges.size(); i++){
      const Real cur_dist = edges[i]->signed_distance(a_x0);
      if(Abs(cur_dist) < Abs(retval)){
	retval = cur_dist;
      }
    }
  }
  
  return retval;
}
