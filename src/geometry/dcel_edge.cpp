/*!
  @file   dcel_edge.cpp
  @brief  Implementation of dcel_edge.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_vert.H"
#include "dcel_edge.H"
#include "dcel_poly.H"
#include "dcel_iterator.H"
#include "PolyGeom.H"

dcel_edge::dcel_edge(){

}

dcel_edge::~dcel_edge(){

}



void dcel_edge::define(const dcel_vert* const a_vert,
		       const dcel_edge* const a_pair,
		       const dcel_edge* const a_next,
		       const dcel_edge* const a_prev,
		       const RealVect         a_normal){
  m_vert   = a_vert;
  m_pair   = a_pair;
  m_next   = a_next;
  m_prev   = a_prev;
  m_normal = a_normal;
}

void dcel_edge::set_poly(const dcel_poly* const a_poly){
  m_poly = a_poly;
}

void dcel_edge::set_normal(const RealVect a_normal){
  m_normal = a_normal;
}

const dcel_vert* dcel_edge::get_vert() const {
  return m_vert;
}

const dcel_vert* dcel_edge::get_other_vert() const {
  return m_pair->get_vert();
}

const dcel_edge* dcel_edge::get_pair() const {
  return m_pair;
}

const dcel_edge* dcel_edge::get_prev() const {
  return m_prev;
}

const dcel_edge* dcel_edge::get_next() const {
  return m_next;
}

const dcel_poly* dcel_edge::get_poly() const{
  return m_poly;
}

RealVect dcel_edge::get_normal() const {
  return m_normal;
}

Real dcel_edge::signed_distance(const RealVect a_x0) const {

  Real retval = 1.234567E89;
  
  // Involved vertices
  const RealVect x1 = m_vert->get_pos();
  const RealVect x2 = this->get_other_vert()->get_pos();

  const RealVect R  = PolyGeom::cross(PolyGeom::cross(x1-a_x0, x2-a_x0), x2-x1);
  const RealVect r  = R/R.vectorLength();
  const RealVect xp = a_x0 + PolyGeom::dot(x1-a_x0,r)*r;
  const Real t      = PolyGeom::dot(xp-x1,x2-x1)/PolyGeom::dot(x2-x1, x2-x1);

  RealVect p;
  RealVect n;
  if(t <= 0.0){ // Closest to x1, vertex normal takes precedence.
    p = x1;
    n = m_vert->get_normal();
  }
  else if (t >= 1.0){ // Closest to x2, call signed distance function for that vertex
    p = x2;
    n = this->get_other_vert()->get_normal();
  }
  else{ // Projection onto line lies on the line segment
    p = xp;
    n = m_normal;
  }

  const Real dot = PolyGeom::dot(n, (a_x0 - p)); // Determine sign from projection. However, I don't think this is a good way of
  const int sgn = (dot > 0.0) - (dot < 0.0);     // doing it since the normal vector might be orthogonal to x0-p. 

  //  retval = (a_x0-xp).vectorLength()*sgn;
  retval = (a_x0-xp).vectorLength();

  return retval;
}
