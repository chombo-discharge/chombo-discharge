/*!
  @file   dcel_edge.cpp
  @brief  Implementation of dcel_edge.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_polygon.H"
#include "dcel_iterator.H"

#include "PolyGeom.H"

using namespace dcel;

edge::edge(){

}

edge::~edge(){

}


void edge::define(const std::shared_ptr<vertex>& a_vert,
		       const std::shared_ptr<edge>& a_pair,
		       const std::shared_ptr<edge>& a_next,
		       const std::shared_ptr<edge>& a_prev,
		       const RealVect                    a_normal){
  this->set_vert(a_vert);
  this->set_pair(a_pair);
  this->set_next(a_next);
  this->set_prev(a_prev);
  this->set_normal(a_normal);
}


void edge::set_poly(const std::shared_ptr<polygon>& a_poly){
  m_poly = a_poly;
}


void edge::set_vert(const std::shared_ptr<vertex>& a_vert){
  m_vert = a_vert;
}


void edge::set_pair(const std::shared_ptr<edge>& a_pair){
  m_pair = a_pair;
}


void edge::set_next(const std::shared_ptr<edge>& a_next){
  m_next = a_next;
}


void edge::set_prev(const std::shared_ptr<edge>& a_prev){
  m_prev = a_prev;
}


void edge::set_normal(const RealVect a_normal){
  m_normal = a_normal;
}


const std::shared_ptr<vertex>& edge::get_vert() const {
  return m_vert;
}


std::shared_ptr<vertex>& edge::get_vert() {
  return m_vert;
}


const std::shared_ptr<vertex>& edge::get_other_vert() const {
  return m_pair->get_vert();
}


std::shared_ptr<vertex>& edge::get_other_vert(){
  return m_pair->get_vert();
}


const std::shared_ptr<edge>& edge::get_pair() const {
  return m_pair;
}


std::shared_ptr<edge>& edge::get_pair(){
  return m_pair;
}


const std::shared_ptr<edge>& edge::get_prev() const {
  return m_prev;
}


std::shared_ptr<edge>& edge::get_prev() {
  return m_prev;
}


const std::shared_ptr<edge>& edge::get_next() const {
  return m_next;
}


std::shared_ptr<edge>& edge::get_next(){
  return m_next;
}


const std::shared_ptr<polygon>& edge::get_poly() const{
  return m_poly;
}


std::shared_ptr<polygon>& edge::get_poly(){
  return m_poly;
}


RealVect edge::normal() const {
  return m_normal;
}


Real edge::signed_distance(const RealVect a_x0) const {
  Real retval = 1.234567E89;
  
  // Involved vertices
  const RealVect x1 = this->get_other_vert()->position();
  const RealVect x2 = m_vert->position();

  const RealVect R  = PolyGeom::cross(x2-x1,PolyGeom::cross(x1-a_x0, x2-a_x0));
  const RealVect r  = R/R.vectorLength();
  const RealVect xp = a_x0 + PolyGeom::dot(x1-a_x0,r)*r;
  const Real t      = PolyGeom::dot(xp-x1,x2-x1)/PolyGeom::dot(x2-x1, x2-x1);


  RealVect p;
  RealVect n;
  if(t < 0.0){ // Closest to x1, vertex normal takes precedence.
    p = x1;
    n = this->get_other_vert()->normal();
  }
  else if (t > 1.0){ // Closest to x2, vertex normal takes precedence
    p = x2;
    n = m_vert->normal();
  }
  else{ // Projection onto line lies on the line segment
    p = xp;
    n = m_normal;
  }

  const Real dot = PolyGeom::dot(n, (a_x0 - p)); // Determine sign from projection. If the point is orthogonal to the normal,
  const int sgn = dot >= 0.0 ? 1 : -1;           // it must (I think) be outside

  retval = (a_x0-p).vectorLength()*sgn;

  return retval;
}

