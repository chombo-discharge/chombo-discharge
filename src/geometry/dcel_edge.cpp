/*!
  @file   dcel_edge.cpp
  @brief  Implementation of dcel_edge.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_vert.H"
#include "dcel_edge.H"
#include "dcel_poly.H"
#include "PolyGeom.H"

dcel_edge::dcel_edge(){

}

dcel_edge::~dcel_edge(){

}

edge_iterator::edge_iterator(){

}

edge_iterator::edge_iterator(const dcel_poly* const a_poly){
  m_polymode  = true;
  m_begin     = a_poly->get_edge();
  m_current   = m_begin;
  m_full_loop = false;

  //  pout() << "begin vert = " << m_begin->get_vert()->get_pos() << endl;
}

edge_iterator::edge_iterator(const dcel_vert* const a_vert){
  m_polymode  = false;
  m_begin     = a_vert->get_edge();
  m_current   = m_begin;
  m_full_loop = false;
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

const dcel_edge* edge_iterator::operator() (){
  return m_current;
}

RealVect dcel_edge::get_normal() const {
  return m_normal;
}

Real dcel_edge::signed_distance(const RealVect a_x0) const {
  // Involved vertices
  const RealVect x1 = m_vert->get_pos();
  const RealVect x2 = this->get_other_vert()->get_pos();

  // Must now find the shortest distance from x0 to the line connecting x1,x2
  const RealVect t = (x2-x1)/(x2-x1).vectorLength();           // Line tangent
  const RealVect d = (x1-a_x0) - t*PolyGeom::dot(x1-a_x0, t); // Component of vector to point that is normal to the line
  const Real ret   = PolyGeom::dot(d, -m_normal);              // Project with normal vector

  return ret;
}

void edge_iterator::reset(){
  m_current = m_begin;
}

void edge_iterator::operator++(){
  if(m_polymode){
    m_current = m_current->get_next();
  }
  else{
    m_current = m_current->get_next()->get_pair();
  }

  // Signal a full loop around the polygon. 
  if(m_current == m_begin){ 
    m_full_loop = true;
  }
}

bool edge_iterator::ok(){
  if(m_current != NULL && !m_full_loop){
    return true;
  }
  else {
    return false;
  }
}
