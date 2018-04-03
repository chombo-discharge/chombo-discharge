/*!
  @file   dcel_edge.cpp
  @brief  Implementation of dcel_edge.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_vert.H"
#include "dcel_edge.H"
#include "dcel_poly.H"

dcel_edge::dcel_edge(){

}

dcel_edge::~dcel_edge(){

}

edge_iterator::edge_iterator(){

}

edge_iterator::edge_iterator(const dcel_poly* const a_poly){
  m_polymode = false;
  m_begin    = a_poly->get_edge();
  m_current  = m_begin;
}

edge_iterator::edge_iterator(const dcel_vert* const a_vert){
  m_polymode = true;
  m_begin    = a_vert->get_edge();
  m_current  = m_begin;
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

  const RealVect x1 = m_vert->get_pos();
  const RealVect x2 = this->get_other_vert()->get_pos();

  // Must now find the shortest distance from x0 to the line connecting x1,x2
  
  return 0.0; 
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
}

bool edge_iterator::ok(){
  if(m_current != NULL && m_current->get_next() != m_begin){
    return true;
  }
  else {
    return false;
  }
}
