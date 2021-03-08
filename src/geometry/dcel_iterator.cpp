/*!
  @file   dcel_iterator.cpp
  @brief  Implementation of dcel_iterator.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_iterator.H"
#include "dcel_polygon.H"
#include "dcel_edge.H"
#include "dcel_vertex.H"

using namespace dcel;

edge_iterator::edge_iterator(){

}

edge_iterator::edge_iterator(polygon& a_poly){
  m_polymode  = true;
  m_begin     = a_poly.getEdge();
  m_current   = m_begin;
  m_full_loop = false;
}

edge_iterator::edge_iterator(const polygon& a_poly){
  m_polymode  = true;
  m_begin     = a_poly.getEdge();
  m_current   = m_begin;
  m_full_loop = false;
}

edge_iterator::edge_iterator(vertex& a_vert){
  m_polymode  = false;
  m_begin     = a_vert.getEdge();
  m_current   = m_begin;
  m_full_loop = false;
}

edge_iterator::edge_iterator(const vertex& a_vert){
  m_polymode  = false;
  m_begin     = a_vert.getEdge();
  m_current   = m_begin;
  m_full_loop = false;
}

std::shared_ptr<edge>& edge_iterator::operator() (){
  return m_current;
}


void edge_iterator::reset(){
  m_current = m_begin;
}


void edge_iterator::operator++(){
  if(m_polymode){
    m_current = m_current->getNextEdge();
  }
  else{
    m_current = m_current->getPreviousEdge()->getPairEdge();
  }

  // Have now done a full loop around the polygon or vertex. 
  if(m_current == m_begin){ 
    m_full_loop = true;
  }
}


bool edge_iterator::ok(){
  if(m_current != nullptr && !m_full_loop){
    return true;
  }
  else {
    return false;
  }
}

