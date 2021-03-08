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

edge_iterator::edge_iterator(polygon& a_poly){
  m_startEdge = a_poly.getEdge();
  m_curEdge   = m_startEdge;
  m_full_loop = false;
  
  m_iterMode = IterationMode::Polygon;
}

edge_iterator::edge_iterator(const polygon& a_poly){
  m_startEdge = a_poly.getEdge();
  m_curEdge   = m_startEdge;
  m_full_loop = false;

  m_iterMode = IterationMode::Polygon;
}

edge_iterator::edge_iterator(vertex& a_vert){
  m_startEdge = a_vert.getEdge();
  m_curEdge   = m_startEdge;
  m_full_loop = false;

  m_iterMode = IterationMode::Vertex;
}

edge_iterator::edge_iterator(const vertex& a_vert){
  m_startEdge = a_vert.getEdge();
  m_curEdge   = m_startEdge;
  m_full_loop = false;

  m_iterMode = IterationMode::Vertex;
}

std::shared_ptr<edge>& edge_iterator::operator() (){
  return m_curEdge;
}

void edge_iterator::reset(){
  m_curEdge = m_startEdge;
}

void edge_iterator::operator++(){
  switch(m_iterMode){
  case IterationMode::Polygon:
    m_curEdge = m_curEdge->getNextEdge();
    break;
  case IterationMode::Vertex:
    m_curEdge = m_curEdge->getPreviousEdge()->getPairEdge();
    break;
  }

  // Check if we have done a full loop.
  if(m_curEdge == m_startEdge){ 
    m_full_loop = true;
  }
}

bool edge_iterator::ok(){
  if(m_curEdge != nullptr && !m_full_loop){
    return true;
  }
  else {
    return false;
  }
}
