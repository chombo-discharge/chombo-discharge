/*!
  @file   dcel_iterator.cpp
  @brief  Implementation of dcel_iterator.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_iterator.H"
#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"

using namespace dcel;

edge_iterator::edge_iterator(face& a_face){
  m_startEdge = a_face.getHalfEdge();
  m_curEdge   = m_startEdge;
  m_fullLoop  = false;
  
  m_iterMode = IterationMode::Face;
}

edge_iterator::edge_iterator(const face& a_face){
  m_startEdge = a_face.getHalfEdge();
  m_curEdge   = m_startEdge;
  m_fullLoop  = false;

  m_iterMode = IterationMode::Face;
}

edge_iterator::edge_iterator(vertex& a_vert){
  m_startEdge = a_vert.getEdge();
  m_curEdge   = m_startEdge;
  m_fullLoop  = false;

  m_iterMode = IterationMode::Vertex;
}

edge_iterator::edge_iterator(const vertex& a_vert){
  m_startEdge = a_vert.getEdge();
  m_curEdge   = m_startEdge;
  m_fullLoop  = false;

  m_iterMode = IterationMode::Vertex;
}

std::shared_ptr<edge>& edge_iterator::operator() (){
  return m_curEdge;
}

void edge_iterator::reset(){
  m_curEdge  = m_startEdge;
  m_fullLoop = false;
}

void edge_iterator::operator++(){
  switch(m_iterMode){
  case IterationMode::Face:
    m_curEdge = m_curEdge->getNextEdge();
    break;
  case IterationMode::Vertex:
    m_curEdge = m_curEdge->getPreviousEdge()->getPairEdge();
    break;
  }

  // Check if we have done a full loop.
  if(m_curEdge == m_startEdge){ 
    m_fullLoop = true;
  }
}

bool edge_iterator::ok(){
  if(m_curEdge != nullptr && !m_fullLoop){
    return true;
  }
  else {
    return false;
  }
}
