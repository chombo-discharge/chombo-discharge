/*!
  @file   dcel_iteratorI.H
  @brief  Implementation of dcel_iterator.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_ITERATORI_H_
#define _DCEL_ITERATORI_H_

#include "dcel_iterator.H"
#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"

namespace dcel {

  inline
  edge_iterator::edge_iterator(face& a_face){
    m_startEdge = a_face.getHalfEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;
  
    m_iterMode = IterationMode::Face;
  }

  inline
  edge_iterator::edge_iterator(const face& a_face){
    m_startEdge = a_face.getHalfEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;

    m_iterMode = IterationMode::Face;
  }

  inline
  edge_iterator::edge_iterator(vertex& a_vert){
    m_startEdge = a_vert.getEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;

    m_iterMode  = IterationMode::Vertex;
  }

  inline
  edge_iterator::edge_iterator(const vertex& a_vert){
    m_startEdge = a_vert.getEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;

    m_iterMode  = IterationMode::Vertex;
  }

  inline
  std::shared_ptr<edge>& edge_iterator::operator() () noexcept{
    return m_curEdge;
  }

  inline
  const std::shared_ptr<edge>& edge_iterator::operator() () const noexcept{
    return m_curEdge;
  }

  inline
  void edge_iterator::reset() noexcept {
    m_curEdge  = m_startEdge;
    m_fullLoop = false;
  }

  inline
  void edge_iterator::operator++() noexcept {
    switch(m_iterMode){
    case IterationMode::Face:
      m_curEdge = m_curEdge->getNextEdge();
      break;
    case IterationMode::Vertex:
      m_curEdge = m_curEdge->getPreviousEdge()->getPairEdge();
      break;
    }

    m_fullLoop = (m_curEdge == m_startEdge);
  }

  inline
  bool edge_iterator::ok() const noexcept{
    return !m_fullLoop;
  }
}

#endif
