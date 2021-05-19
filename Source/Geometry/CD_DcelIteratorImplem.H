/*!
  @file   dcel_iteratorI.H
  @brief  Implementation of dcel_iterator.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_ITERATORI_H_
#define _DCEL_ITERATORI_H_

#include <CD_DcelVertex.H>
#include <CD_DcelEdge.H>
#include <CD_DcelFace.H>
#include <CD_DcelIterator.H>

#include "CD_NamespaceHeader.H"
namespace Dcel {

  template <class T>
  inline
  edgeIteratorT<T>::edgeIteratorT(Face& a_face){
    m_startEdge = a_face.getHalfEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;
  
    m_iterMode = IterationMode::Faces;
  }

  template <class T>
  inline
  edgeIteratorT<T>::edgeIteratorT(const Face& a_face){
    m_startEdge = a_face.getHalfEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;

    m_iterMode = IterationMode::Faces;
  }

  template <class T>
  inline
  edgeIteratorT<T>::edgeIteratorT(Vertex& a_vert){
    m_startEdge = a_vert.getOutgoingEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;

    m_iterMode  = IterationMode::Vertices;
  }

  template <class T>
  inline
  edgeIteratorT<T>::edgeIteratorT(const Vertex& a_vert){
    m_startEdge = a_vert.getOutgoingEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;

    m_iterMode  = IterationMode::Vertices;
  }

  template <class T>
  inline
  std::shared_ptr<EdgeT<T> >& edgeIteratorT<T>::operator() () noexcept {
    return (m_curEdge);
  }

  template <class T>
  inline
  const std::shared_ptr<EdgeT<T> >& edgeIteratorT<T>::operator() () const noexcept {
    return (m_curEdge);
  }

  template <class T>
  inline
  void edgeIteratorT<T>::reset() noexcept {
    m_curEdge  = m_startEdge;
    m_fullLoop = false;
  }

  template <class T>
  inline
  void edgeIteratorT<T>::operator++() noexcept {
    switch(m_iterMode){
    case IterationMode::Faces:
      m_curEdge = m_curEdge->getNextEdge();
      break;
    case IterationMode::Vertices:
      // For vertices, we want to compute the 
      m_curEdge = m_curEdge->getPreviousEdge()->getPairEdge();
      break;
    }

    m_fullLoop = (m_curEdge == m_startEdge);
  }

  template <class T>
  inline
  bool edgeIteratorT<T>::ok() const noexcept {
    return !m_fullLoop && m_curEdge;
  }
}

#include "CD_NamespaceFooter.H"

#endif
