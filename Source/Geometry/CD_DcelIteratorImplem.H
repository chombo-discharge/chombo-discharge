/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DcelIteratorImplem.H
  @brief  Implementation of CD_DcelIterator.H
  @author Robert Marskar

*/

#ifndef CD_DcelIteratorImplem_H
#define CD_DcelIteratorImplem_H

// Our includes
#include <CD_DcelVertex.H>
#include <CD_DcelEdge.H>
#include <CD_DcelFace.H>
#include <CD_DcelIterator.H>
#include <CD_NamespaceHeader.H>

namespace Dcel {

  template <class T>
  inline
  EdgeIteratorT<T>::EdgeIteratorT(Face& a_face){
    m_startEdge = a_face.getHalfEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;
    m_iterMode  = IterationMode::Faces;
  }

  template <class T>
  inline
  EdgeIteratorT<T>::EdgeIteratorT(const Face& a_face){
    m_startEdge = a_face.getHalfEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;
    m_iterMode  = IterationMode::Faces;
  }

  template <class T>
  inline
  EdgeIteratorT<T>::EdgeIteratorT(Vertex& a_vert){
    m_startEdge = a_vert.getOutgoingEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;
    m_iterMode  = IterationMode::Vertices;
  }

  template <class T>
  inline
  EdgeIteratorT<T>::EdgeIteratorT(const Vertex& a_vert){
    m_startEdge = a_vert.getOutgoingEdge();
    m_curEdge   = m_startEdge;
    m_fullLoop  = false;
    m_iterMode  = IterationMode::Vertices;
  }

  template <class T>
  inline
  std::shared_ptr<EdgeT<T> >& EdgeIteratorT<T>::operator() () noexcept {
    return (m_curEdge);
  }

  template <class T>
  inline
  const std::shared_ptr<EdgeT<T> >& EdgeIteratorT<T>::operator() () const noexcept {
    return (m_curEdge);
  }

  template <class T>
  inline
  void EdgeIteratorT<T>::reset() noexcept {
    m_curEdge  = m_startEdge;
    m_fullLoop = false;
  }

  template <class T>
  inline
  void EdgeIteratorT<T>::operator++() noexcept {
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
  bool EdgeIteratorT<T>::ok() const noexcept {
    return !m_fullLoop && m_curEdge;
  }
}

#include <CD_NamespaceFooter.H>

#endif
