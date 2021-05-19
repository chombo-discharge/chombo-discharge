/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DcelEdgeImplem.H
  @brief  Implementation of CD_DcelEdge.H
  @author Robert Marskar
*/

#ifndef CD_DcelEdgeImplem_H
#define CD_DcelEdgeImplem_H

// Our includes
#include <CD_DcelVertex.H>
#include <CD_DcelEdge.H>
#include <CD_DcelFace.H>
#include <CD_DcelIterator.H>
#include <CD_NamespaceHeader.H>

namespace Dcel {

  template <class T>
  inline
  edgeT<T>::edgeT(){
    m_vertex       = nullptr;
    m_pairEdge     = nullptr;
    m_nextEdge     = nullptr;
    m_previousEdge = nullptr;
    m_normal       = Vec3::zero();
    m_x2x1         = Vec3::zero();
    m_invLen2      = 0.0; 
  }

  template <class T>
  inline
  edgeT<T>::edgeT(const VertexPtr& a_vertex) : edgeT<T>() {
    m_vertex = a_vertex;
  }

  template <class T>
  inline
  edgeT<T>::edgeT(const edge& a_otherEdge) : edgeT<T>() {
    m_vertex       = a_otherEdge.m_vertex;
    m_pairEdge     = a_otherEdge.m_pairEdge;
    m_nextEdge     = a_otherEdge.m_nextEdge;
    m_previousEdge = a_otherEdge.m_previousEdge;
    m_normal       = a_otherEdge.m_normal;
    m_x2x1         = a_otherEdge.m_x2x1;
    m_invLen2      = a_otherEdge.m_invLen2;
  }

  template <class T>
  inline
  edgeT<T>::~edgeT(){
  }

  template <class T>
  inline
  void edgeT<T>::define(const VertexPtr& a_vertex,
			const edgePtr&   a_pairEdge,
			const edgePtr&   a_nextEdge,
			const edgePtr&   a_previousEdge,
			const Vec3       a_normal) noexcept {
    m_vertex       = a_vertex;
    m_pairEdge     = a_pairEdge;
    m_nextEdge     = a_nextEdge;
    m_previousEdge = a_previousEdge;
    m_normal       = a_normal;
  }

  template <class T>
  inline
  void edgeT<T>::setVertex(const VertexPtr& a_vertex) noexcept {
    m_vertex = a_vertex;
  }

  template <class T>
  inline
  void edgeT<T>::setPairEdge(const edgePtr& a_pairEdge) noexcept {
    m_pairEdge = a_pairEdge;
  }

  template <class T>
  inline
  void edgeT<T>::setNextEdge(const edgePtr& a_nextEdge) noexcept {
    m_nextEdge = a_nextEdge;
  }

  template <class T>
  inline
  void edgeT<T>::setPreviousEdge(const edgePtr& a_previousEdge) noexcept {
    m_previousEdge = a_previousEdge;
  }

  template <class T>
  inline
  void edgeT<T>::setFace(const FacePtr& a_face) noexcept {
    m_face = a_face;
  }

  template <class T>
  inline
  void edgeT<T>::normalizeNormalVector() noexcept {
    m_normal = m_normal/m_normal.length();
  }

  template <class T>
  inline
  void edgeT<T>::computeEdgeLength() noexcept {
    const auto& x1 = this->getVertex()->getPosition();
    const auto& x2 = this->getOtherVertex()->getPosition();

    m_x2x1 = x2-x1;

    const auto len2 = m_x2x1.dot(m_x2x1);

    m_invLen2 = 1./len2;
  }

  template <class T>
  inline
  void edgeT<T>::computeNormal() noexcept {

    m_normal = m_face->getNormal();

    if(m_pairEdge){
      m_normal += m_pairEdge->getFace()->getNormal();
    }

    this->normalizeNormalVector();
  }

  template <class T>
  inline
  void edgeT<T>::reconcile() noexcept {
    this->computeNormal();
    this->computeEdgeLength();
  }

  template <class T>
  inline
  std::shared_ptr<VertexT<T> >& edgeT<T>::getVertex() noexcept {
    return (m_vertex);
  }

  template <class T>
  inline
  const std::shared_ptr<VertexT<T> >& edgeT<T>::getVertex() const noexcept {
    return (m_vertex);
  }

  template <class T>
  inline
  std::shared_ptr<VertexT<T> >& edgeT<T>::getOtherVertex() noexcept {
    return (m_nextEdge->getVertex());
  }

  template <class T>
  inline
  const std::shared_ptr<VertexT<T> >& edgeT<T>::getOtherVertex() const noexcept {
    return (m_nextEdge->getVertex());
  }

  template <class T>
  inline
  std::shared_ptr<edgeT<T> >& edgeT<T>::getPairEdge() noexcept {
    return (m_pairEdge);
  }

  template <class T>
  inline
  const std::shared_ptr<edgeT<T> >& edgeT<T>::getPairEdge() const noexcept {
    return (m_pairEdge);
  }

  template <class T>
  inline
  std::shared_ptr<edgeT<T> >& edgeT<T>::getPreviousEdge() noexcept {
    return (m_previousEdge);
  }

  template <class T>
  inline
  const std::shared_ptr<edgeT<T> >& edgeT<T>::getPreviousEdge() const noexcept {
    return (m_previousEdge);
  }

  template <class T>
  inline
  std::shared_ptr<edgeT<T> >& edgeT<T>::getNextEdge() noexcept {
    return (m_nextEdge);
  }

  template <class T>
  inline
  const std::shared_ptr<edgeT<T> >& edgeT<T>::getNextEdge() const noexcept {
    return (m_nextEdge);
  }

  template <class T>
  inline
  Vec3T<T>& edgeT<T>::getNormal() noexcept {
    return (m_normal);
  }

  template <class T>
  inline
  const Vec3T<T>& edgeT<T>::getNormal() const noexcept {
    return (m_normal);
  }

  template <class T>
  inline
  std::shared_ptr<FaceT<T> >& edgeT<T>::getFace() noexcept {
    return (m_face);
  }

  template <class T>
  inline
  const std::shared_ptr<FaceT<T> >& edgeT<T>::getFace() const noexcept {
    return (m_face);
  }

  template <class T>
  inline
  T edgeT<T>::projectPointToEdge(const Vec3& a_x0) const noexcept {
    const auto p = a_x0 - m_vertex->getPosition();

    return p.dot(m_x2x1)*m_invLen2;
  }

  template <class T>
  inline
  T edgeT<T>::signedDistance(const Vec3& a_x0) const noexcept {
    const T t = this->projectPointToEdge(a_x0);

    T retval;
    if(t <= 0.0) {
      retval = this->getVertex()->signedDistance(a_x0);
    }
    else if(t >= 1.0){
      retval = this->getOtherVertex()->signedDistance(a_x0);
    }
    else{
      const Vec3 linePoint = m_vertex->getPosition() + t*m_x2x1;
      const Vec3 delta     = a_x0 - linePoint;
      const T dot          = m_normal.dot(delta);

      const int sgn = (dot > 0.0) ? 1 : -1;

      retval = sgn*delta.length();
    }

    return retval;
  }

  template <class T>
  inline
  T edgeT<T>::unsignedDistance2(const Vec3& a_x0) const noexcept {
    T t = this->projectPointToEdge(a_x0);

    constexpr T zero = 0.0;
    constexpr T one  = 1.0;
    
    t = std::min(std::max(zero,t), one);

    const Vec3T<T> linePoint = m_vertex->getPosition() + t*m_x2x1;
    const Vec3T<T> delta     = a_x0 - linePoint;
    
    return delta.dot(delta);
  }
}

#include <CD_NamespaceFooter.H>

#endif
