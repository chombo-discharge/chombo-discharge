/*!
  @file   dcel_edgeI.cpp
  @brief  Implementation of dcel_edge.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_EDGEI_H_
#define _DCEL_EDGEI_H_

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"
#include "dcel_iterator.H"

namespace dcel {

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
  edgeT<T>::edgeT(const vertexPtr& a_vertex) : edgeT<T>() {
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
  decltype(auto) edgeT<T>::define(const vertexPtr& a_vertex,
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
  decltype(auto) edgeT<T>::setVertex(const vertexPtr& a_vertex) noexcept {
    m_vertex = a_vertex;
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::setPairEdge(const edgePtr& a_pairEdge) noexcept {
    m_pairEdge = a_pairEdge;
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::setNextEdge(const edgePtr& a_nextEdge) noexcept {
    m_nextEdge = a_nextEdge;
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::setPreviousEdge(const edgePtr& a_previousEdge) noexcept {
    m_previousEdge = a_previousEdge;
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::setFace(const facePtr& a_face) noexcept {
    m_face = a_face;
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::normalizeNormalVector() noexcept {
    m_normal = m_normal/m_normal.length();
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::computeEdgeLength() noexcept {
    const auto& x1 = this->getVertex()->getPosition();
    const auto& x2 = this->getOtherVertex()->getPosition();

    m_x2x1 = x2-x1;

    const auto len2 = m_x2x1.dot(m_x2x1);

    m_invLen2 = 1./len2;
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::computeNormal() noexcept {
    const auto& pairEdge = this->getPairEdge();

    const auto& F     = this->getFace();
    const auto& pairF = pairEdge->getFace();
    
    const auto& n1 = F->getNormal();
    const auto& n2 = pairF->getNormal();

    m_normal = n1 + n2;

    this->normalizeNormalVector();
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::reconcile() noexcept {
    this->computeNormal();
    this->computeEdgeLength();
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getVertex() noexcept {
    return (m_vertex);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getVertex() const noexcept {
    return (m_vertex);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getOtherVertex() noexcept {
    return (m_pairEdge->getVertex());
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getOtherVertex() const noexcept {
    return (m_pairEdge->getVertex());
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getPairEdge() noexcept {
    return (m_pairEdge);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getPairEdge() const noexcept {
    return (m_pairEdge);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getPreviousEdge() noexcept {
    return (m_previousEdge);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getPreviousEdge() const noexcept {
    return (m_previousEdge);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getNextEdge() noexcept {
    return (m_nextEdge);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getNextEdge() const noexcept {
    return (m_nextEdge);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getNormal() noexcept {
    return (m_normal);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getNormal() const noexcept {
    return (m_normal);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getFace() noexcept {
    return (m_face);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getFace() const noexcept {
    return (m_face);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::getFlippedEdge() const noexcept {

    edgeT<T> flippedEdge = edgeT<T>();

    flippedEdge.m_vertex       = this->getOtherVertex();
    flippedEdge.m_pairEdge     = m_pairEdge;
    flippedEdge.m_previousEdge = m_nextEdge;
    flippedEdge.m_nextEdge     = m_previousEdge;
    flippedEdge.m_normal       = m_normal;
    flippedEdge.m_x2x1         = -m_x2x1;
    flippedEdge.m_invLen2      = m_invLen2;
    flippedEdge.m_face         = m_face;

    return std::make_shared<edgeT<T> >(flippedEdge);
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::projectPointToEdge(const Vec3& a_x0) const noexcept {
    const auto p = a_x0 - m_vertex->getPosition();

    return p.dot(m_x2x1)*m_invLen2;
  }

  template <class T>
  inline
  decltype(auto) edgeT<T>::signedDistance(const Vec3& a_x0) const noexcept {
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
  decltype(auto) edgeT<T>::unsignedDistance2(const Vec3& a_x0) const noexcept {
    T t = this->projectPointToEdge(a_x0);

    constexpr T zero = 0.0;
    constexpr T one  = 1.0;
    
    t = std::min(std::max(zero,t), one);

    const Vec3T<T> linePoint = m_vertex->getPosition() + t*m_x2x1;
    const Vec3T<T> delta     = a_x0 - linePoint;
    
    return delta.dot(delta);
  }
}

#endif
