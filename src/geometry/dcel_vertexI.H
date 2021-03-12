/*!
  @file   dcel_vertexI.H
  @brief  Implementaion of dcel_vertex.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_VERTEXI_
#define _DCEL_VERTEXI_

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"
#include "dcel_iterator.H"

namespace dcel {

  template <class T>
  inline
  vertexT<T>::vertexT(){
    m_pos    = Vec3::zero();
    m_normal = Vec3::zero();

    m_faceCache.resize(0, nullptr);
  }

  template <class T>
  inline
  vertexT<T>::vertexT(const Vec3& a_pos) : vertexT() {
    m_pos    = a_pos;
  }

  template <class T>
  inline
  vertexT<T>::vertexT(const Vec3& a_pos, const Vec3& a_normal) : vertexT() {
    m_pos    = a_pos;
    m_normal = a_normal;
  }

  template <class T>
  inline
  vertexT<T>::vertexT(const vertexT<T>& a_otherVertex){
    m_pos    = a_otherVertex.m_pos;
    m_normal = a_otherVertex.m_m_normal;
    m_edge   = a_otherVertex.m_edge;
  }

  template <class T>
  inline
  vertexT<T>::~vertexT(){

  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::define(const Vec3& a_pos, const std::shared_ptr<edge>& a_edge, const Vec3& a_normal) noexcept {
    m_pos    = a_pos;
    m_edge   = a_edge;
    m_normal = a_normal;
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::setPosition(const Vec3& a_pos) noexcept {
    m_pos = a_pos;
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::setEdge(const std::shared_ptr<edge>& a_edge) noexcept {
    m_edge = a_edge;
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::setNormal(const Vec3& a_normal) noexcept {
    m_normal = a_normal;
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::addFaceToCache(const std::shared_ptr<face>& a_face) noexcept {
    m_faceCache.emplace_back(a_face);
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::clearFaceCache() noexcept {
    m_faceCache.resize(0);
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::normalizeNormalVector() noexcept {
    m_normal = m_normal/m_normal.length();
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::computeVertexNormalAverage() noexcept {
    for (edgeIterator iter(*this); iter.ok(); ++iter){
      m_normal += iter()->getNormal();
    }

    this->normalizeNormalVector();
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::computeVertexNormalAngleWeighted() noexcept {
    m_normal = Vec3::zero();

    for (edgeIterator iter(*this); iter.ok(); ++iter){
      const auto& outgoingEdge = iter();

      const Vec3& x0 = outgoingEdge->getVertex()->getPosition();
      const Vec3& x1 = outgoingEdge->getPreviousEdge()->getVertex()->getPosition();
      const Vec3& x2 = outgoingEdge->getNextEdge()->getVertex()->getPosition();

      Vec3 v1 = x1-x0;
      Vec3 v2 = x2-x0;

      v1 = v1/v1.length();
      v2 = v2/v2.length();

      const Vec3 norm = outgoingEdge->getFace()->getNormal();

      const T alpha = acos(v1.dot(v2));

      m_normal += alpha*norm;
    }

    this->normalizeNormalVector();
  }
  
  template <class T>
  inline
  decltype(auto) vertexT<T>::getPosition() noexcept {
    return (m_pos);
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::getPosition() const noexcept {
    return (m_pos);
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::getNormal() noexcept {
    return (m_normal);
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::getNormal() const noexcept {
    return (m_normal);
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::getEdge() noexcept {
    return (m_edge);
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::getEdge() const noexcept {
    return (m_edge);
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::getFaceCache() noexcept {
    return (m_faceCache);
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::getFaceCache() const noexcept {
    return (m_faceCache);
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::signedDistance(const Vec3& a_x0) const noexcept {
    const auto delta = a_x0 - m_pos;
    const T dist     = delta.length();
    const T dot      = m_normal.dot(delta);
    const int sign   = (dot > 0.) ? 1 : -1;
  
    return dist*sign;
  }

  template <class T>
  inline
  decltype(auto) vertexT<T>::unsignedDistance2(const Vec3& a_x0) const noexcept {
    const auto d = a_x0 - m_pos;

    return d.dot(d);
  }
}

#endif
