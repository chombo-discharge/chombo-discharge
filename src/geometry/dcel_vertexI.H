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
    m_pos    = Vec3::Zero;
    m_normal = Vec3::Zero;

    m_faceCache.resize(0);
  }

  template <class T>
  inline
  vertexT<T>::vertexT(const Vec3& a_pos){
    m_pos    = a_pos;
    m_normal = Vec3::Zero;
  }

  template <class T>
  inline
  vertexT<T>::vertexT(const Vec3& a_pos, const Vec3& a_normal){
    m_pos    = a_pos;
    m_normal = a_normal;
  }

  template <class T>
  inline
  vertexT<T>::vertexT(const vertexT<T>& a_otherVertex){
    this->define(a_otherVertex.getPosition(),
		 a_otherVertex.getEdge(),
		 a_otherVertex.getNormal());
  }

  template <class T>
  inline
  vertexT<T>::~vertexT(){

  }

  template <class T>
  inline
  auto vertexT<T>::define(const Vec3& a_pos, const std::shared_ptr<edge>& a_edge, const Vec3 a_normal) noexcept -> void {
    m_pos = a_pos;
    m_edge = a_edge;
    m_normal = a_normal;
  }

  template <class T>
  inline
  auto vertexT<T>::setPosition(const Vec3& a_pos) noexcept -> void {
    m_pos = a_pos;
  }

  template <class T>
  inline
  auto vertexT<T>::setEdge(const std::shared_ptr<edge>& a_edge) noexcept -> void {
    m_edge = a_edge;
  }

  template <class T>
  inline
  auto vertexT<T>::setNormal(const Vec3& a_normal) noexcept -> void {
    m_normal = a_normal;
  }

  template <class T>
  inline
  auto vertexT<T>::addFaceToCache(const std::shared_ptr<face>& a_face) noexcept -> void {
    m_faceCache.push_back(a_face);
  }

  template <class T>
  inline
  auto vertexT<T>::clearFaceCache() noexcept -> void {
    m_faceCache.resize(0);
  }

  template <class T>
  inline
  auto vertexT<T>::normalizeNormalVector() noexcept -> void {
    m_normal = m_normal/m_normal.length();
  }

  template <class T>
  inline
  auto vertexT<T>::computeVertexNormalAverage() noexcept -> void {
    std::vector<std::shared_ptr<face> > faces;
    for (edge_iterator iter(*this); iter.ok(); ++iter){
      faces.emplace_back(iter()->getFace());
    }

    m_normal = Vec3::Zero;
  
    for (const auto& f : faces){
      m_normal += f->getNormal();
    }

    this->normalizeNormalVector();
  }

  template <class T>
  inline
  auto vertexT<T>::computeVertexNormalAngleWeighted() noexcept -> void {
    m_normal = Vec3::Zero;

    for (edge_iterator iter(*this); iter.ok(); ++iter){
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
  auto vertexT<T>::getPosition() noexcept -> Vec3& {
    return m_pos;
  }

  template <class T>
  inline
  auto vertexT<T>::getPosition() const noexcept -> const Vec3& {
    return m_pos;
  }

  template <class T>
  inline
  auto vertexT<T>::getNormal() noexcept -> Vec3& {
    return m_normal;
  }

  template <class T>
  inline
  auto vertexT<T>::getNormal() const noexcept -> const Vec3& {
    return m_normal;
  }

  template <class T>
  inline
  auto vertexT<T>::getEdge() noexcept -> std::shared_ptr<edge>& {
    return m_edge;
  }

  template <class T>
  inline
  auto vertexT<T>::getEdge() const noexcept -> const std::shared_ptr<edge>& {
    return m_edge;
  }

  template <class T>
  inline
  auto vertexT<T>::getFaceCache() const noexcept -> const std::vector<std::shared_ptr<face> >& {
    return m_faceCache;
  }

  template <class T>
  inline
  auto vertexT<T>::getFaceCache() noexcept -> std::vector<std::shared_ptr<face> >& {
    return m_faceCache;
  }

  template <class T>
  inline
  auto vertexT<T>::signedDistance(const Vec3& a_x0) const noexcept -> T{
    const Vec3 delta = a_x0 - m_pos;
    const T dist        = delta.length();
    const T dot         = m_normal.dot(delta);
    const int sign           = (dot > 0.) ? 1 : -1;
  
    return dist*sign;
  }

  template <class T>
  inline
  auto vertexT<T>::unsignedDistance2(const Vec3& a_x0) const noexcept -> T {
    const Vec3 d = a_x0 - m_pos;

    return d.dot(d);
  }
}

#endif
