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

  inline
  vertex::vertex(){
    m_pos    = Vec3::Zero;
    m_normal = Vec3::Zero;

    m_faceCache.resize(0);
  }

  inline
  vertex::vertex(const Vec3& a_pos){
    m_pos    = a_pos;
    m_normal = Vec3::Zero;
  }

  inline
  vertex::vertex(const Vec3& a_pos, const Vec3& a_normal){
    m_pos    = a_pos;
    m_normal = a_normal;
  }

  inline
  vertex::vertex(const vertex& a_otherVertex){
    this->define(a_otherVertex.getPosition(),
		 a_otherVertex.getEdge(),
		 a_otherVertex.getNormal());
  }

  inline
  vertex::~vertex(){

  }

  inline
  void vertex::define(const Vec3& a_pos, const std::shared_ptr<edge>& a_edge, const Vec3 a_normal) noexcept {
    m_pos = a_pos;
    m_edge = a_edge;
    m_normal = a_normal;
  }

  inline
  void vertex::setPosition(const Vec3& a_pos) noexcept {
    m_pos = a_pos;
  }

  inline
  void vertex::setEdge(const std::shared_ptr<edge>& a_edge) noexcept {
    m_edge = a_edge;
  }

  inline
  void vertex::setNormal(const Vec3& a_normal) noexcept {
    m_normal = a_normal;
  }

  inline
  void vertex::addFaceToCache(const std::shared_ptr<face>& a_face) noexcept {
    m_faceCache.push_back(a_face);
  }

  inline
  void vertex::clearFaceCache() noexcept {
    m_faceCache.resize(0);
  }

  inline
  void vertex::normalizeNormalVector() noexcept {
    m_normal = m_normal/m_normal.length();
  }

  inline
  void vertex::computeVertexNormalAverage() noexcept {
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

  inline
  void vertex::computeVertexNormalAngleWeighted() noexcept {
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

      const double alpha = acos(v1.dot(v2));

      m_normal += alpha*norm;
    }

    this->normalizeNormalVector();
  }

  inline
  Vec3& vertex::getPosition() noexcept {
    return m_pos;
  }

  inline
  const Vec3& vertex::getPosition() const noexcept {
    return m_pos;
  }

  inline
  Vec3& vertex::getNormal() noexcept {
    return m_normal;
  }

  inline
  const Vec3& vertex::getNormal() const noexcept {
    return m_normal;
  }

  inline
  std::shared_ptr<edge>& vertex::getEdge() noexcept {
    return m_edge;
  }

  inline
  const std::shared_ptr<edge>& vertex::getEdge() const noexcept {
    return m_edge;
  }

  inline
  const std::vector<std::shared_ptr<face> >& vertex::getFaceCache() const noexcept{
    return m_faceCache;
  }

  inline
  std::vector<std::shared_ptr<face> >& vertex::getFaceCache() noexcept {
    return m_faceCache;
  }

  inline
  double vertex::signedDistance(const Vec3& a_x0) const noexcept {
    const Vec3 delta = a_x0 - m_pos;
    const double dist        = delta.length();
    const double dot         = m_normal.dot(delta);
    const int sign           = (dot > 0.) ? 1 : -1;
  
    return dist*sign;
  }

  inline
  double vertex::unsignedDistance2(const Vec3& a_x0) const noexcept {
    const Vec3 d = a_x0 - m_pos;

    return d.dot(d);
  }
}

#endif
