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

  inline
  edge::edge(){
  }

  inline
  edge::edge(const std::shared_ptr<vertex>& a_vertex){
    this->setVertex(a_vertex);
    this->setPairEdge(nullptr);
    this->setNextEdge(nullptr);
    this->setPreviousEdge(nullptr);
  }

  inline
  edge::edge(const edge& a_otherEdge){
    this->define(a_otherEdge.getVertex(),
		 a_otherEdge.getPairEdge(),
		 a_otherEdge.getNextEdge(),
		 a_otherEdge.getPreviousEdge(),
		 a_otherEdge.getNormal());

    m_len2 = a_otherEdge.m_len2;
    m_x2x1 = a_otherEdge.m_x2x1;
  }

  inline
  edge::~edge(){
  }

  inline
  void edge::define(const std::shared_ptr<vertex>& a_vertex,
		    const std::shared_ptr<edge>&   a_pairEdge,
		    const std::shared_ptr<edge>&   a_nextEdge,
		    const std::shared_ptr<edge>&   a_previousEdge,
		    const Vec3<double>                 a_normal){
    this->setVertex(a_vertex);
    this->setPairEdge(a_pairEdge);
    this->setNextEdge(a_nextEdge);
    this->setPreviousEdge(a_previousEdge);
  }

  inline
  void edge::setVertex(const std::shared_ptr<vertex>& a_vertex) noexcept {
    m_vertex = a_vertex;
  }

  inline
  void edge::setPairEdge(const std::shared_ptr<edge>& a_pairEdge) noexcept {
    m_pairEdge = a_pairEdge;
  }

  inline
  void edge::setNextEdge(const std::shared_ptr<edge>& a_nextEdge) noexcept {
    m_nextEdge = a_nextEdge;
  }

  inline
  void edge::setPreviousEdge(const std::shared_ptr<edge>& a_previousEdge) noexcept {
    m_previousEdge = a_previousEdge;
  }

  inline
  void edge::setFace(const std::shared_ptr<face>& a_face) noexcept {
    m_face = a_face;
  }

  inline
  void edge::normalizeNormalVector() noexcept {
    m_normal = m_normal/m_normal.length();
  }

  inline
  void edge::computeEdgeLength() noexcept {
    const Vec3<double>& x1 = this->getVertex()->getPosition();
    const Vec3<double>& x2 = this->getOtherVertex()->getPosition();

    m_x2x1 = x2-x1;
    m_len2 = m_x2x1.dot(m_x2x1);
  }

  inline
  void edge::computeNormal() noexcept {
    const std::shared_ptr<edge>& pairEdge = this->getPairEdge();

    const std::shared_ptr<face>& F     = this->getFace();
    const std::shared_ptr<face>& pairF = pairEdge->getFace();
    
    const Vec3<double>& n1 = F->getNormal();
    const Vec3<double>& n2 = pairF->getNormal();

    m_normal = n1 + n2;

    this->normalizeNormalVector();
  }

  inline
  void edge::reconcile() noexcept {
    this->computeNormal();
    this->computeEdgeLength();
  }

  inline
  std::shared_ptr<vertex>& edge::getVertex() noexcept {
    return m_vertex;
  }

  inline
  const std::shared_ptr<vertex>& edge::getVertex() const noexcept {
    return m_vertex;
  }

  inline
  std::shared_ptr<vertex>& edge::getOtherVertex() noexcept{
    return m_pairEdge->getVertex();
  }

  inline
  const std::shared_ptr<vertex>& edge::getOtherVertex() const noexcept {
    return m_pairEdge->getVertex();
  }

  inline
  std::shared_ptr<edge>& edge::getPairEdge() noexcept {
    return m_pairEdge;
  }

  inline
  const std::shared_ptr<edge>& edge::getPairEdge() const noexcept {
    return m_pairEdge;
  }

  inline
  std::shared_ptr<edge>& edge::getPreviousEdge() noexcept {
    return m_previousEdge;
  }

  inline
  const std::shared_ptr<edge>& edge::getPreviousEdge() const noexcept {
    return m_previousEdge;
  }

  inline
  std::shared_ptr<edge>& edge::getNextEdge() noexcept {
    return m_nextEdge;
  }

  inline
  const std::shared_ptr<edge>& edge::getNextEdge() const noexcept {
    return m_nextEdge;
  }

  inline
  Vec3<double>& edge::getNormal() noexcept {
    return m_normal;
  }

  inline
  const Vec3<double>& edge::getNormal() const noexcept {
    return m_normal;
  }

  inline
  std::shared_ptr<face>& edge::getFace() noexcept {
    return m_face;
  }

  inline
  const std::shared_ptr<face>& edge::getFace() const noexcept {
    return m_face;
  }

  inline
  double edge::projectPointToEdge(const Vec3<double>& a_x0) const noexcept {
    const Vec3<double> p = a_x0 - m_vertex->getPosition();

    return p.dot(m_x2x1)/m_len2;
  }

  inline
  double edge::signedDistance(const Vec3<double>& a_x0) const noexcept {
    const double t = this->projectPointToEdge(a_x0);

    double retval;
  
    if(t <= 0.0) {
      retval = this->getVertex()->signedDistance(a_x0);
    }
    else if(t >= 1.0){
      retval = this->getOtherVertex()->signedDistance(a_x0);
    }
    else{
      const Vec3<double> linePoint = m_vertex->getPosition() + t*m_x2x1;
      const Vec3<double> delta     = a_x0 - linePoint;
      const double dot             = m_normal.dot(delta);

      // I don't know why, but currently classifying this as "inside" works, this shouldn't matter because the polygon
      // test should trigger in that case. 
      const int sgn = (dot > 0.0) ? 1 : -1;

      retval = sgn*delta.length();
    }

    return retval;
  }

  inline
  double edge::unsignedDistance2(const Vec3<double>& a_x0) const noexcept {
    double t = this->projectPointToEdge(a_x0);
    t = std::min(std::max(0.,t), 1.);

    const Vec3<double> linePoint = m_vertex->getPosition() + t*m_x2x1;

    const Vec3<double> d = a_x0 - linePoint;
    
    return d.dot(d);
  }
}

#endif
