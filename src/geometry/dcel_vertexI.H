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
    m_pos    = RealVect::Zero;
    m_normal = RealVect::Zero;

    m_faceCache.resize(0);
  }

  inline
  vertex::vertex(const RealVect& a_pos){
    m_pos    = a_pos;
    m_normal = RealVect::Zero;
  }

  inline
  vertex::vertex(const RealVect& a_pos, const RealVect& a_normal){
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
  void vertex::define(const RealVect& a_pos, const std::shared_ptr<edge>& a_edge, const RealVect a_normal) noexcept {
    m_pos = a_pos;
    m_edge = a_edge;
    m_normal = a_normal;
  }

  inline
  void vertex::setPosition(const RealVect& a_pos) noexcept {
    m_pos = a_pos;
  }

  inline
  void vertex::setEdge(const std::shared_ptr<edge>& a_edge) noexcept {
    m_edge = a_edge;
  }

  inline
  void vertex::setNormal(const RealVect& a_normal) noexcept {
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
    m_normal = m_normal/m_normal.vectorLength();
  }

  inline
  void vertex::computeVertexNormalAverage() noexcept {
    std::vector<std::shared_ptr<face> > faces;
    for (edge_iterator iter(*this); iter.ok(); ++iter){
      faces.emplace_back(iter()->getFace());
    }

    m_normal = RealVect::Zero;
  
    for (const auto& f : faces){
      m_normal += f->getNormal();
    }

    this->normalizeNormalVector();
  }

  inline
  void vertex::computeVertexNormalAngleWeighted() noexcept {
    m_normal = RealVect::Zero;

    for (edge_iterator iter(*this); iter.ok(); ++iter){
      const auto& outgoingEdge = iter();

      // Edges circulate around the face. Should be 
      const RealVect& x0 = outgoingEdge->getVertex()->getPosition();
      const RealVect& x1 = outgoingEdge->getPreviousEdge()->getVertex()->getPosition();
      const RealVect& x2 = outgoingEdge->getNextEdge()->getVertex()->getPosition();

      RealVect v1 = x1-x0;
      RealVect v2 = x2-x0;

      v1 = v1/v1.vectorLength();
      v2 = v2/v2.vectorLength();

      //    const RealVect norm = PolyGeom::cross(v1, v2);

      // We could use the face normal, but I'm decoupling the way these are computed....
      const RealVect norm = outgoingEdge->getFace()->getNormal();

      const Real alpha = acos(v1.dotProduct(v2));

      m_normal += alpha*norm;
    }

    this->normalizeNormalVector();
  }

  inline
  RealVect& vertex::getPosition() noexcept {
    return m_pos;
  }

  inline
  const RealVect& vertex::getPosition() const noexcept {
    return m_pos;
  }

  inline
  RealVect& vertex::getNormal() noexcept {
    return m_normal;
  }

  inline
  const RealVect& vertex::getNormal() const noexcept {
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
  Real vertex::signedDistance(const RealVect& a_x0) const noexcept {
    const RealVect delta = a_x0 - m_pos;
    const Real dist      = delta.vectorLength();
    const Real dot       = m_normal.dotProduct(delta);
    const int sign       = (dot > 0.) ? 1 : -1;
  
    return dist*sign;
  }

  inline
  Real vertex::unsignedDistance2(const RealVect& a_x0) const noexcept {
    const RealVect d = a_x0 - m_pos;

    return d.dotProduct(d);
  }
}

#endif
