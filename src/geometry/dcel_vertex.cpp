/*!
  @file   vertex.cpp
  @brief  Implementation of vertex.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"
#include "dcel_iterator.H"

using namespace dcel;

void vertex::define(const RealVect& a_pos, const std::shared_ptr<edge>& a_edge, const RealVect a_normal) noexcept {
  this->setPosition(a_pos);
  this->setEdge(a_edge);
  this->setNormal(a_normal);
}

void vertex::setPosition(const RealVect& a_pos) noexcept {
  m_pos = a_pos;
}

void vertex::setEdge(const std::shared_ptr<edge>& a_edge) noexcept {
  m_edge = a_edge;
}

void vertex::setNormal(const RealVect& a_normal) noexcept {
  m_normal = a_normal;
}

void vertex::addFaceToCache(const std::shared_ptr<face>& a_face) noexcept {
  m_faceCache.push_back(a_face);
}

void vertex::clearFaceCache() noexcept {
  m_faceCache.resize(0);
}

void vertex::normalizeNormalVector() noexcept {
  m_normal = m_normal/m_normal.vectorLength();
}

void vertex::computeVertexNormalAverage() noexcept {
  const auto& faces = this->getFaces();

  m_normal = RealVect::Zero;
  
  for (const auto& f : faces){
    m_normal += f->getNormal();
  }

  this->normalizeNormalVector();
}

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

std::vector<std::shared_ptr<face> > vertex::getFaces() noexcept {
  std::vector<std::shared_ptr<face> > faces;
  for (edge_iterator iter(*this); iter.ok(); ++iter){
    faces.push_back(iter()->getFace());
  }

  return faces;
}

const std::vector<std::shared_ptr<face> >& vertex::getFaceCache() const noexcept{
  return m_faceCache;
}

std::vector<std::shared_ptr<face> >& vertex::getFaceCache() noexcept {
  return m_faceCache;
}

Real vertex::signedDistance(const RealVect& a_x0) const noexcept {
  const RealVect delta = a_x0 - m_pos;
  const Real dist      = delta.vectorLength();
  const Real dot       = m_normal.dotProduct(delta);
  const int sign       = (dot > 0.) ? 1 : -1;
  
  return dist*sign;
}

Real vertex::unsignedDistance2(const RealVect& a_x0) const noexcept {
  const RealVect d = a_x0 - m_pos;

  return d.dotProduct(d);
}
