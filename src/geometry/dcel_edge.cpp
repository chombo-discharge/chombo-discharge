/*!
  @file   dcel_edge.cpp
  @brief  Implementation of dcel_edge.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_polygon.H"
#include "dcel_iterator.H"

#include "PolyGeom.H"

using namespace dcel;

edge::edge(){
}

edge::edge(const std::shared_ptr<vertex>& a_vertex){
  this->setVertex(a_vertex);
  this->setPairEdge(nullptr);
  this->setNextEdge(nullptr);
  this->setPreviousEdge(nullptr);
  this->setNormal(RealVect::Zero);
}

edge::edge(const edge& a_otherEdge){
  this->define(a_otherEdge.getVertex(),
	       a_otherEdge.getPairEdge(),
	       a_otherEdge.getNextEdge(),
	       a_otherEdge.getPreviousEdge(),
	       a_otherEdge.getNormal());

  m_len2 = a_otherEdge.m_len2;
  m_x2x1 = a_otherEdge.m_x2x1;
}

edge::~edge(){
}

void edge::define(const std::shared_ptr<vertex>& a_vertex,
		  const std::shared_ptr<edge>&   a_pairEdge,
		  const std::shared_ptr<edge>&   a_nextEdge,
		  const std::shared_ptr<edge>&   a_previousEdge,
		  const RealVect                 a_normal){
  this->setVertex(a_vertex);
  this->setPairEdge(a_pairEdge);
  this->setNextEdge(a_nextEdge);
  this->setPreviousEdge(a_previousEdge);
  this->setNormal(a_normal);
}

void edge::setVertex(const std::shared_ptr<vertex>& a_vertex) noexcept {
  m_vertex = a_vertex;
}

void edge::setPairEdge(const std::shared_ptr<edge>& a_pairEdge) noexcept {
  m_pairEdge = a_pairEdge;
}

void edge::setNextEdge(const std::shared_ptr<edge>& a_nextEdge) noexcept {
  m_nextEdge = a_nextEdge;
}

void edge::setPreviousEdge(const std::shared_ptr<edge>& a_previousEdge) noexcept {
  m_previousEdge = a_previousEdge;
}

void edge::setNormal(const RealVect a_normal) noexcept {
  m_normal = a_normal;
}

void edge::setPolygon(const std::shared_ptr<polygon>& a_polygon) noexcept {
  m_polygon = a_polygon;
}

void edge::normalizeNormalVector() noexcept {
  m_normal = m_normal/m_normal.vectorLength();
}

void edge::computeEdgeLength() noexcept {
  const auto& x1 = this->getVertex()->getPosition();
  const auto& x2 = this->getOtherVertex()->getPosition();

  m_x2x1 = x2-x1;
  m_len2 = m_x2x1.dotProduct(m_x2x1);
}

std::shared_ptr<vertex>& edge::getVertex() noexcept {
  return m_vertex;
}

const std::shared_ptr<vertex>& edge::getVertex() const noexcept {
  return m_vertex;
}

std::shared_ptr<vertex>& edge::getOtherVertex() noexcept{
  return m_pairEdge->getVertex();
}

const std::shared_ptr<vertex>& edge::getOtherVertex() const noexcept {
  return m_pairEdge->getVertex();
}

std::shared_ptr<edge>& edge::getPairEdge() noexcept {
  return m_pairEdge;
}

const std::shared_ptr<edge>& edge::getPairEdge() const noexcept {
  return m_pairEdge;
}

std::shared_ptr<edge>& edge::getPreviousEdge() noexcept {
  return m_previousEdge;
}

const std::shared_ptr<edge>& edge::getPreviousEdge() const noexcept {
  return m_previousEdge;
}

std::shared_ptr<edge>& edge::getNextEdge() noexcept {
  return m_nextEdge;
}

const std::shared_ptr<edge>& edge::getNextEdge() const noexcept {
  return m_nextEdge;
}

RealVect& edge::getNormal() noexcept {
  return m_normal;
}

const RealVect& edge::getNormal() const noexcept {
  return m_normal;
}

std::shared_ptr<polygon>& edge::getPolygon() noexcept {
  return m_polygon;
}

const std::shared_ptr<polygon>& edge::getPolygon() const noexcept {
  return m_polygon;
}

Real edge::signedDistance(const RealVect& a_x0) const noexcept {

  const Real t = this->projectPointToEdge(a_x0);

  Real retval;
  if(t < 0.0) {
    retval = this->getVertex()->signedDistance(a_x0);
  }
  else if(t > 1.0){
    retval = this->getOtherVertex()->signedDistance(a_x0);
  }
  else{
    const RealVect linePoint = m_vertex->getPosition() + t*m_x2x1;
    const RealVect delta     = a_x0 - linePoint;
    const Real dot           = m_normal.dotProduct(delta);

    const int sgn = (dot >= 0.0) ? 1 : -1;

    retval = sgn*delta.vectorLength();
  }

  return retval;
}

Real edge::unsignedDistance2(const RealVect& a_x0) const noexcept {
  Real t = this->projectPointToEdge(a_x0);
  t = std::min(std::max(0.,t), 1.);

  const RealVect linePoint = m_vertex->getPosition() + t*m_x2x1;

  const RealVect d = a_x0 - linePoint;
    
  return d.dotProduct(d);
}

Real edge::projectPointToEdge(const RealVect& a_x0) const noexcept {
  const RealVect p = a_x0 - m_vertex->getPosition();

  return p.dotProduct(m_x2x1)/m_len2;
}
