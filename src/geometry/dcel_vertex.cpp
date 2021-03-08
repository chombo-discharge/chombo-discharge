/*!
  @file   vertex.cpp
  @brief  Implementation of vertex.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_polygon.H"
#include "dcel_iterator.H"

using namespace dcel;

vertex::vertex(){
  m_pos    = RealVect::Zero;
  m_normal = RealVect::Zero;

  m_polycache.resize(0);
}

vertex::vertex(const RealVect& a_pos){
  m_pos    = a_pos;
  m_normal = RealVect::Zero;
}

vertex::vertex(const RealVect& a_pos, const RealVect& a_normal){
  m_pos    = a_pos;
  m_normal = a_normal;
}

vertex::vertex(const vertex& a_otherVertex){
  this->define(a_otherVertex.getPosition(),
	       a_otherVertex.getEdge(),
	       a_otherVertex.getNormal());
}

vertex::~vertex(){

}

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

void vertex::addPolygonToCache(const std::shared_ptr<polygon>& a_poly) noexcept {
  m_polycache.push_back(a_poly);
}

void vertex::clearPolygonCache() noexcept {
  m_polycache.resize(0);
}

void vertex::normalizeNormalVector() noexcept {
  m_normal = m_normal/m_normal.vectorLength();
}

RealVect& vertex::getPosition() noexcept {
  return m_pos;
}

const RealVect& vertex::getPosition() const noexcept {
  return m_pos;
}

std::shared_ptr<edge>& vertex::getEdge() noexcept {
  return m_edge;
}

const std::shared_ptr<edge>& vertex::getEdge() const noexcept {
  return m_edge;
}

RealVect& vertex::getNormal() noexcept {
  return m_normal;
}

const RealVect& vertex::getNormal() const noexcept {
  return m_normal;
}

std::vector<std::shared_ptr<polygon> > vertex::getPolygons() noexcept {
  std::vector<std::shared_ptr<polygon> > polygons;
  for (edge_iterator iter(*this); iter.ok(); ++iter){
    polygons.push_back(iter()->getPolygon());
  }

  return polygons;
}

const std::vector<std::shared_ptr<polygon> >& vertex::getPolycache() const noexcept{
  return m_polycache;
}

std::vector<std::shared_ptr<polygon> >& vertex::getPolycache() noexcept {
  return m_polycache;
}



Real vertex::signedDistance(const RealVect a_x0) const noexcept {
  const RealVect delta = a_x0 - m_pos;
  const Real dist      = delta.vectorLength();
  const Real dot       = m_normal.dotProduct(delta);
  const int sign       = (dot > 0) ? 1 : -1;
  
  return dist*sign;
}
