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

vertex::vertex(const vertex& a_otherVertex){
  this->define(a_otherVertex.getPosition(),
	       a_otherVertex.getEdge(),
	       a_otherVertex.getNormal());
}

vertex::~vertex(){

}

void vertex::define(const RealVect& a_pos, const std::shared_ptr<edge>& a_edge, const RealVect a_normal){
  this->setPosition(a_pos);
  this->setEdge(a_edge);
  this->setNormal(a_normal);
}

void vertex::setPosition(const RealVect& a_pos){
  m_pos = a_pos;
}

void vertex::setEdge(const std::shared_ptr<edge>& a_edge){
  m_edge = a_edge;
}

void vertex::setNormal(const RealVect& a_normal){
  m_normal = a_normal;
}

void vertex::addPolygon(const std::shared_ptr<polygon>& a_poly){
  m_polycache.push_back(a_poly);
}

std::vector<std::shared_ptr<polygon> > vertex::getPolygons() {
  std::vector<std::shared_ptr<polygon> > polygons;
  for (edge_iterator iter(*this); iter.ok(); ++iter){
    polygons.push_back(iter()->getPolygon());
  }

  return polygons;
}

const std::vector<std::shared_ptr<polygon> >& vertex::getPolycache() const{
  return m_polycache;
}

std::vector<std::shared_ptr<polygon> >& vertex::getPolycache(){
  return m_polycache;
}

std::shared_ptr<edge>& vertex::getEdge(){
  return m_edge;
}

const std::shared_ptr<edge>& vertex::getEdge() const{
  return m_edge;
}

RealVect& vertex::getPosition() {
  return m_pos;
}

const RealVect& vertex::getPosition() const {
  return m_pos;
}

RealVect& vertex::getNormal() {
  return m_normal;
}

const RealVect& vertex::getNormal() const {
  return m_normal;
}

Real vertex::signedDistance(const RealVect a_x0) const {
  const RealVect delta = a_x0 - m_pos;
  const Real dist      = delta.vectorLength();
  const Real dot       = m_normal.dotProduct(delta);
  const int sign       = (dot > 0) ? 1 : -1;
  
  return dist*sign;
}
