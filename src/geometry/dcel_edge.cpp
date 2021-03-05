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

edge::edge(const edge& a_otherEdge){
  this->define(a_otherEdge.getVertex(),
	       a_otherEdge.getPairEdge(),
	       a_otherEdge.getNextEdge(),
	       a_otherEdge.getPreviousEdge(),
	       a_otherEdge.getNormal());
}

edge::~edge(){
}

void edge::define(const std::shared_ptr<vertex>&    a_vertex,
		  const std::shared_ptr<edge>& a_pairEdge,
		  const std::shared_ptr<edge>& a_nextEdge,
		  const std::shared_ptr<edge>& a_previousEdge,
		  const RealVect               a_normal){
  this->setVertex(a_vertex);
  this->setPairEdge(a_pairEdge);
  this->setNextEdge(a_nextEdge);
  this->setPreviousEdge(a_previousEdge);
  this->setNormal(a_normal);
}

void edge::setPolygon(const std::shared_ptr<polygon>& a_polygon){
  m_polygon = a_polygon;
}

void edge::setVertex(const std::shared_ptr<vertex>& a_vertex){
  m_vertex = a_vertex;
}

void edge::setPairEdge(const std::shared_ptr<edge>& a_pairEdge){
  m_pairEdge = a_pairEdge;
}

void edge::setNextEdge(const std::shared_ptr<edge>& a_nextEdge){
  m_nextEdge = a_nextEdge;
}

void edge::setPreviousEdge(const std::shared_ptr<edge>& a_previousEdge){
  m_previousEdge = a_previousEdge;
}

void edge::setNormal(const RealVect a_normal){
  m_normal = a_normal;
}

std::shared_ptr<vertex>& edge::getVertex() {
  return m_vertex;
}

const std::shared_ptr<vertex>& edge::getVertex() const {
  return m_vertex;
}

std::shared_ptr<vertex>& edge::getOtherVertex(){
  return m_pairEdge->getVertex();
}

const std::shared_ptr<vertex>& edge::getOtherVertex() const {
  return m_pairEdge->getVertex();
}

std::shared_ptr<edge>& edge::getPairEdge(){
  return m_pairEdge;
}

const std::shared_ptr<edge>& edge::getPairEdge() const {
  return m_pairEdge;
}

std::shared_ptr<edge>& edge::getPreviousEdge() {
  return m_previousEdge;
}

const std::shared_ptr<edge>& edge::getPreviousEdge() const {
  return m_previousEdge;
}

std::shared_ptr<edge>& edge::getNextEdge(){
  return m_nextEdge;
}

const std::shared_ptr<edge>& edge::getNextEdge() const {
  return m_nextEdge;
}

std::shared_ptr<polygon>& edge::getPolygon(){
  return m_polygon;
}

const std::shared_ptr<polygon>& edge::getPolygon() const{
  return m_polygon;
}

RealVect& edge::getNormal() {
  return m_normal;
}

const RealVect& edge::getNormal() const {
  return m_normal;
}

Real edge::signedDistance(const RealVect a_x0) const {
  Real retval = 1.234567E89;
  
  // Involved vertices
  const RealVect x1 = this->getOtherVertex()->getPosition();
  const RealVect x2 = m_vertex->getPosition();

  const RealVect R  = PolyGeom::cross(x2-x1,PolyGeom::cross(x1-a_x0, x2-a_x0));
  const RealVect r  = R/R.vectorLength();
  const RealVect xp = a_x0 + PolyGeom::dot(x1-a_x0,r)*r;
  const Real t      = PolyGeom::dot(xp-x1,x2-x1)/PolyGeom::dot(x2-x1, x2-x1);


  RealVect p;
  RealVect n;
  if(t < 0.0){ // Closest to x1, vertex normal takes precedence.
    p = x1;
    n = this->getOtherVertex()->getNormal();
  }
  else if (t > 1.0){ // Closest to x2, vertex normal takes precedence
    p = x2;
    n = m_vertex->getNormal();
  }
  else{ // Projection onto line lies on the line segment
    p = xp;
    n = m_normal;
  }

  const Real dot = PolyGeom::dot(n, (a_x0 - p)); // Determine sign from projection. If the point is orthogonal to the normal,
  const int sgn = dot >= 0.0 ? 1 : -1;           // it must (I think) be outside

  retval = (a_x0-p).vectorLength()*sgn;

  return retval;
}

