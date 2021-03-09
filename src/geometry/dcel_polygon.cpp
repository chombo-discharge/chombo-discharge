/*!
  @file   dcel_poly.cpp
  @brief  Implementation of dcel_poly.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_polygon.H"
#include "dcel_iterator.H"

#include <PolyGeom.H>
  
#define EPSILON 1.E-8
#define TWOPI 6.283185307179586476925287

using namespace dcel;

polygon::polygon(){
  m_normal = RealVect::Zero;
}

polygon::polygon(const std::shared_ptr<edge>& a_edge){
  this->setHalfEdge(a_edge);
}

polygon::polygon(const polygon& a_otherPolygon){
  this->define(a_otherPolygon.getNormal(),
	       a_otherPolygon.getHalfEdge());
}

polygon::~polygon(){

}

void polygon::define(const RealVect& a_normal, const std::shared_ptr<edge>& a_edge) noexcept {
  this->setNormal(a_normal);
  this->setHalfEdge(a_edge);
}

void polygon::setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept {
  m_halfEdge = a_halfEdge;
}

void polygon::setNormal(const RealVect& a_normal) noexcept {
  m_normal = a_normal;
}

void polygon::normalizeNormalVector() noexcept {
  m_normal *= 1./m_normal.vectorLength();
}

void polygon::computeArea() noexcept {
  const std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();

  m_area = 0.0;

  for (int i = 0; i < vertices.size() - 1; i++){
    const RealVect& v1 = vertices[i]->getPosition();
    const RealVect& v2 = vertices[i+1]->getPosition();
    m_area += m_normal.dotProduct(PolyGeom::cross(v2,v1));
  }

  m_area = 0.5*std::abs(m_area);
}

void polygon::computeCentroid() noexcept {
  const std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();

  m_centroid = RealVect::Zero;
  
  for (const auto& v : vertices){
    m_centroid += v->getPosition();
  }
  
  m_centroid = m_centroid/vertices.size();
}

void polygon::computeNormal(const bool a_outwardNormal) noexcept {
  std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();  

  // Go through all vertices because some vertices may (correctly) lie on a line (but all of them shouldn't).

  const int n = vertices.size();
  
  for (int i = 0; i < n; i++){
    const RealVect& x0 = vertices[i]      ->getPosition();
    const RealVect& x1 = vertices[(i+1)%n]->getPosition();
    const RealVect& x2 = vertices[(i+2)%n]->getPosition();

    m_normal = PolyGeom::cross(x2-x1, x2-x0);
    
    if(m_normal.vectorLength() > 0.0) break;
  }

  this->normalizeNormalVector();

  if(!a_outwardNormal){    // If normal points inwards, make it point outwards
    m_normal = -m_normal;
  }
}

void polygon::computeBoundingSphere() noexcept {
  m_boundingSphere.define(this->getAllVertexCoordinates(), BoundingSphere::Algorithm::Ritter);
}

void polygon::computeBoundingBox() noexcept {
  m_boundingBox.define(this->getAllVertexCoordinates());
}

const std::shared_ptr<edge>& polygon::getHalfEdge() const noexcept{
  return m_halfEdge;
}

std::shared_ptr<edge>& polygon::getHalfEdge() noexcept {
  return m_halfEdge;
}

const std::vector<RealVect> polygon::getAllVertexCoordinates() const noexcept {
  std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();

  std::vector<RealVect> pos;
  
  for (const auto& v : vertices){
    pos.emplace_back(v->getPosition());
  }
  
  return pos;
}

const std::vector<std::shared_ptr<vertex> > polygon::getVertices() const noexcept {
  std::vector<std::shared_ptr<vertex> > vertices;

  for (edge_iterator iter(*this); iter.ok(); ++iter){
    std::shared_ptr<edge>& edge = iter();
    vertices.push_back(edge->getVertex());
  }

  return vertices;
}

const std::vector<std::shared_ptr<edge> > polygon::getEdges() const noexcept {
  std::vector<std::shared_ptr<edge> > edges;

  for (edge_iterator iter(*this); iter.ok(); ++iter){
    edges.push_back(iter());
  }

  return edges;
}

RealVect& polygon::getNormal() noexcept {
  return m_normal;
}

const RealVect& polygon::getNormal() const noexcept {
  return m_normal;
}

RealVect& polygon::getCentroid() noexcept {
  return m_centroid;
}

const RealVect& polygon::getCentroid() const noexcept {
  return m_centroid;
}

Real& polygon::getArea() noexcept {
  return m_area;
}

const Real& polygon::getArea() const noexcept {
  return m_area;
}

RealVect& polygon::getBoundingBoxLo() noexcept {
  return m_boundingBox.getLowCorner();
}

const RealVect& polygon::getBoundingBoxLo() const noexcept {
  return m_boundingBox.getLowCorner();
}

RealVect& polygon::getBoundingBoxHi() noexcept {
  return m_boundingBox.getHighCorner();
}

const RealVect& polygon::getBoundingBoxHi() const noexcept {
  return m_boundingBox.getHighCorner();
}

Real polygon::signedDistance(const RealVect& a_x0) const noexcept {
  std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();
  
  Real retval = 1.234567E89;

  // Compute projection of x0 on the polygon plane
  const RealVect x1          = vertices.front()->getPosition();
  const Real normalComponent = PolyGeom::dot(a_x0-x1, m_normal);
  const RealVect xp          = a_x0 - normalComponent*m_normal;

  // Use angle rule to check if projected point lies inside the polygon. Very expensive because of the acos(cosTheta). 
  Real anglesum = 0.0;
  const int n = vertices.size();
  for(int i = 0; i < n; i++){
    const RealVect p1 = vertices[i]      ->getPosition() - xp;
    const RealVect p2 = vertices[(i+1)%n]->getPosition() - xp;

    const Real m1 = p1.vectorLength();
    const Real m2 = p2.vectorLength();

    if(m1*m2 <= EPSILON) {// Projected point hits a vertex, return early. 
      anglesum = 0.;
      break;
    }
    else {
      const Real cosTheta = PolyGeom::dot(p1, p2)/(m1*m2);
      anglesum += acos(cosTheta);
    }
  }

  // Projected point is inside if angles sum to 2*pi
  if(std::abs(std::abs(anglesum) - TWOPI) < EPSILON){ // Ok, the projection onto the polygon plane places the point "inside" the planer polygon
    retval = normalComponent;
  }
  else{ // The projected point lies outside the triangle. Check distance to edges/vertices
    const std::vector<std::shared_ptr<edge> > edges = this->getEdges();

    for (const auto& e : edges){
      const Real curDist = e->signedDistance(a_x0);
      retval = (std::abs(curDist) < std::abs(retval)) ? curDist : retval;
    }
  }

  return retval;
}
