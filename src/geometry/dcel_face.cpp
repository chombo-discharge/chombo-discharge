/*!
  @file   dcel_face.cpp
  @brief  Implementation of dcel_face.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"
#include "dcel_poly2.H"
#include "dcel_iterator.H"

#include <PolyGeom.H>

using namespace dcel;

face::face(){
  m_normal = RealVect::Zero;
  m_poly2Algorithm = InsideOutsideAlgorithm::CrossingNumber;
}

face::face(const std::shared_ptr<edge>& a_edge) : face() {
  this->setHalfEdge(a_edge);
}

face::face(const face& a_otherFace) : face() {
  this->define(a_otherFace.getNormal(),
	       a_otherFace.getHalfEdge());
}

face::~face(){

}

void face::define(const RealVect& a_normal, const std::shared_ptr<edge>& a_edge) noexcept {
  this->setNormal(a_normal);
  this->setHalfEdge(a_edge);
}

void face::setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept {
  m_halfEdge = a_halfEdge;
}

void face::setNormal(const RealVect& a_normal) noexcept {
  m_normal = a_normal;
}

void face::normalizeNormalVector() noexcept {
  m_normal = m_normal/m_normal.vectorLength();
}

void face::setInsideOutsideAlgorithm(const InsideOutsideAlgorithm& a_algorithm) noexcept {
  m_poly2Algorithm = a_algorithm;
}

void face::computeArea() noexcept {
  m_area = 0.0;

  for (int i = 0; i < m_vertices.size() - 1; i++){
    const RealVect& v1 = m_vertices[i]->getPosition();
    const RealVect& v2 = m_vertices[i+1]->getPosition();
    m_area += m_normal.dotProduct(PolyGeom::cross(v2,v1));
  }

  m_area = 0.5*std::abs(m_area);
}

void face::computeCentroid() noexcept {
  m_centroid = RealVect::Zero;

  const int N = m_vertices.size();
  
  for (const auto& v : m_vertices){
    m_centroid += v->getPosition();
  }
  
  m_centroid = m_centroid/N; 
}

void face::computeNormal() noexcept {
  // Go through all vertices because some vertices may (correctly) lie on a line (but all of them shouldn't).

  const int n = m_vertices.size();
  
  for (int i = 0; i < n; i++){
    const RealVect& x0 = m_vertices[i]      ->getPosition();
    const RealVect& x1 = m_vertices[(i+1)%n]->getPosition();
    const RealVect& x2 = m_vertices[(i+2)%n]->getPosition();

    m_normal = PolyGeom::cross(x2-x1, x2-x0);
    
    if(m_normal.vectorLength() > 0.0) break;
  }

  this->normalizeNormalVector();
}

void face::computeBoundingSphere() noexcept {
  m_boundingSphere.define(this->getAllVertexCoordinates(), BoundingSphere::Algorithm::Ritter);
}

void face::computeBoundingBox() noexcept {
  m_boundingBox.define(this->getAllVertexCoordinates());
}

void face::computeVerticesAndEdges() noexcept {
  m_edges    = this->gatherEdges();
  m_vertices = this->gatherVertices();
}

void face::computePolygon2D() noexcept {
  m_poly2 = std::make_shared<Polygon2D>(*this);
}

const std::shared_ptr<edge>& face::getHalfEdge() const noexcept{
  return m_halfEdge;
}

std::shared_ptr<edge>& face::getHalfEdge() noexcept {
  return m_halfEdge;
}

const std::vector<RealVect> face::getAllVertexCoordinates() const noexcept {
  std::vector<RealVect> pos;
  
  for (const auto& v : m_vertices){
    pos.emplace_back(v->getPosition());
  }
  
  return pos;
}

std::vector<std::shared_ptr<vertex> >& face::getVertices() noexcept {
  return m_vertices;
}

const std::vector<std::shared_ptr<vertex> >& face::getVertices() const noexcept {
  return m_vertices;
}

std::vector<std::shared_ptr<edge> >& face::getEdges() noexcept {
  return m_edges;
}

const std::vector<std::shared_ptr<edge> >& face::getEdges() const noexcept {
  return m_edges;
}

const std::vector<std::shared_ptr<vertex> > face::gatherVertices() const noexcept {
  std::vector<std::shared_ptr<vertex> > vertices;

  for (edge_iterator iter(*this); iter.ok(); ++iter){
    std::shared_ptr<edge>& edge = iter();
    vertices.emplace_back(edge->getVertex());
  }

  return vertices;
}

const std::vector<std::shared_ptr<edge> > face::gatherEdges() const noexcept {
  std::vector<std::shared_ptr<edge> > edges;

  for (edge_iterator iter(*this); iter.ok(); ++iter){
    edges.emplace_back(iter());
  }

  return edges;
}

RealVect& face::getNormal() noexcept {
  return m_normal;
}

const RealVect& face::getNormal() const noexcept {
  return m_normal;
}

RealVect& face::getCentroid() noexcept {
  return m_centroid;
}

const RealVect& face::getCentroid() const noexcept {
  return m_centroid;
}

Real& face::getArea() noexcept {
  return m_area;
}

const Real& face::getArea() const noexcept {
  return m_area;
}

RealVect& face::getBoundingBoxLo() noexcept {
  return m_boundingBox.getLowCorner();
}

const RealVect& face::getBoundingBoxLo() const noexcept {
  return m_boundingBox.getLowCorner();
}

RealVect& face::getBoundingBoxHi() noexcept {
  return m_boundingBox.getHighCorner();
}

const RealVect& face::getBoundingBoxHi() const noexcept {
  return m_boundingBox.getHighCorner();
}

Real face::signedDistance(const RealVect& a_x0) const noexcept {
  Real retval = 1.234567E89;

  const bool inside = this->isPointInsideFace(a_x0);

  if(inside){ 
    retval = m_normal.dotProduct(a_x0 - m_centroid);
  }

  // Now check the edges. 
  for (const auto& e : m_edges){
    const Real curDist = e->signedDistance(a_x0);
    
    if(std::abs(curDist) <= std::abs(retval)){ // <= because edge normals are more important than polygon normals. 
      retval = curDist;
    }
  }

  return retval;
}

Real face::unsignedDistance2(const RealVect& a_x0) const noexcept {
  Real retval = 1.234567E89;

  const bool inside = this->isPointInsideFace(a_x0);

  if(inside){ 
    retval  = m_normal.dotProduct(a_x0 - m_centroid);
    retval *= retval;
  }
  else{ // The projected point lies outside the triangle. Check distance to edges/vertices
    for (const auto& e : m_edges){
      const Real curDist = e->unsignedDistance2(a_x0);
      retval = std::min(retval, curDist);
    }
  }

  return retval;
}

RealVect face::projectPointIntoFacePlane(const RealVect& a_p) const noexcept {
  const RealVect& planePoint     = m_vertices.front()->getPosition();
  const RealVect normalComponent = m_normal.dotProduct(a_p - planePoint) * m_normal;
  const RealVect projectedPoint  = a_p - normalComponent;

  return projectedPoint;
}

bool face::isPointInsideFace(const RealVect& a_p) const noexcept {
  const RealVect projectedPoint = this->projectPointIntoFacePlane(a_p);

  return m_poly2->isPointInside(projectedPoint, m_poly2Algorithm);
}
