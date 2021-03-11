/*!
  @file   dcel_faceI.H
  @brief  Implementation of dcel_face.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_FACEI_H_
#define _DCEL_FACEI_H_

#include "dcel_face.H"
#include "dcel_iterator.H"

namespace dcel {

  inline
  face::face(){
    m_halfEdge       = nullptr;
    m_normal         = Vec3<double>::Zero;
    m_poly2Algorithm = InsideOutsideAlgorithm::CrossingNumber;
  }

  inline
  face::face(const std::shared_ptr<edge>& a_edge) : face() {
    m_halfEdge = a_edge;
  }

  inline
  face::face(const face& a_otherFace) : face() {
    m_normal   = a_otherFace.getNormal();
    m_halfEdge = a_otherFace.getHalfEdge();
  }

  inline
  face::~face(){
  }

  inline
  void face::define(const Vec3<double>& a_normal, const std::shared_ptr<edge>& a_edge) noexcept {
    m_normal   = a_normal;
    m_halfEdge = a_edge;
  }

  inline
  void face::reconcile() noexcept {
    this->computeNormal();
    this->normalizeNormalVector();
    this->computeCentroid();
    this->computeArea();
    this->computePolygon2D();
    this->computeBoundingBox();
    this->computeBoundingSphere();
  }

  inline
  void face::setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept {
    m_halfEdge = a_halfEdge;
  }

  inline
  void face::normalizeNormalVector() noexcept {
    m_normal = m_normal/m_normal.length();
  }

  inline
  void face::setInsideOutsideAlgorithm(const InsideOutsideAlgorithm& a_algorithm) noexcept {
    m_poly2Algorithm = a_algorithm;
  }

  inline
  void face::computeArea() noexcept {
    m_area = 0.0;

    const auto vertices = this->gatherVertices();

    for (int i = 0; i < vertices.size() - 1; i++){
      const Vec3<double>& v1 = vertices[i]  ->getPosition();
      const Vec3<double>& v2 = vertices[i+1]->getPosition();
      m_area += m_normal.dot(v2.cross(v1));
    }

    m_area = 0.5*std::abs(m_area);
  }

  inline
  void face::computeCentroid() noexcept {
    m_centroid = Vec3<double>::Zero;

    const auto vertices = this->gatherVertices();
  
    for (const auto& v : vertices){
      m_centroid += v->getPosition();
    }
  
    m_centroid = m_centroid/vertices.size();
  }

  inline
  void face::computeNormal() noexcept {
    const auto vertices = this->gatherVertices();
    
    const int N = vertices.size();
  
    for (int i = 0; i < N; i++){
      const Vec3<double>& x0 = vertices[i]      ->getPosition();
      const Vec3<double>& x1 = vertices[(i+1)%N]->getPosition();
      const Vec3<double>& x2 = vertices[(i+2)%N]->getPosition();

      m_normal = (x2-x1).cross(x2-x0);
    
      if(m_normal.length() > 0.0) break; // Found one. 
    }

    this->normalizeNormalVector();
  }

  inline
  void face::computeBoundingBox() noexcept {
    m_boundingBox.define(this->getAllVertexCoordinates());
  }

  inline
  void face::computeBoundingSphere() noexcept {
    m_boundingSphere.define(this->getAllVertexCoordinates(), BoundingSphere::Algorithm::Ritter);
  }

  inline
  void face::computePolygon2D() noexcept {
    m_poly2 = std::make_shared<Polygon2D<double> >(m_normal, this->getAllVertexCoordinates());
  }

  inline
  Vec3<double>& face::getCentroid() noexcept {
    return m_centroid;
  }

  inline
  const Vec3<double>& face::getCentroid() const noexcept {
    return m_centroid;
  }

  inline
  Vec3<double>& face::getNormal() noexcept {
    return m_normal;
  }

  inline
  const Vec3<double>& face::getNormal() const noexcept {
    return m_normal;
  }

  inline
  Real& face::getArea() noexcept {
    return m_area;
  }

  inline
  const Real& face::getArea() const noexcept {
    return m_area;
  }

  inline
  Vec3<double>& face::getBoundingBoxLo() noexcept {
    return m_boundingBox.getLowCorner();
  }

  inline
  const Vec3<double>& face::getBoundingBoxLo() const noexcept {
    return m_boundingBox.getLowCorner();
  }

  inline
  Vec3<double>& face::getBoundingBoxHi() noexcept {
    return m_boundingBox.getHighCorner();
  }

  inline
  const Vec3<double>& face::getBoundingBoxHi() const noexcept {
    return m_boundingBox.getHighCorner();
  }

  inline
  std::shared_ptr<edge>& face::getHalfEdge() noexcept {
    return m_halfEdge;
  }

  inline
  const std::shared_ptr<edge>& face::getHalfEdge() const noexcept{
    return m_halfEdge;
  }

  inline
  const std::vector<std::shared_ptr<vertex> > face::gatherVertices() const noexcept {
    std::vector<std::shared_ptr<vertex> > vertices;

    for (edge_iterator iter(*this); iter.ok(); ++iter){
      std::shared_ptr<edge>& edge = iter();
      vertices.emplace_back(edge->getVertex());
    }

    return vertices;
  }

  inline
  const std::vector<Vec3<double> > face::getAllVertexCoordinates() const noexcept {
    std::vector<Vec3<double> > ret;

    for (edge_iterator iter(*this); iter.ok(); ++iter){
      std::shared_ptr<edge>& edge = iter();
      ret.emplace_back(edge->getVertex()->getPosition());
    }

    return ret;
  }

  inline
  Vec3<double> face::projectPointIntoFacePlane(const Vec3<double>& a_p) const noexcept {
    const Vec3<double> normalComponent = m_normal.dot(a_p - m_centroid) * m_normal;
    const Vec3<double> projectedPoint  = a_p - normalComponent;

    return projectedPoint;
  }

  inline
  bool face::isPointInsideFace(const Vec3<double>& a_p) const noexcept {
    const Vec3<double> p = this->projectPointIntoFacePlane(a_p);

    return m_poly2->isPointInside(p, m_poly2Algorithm);
  }

  inline
  Real face::signedDistance(const Vec3<double>& a_x0) const noexcept {
    Real retval = 1.234567E89;

    const bool inside = this->isPointInsideFace(a_x0);

    if(inside){ 
      retval = m_normal.dot(a_x0 - m_centroid);
    }
    else {
      for (edge_iterator edgeIt(*this); edgeIt.ok(); ++edgeIt){
	const auto& e = edgeIt();
	
	const Real curDist = e->signedDistance(a_x0);
    
	if(std::abs(curDist) < std::abs(retval)){ // <= because edge normals are more important than polygon normals. 
	  retval = curDist;
	}
      }
    }

    return retval;
  }

  inline
  Real face::unsignedDistance2(const Vec3<double>& a_x0) const noexcept {
    Real retval = 1.234567E89;

    const bool inside = this->isPointInsideFace(a_x0);

    if(inside){ 
      retval  = m_normal.dot(a_x0 - m_centroid);
      retval *= retval;
    }
    else{ // The projected point lies outside the triangle. Check distance to edges/vertices
      for (edge_iterator edgeIt(*this); edgeIt.ok(); ++edgeIt){
	const auto& e = edgeIt();
	
	const Real curDist = e->unsignedDistance2(a_x0);
	retval = std::min(retval, curDist);
      }
    }

    return retval;
  }
}
#endif
