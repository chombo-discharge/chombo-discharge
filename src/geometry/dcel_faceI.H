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
    m_normal         = Vec3::Zero;
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
  auto face::define(const Vec3& a_normal, const std::shared_ptr<edge>& a_edge) noexcept -> void{
    m_normal   = a_normal;
    m_halfEdge = a_edge;
  }

  inline
  auto face::reconcile() noexcept -> void {
    this->computeNormal();
    this->normalizeNormalVector();
    this->computeCentroid();
    this->computeArea();
    this->computePolygon2D();
    this->computeBoundingBox();
    this->computeBoundingSphere();
  }

  inline
  auto face::setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept -> void{
    m_halfEdge = a_halfEdge;
  }

  inline
  auto face::normalizeNormalVector() noexcept -> void{
    m_normal = m_normal/m_normal.length();
  }

  inline
  auto face::setInsideOutsideAlgorithm(const InsideOutsideAlgorithm& a_algorithm) noexcept -> void{
    m_poly2Algorithm = a_algorithm;
  }

  inline
  auto face::computeArea() noexcept -> void{
    m_area = 0.0;

    const auto vertices = this->gatherVertices();

    for (int i = 0; i < vertices.size() - 1; i++){
      const Vec3& v1 = vertices[i]  ->getPosition();
      const Vec3& v2 = vertices[i+1]->getPosition();
      m_area += m_normal.dot(v2.cross(v1));
    }

    m_area = 0.5*std::abs(m_area);
  }

  inline
  auto face::computeCentroid() noexcept -> void{
    m_centroid = Vec3::Zero;

    const auto vertices = this->gatherVertices();
  
    for (const auto& v : vertices){
      m_centroid += v->getPosition();
    }
  
    m_centroid = m_centroid/vertices.size();
  }

  inline
  auto face::computeNormal() noexcept -> void{
    const auto vertices = this->gatherVertices();
    
    const int N = vertices.size();
  
    for (int i = 0; i < N; i++){
      const Vec3& x0 = vertices[i]      ->getPosition();
      const Vec3& x1 = vertices[(i+1)%N]->getPosition();
      const Vec3& x2 = vertices[(i+2)%N]->getPosition();

      m_normal = (x2-x1).cross(x2-x0);
    
      if(m_normal.length() > 0.0) break; // Found one. 
    }

    this->normalizeNormalVector();
  }

  inline
  auto face::computeBoundingBox() noexcept -> void{
    m_boundingBox.define(this->getAllVertexCoordinates());
  }

  inline
  auto face::computeBoundingSphere() noexcept -> void{
    m_boundingSphere.define(this->getAllVertexCoordinates(), BoundingSphere::Algorithm::Ritter);
  }

  inline
  auto face::computePolygon2D() noexcept -> void{
    m_poly2 = std::make_shared<Polygon2D<double> >(m_normal, this->getAllVertexCoordinates());
  }

  inline
  auto face::getCentroid() noexcept -> Vec3& {
    return m_centroid;
  }

  inline
  auto face::getCentroid() const noexcept -> const Vec3& {
    return m_centroid;
  }

  inline
  auto face::getNormal() noexcept -> Vec3&{
    return m_normal;
  }

  inline
  auto face::getNormal() const noexcept -> const Vec3&{
    return m_normal;
  }

  inline
  auto face::getArea() noexcept -> double&{
    return m_area;
  }

  inline
  auto face::getArea() const noexcept -> const double& {
    return m_area;
  }

  inline
  auto face::getBoundingBoxLo() noexcept -> Vec3& {
    return m_boundingBox.getLowCorner();
  }

  inline
  auto face::getBoundingBoxLo() const noexcept -> const Vec3& {
    return m_boundingBox.getLowCorner();
  }

  inline
  auto face::getBoundingBoxHi() noexcept -> Vec3& {
    return m_boundingBox.getHighCorner();
  }

  inline
  auto face::getBoundingBoxHi() const noexcept -> const Vec3& {
    return m_boundingBox.getHighCorner();
  }

  inline
  auto face::getHalfEdge() noexcept -> std::shared_ptr<edge>& {
    return m_halfEdge;
  }

  inline
  auto face::getHalfEdge() const noexcept -> const std::shared_ptr<edge>& {
    return m_halfEdge;
  }

  inline
  auto face::gatherVertices() const noexcept -> std::vector<std::shared_ptr<vertex> > {
    std::vector<std::shared_ptr<vertex> > vertices;

    for (edge_iterator iter(*this); iter.ok(); ++iter){
      std::shared_ptr<edge>& edge = iter();
      vertices.emplace_back(edge->getVertex());
    }

    return vertices;
  }

  inline
  auto face::getAllVertexCoordinates() const noexcept -> std::vector<Vec3> {
    std::vector<Vec3> ret;

    for (edge_iterator iter(*this); iter.ok(); ++iter){
      std::shared_ptr<edge>& edge = iter();
      ret.emplace_back(edge->getVertex()->getPosition());
    }

    return ret;
  }

  inline
  auto face::projectPointIntoFacePlane(const Vec3& a_p) const noexcept -> Vec3 {
    const Vec3 normalComponent = m_normal.dot(a_p - m_centroid) * m_normal;
    const Vec3 projectedPoint  = a_p - normalComponent;

    return projectedPoint;
  }

  inline
  auto face::isPointInsideFace(const Vec3& a_p) const noexcept -> bool {
    const Vec3 p = this->projectPointIntoFacePlane(a_p);

    return m_poly2->isPointInside(p, m_poly2Algorithm);
  }

  inline
  auto face::signedDistance(const Vec3& a_x0) const noexcept -> double{
    double retval = 1.234567E89;

    const bool inside = this->isPointInsideFace(a_x0);

    if(inside){ 
      retval = m_normal.dot(a_x0 - m_centroid);
    }
    else {
      for (edge_iterator edgeIt(*this); edgeIt.ok(); ++edgeIt){
	const auto& e = edgeIt();
	
	const double curDist = e->signedDistance(a_x0);
    
	if(std::abs(curDist) < std::abs(retval)){ // <= because edge normals are more important than polygon normals. 
	  retval = curDist;
	}
      }
    }

    return retval;
  }

  inline
  auto face::unsignedDistance2(const Vec3& a_x0) const noexcept -> double {
    double retval = 1.234567E89;

    const bool inside = this->isPointInsideFace(a_x0);

    if(inside){ 
      retval  = m_normal.dot(a_x0 - m_centroid);
      retval *= retval;
    }
    else{ // The projected point lies outside the triangle. Check distance to edges/vertices
      for (edge_iterator edgeIt(*this); edgeIt.ok(); ++edgeIt){
	const auto& e = edgeIt();
	
	const double curDist = e->unsignedDistance2(a_x0);
	retval = std::min(retval, curDist);
      }
    }

    return retval;
  }
}
#endif
