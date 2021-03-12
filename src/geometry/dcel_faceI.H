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

  template <class T>
  inline
  faceT<T>::faceT(){
    m_halfEdge       = nullptr;
    m_normal         = Vec3::Zero;
    m_poly2Algorithm = InsideOutsideAlgorithm::CrossingNumber;
  }

  template <class T>
  inline
  faceT<T>::faceT(const std::shared_ptr<edge>& a_edge) : face() {
    m_halfEdge = a_edge;
  }

  template <class T>
  inline
  faceT<T>::faceT(const face& a_otherFace) : face() {
    m_normal   = a_otherFace.getNormal();
    m_halfEdge = a_otherFace.getHalfEdge();
  }

  template <class T>
  inline
  faceT<T>::~faceT(){
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::define(const Vec3& a_normal, const std::shared_ptr<edge>& a_edge) noexcept {
    m_normal   = a_normal;
    m_halfEdge = a_edge;
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::reconcile() noexcept {
    this->computeNormal();
    this->normalizeNormalVector();
    this->computeCentroid();
    this->computeArea();
    this->computePolygon2D();
    this->computeBoundingBox();
    this->computeBoundingSphere();
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept {
    m_halfEdge = a_halfEdge;
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::normalizeNormalVector() noexcept {
    m_normal = m_normal/m_normal.length();
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::setInsideOutsideAlgorithm(const InsideOutsideAlgorithm& a_algorithm) noexcept {
    m_poly2Algorithm = a_algorithm;
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::computeArea() noexcept {
    m_area = 0.0;

    const auto vertices = this->gatherVertices();

    for (int i = 0; i < vertices.size() - 1; i++){
      const auto& v1 = vertices[i]  ->getPosition();
      const auto& v2 = vertices[i+1]->getPosition();
      m_area += m_normal.dot(v2.cross(v1));
    }

    m_area = 0.5*std::abs(m_area);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::computeCentroid() noexcept {
    m_centroid = Vec3::Zero;

    const auto vertices = this->gatherVertices();
  
    for (const auto& v : vertices){
      m_centroid += v->getPosition();
    }
  
    m_centroid = m_centroid/vertices.size();
  }

  template <class T>  
  inline
  decltype(auto) faceT<T>::computeNormal() noexcept {
    const auto vertices = this->gatherVertices();
    
    const int N = vertices.size();
  
    for (int i = 0; i < N; i++){
      const auto& x0 = vertices[i]      ->getPosition();
      const auto& x1 = vertices[(i+1)%N]->getPosition();
      const auto& x2 = vertices[(i+2)%N]->getPosition();

      m_normal = (x2-x1).cross(x2-x0);
    
      if(m_normal.length() > 0.0) break; // Found one. 
    }

    this->normalizeNormalVector();
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::computeBoundingBox() noexcept{
    m_boundingBox.define(this->getAllVertexCoordinates());
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::computeBoundingSphere() noexcept {
    m_boundingSphere.define(this->getAllVertexCoordinates(), BoundingVolumeAlgorithm::Ritter);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::computePolygon2D() noexcept{
    m_poly2 = std::make_shared<Polygon2D<T> >(m_normal, this->getAllVertexCoordinates());
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getCentroid() noexcept {
    return (m_centroid);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getCentroid() const noexcept {
    return (m_centroid);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getNormal() noexcept {
    return (m_normal);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getNormal() const noexcept {
    return (m_normal);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getArea() noexcept {
    return (m_area);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getArea() const noexcept {
    return (m_area);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getBoundingBoxLo() noexcept {
    return (m_boundingBox.getLowCorner());
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getBoundingBoxLo() const noexcept {
    return (m_boundingBox.getLowCorner());
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getBoundingBoxHi() noexcept {
    return (m_boundingBox.getHighCorner());
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getBoundingBoxHi() const noexcept {
    return (m_boundingBox.getHighCorner());
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getHalfEdge() noexcept {
    return (m_halfEdge);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getHalfEdge() const noexcept {
    return (m_halfEdge);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::gatherVertices() const noexcept {
    std::vector<std::shared_ptr<vertex> > vertices;

    for (edgeIterator iter(*this); iter.ok(); ++iter){
      std::shared_ptr<edge>& edge = iter();
      vertices.emplace_back(edge->getVertex());
    }

    return vertices;
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::getAllVertexCoordinates() const noexcept {
    std::vector<Vec3> ret;

    for (edgeIterator iter(*this); iter.ok(); ++iter){
      std::shared_ptr<edge>& edge = iter();
      ret.emplace_back(edge->getVertex()->getPosition());
    }

    return ret;
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::projectPointIntoFacePlane(const Vec3& a_p) const noexcept {
    const auto normalComponent = m_normal.dot(a_p - m_centroid) * m_normal;
    const auto projectedPoint  = a_p - normalComponent;

    return projectedPoint;
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::isPointInsideFace(const Vec3& a_p) const noexcept {
    const Vec3 p = this->projectPointIntoFacePlane(a_p);

    return m_poly2->isPointInside(p, m_poly2Algorithm);
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::signedDistance(const Vec3& a_x0) const noexcept {
    T retval = std::numeric_limits<T>::max();

    const bool inside = this->isPointInsideFace(a_x0);

    if(inside){ 
      retval = m_normal.dot(a_x0 - m_centroid);
    }
    else {
      for (edgeIterator edgeIt(*this); edgeIt.ok(); ++edgeIt){
	const auto& e = edgeIt();
	
	const T curDist = e->signedDistance(a_x0);
    
	if(std::abs(curDist) < std::abs(retval)){ // <= because edge normals are more important than polygon normals. 
	  retval = curDist;
	}
      }
    }

    return retval;
  }

  template <class T>
  inline
  decltype(auto) faceT<T>::unsignedDistance2(const Vec3& a_x0) const noexcept {
    T retval = std::numeric_limits<T>::max();

    const bool inside = this->isPointInsideFace(a_x0);

    if(inside){ 
      retval  = m_normal.dot(a_x0 - m_centroid);
      retval *= retval;
    }
    else{ // The projected point lies outside the triangle. Check distance to edges/vertices
      for (edgeIterator edgeIt(*this); edgeIt.ok(); ++edgeIt){
	const auto& e = edgeIt();
	
	const T curDist = e->unsignedDistance2(a_x0);
	retval = std::min(retval, curDist);
      }
    }

    return retval;
  }
}
#endif
