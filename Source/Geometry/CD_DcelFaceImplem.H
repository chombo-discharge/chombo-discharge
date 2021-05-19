/*!
  @file   dcel_faceI.H
  @brief  Implementation of dcel_face.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_FACEI_H_
#define _DCEL_FACEI_H_

#include <CD_DcelFace.H>
#include <CD_DcelIterator.H>

#include "CD_NamespaceHeader.H"
namespace Dcel {

  template <class T>
  inline
  faceT<T>::faceT(){
    m_halfEdge       = nullptr;
    m_normal         = Vec3::zero();
    m_poly2Algorithm = InsideOutsideAlgorithm::CrossingNumber;
  }

  template <class T>
  inline
  faceT<T>::faceT(const edgePtr& a_edge) : face() {
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
  void faceT<T>::define(const Vec3& a_normal, const edgePtr& a_edge) noexcept {
    m_normal   = a_normal;
    m_halfEdge = a_edge;
  }

  template <class T>
  inline
  void faceT<T>::reconcile() noexcept {
    this->computeNormal();
    this->normalizeNormalVector();
    this->computeCentroid();
    this->computeArea();
    this->computePolygon2D();
    this->computeAndStoreEdges();
  }

  template <class T>
  inline
  void faceT<T>::computeAndStoreEdges() noexcept {
    m_edges.resize(0);

    for (edgeIterator edgeIt(*this); edgeIt.ok(); ++edgeIt){
      m_edges.emplace_back(edgeIt());
    }
  }

  template <class T>
  inline
  void faceT<T>::setHalfEdge(const edgePtr& a_halfEdge) noexcept {
    m_halfEdge = a_halfEdge;
  }

  template <class T>
  inline
  void faceT<T>::normalizeNormalVector() noexcept {
    m_normal = m_normal/m_normal.length();
  }

  template <class T>
  inline
  void faceT<T>::setInsideOutsideAlgorithm(const InsideOutsideAlgorithm& a_algorithm) noexcept {
    m_poly2Algorithm = a_algorithm;
  }

  template <class T>
  inline
  void faceT<T>::computeArea() noexcept {
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
  void faceT<T>::computeCentroid() noexcept {
    m_centroid = Vec3::zero();

    const auto vertices = this->gatherVertices();
  
    for (const auto& v : vertices){
      m_centroid += v->getPosition();
    }
  
    m_centroid = m_centroid/vertices.size();
  }

  template <class T>  
  inline
  void faceT<T>::computeNormal() noexcept {
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
  void faceT<T>::computePolygon2D() noexcept{
    m_poly2 = std::make_shared<Polygon2D<T> >(m_normal, this->getAllVertexCoordinates());
  }

  template <class T>
  inline
  T& faceT<T>::getCentroid(const int a_dir) noexcept {
    return m_centroid[a_dir];
  }

  template <class T>
  inline
  const T& faceT<T>::getCentroid(const int a_dir) const noexcept {
    return m_centroid[a_dir];
  }

  template <class T>
  inline
  Vec3T<T>& faceT<T>::getCentroid() noexcept {
    return (m_centroid);
  }

  template <class T>
  inline
  const Vec3T<T>& faceT<T>::getCentroid() const noexcept {
    return (m_centroid);
  }

  template <class T>
  inline
  Vec3T<T>& faceT<T>::getNormal() noexcept {
    return (m_normal);
  }

  template <class T>
  inline
  const Vec3T<T>& faceT<T>::getNormal() const noexcept {
    return (m_normal);
  }

  template <class T>
  inline
  T faceT<T>::getArea() noexcept {
    return (m_area);
  }

  template <class T>
  inline
  T faceT<T>::getArea() const noexcept {
    return (m_area);
  }

  template <class T>
  inline
  std::shared_ptr<edgeT<T> >& faceT<T>::getHalfEdge() noexcept {
    return (m_halfEdge);
  }

  template <class T>
  inline
  const std::shared_ptr<edgeT<T> >& faceT<T>::getHalfEdge() const noexcept {
    return (m_halfEdge);
  }

  template <class T>
  inline
  std::vector<std::shared_ptr<vertexT<T> > > faceT<T>::gatherVertices() const noexcept {
    std::vector<vertexPtr> vertices;

    for (edgeIterator iter(*this); iter.ok(); ++iter){
      edgePtr& edge = iter();
      vertices.emplace_back(edge->getVertex());
    }

    return vertices;
  }

  template <class T>
  inline
  std::vector<Vec3T<T> > faceT<T>::getAllVertexCoordinates() const noexcept {
    std::vector<Vec3> ret;

    for (edgeIterator iter(*this); iter.ok(); ++iter){
      edgePtr& edge = iter();
      ret.emplace_back(edge->getVertex()->getPosition());
    }

    return ret;
  }

  template <class T>
  inline
  Vec3T<T> faceT<T>::getSmallestCoordinate() const noexcept {
    const auto coords = this->getAllVertexCoordinates();

    auto minCoord = coords.front();

    for (const auto& c : coords){
      minCoord = min(minCoord, c);
    }

    return minCoord;
  }

  template <class T>
  inline
  Vec3T<T> faceT<T>::getHighestCoordinate() const noexcept {
    const auto coords = this->getAllVertexCoordinates();

    auto maxCoord = coords.front();

    for (const auto& c : coords){
      maxCoord = max(maxCoord, c);
    }

    return maxCoord;
  }

  template <class T>
  inline
  Vec3T<T> faceT<T>::projectPointIntoFacePlane(const Vec3& a_p) const noexcept {
    return a_p - m_normal * (m_normal.dot(a_p-m_centroid));
  }

  template <class T>
  inline
  bool faceT<T>::isPointInsideFace(const Vec3& a_p) const noexcept {
    const Vec3 p = this->projectPointIntoFacePlane(a_p);

    return m_poly2->isPointInside(p, m_poly2Algorithm);
  }

  template <class T>
  inline
  T faceT<T>::signedDistance(const Vec3& a_x0) const noexcept {
    T retval = std::numeric_limits<T>::infinity();

    const bool inside = this->isPointInsideFace(a_x0);

    if(inside){ 
      retval = m_normal.dot(a_x0 - m_centroid);
    }
    else {
      for (const auto& e : m_edges){
	const T curDist = e->signedDistance(a_x0);

	retval = (curDist*curDist < retval*retval) ? curDist : retval;
      }
    }

    return retval;
  }

  template <class T>
  inline
  T faceT<T>::unsignedDistance2(const Vec3& a_x0) const noexcept {
    T retval = std::numeric_limits<T>::infinity();

    const bool inside = this->isPointInsideFace(a_x0);

    if(inside){
      const T normDist = m_normal.dot(a_x0-m_centroid);
      
      retval = normDist*normDist;
    }
    else{
      for (const auto& e : m_edges){
	const T curDist2 = e->unsignedDistance2(a_x0);

	retval = (curDist2 < retval) ? curDist2 : retval;
      }
    }

    return retval;
  }
}
#include "CD_NamespaceFooter.H"


#endif
