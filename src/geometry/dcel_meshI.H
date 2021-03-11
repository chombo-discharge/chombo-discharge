/*!
  @file   dcel_meshI.H
  @brief  Implementation of dcel_mesh.H
  @author Robert Marskar
  @date   March 2021.
*/

#ifndef _DCEL_MESHI_H_
#define _DCEL_MESHI_H_

#include "dcel_mesh.H"
#include "dcel_iterator.H"
#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"

namespace dcel {

  template <class T>
  inline
  mesh<T>::mesh(){
    m_use_tree   = false;
  }

  template <class T>
  inline
  mesh<T>::mesh(std::vector<std::shared_ptr<face> >&   a_faces,
	     std::vector<std::shared_ptr<edge> >&   a_edges,
	     std::vector<std::shared_ptr<vertex> >& a_vertices){
    m_use_tree   = false;
  
    this->define(a_faces, a_edges, a_vertices);
  }

  template <class T>
  inline
  mesh<T>::~mesh(){

  }

  template <class T>
  inline
  void mesh<T>::define(std::vector<std::shared_ptr<face> >&   a_faces,
		    std::vector<std::shared_ptr<edge> >&   a_edges,
		    std::vector<std::shared_ptr<vertex> >& a_vertices) noexcept{
    m_faces = a_faces;
    m_edges    = a_edges;
    m_vertices = a_vertices;
  }

  template <class T>
  inline
  auto mesh<T>::sanityCheck() const noexcept -> void {
    for (const auto& e : m_edges){
      const auto& nextEdge  = e->getNextEdge();
      const auto& prevEdge  = e->getPreviousEdge();
      const auto& pairEdge  = e->getPairEdge();
      const auto& curVertex = e->getVertex();
      const auto& curPoly   = e->getFace();

      // Check basic points for current edge. 
      if(e == nullptr) {
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh<T>::sanityCheck  - edge is nullptr\n";
      }
      else if(pairEdge == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh<T>::sanityCheck  - pair edge is nullptr, your geometry probably isn't watertight.\n";
      }
      else if(nextEdge == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh<T>::sanityCheck  - next edge is nullptr, something has gone wrong with edge generation.\n";
      }
      else if(prevEdge == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh<T>::sanityCheck  - previous edge is nullptr, something has gone wrong with edge generation.\n";
      }
      else if(curVertex == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh<T>::sanityCheck  - vertex is nullptr, something has gone wrong with edge generation.\n";
      }
      else if(curPoly == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh<T>::sanityCheck  - face is nullptr, something has gone wrong with edge generation.\n";
      }

      // Check that the next edge's previous edge is this edge. 
      if(prevEdge->getNextEdge() != e){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh<T>::sanityCheck  - this->getPreviousEdge()->getNextEdge() is not the current edge, but it should be.\n";
      }
      else if(nextEdge->getPreviousEdge() != e){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh<T>::sanityCheck  - this->getNextEdge()->getPreviousEdge() is not the current edge, but it should be.\n";
      }
    }

    // Vertex check
    for (const auto& v : m_vertices){
      if(v == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh<T>::sanityCheck  - got a nullptr vertex, something has gone wrong with vertex generation.\n";
      }
      else if(v->getEdge() == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh<T>::sanityCheck  - vertex has a nullptr edge, something has gone wrong with vertex generation.\n";
      }
    }
  }

  template <class T>
  inline
  auto mesh<T>::setSearchAlgorithm(const SearchAlgorithm a_algorithm) noexcept -> void {
    m_algorithm = a_algorithm;
  }

  template <class T>
  inline
  auto mesh<T>::setInsideOutsideAlgorithm(const InsideOutsideAlgorithm a_algorithm) noexcept -> void {
    for (auto& f : m_faces){
      f->setInsideOutsideAlgorithm(a_algorithm);
    }
  }

  template <class T>
  inline
  auto mesh<T>::reconcile(VertexNormalWeight a_weight) noexcept -> void {
    this->reconcileFaces();
    this->reconcileEdges();
    this->reconcileVertices(a_weight);
  }

  template <class T>
  inline
  auto mesh<T>::reconcileFaces() noexcept -> void {
    // Reconcile faces; compute face area and provide edges explicit access
    // to their faces
    for (auto& f : m_faces){
      f->reconcile();
    }
  }

  template <class T>
  inline
  auto mesh<T>::reconcileEdges() noexcept -> void {
    for (auto& e : m_edges){
      e->reconcile();
    }
  }

  template <class T>
  inline
  auto mesh<T>::reconcileVertices(VertexNormalWeight a_weight) noexcept -> void {
    for (auto& v : m_vertices){
      switch(a_weight) {
      case VertexNormalWeight::None:
	v->computeVertexNormalAverage();
	break;
      case VertexNormalWeight::Angle:
	v->computeVertexNormalAngleWeighted();
	break;
      default:
	std::cerr << "In file dcel_mesh function dcel::mesh<T>::reconcileVertices(VertexNormalWeighting) - unsupported algorithm requested\n";
      }
    }
  }

  template <class T>
  inline
  auto mesh<T>::computeBoundingSphere() noexcept -> void {
    m_boundingSphere.define(this->getAllVertexCoordinates(), BoundingVolumeAlgorithm::Ritter);
  }

  template <class T>
  inline
  auto mesh<T>::computeBoundingBox() noexcept -> void {
    m_boundingBox.define(this->getAllVertexCoordinates());
  }

  template <class T>
  inline
  auto mesh<T>::getVertices() noexcept -> std::vector<std::shared_ptr<vertex> >& {
    return m_vertices;
  }

  template <class T>
  inline
  auto mesh<T>::getVertices() const noexcept -> const std::vector<std::shared_ptr<vertex> >& {
    return m_vertices;
  }

  template <class T>  
  inline
  auto mesh<T>::getEdges() noexcept -> std::vector<std::shared_ptr<edge> >& {
    return m_edges;
  }

  template <class T>  
  inline
  auto mesh<T>::getEdges() const noexcept -> const std::vector<std::shared_ptr<edge> >& {
    return m_edges;
  }

  template <class T>
  inline
  auto mesh<T>::getFaces() noexcept -> std::vector<std::shared_ptr<face> >& {
    return m_faces;
  }

  template <class T>  
  inline
  auto mesh<T>::getFaces() const noexcept -> const std::vector<std::shared_ptr<face> >& {
    return m_faces;
  }

  template <class T>  
  inline
  auto mesh<T>::getAllVertexCoordinates() const noexcept -> std::vector<Vec3 >{
    std::vector<Vec3 > vertexCoordinates;
    for (const auto& v : m_vertices){
      vertexCoordinates.emplace_back(v->getPosition());
    }

    return vertexCoordinates;
  }

  template <class T>  
  inline
  auto mesh<T>::signedDistance(const Vec3& a_point) const noexcept -> double {
    return this->signedDistance(a_point, m_algorithm);
  }

  template <class T>  
  inline
  auto mesh<T>::signedDistance(const Vec3& a_point, SearchAlgorithm a_algorithm) const noexcept -> double {
    double minDist;
  
    switch(a_algorithm){
    case SearchAlgorithm::Direct:
      minDist = this->DirectSignedDistance(a_point);
      break;
    case SearchAlgorithm::Direct2:
      minDist = this->DirectSignedDistance2(a_point);
      break;
    default:
      std::cerr << "Error in file dcel_mesh mesh<T>::signedDistance unsupported algorithm requested\n";
      break;
    }

    return minDist;
  }

  template <class T>  
  inline
  auto mesh<T>::DirectSignedDistance(const Vec3& a_point) const noexcept -> double {
    double minDist = m_faces.front()->signedDistance(a_point);
    
    for (const auto& f : m_faces){
      const double curDist = f->signedDistance(a_point);

      if(std::abs(curDist) < std::abs(minDist)){
	minDist = curDist;
      }
    }

    return minDist;
  }

  template <class T>  
  inline
  auto mesh<T>::DirectSignedDistance2(const Vec3& a_point) const noexcept -> double {
    std::shared_ptr<face> closest = m_faces.front();
    double minDist2 = closest->unsignedDistance2(a_point);

    for (const auto& f : m_faces){
      const double curDist2 = f->unsignedDistance2(a_point);

      if(curDist2 < minDist2) {
	closest = f;
	minDist2 = curDist2;
      }
    }

    return closest->signedDistance(a_point);
  }
}

#endif
