/*!
  @file   dcel_meshI.H
  @brief  Implementation of dcel_mesh.H
  @author Robert Marskar
  @date   March 2021
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
  meshT<T>::meshT(){

  }

  template <class T>
  inline
  meshT<T>::meshT(std::vector<std::shared_ptr<face> >&   a_faces,
		  std::vector<std::shared_ptr<edge> >&   a_edges,
		  std::vector<std::shared_ptr<vertex> >& a_vertices){
    this->define(a_faces, a_edges, a_vertices);
  }

  template <class T>
  inline
  meshT<T>::~meshT(){

  }

  template <class T>
  inline
  decltype(auto) meshT<T>::define(std::vector<std::shared_ptr<face> >&   a_faces,
				  std::vector<std::shared_ptr<edge> >&   a_edges,
				  std::vector<std::shared_ptr<vertex> >& a_vertices) noexcept{
    m_faces    = a_faces;
    m_edges    = a_edges;
    m_vertices = a_vertices;
  }

  template <class T>
  inline
  decltype(auto) meshT<T>::sanityCheck() const noexcept {
    for (const auto& e : m_edges){
      const auto& nextEdge  = e->getNextEdge();
      const auto& prevEdge  = e->getPreviousEdge();
      const auto& pairEdge  = e->getPairEdge();
      const auto& curVertex = e->getVertex();
      const auto& curPoly   = e->getFace();

      // Check basic points for current edge. 
      if(e == nullptr) {
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::meshT<T>::sanityCheck  - edge is nullptr\n";
      }
      else if(pairEdge == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::meshT<T>::sanityCheck  - pair edge is nullptr, your geometry probably isn't watertight.\n";
      }
      else if(nextEdge == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::meshT<T>::sanityCheck  - next edge is nullptr, something has gone wrong with edge generation.\n";
      }
      else if(prevEdge == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::meshT<T>::sanityCheck  - previous edge is nullptr, something has gone wrong with edge generation.\n";
      }
      else if(curVertex == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::meshT<T>::sanityCheck  - vertex is nullptr, something has gone wrong with edge generation.\n";
      }
      else if(curPoly == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::meshT<T>::sanityCheck  - face is nullptr, something has gone wrong with edge generation.\n";
      }

      // Check that the next edge's previous edge is this edge. 
      if(prevEdge->getNextEdge() != e){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::meshT<T>::sanityCheck  - this->getPreviousEdge()->getNextEdge() is not the current edge, but it should be.\n";
      }
      else if(nextEdge->getPreviousEdge() != e){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::meshT<T>::sanityCheck  - this->getNextEdge()->getPreviousEdge() is not the current edge, but it should be.\n";
      }
    }

    // Vertex check
    for (const auto& v : m_vertices){
      if(v == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::meshT<T>::sanityCheck  - got a nullptr vertex, something has gone wrong with vertex generation.\n";
      }
      else if(v->getEdge() == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::meshT<T>::sanityCheck  - vertex has a nullptr edge, something has gone wrong with vertex generation.\n";
      }
    }
  }

  template <class T>
  inline
  decltype(auto) meshT<T>::setSearchAlgorithm(const SearchAlgorithm a_algorithm) noexcept {
    m_algorithm = a_algorithm;
  }

  template <class T>
  inline
  decltype(auto) meshT<T>::setInsideOutsideAlgorithm(const InsideOutsideAlgorithm a_algorithm) noexcept {
    for (auto& f : m_faces){
      f->setInsideOutsideAlgorithm(a_algorithm);
    }
  }

  template <class T>
  inline
  decltype(auto) meshT<T>::reconcile(VertexNormalWeight a_weight) noexcept {
    this->reconcileFaces();
    this->reconcileEdges();
    this->reconcileVertices(a_weight);
  }

  template <class T>
  inline
  decltype(auto) meshT<T>::reconcileFaces() noexcept {
    for (auto& f : m_faces){
      f->reconcile();
    }
  }

  template <class T>
  inline
  decltype(auto) meshT<T>::reconcileEdges() noexcept {
    for (auto& e : m_edges){
      e->reconcile();
    }
  }

  template <class T>
  inline
  decltype(auto) meshT<T>::reconcileVertices(VertexNormalWeight a_weight) noexcept {
    for (auto& v : m_vertices){
      switch(a_weight) {
      case VertexNormalWeight::None:
	v->computeVertexNormalAverage();
	break;
      case VertexNormalWeight::Angle:
	v->computeVertexNormalAngleWeighted();
	break;
      default:
	std::cerr << "In file dcel_mesh function dcel::meshT<T>::reconcileVertices(VertexNormalWeighting) - unsupported algorithm requested\n";
      }
    }
  }

  template <class T>
  inline
  decltype(auto) meshT<T>::getVertices() noexcept {
    return (m_vertices);
  }

  template <class T>
  inline
  decltype(auto) meshT<T>::getVertices() const noexcept {
    return (m_vertices);
  }

  template <class T>  
  inline
  decltype(auto) meshT<T>::getEdges() noexcept {
    return (m_edges);
  }

  template <class T>  
  inline
  decltype(auto) meshT<T>::getEdges() const noexcept {
    return (m_edges);
  }

  template <class T>
  inline
  decltype(auto) meshT<T>::getFaces() noexcept {
    return (m_faces);
  }

  template <class T>  
  inline
  decltype(auto) meshT<T>::getFaces() const noexcept {
    return (m_faces);
  }

  template <class T>  
  inline
  decltype(auto) meshT<T>::getAllVertexCoordinates() const noexcept {
    std::vector<Vec3 > vertexCoordinates;
    for (const auto& v : m_vertices){
      vertexCoordinates.emplace_back(v->getPosition());
    }

    return vertexCoordinates;
  }



  template <class T>  
  inline
  decltype(auto) meshT<T>::signedDistance(const Vec3& a_point) const noexcept {
    return this->signedDistance(a_point, m_algorithm);
  }

  template <class T>  
  inline
  decltype(auto) meshT<T>::signedDistance(const Vec3& a_point, SearchAlgorithm a_algorithm) const noexcept {
    T minDist;
  
    switch(a_algorithm){
    case SearchAlgorithm::Direct:
      minDist = this->DirectSignedDistance(a_point);
      break;
    case SearchAlgorithm::Direct2:
      minDist = this->DirectSignedDistance2(a_point);
      break;
    default:
      std::cerr << "Error in file dcel_mesh meshT<T>::signedDistance unsupported algorithm requested\n";
      break;
    }

    return minDist;
  }

  template <class T>  
  inline
  decltype(auto) meshT<T>::DirectSignedDistance(const Vec3& a_point) const noexcept {
    T minDist = m_faces.front()->signedDistance(a_point);
    
    for (const auto& f : m_faces){
      const T curDist = f->signedDistance(a_point);

      if(std::abs(curDist) < std::abs(minDist)){
	minDist = curDist;
      }
    }

    return minDist;
  }

  template <class T>  
  inline
  decltype(auto) meshT<T>::DirectSignedDistance2(const Vec3& a_point) const noexcept {
    std::shared_ptr<face> closest = m_faces.front();
    T minDist2 = closest->unsignedDistance2(a_point);

    for (const auto& f : m_faces){
      const T curDist2 = f->unsignedDistance2(a_point);

      if(curDist2 < minDist2) {
	closest = f;
	minDist2 = curDist2;
      }
    }

    return closest->signedDistance(a_point);
  }
}

#endif
