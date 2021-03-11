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

  inline
  mesh::mesh(){
    m_use_tree   = false;
  }

  inline
  mesh::mesh(std::vector<std::shared_ptr<face> >&   a_faces,
	     std::vector<std::shared_ptr<edge> >&   a_edges,
	     std::vector<std::shared_ptr<vertex> >& a_vertices){
    m_use_tree   = false;
  
    this->define(a_faces, a_edges, a_vertices);
  }

  inline
  mesh::~mesh(){

  }

  inline
  void mesh::define(std::vector<std::shared_ptr<face> >&   a_faces,
		    std::vector<std::shared_ptr<edge> >&   a_edges,
		    std::vector<std::shared_ptr<vertex> >& a_vertices) noexcept{
    m_faces = a_faces;
    m_edges    = a_edges;
    m_vertices = a_vertices;
  }

  inline
  auto mesh::sanityCheck() const noexcept -> void {
    for (const auto& e : m_edges){
      const auto& nextEdge  = e->getNextEdge();
      const auto& prevEdge  = e->getPreviousEdge();
      const auto& pairEdge  = e->getPairEdge();
      const auto& curVertex = e->getVertex();
      const auto& curPoly   = e->getFace();

      // Check basic points for current edge. 
      if(e == nullptr) {
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh::sanityCheck  - edge is nullptr\n";
      }
      else if(pairEdge == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh::sanityCheck  - pair edge is nullptr, your geometry probably isn't watertight.\n";
      }
      else if(nextEdge == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh::sanityCheck  - next edge is nullptr, something has gone wrong with edge generation.\n";
      }
      else if(prevEdge == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh::sanityCheck  - previous edge is nullptr, something has gone wrong with edge generation.\n";
      }
      else if(curVertex == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh::sanityCheck  - vertex is nullptr, something has gone wrong with edge generation.\n";
      }
      else if(curPoly == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh::sanityCheck  - face is nullptr, something has gone wrong with edge generation.\n";
      }

      // Check that the next edge's previous edge is this edge. 
      if(prevEdge->getNextEdge() != e){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh::sanityCheck  - this->getPreviousEdge()->getNextEdge() is not the current edge, but it should be.\n";
      }
      else if(nextEdge->getPreviousEdge() != e){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh::sanityCheck  - this->getNextEdge()->getPreviousEdge() is not the current edge, but it should be.\n";
      }
    }

    // Vertex check
    for (const auto& v : m_vertices){
      if(v == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh::sanityCheck  - got a nullptr vertex, something has gone wrong with vertex generation.\n";
      }
      else if(v->getEdge() == nullptr){
	std::cerr << "In file 'dcel_mesh.cpp' function dcel::mesh::sanityCheck  - vertex has a nullptr edge, something has gone wrong with vertex generation.\n";
      }
    }
  }

  inline
  auto mesh::setSearchAlgorithm(const SearchAlgorithm a_algorithm) noexcept -> void {
    m_algorithm = a_algorithm;
  }

  inline
  auto mesh::setInsideOutsideAlgorithm(const InsideOutsideAlgorithm a_algorithm) noexcept -> void {
    for (auto& f : m_faces){
      f->setInsideOutsideAlgorithm(a_algorithm);
    }
  }

  inline
  auto mesh::reconcile(VertexNormalWeight a_weight) noexcept -> void {
    this->reconcileFaces();
    this->reconcileEdges();
    this->reconcileVertices(a_weight);
  }

  inline
  auto mesh::reconcileFaces() noexcept -> void {
    // Reconcile faces; compute face area and provide edges explicit access
    // to their faces
    for (auto& f : m_faces){
      f->reconcile();
    }
  }

  inline
  auto mesh::reconcileEdges() noexcept -> void {
    for (auto& e : m_edges){
      e->reconcile();
    }
  }

  inline
  auto mesh::reconcileVertices(VertexNormalWeight a_weight) noexcept -> void {
    for (auto& v : m_vertices){
      switch(a_weight) {
      case VertexNormalWeight::None:
	v->computeVertexNormalAverage();
	break;
      case VertexNormalWeight::Angle:
	v->computeVertexNormalAngleWeighted();
	break;
      default:
	std::cerr << "In file dcel_mesh function dcel::mesh::reconcileVertices(VertexNormalWeighting) - unsupported algorithm requested\n";
      }
    }
  }

  inline
  auto mesh::computeBoundingSphere() noexcept -> void {
    m_boundingSphere.define(this->getAllVertexCoordinates(), BoundingVolumeAlgorithm::Ritter);
  }

  inline
  auto mesh::computeBoundingBox() noexcept -> void {
    m_boundingBox.define(this->getAllVertexCoordinates());
  }

  inline
  auto mesh::getVertices() noexcept -> std::vector<std::shared_ptr<vertex> >& {
    return m_vertices;
  }

  inline
  auto mesh::getVertices() const noexcept -> const std::vector<std::shared_ptr<vertex> >& {
    return m_vertices;
  }

  inline
  auto mesh::getEdges() noexcept -> std::vector<std::shared_ptr<edge> >& {
    return m_edges;
  }

  inline
  auto mesh::getEdges() const noexcept -> const std::vector<std::shared_ptr<edge> >& {
    return m_edges;
  }

  inline
  auto mesh::getFaces() noexcept -> std::vector<std::shared_ptr<face> >& {
    return m_faces;
  }

  inline
  auto mesh::getFaces() const noexcept -> const std::vector<std::shared_ptr<face> >& {
    return m_faces;
  }

  inline
  auto mesh::getAllVertexCoordinates() const noexcept -> std::vector<Vec3 >{
    std::vector<Vec3 > vertexCoordinates;
    for (const auto& v : m_vertices){
      vertexCoordinates.emplace_back(v->getPosition());
    }

    return vertexCoordinates;
  }

  inline
  auto mesh::signedDistance(const Vec3& a_point) const noexcept -> double {
    return this->signedDistance(a_point, m_algorithm);
  }

  inline
  auto mesh::signedDistance(const Vec3& a_point, SearchAlgorithm a_algorithm) const noexcept -> double {
    double minDist;
  
    switch(a_algorithm){
    case SearchAlgorithm::Direct:
      minDist = this->DirectSignedDistance(a_point);
      break;
    case SearchAlgorithm::Direct2:
      minDist = this->DirectSignedDistance2(a_point);
      break;
    default:
      std::cerr << "Error in file dcel_mesh mesh::signedDistance unsupported algorithm requested\n";
      break;
    }

    return minDist;
  }

  inline
  auto mesh::DirectSignedDistance(const Vec3& a_point) const noexcept -> double {
    double minDist = m_faces.front()->signedDistance(a_point);
    
    for (const auto& f : m_faces){
      const double curDist = f->signedDistance(a_point);

      if(std::abs(curDist) < std::abs(minDist)){
	minDist = curDist;
      }
    }

    return minDist;
  }

  inline
  auto mesh::DirectSignedDistance2(const Vec3& a_point) const noexcept -> double {
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
