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
  void mesh::sanityCheck() const {
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
  void mesh::setSearchAlgorithm(const SearchAlgorithm a_algorithm) noexcept {
    m_algorithm = a_algorithm;
  }

  inline
  void mesh::setInsideOutsideAlgorithm(const InsideOutsideAlgorithm a_algorithm) noexcept {
    for (auto& f : m_faces){
      f->setInsideOutsideAlgorithm(a_algorithm);
    }
  }

  inline
  void mesh::reconcile(VertexNormalWeight a_weight) noexcept {
    this->reconcileFaces();
    this->reconcileEdges();
    this->reconcileVertices(a_weight);
  }

  inline
  void mesh::reconcileFaces() noexcept{
    // Reconcile faces; compute face area and provide edges explicit access
    // to their faces
    for (auto& f : m_faces){
      f->computeVerticesAndEdges();
      f->computeNormal();
      f->normalizeNormalVector();
      f->computeCentroid();
      f->computeArea();
      f->computePolygon2D();
      f->computeBoundingBox();
      f->computeBoundingSphere();
    }
  }

  inline
  void mesh::reconcileEdges() noexcept {
    for (auto& e : m_edges){
      e->computeNormal();
      e->computeEdgeLength();
    }
  }

  inline
  void mesh::reconcileVertices(VertexNormalWeight a_weight) noexcept {
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
  void mesh::computeBoundingSphere() noexcept {
    m_boundingSphere.define(this->getAllVertexCoordinates(), BoundingSphere::Algorithm::Ritter);
  }

  inline
  void mesh::computeBoundingBox() noexcept {
    m_boundingBox.define(this->getAllVertexCoordinates());
  }

  inline
  std::vector<std::shared_ptr<vertex> >& mesh::getVertices() noexcept {
    return m_vertices;
  }

  inline
  const std::vector<std::shared_ptr<vertex> >& mesh::getVertices() const noexcept {
    return m_vertices;
  }

  inline
  std::vector<std::shared_ptr<edge> >& mesh::getEdges() noexcept {
    return m_edges;
  }

  inline
  const std::vector<std::shared_ptr<edge> >& mesh::getEdges() const noexcept {
    return m_edges;
  }

  inline
  std::vector<std::shared_ptr<face> >& mesh::getFaces() noexcept {
    return m_faces;
  }

  inline
  const std::vector<std::shared_ptr<face> >& mesh::getFaces() const noexcept {
    return m_faces;
  }

  inline
  std::vector<RealVect> mesh::getAllVertexCoordinates() const noexcept{
    std::vector<RealVect> vertexCoordinates;
    for (const auto& v : m_vertices){
      vertexCoordinates.emplace_back(v->getPosition());
    }

    return vertexCoordinates;
  }

  inline
  Real mesh::signedDistance(const RealVect& a_point) const noexcept {
    return this->signedDistance(a_point, m_algorithm);
  }

  inline
  Real mesh::signedDistance(const RealVect& a_point, SearchAlgorithm a_algorithm) const noexcept {
    Real minDist;
  
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
  Real mesh::DirectSignedDistance(const RealVect& a_point) const noexcept {
    Real minDist = m_faces.front()->signedDistance(a_point);
    
    for (const auto& f : m_faces){
      const Real curDist = f->signedDistance(a_point);

      if(std::abs(curDist) < std::abs(minDist)){
	minDist = curDist;
      }
    }

    return minDist;
  }

  inline
  Real mesh::DirectSignedDistance2(const RealVect& a_point) const noexcept {
    std::shared_ptr<face> closest = m_faces.front();
    Real minDist2 = closest->unsignedDistance2(a_point);

    for (const auto& f : m_faces){
      const Real curDist2 = f->unsignedDistance2(a_point);

      if(curDist2 < minDist2) {
	closest = f;
	minDist2 = curDist2;
      }
    }

    return closest->signedDistance(a_point);
  }
}

#endif
