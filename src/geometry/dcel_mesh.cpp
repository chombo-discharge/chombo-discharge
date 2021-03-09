/*!
  @file   dcel_mesh.cpp
  @brief  Implementation of dcel_mesh.H
  @author Robert Marskar
  @date   March 2021.
*/

#include "dcel_mesh.H"
#include "dcel_iterator.H"
#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"
#include <PolyGeom.H>

using namespace dcel;

mesh::mesh(){
  m_use_tree   = false;
}

mesh::mesh(std::vector<std::shared_ptr<face> >& a_faces,
	   std::vector<std::shared_ptr<edge> >&    a_edges,
	   std::vector<std::shared_ptr<vertex> >&  a_vertices){
  m_use_tree   = false;
  
  this->define(a_faces, a_edges, a_vertices);
}

mesh::~mesh(){

}

void mesh::define(std::vector<std::shared_ptr<face> >& a_faces,
		  std::vector<std::shared_ptr<edge> >&    a_edges,
		  std::vector<std::shared_ptr<vertex> >&  a_vertices) noexcept{
  m_faces = a_faces;
  m_edges    = a_edges;
  m_vertices = a_vertices;
}

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

void mesh::setAlgorithm(SearchAlgorithm a_algorithm) noexcept {
  m_algorithm = a_algorithm;
}

std::vector<std::shared_ptr<vertex> >& mesh::getVertices() noexcept {
  return m_vertices;
}

const std::vector<std::shared_ptr<vertex> >& mesh::getVertices() const noexcept {
  return m_vertices;
}

std::vector<std::shared_ptr<edge> >& mesh::getEdges() noexcept {
  return m_edges;
}

const std::vector<std::shared_ptr<edge> >& mesh::getEdges() const noexcept {
  return m_edges;
}

std::vector<std::shared_ptr<face> >& mesh::getFaces() noexcept {
  return m_faces;
}

const std::vector<std::shared_ptr<face> >& mesh::getFaces() const noexcept {
  return m_faces;
}

std::vector<RealVect> mesh::getAllVertexCoordinates() const noexcept{
  std::vector<RealVect> vertexCoordinates;
  for (const auto& v : m_vertices){
    vertexCoordinates.emplace_back(v->getPosition());
  }

  return vertexCoordinates;
}

void mesh::computeBoundingSphere() noexcept {
  m_boundingSphere.define(this->getAllVertexCoordinates(), BoundingSphere::Algorithm::Ritter);
}

void mesh::computeBoundingBox() noexcept {
  m_boundingBox.define(this->getAllVertexCoordinates());
}

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

void mesh::computeVertexNormals(VertexNormalWeight a_weight) noexcept {
  for (auto& v : m_vertices){
    if (v == nullptr) std::cerr << "In file dcel_mesh function dcel::mesh::computeVertexNormals(VertexNormalWeighting) - vertex is 'nullptr'\n";

    switch(a_weight) {
    case VertexNormalWeight::None:
      this->computeVertexNormalAverage(v);
      break;
    case VertexNormalWeight::Angle:
      this->computeVertexNormalAngleWeighted(v);
      break;
    default:
      std::cerr << "In file dcel_mesh function dcel::mesh::computeVertexNormal(VertexNormalWeighting) - unsupported algorithm requested\n";
    }
  }
}

void mesh::computeEdgeNormals() noexcept {
  for (auto& e : m_edges){
    const std::shared_ptr<edge>& pairEdge = e->getPairEdge();

    const std::shared_ptr<face>& F     = e->getFace();
    const std::shared_ptr<face>& pairF = pairEdge->getFace();
    
    const RealVect& n1 = F->getNormal();
    const RealVect& n2 = pairF->getNormal();

    RealVect& normal = e->getNormal();
    
    normal = n1 + n2;

    e->normalizeNormalVector();

    e->computeEdgeLength();
  }
}

void mesh::computeVerticesAndEdges() noexcept {
  for (auto& f : m_faces){
    f->computeVerticesAndEdges();
    f->computePolygon2D();
  }
}

void mesh::buildKdTree(const int a_max_depth, const int a_max_elements) noexcept {
  m_tree = std::make_shared<kd_tree<face> > (m_faces, a_max_depth, a_max_elements);
}

Real mesh::signedDistance(const RealVect& a_point) const noexcept {
  return this->signedDistance(a_point, m_algorithm);
}

Real mesh::DirectSignedDistance(const RealVect& a_point) const noexcept {
  Real minDist = m_faces.front()->signedDistance(a_point);
    
  for (const auto& f : m_faces){
    const Real curDist = f->signedDistance(a_point);

    if(std::abs(curDist) < std::abs(minDist)) minDist = curDist;
  }

  return minDist;
}

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

Real mesh::KdTreeSignedDistance(const RealVect& a_point) const noexcept {

  Real minDist;

  if(m_tree->isDefined()){
    std::vector<std::shared_ptr<face> > candidates = m_tree->find_closest(a_point);

    minDist = candidates.front()->signedDistance(a_point);
    
    for (const auto& p : candidates){
      const Real curDist = p->signedDistance(a_point);

      if(std::abs(curDist) < std::abs(minDist)) minDist = curDist;
    }
  }
  else{
    std::cerr << "In file dcel_mesh function mesh::KdTreeSignedDistance - tree is not defined!\n";
  }

  return minDist;
}

Real mesh::signedDistance(const RealVect& a_point, SearchAlgorithm a_algorithm) const noexcept {
  Real minDist;
  
  switch(a_algorithm){
  case SearchAlgorithm::Direct:
    minDist = this->DirectSignedDistance(a_point);
    break;
  case SearchAlgorithm::Direct2:
    minDist = this->DirectSignedDistance2(a_point);
    break;
  case SearchAlgorithm::KdTree:
    minDist = this->KdTreeSignedDistance(a_point);
    break;
  default:
    std::cerr << "Error in file dcel_mesh mesh::signedDistance unsupported algorithm requested\n";
    break;
  }

  return minDist;
}

void mesh::computeVertexNormalAverage(std::shared_ptr<vertex>& a_vert) noexcept {
  const auto& faces = a_vert->getFaces();

  RealVect& normal = a_vert->getNormal();

  normal = RealVect::Zero;
  
  for (const auto& p : faces){
    normal += p->getNormal();
  }

  a_vert->normalizeNormalVector();
}

void mesh::computeVertexNormalAngleWeighted(std::shared_ptr<vertex>& a_vert) noexcept {

  RealVect& normal = a_vert->getNormal();
  
  normal = RealVect::Zero;

  for (edge_iterator iter(*a_vert); iter.ok(); ++iter){
    const auto& outgoingEdge = iter();

    // Edges circulate around the face. Should be 
    const RealVect& x0 = outgoingEdge->getVertex()->getPosition();
    const RealVect& x1 = outgoingEdge->getPreviousEdge()->getVertex()->getPosition();
    const RealVect& x2 = outgoingEdge->getNextEdge()->getVertex()->getPosition();

    RealVect v1 = x1-x0;
    RealVect v2 = x2-x0;

    v1 = v1/v1.vectorLength();
    v2 = v2/v2.vectorLength();

    const RealVect norm = PolyGeom::cross(v1, v2);

    // We could use the face normal, but I'm decoupling the way these are computed....
    /// const RealVect norm = outgoingEdge->getFace()->getNormal();

    const Real alpha = acos(PolyGeom::dot(v1, v2));

    normal += alpha*norm;
  }

  a_vert->normalizeNormalVector();
}


