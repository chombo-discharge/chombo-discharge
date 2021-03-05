/*!
  @file   dcel_meshI.H
  @brief  Implementation of dcel_mesh.H
  @author Robert Marskar
  @date   Apr. 2018
  @todo   Bugs in compute_vertex_normals. 
*/

#include "dcel_mesh.H"
#include "dcel_iterator.H"
#include "dcel_vertex.H"
#include "dcel_polygon.H"
#include "dcel_edge.H"
#include <PolyGeom.H>

#include <chrono>

using namespace dcel;

bool mesh::s_angle_weighted = false; // Use angle-weighted vertex normal vectors. This currently breaks. 

mesh::mesh(){
  m_reconciled = false;
  m_use_tree   = false;
}

mesh::mesh(std::vector<std::shared_ptr<polygon> >& a_polygons,
	   std::vector<std::shared_ptr<edge> >&    a_edges,
	   std::vector<std::shared_ptr<vertex> >&  a_vertices){
  
  m_reconciled = false;
  m_use_tree   = false;
  
  this->define(a_polygons, a_edges, a_vertices);
}

mesh::~mesh(){

}

void mesh::define(std::vector<std::shared_ptr<polygon> >& a_polygons,
		  std::vector<std::shared_ptr<edge> >&    a_edges,
		  std::vector<std::shared_ptr<vertex> >&  a_vertices) noexcept{
  m_polygons = a_polygons;
  m_edges    = a_edges;
  m_vertices = a_vertices;
}

bool mesh::sanityCheck() const {
  for (int i = 0; i < m_edges.size(); i++){
    if(m_edges[i] == nullptr){
      std::cerr << "mesh::sanity_check - edge is NULL\n";
    }
    else{
      const std::shared_ptr<edge>& edge = m_edges[i];
    
      if(edge->getPairEdge()== nullptr){
	std::cerr << "mesh::sanity_check - pair edge is NULL, your geometry probably isn't watertight.\n";
      }
      else if(edge->getNextEdge()== nullptr){
	std::cerr << "mesh::sanity_check - next edge is NULL, something has gone wrong with edge generation.\n";
      }
      else if(edge->getPreviousEdge()== nullptr){
	std::cerr << "mesh::sanity_check - prev edge is NULL, something has gone wrong with edge generation.\n";
      }
      else if(edge->getVertex()== nullptr){
	std::cerr << "mesh::sanity_check - vertex is NULL, something has gone wrong with edge generation.\n";
      }
    }
  }

  for (int i = 0; i < m_vertices.size(); i++){
    if(m_vertices[i]== nullptr){
      std::cerr << "mesh::sanity_check - m_vertices[i] is NULL, something has gone wrong with vertex generation.\n";
    }
    else{
      if(m_vertices[i]->getEdge()== nullptr){
	CH_assert(m_vertices[i]->getPolycache().size() == 0);
	pout() << "mesh::sanity_check - vertex edge is NULL, you may have an unreferenced vertex." << endl;
      }
    }
  }

  return true;
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

std::vector<std::shared_ptr<polygon> >& mesh::getPolygons() noexcept {
  return m_polygons;
}

const std::vector<std::shared_ptr<polygon> >& mesh::getPolygons() const noexcept {
  return m_polygons;
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

void mesh::reconcilePolygons(const bool a_outwardNormal, const bool a_recompute_vnormal) noexcept{

  /*!
    @brief Reconcile polygon edges. This gives each edge a reference to the polygon they circulate, and also computes the 
    polygon area
  */

  // Reconcile polygons; compute polygon area and provide edges explicit access
  // to their polygons
  for (int i = 0; i < m_polygons.size(); i++){
    std::shared_ptr<polygon>& poly = m_polygons[i];

    // Every edge gets a reference to this polygon
    for (edge_iterator iter(*poly); iter.ok(); ++iter){
      std::shared_ptr<edge>& edge = iter();
      edge->setPolygon(poly);
    }
    
    poly->computeNormal(a_outwardNormal);
    poly->computeCentroid();
    poly->computeArea();
    poly->normalizeNormalVector();
    poly->computeBoundingBox();
  }

  if(a_recompute_vnormal){   // Compute pseudonormals for vertices 
    std::cerr << "mesh::reconcile_polygons - there is probably a bug in the vertex normal computation somewhere\n";
    this->computeVertexNormals();
  }
  
  this->computeEdgeNormals();
  this->computeBoundingSphere();

  m_reconciled = true;
}

void mesh::computeVertexNormals() noexcept {
#define debug_func 1

#if debug_func
  pout() << "starting computation" << endl;
#endif
  for (int i = 0; i < m_vertices.size(); i++){
    if(!(m_vertices[i]->getEdge()== nullptr)){
#if 1 // This doesn't work, why?!?
      const std::vector<std::shared_ptr<polygon> > polygons = m_vertices[i]->getPolygons();
#else
      const std::vector<std::shared_ptr<polygon> > polygons = m_vertices[i]->getPolycache();
#endif

      // Mean or area weighted
      if(!s_angle_weighted){
	RealVect normal = RealVect::Zero;
	for (int j = 0; j < polygons.size(); j++){
	  //normal += polygons[j]->get_area()*polygons[j]->getNormal(); // Area weighted
	  normal += polygons[j]->getNormal(); // Mean
	}

	// Set normal
	if(normal.vectorLength() > 0.0){
	  normal *= 1./normal.vectorLength();
	  m_vertices[i]->setNormal(normal);
	}
	else{
	  normal = polygons[1]->getNormal();
	  m_vertices[i]->setNormal(normal);
	}
      }
      else { // Angle-weighted normal vector
	RealVect normal = RealVect::Zero;
#if debug_func
	int num = 0;
#endif
	for (edge_iterator iter(*m_vertices[i]); iter.ok(); ++iter){ 
	  const std::shared_ptr<edge>& outgoing = iter();
	  const std::shared_ptr<edge>& incoming = outgoing->getPreviousEdge();

	  const RealVect origin = incoming->getVertex()->getPosition();
	  const RealVect x2     = outgoing->getVertex()->getPosition();
	  const RealVect x1     = incoming->getOtherVertex()->getPosition();
	  const Real len1       = (x1-origin).vectorLength();
	  const Real len2       = (x2-origin).vectorLength();

	  const RealVect norm = PolyGeom::cross(x2-origin, x1-origin)/(len1*len2);

#if debug_func
	  CH_assert(len1 > 0.0);
	  CH_assert(len2 > 0.0);
	  CH_assert(norm.vectorLength() > 0.0);
#endif
	  //	const Real alpha = asin(norm.vectorLength());

	  const Real alpha = acos(PolyGeom::dot(x2-origin, x1-origin)/len1*len2);

	  normal += alpha*norm/norm.vectorLength();
#if debug_func
	  num++;
	  pout() << num << endl;

	  if(num > 20){
	    pout() << "problem vertex = " << m_vertices[i]->getPosition() << endl;
	    std::cerr << "dcel_compute_vertex_normals - stop\n";
	  }
#endif
	}
	normal *= 1./normal.vectorLength();

	m_vertices[i]->setNormal(normal);
      }
    }

#if debug_func
    pout() << "done computing vertex vectors" << endl;
#endif
  }
}

void mesh::computeEdgeNormals() noexcept {
  for (int i = 0; i < m_edges.size(); i++){

    std::shared_ptr<edge>& cur_edge        = m_edges[i];
    const std::shared_ptr<edge>& pair_edge = cur_edge->getPairEdge();

    const std::shared_ptr<polygon>& poly      = cur_edge->getPolygon();
    const std::shared_ptr<polygon>& pair_poly = pair_edge->getPolygon();


    const RealVect n1 = poly->getNormal();
    const RealVect n2 = pair_poly->getNormal();
    
    RealVect normal = (n1 + n2);
    normal = normal/normal.vectorLength();
    
    cur_edge->setNormal(normal);
  }
}

void mesh::buildKdTree(const int a_max_depth, const int a_max_elements) noexcept {
  m_tree = std::shared_ptr<kd_tree<polygon> > (new kd_tree<polygon>(m_polygons, a_max_depth, a_max_elements));
}

Real mesh::signedDistance(const RealVect& a_point) const noexcept {
  return this->signedDistance(a_point, m_algorithm);
}

Real mesh::DirectSignedDistance(const RealVect& a_point) const noexcept {
  Real minDist = m_polygons.front()->signedDistance(a_point);
    
  for (const auto& poly : m_polygons){
    const Real curDist = poly->signedDistance(a_point);

    if(std::abs(curDist) < std::abs(minDist)) minDist = curDist;
  }

  return minDist;
}

Real mesh::KdTreeSignedDistance(const RealVect& a_point) const noexcept {

  Real minDist;

  if(m_tree->isDefined()){
    std::vector<std::shared_ptr<polygon> > candidates = m_tree->find_closest(a_point);

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
  case SearchAlgorithm::KdTree:
    minDist = this->KdTreeSignedDistance(a_point);
    break;
  default:
    std::cerr << "Error in file dcel_mesh mesh::signedDistance unsupported algorithm requested\n";
    break;
  }

  return minDist;
}

void mesh::computeVertexNormals(VertexNormalComputation a_comp) noexcept {
  for (auto& v : m_vertices){
    if (v == nullptr) std::cerr << "In file dcel_mesh function dcel::mesh::computeVertexNormals(VertexNormalComputation) - vertex is 'nullptr'\n";

    switch(a_comp) {
    case VertexNormalComputation::Average:
      this->computeVertexNormalAverage(v);
      break;
    case VertexNormalComputation::AngleWeighted:
      this->computeVertexNormalAngleWeighted(v);
      break;
    default:
      std::cerr << "In file dcel_mesh function dcel::mesh::computeVertexNormal(VertexNormalComputation) - unsupported algorithm requested\n";
    }
  }
}

void mesh::computeVertexNormalAverage(std::shared_ptr<vertex>& a_vert) noexcept {
#if 1 // This doesn't work, why?!?
  auto polygons = a_vert->getPolygons();
#else
  auto polygons = a_vert->getPolycache();
#endif

  auto& normal = a_vert->getNormal();

  normal = RealVect::Zero;
  for (const auto& p : polygons){
    normal += p->getNormal();
  }

  a_vert->normalizeNormalVector();
}

void mesh::computeVertexNormalAngleWeighted(std::shared_ptr<vertex>& a_vert) noexcept {

  RealVect& normal = a_vert->getNormal();
  
  normal = RealVect::Zero;
  
  for (edge_iterator iter(*a_vert); iter.ok(); ++iter){ 
    const std::shared_ptr<edge>& outgoing = iter();
    const std::shared_ptr<edge>& incoming = outgoing->getPreviousEdge();

    const RealVect origin = incoming->getVertex()->getPosition();
    const RealVect x2     = outgoing->getVertex()->getPosition();
    const RealVect x1     = incoming->getOtherVertex()->getPosition();
    const Real len1       = (x1-origin).vectorLength();
    const Real len2       = (x2-origin).vectorLength();

    const RealVect norm = PolyGeom::cross(x2-origin, x1-origin)/(len1*len2);

    //	const Real alpha = asin(norm.vectorLength());

    const Real alpha = acos(PolyGeom::dot(x2-origin, x1-origin)/len1*len2);

    normal += alpha*norm/norm.vectorLength();

  }

  a_vert->normalizeNormalVector();
}


