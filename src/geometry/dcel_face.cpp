/*!
  @file   dcel_face.cpp
  @brief  Implementation of dcel_face.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"
#include "dcel_iterator.H"

#include <PolyGeom.H>

using namespace dcel;

Polygon2D::Point2D::Point2D(const Real a_x, const Real a_y){
  x = a_x;
  y = a_y;
}

Polygon2D::Polygon2D(const face& a_face){
  this->define(a_face);
}

Polygon2D::Point2D Polygon2D::projectPoint(const RealVect& a_point) const noexcept {
  return Point2D(a_point[m_xDir], a_point[m_yDir]);
}

void Polygon2D::define(const face& a_face) noexcept {
  m_points.resize(0);

  m_ignoreDir = 0;

  const RealVect& normal = a_face.getNormal();
  
  for (int dir = 0; dir < SpaceDim; dir++){
    m_ignoreDir = (normal[dir] > normal[m_ignoreDir]) ? dir : m_ignoreDir;
  }

  m_xDir = 3;
  m_yDir = -1;

  for (int dir = 0; dir < SpaceDim; dir++){
    if(dir != m_ignoreDir){
      m_xDir = std::min(m_xDir, dir);
      m_yDir = std::max(m_yDir, dir);
    }
  }

  // Ignore coordinate with biggest normal component
  for (const auto& v : a_face.getVertices()){
    const RealVect& p = v->getPosition();
    
    m_points.emplace_back(projectPoint(p));
  }
}

int Polygon2D::isLeft(const Point2D& P0, const Point2D& P1, const Point2D& P2) const noexcept {
  return ( (P1.x - P0.x) * (P2.y - P0.y) - (P2.x -  P0.x) * (P1.y - P0.y) );
}

int Polygon2D::wn_PnPoly(const Point2D& P) const noexcept {
  int wn = 0;    // the  winding number counter

  const int N = m_points.size();

  // loop through all edges of the polygon
  for (int i = 0; i < N; i++) {   // edge from V[i] to  V[i+1]

    const Point2D& P1 = m_points[i];
    const Point2D& P2 = m_points[(i+1)%N];
    
    if (P1.y <= P.y) {          // start y <= P.y
      if (P2.y  > P.y)      // an upward crossing
	if (this->isLeft( P1, P2, P) > 0)  // P left of  edge
	  ++wn;            // have  a valid up intersect
    }
    else {                        // start y > P.y (no test needed)
      if (P2.y  <= P.y)     // a downward crossing
	if (this->isLeft( P1, P2, P) < 0)  // P right of  edge
	  --wn;            // have  a valid down intersect
    }
  }
  
  return wn;
}

int Polygon2D::wn_PnPoly(const RealVect& a_point) const noexcept {
  const Point2D point2 = this->projectPoint(a_point);
  const int WN         = this->wn_PnPoly(point2);

  return WN;
}

bool Polygon2D::isPointInsidePolygon(const RealVect& a_point) const noexcept {
  const int wn = this->wn_PnPoly(a_point);

  return wn != 0;
}

face::face(){
  m_normal = RealVect::Zero;
}

face::face(const std::shared_ptr<edge>& a_edge){
  this->setHalfEdge(a_edge);
}

face::face(const face& a_otherFace){
  this->define(a_otherFace.getNormal(),
	       a_otherFace.getHalfEdge());
}

face::~face(){

}

void face::define(const RealVect& a_normal, const std::shared_ptr<edge>& a_edge) noexcept {
  this->setNormal(a_normal);
  this->setHalfEdge(a_edge);
}

void face::setHalfEdge(const std::shared_ptr<edge>& a_halfEdge) noexcept {
  m_halfEdge = a_halfEdge;
}

void face::setNormal(const RealVect& a_normal) noexcept {
  m_normal = a_normal;
}

void face::normalizeNormalVector() noexcept {
  m_normal *= 1./m_normal.vectorLength();
}

void face::computeArea() noexcept {
  m_area = 0.0;

  for (int i = 0; i < m_vertices.size() - 1; i++){
    const RealVect& v1 = m_vertices[i]->getPosition();
    const RealVect& v2 = m_vertices[i+1]->getPosition();
    m_area += m_normal.dotProduct(PolyGeom::cross(v2,v1));
  }

  m_area = 0.5*std::abs(m_area);
}

void face::computeCentroid() noexcept {
  m_centroid = RealVect::Zero;
  
  for (const auto& v : m_vertices){
    m_centroid += v->getPosition();
  }
  
  m_centroid = m_centroid/m_vertices.size();
}

void face::computeNormal(const bool a_flipNormal) noexcept {
  // Go through all vertices because some vertices may (correctly) lie on a line (but all of them shouldn't).

  const int n = m_vertices.size();
  
  for (int i = 0; i < n; i++){
    const RealVect& x0 = m_vertices[i]      ->getPosition();
    const RealVect& x1 = m_vertices[(i+1)%n]->getPosition();
    const RealVect& x2 = m_vertices[(i+2)%n]->getPosition();

    m_normal = PolyGeom::cross(x2-x1, x2-x0);
    
    if(m_normal.vectorLength() > 0.0) break;
  }

  this->normalizeNormalVector();

  if(a_flipNormal){    // If normal points inwards, make it point outwards
    m_normal = -m_normal;
  }
}

void face::computeBoundingSphere() noexcept {
  m_boundingSphere.define(this->getAllVertexCoordinates(), BoundingSphere::Algorithm::Ritter);
}

void face::computeBoundingBox() noexcept {
  m_boundingBox.define(this->getAllVertexCoordinates());
}

void face::computeVerticesAndEdges() noexcept {
  m_edges    = this->gatherEdges();
  m_vertices = this->gatherVertices();
}

void face::computePolygon2D() noexcept {
  m_poly2 = std::make_shared<Polygon2D>(*this);
}

const std::shared_ptr<edge>& face::getHalfEdge() const noexcept{
  return m_halfEdge;
}

std::shared_ptr<edge>& face::getHalfEdge() noexcept {
  return m_halfEdge;
}

const std::vector<RealVect> face::getAllVertexCoordinates() const noexcept {
  std::vector<RealVect> pos;
  
  for (const auto& v : m_vertices){
    pos.emplace_back(v->getPosition());
  }
  
  return pos;
}

std::vector<std::shared_ptr<vertex> >& face::getVertices() noexcept {
  return m_vertices;
}

const std::vector<std::shared_ptr<vertex> >& face::getVertices() const noexcept {
  return m_vertices;
}

std::vector<std::shared_ptr<edge> >& face::getEdges() noexcept {
  return m_edges;
}

const std::vector<std::shared_ptr<edge> >& face::getEdges() const noexcept {
  return m_edges;
}

const std::vector<std::shared_ptr<vertex> > face::gatherVertices() const noexcept {
  std::vector<std::shared_ptr<vertex> > vertices;

  for (edge_iterator iter(*this); iter.ok(); ++iter){
    std::shared_ptr<edge>& edge = iter();
    vertices.emplace_back(edge->getVertex());
  }

  return vertices;
}

const std::vector<std::shared_ptr<edge> > face::gatherEdges() const noexcept {
  std::vector<std::shared_ptr<edge> > edges;

  for (edge_iterator iter(*this); iter.ok(); ++iter){
    edges.emplace_back(iter());
  }

  return edges;
}

RealVect& face::getNormal() noexcept {
  return m_normal;
}

const RealVect& face::getNormal() const noexcept {
  return m_normal;
}

RealVect& face::getCentroid() noexcept {
  return m_centroid;
}

const RealVect& face::getCentroid() const noexcept {
  return m_centroid;
}

Real& face::getArea() noexcept {
  return m_area;
}

const Real& face::getArea() const noexcept {
  return m_area;
}

RealVect& face::getBoundingBoxLo() noexcept {
  return m_boundingBox.getLowCorner();
}

const RealVect& face::getBoundingBoxLo() const noexcept {
  return m_boundingBox.getLowCorner();
}

RealVect& face::getBoundingBoxHi() noexcept {
  return m_boundingBox.getHighCorner();
}

const RealVect& face::getBoundingBoxHi() const noexcept {
  return m_boundingBox.getHighCorner();
}

Real face::signedDistance(const RealVect& a_x0) const noexcept {
  Real retval = 1.234567E89;

  // Compute projection of x0 on the face plane
#if 1 // Original code
  const bool inside = this->isPointInsideFaceAngleSum(a_x0);
#else
  const bool inside = this->isPointInsideFaceWindingNumber(a_x0);
#endif

  // Projected point is inside if angles sum to 2*pi
  if(inside){ // Ok, the projection onto the face plane places the point "inside" the planer face
    const RealVect& x1         = m_vertices.front()->getPosition();
    const Real normalComponent = PolyGeom::dot(a_x0-x1, m_normal);
    retval = normalComponent;
  }
  else{ // The projected point lies outside the triangle. Check distance to edges/vertices
    for (const auto& e : m_edges){
      const Real curDist = e->signedDistance(a_x0);
      retval = (std::abs(curDist) < std::abs(retval)) ? curDist : retval;
    }
  }

  return retval;
}

Real face::unsignedDistance2(const RealVect& a_x0) const noexcept {
  std::cerr << "In file 'dcel_face.cpp' function dcel::face::unsignedDistance2 - not implemented!\n";

  return 0.0;
}

RealVect face::projectPointIntoFacePlane(const RealVect& a_p) const noexcept {
  const RealVect& planePoint     = m_vertices.front()->getPosition();
  const RealVect normalComponent = m_normal.dotProduct(a_p - planePoint) * m_normal;
  const RealVect projectedPoint  = a_p - normalComponent;

  return projectedPoint;
}

bool face::isPointInsideFaceAngleSum(const RealVect& a_p) const noexcept {
  bool retval;
  
  const RealVect projectedPoint = this->projectPointIntoFacePlane(a_p);

  Real sum = 0.0;

  constexpr Real thresh = 1.E-6;

  const int N = m_vertices.size();
  
  for (int i = 0; i < N; i++){
    const RealVect p1 = m_vertices[i]      ->getPosition() - projectedPoint;
    const RealVect p2 = m_vertices[(i+1)%N]->getPosition() - projectedPoint;

    const Real m1 = p1.vectorLength();
    const Real m2 = p2.vectorLength();

    if(m1*m2 < thresh){
      return false;
    }

    const Real cosTheta = p1.dotProduct(p2)/(m1*m2);

    sum += acos(cosTheta);
  }

  sum = std::abs(sum)/(2.0*M_PI) - 1.0;

  return std::abs(sum) < thresh;
}

bool face::isPointInsideFaceWindingNumber(const RealVect& a_p) const noexcept {
  return m_poly2->isPointInsidePolygon(a_p);
}
