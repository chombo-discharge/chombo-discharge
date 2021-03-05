/*!
  @file   dcel_poly.cpp
  @brief  Implementation of dcel_poly.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_polygon.H"
#include "dcel_iterator.H"

#include <PolyGeom.H>
  
#define EPSILON 1.E-8
#define TWOPI 6.283185307179586476925287

using namespace dcel;

polygon::polygon(){
  m_normal = RealVect::Zero;
}

polygon::polygon(const polygon& a_otherPolygon){
  this->define(a_otherPolygon.getNormal(),
	       a_otherPolygon.getEdge());
}

polygon::~polygon(){

}

void polygon::define(const RealVect& a_normal, const std::shared_ptr<edge>& a_edge){
  this->setNormal(a_normal);
  this->setEdge(a_edge);
}

void polygon::setEdge(const std::shared_ptr<edge>& a_edge){
  m_edge = a_edge;
}

void polygon::setNormal(const RealVect& a_normal){
  m_normal = a_normal;
}





void polygon::normalizeNormalVector(){
  m_normal *= 1./m_normal.vectorLength();
}


void polygon::computeArea() {
  const std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();

  Real area = 0.0;

  for (int i = 0; i < vertices.size() - 1; i++){
    const RealVect& v1 = vertices[i]->getPosition();
    const RealVect& v2 = vertices[i+1]->getPosition();
    area += PolyGeom::dot(PolyGeom::cross(v2,v1), m_normal);
  }

  m_area = Abs(0.5*area);
}


void polygon::computeCentroid() {
  m_centroid = RealVect::Zero;
  
  const std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();

  for (const auto& v : vertices){
    m_centroid += v->getPosition();
  }
  
  m_centroid = m_centroid/vertices.size();
}

void polygon::computeNormal(const bool a_outwardNormal){
  
  // TLDR: We assume that the normal is defined by right-hand rule where the rotation direction is along the half edges
  
  bool found_normal = false;
  
  std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();


  // Funky code - I guess we do this since some cross products don't exist...?
  const int n = vertices.size();
  for (int i = 0; i < n; i++){
    const RealVect& x0 = vertices[i]      ->getPosition();
    const RealVect& x1 = vertices[(i+1)%n]->getPosition();
    const RealVect& x2 = vertices[(i+2)%n]->getPosition();

    m_normal = PolyGeom::cross(x2-x1, x2-x0);
    
    if(m_normal.vectorLength() > 0.0){
      found_normal = true;
      break;
    }
  }

  this->normalizeNormalVector();

  
  if(!found_normal){
    pout() << "polygon::compute_normal - vertex vectors:" << endl;
    for (int i = 0; i < vertices.size(); i++){
      pout() << "\t" << vertices[i]->getPosition() << endl;
    }
    pout() << "polygon::compute_normal - From this I computed n = " << m_normal << endl;
    pout() << "polygon::compute_normal - Aborting..." << endl;
    MayDay::Warning("polygon::compute_normal - Cannot compute normal vector. The polygon is probably degenerate");
  }

#if 0
  const std::shared_ptr<vertex>& v0 = m_edge->getPreviousEdge()->getVertex();
  const std::shared_ptr<vertex>& v1 = m_edge->getVertex();
  const std::shared_ptr<vertex>& v2 = m_edge->getNextEdge()->getVertex();
  
  const RealVect& x0 = v0->getPosition();
  const RealVect& x1 = v1->getPosition();
  const RealVect& x2 = v2->getPosition();
  
  m_normal = PolyGeom::cross(x2-x1,x1-x0);
  if(m_normal.vectorLength() < 1.E-40){
    MayDay::Abort("polygon::compute_normal - vertices lie on a line. Cannot compute normal vector");
  }
  else{
    m_normal = m_normal/m_normal.vectorLength();
  }
#endif

  if(!a_outwardNormal){    // If normal points inwards, make it point outwards
    m_normal = -m_normal;
  }
}

void polygon::computeBoundingSphere(){
  const std::vector<RealVect> points = this->getAllVertexCoordinates();

  m_boundingSphere.define(points, BoundingSphere::Algorithm::Ritter);
}

void polygon::computeBoundingBox(){
  std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();
  std::vector<RealVect> coords;

#if 1 // New code
  for (const auto& v : vertices){
    coords.emplace_back(v->getPosition());
  }
#else
  // original code
  for (int i = 0; i < vertices.size(); i++){
    coords.push_back(vertices[i]->getPosition());
  }
#endif

  m_lo =  1.23456E89*RealVect::Unit;
  m_hi = -1.23456E89*RealVect::Unit;

  for (int i = 0; i < coords.size(); i++){
    for (int dir = 0; dir < SpaceDim; dir++){
      if(coords[i][dir] < m_lo[dir]){
	m_lo[dir] = coords[i][dir];
      }
      if(coords[i][dir] > m_hi[dir]){
	m_hi[dir] = coords[i][dir];
      }
    }
  }

  // Define ritter sphere
  //  m_ritter.define(coords);

#if 1 // Debug test
  for (int i = 0; i < vertices.size(); i++){
    const RealVect pos = vertices[i]->getPosition();
    for (int dir = 0; dir < SpaceDim; dir++){
      if(pos[dir] < m_lo[dir] || pos[dir] > m_hi[dir]){
	pout() << "pos = " << pos << "\t Lo = " << m_lo << "\t Hi = " << m_hi << endl;
	MayDay::Abort("polygon::compute_bbox - Box does not contain vertices");
      }
    }
  }
#endif

#if 0 // Disabled because I want a tight-fitting box
  Real widest = 0;
  for (int dir = 0; dir < SpaceDim; dir++){
    const Real cur = m_hi[dir] - m_lo[dir];
    widest = (cur > widest) ? cur : widest;
  }


  // Grow box by 5% in each direction
  m_hi += 5.E-2*widest*RealVect::Unit;
  m_lo -= 5.E-2*widest*RealVect::Unit;
#endif
}


const std::shared_ptr<edge>& polygon::getEdge() const{
  return m_edge;
}

std::shared_ptr<edge>& polygon::getEdge(){
  return m_edge;
}

std::vector<RealVect> polygon::getAllVertexCoordinates(){
  std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();

  std::vector<RealVect> pos;
  
  for (const auto& v : vertices){
    pos.emplace_back(v->getPosition());
  }
  
  return pos;
}

std::vector<std::shared_ptr<vertex> > polygon::getVertices(){
  std::vector<std::shared_ptr<vertex> > vertices;

  for (edge_iterator iter(*this); iter.ok(); ++iter){
    std::shared_ptr<edge>& edge = iter();
    vertices.push_back(edge->getVertex());
  }

  return vertices;
}

std::vector<std::shared_ptr<edge> > polygon::getEdges(){
  std::vector<std::shared_ptr<edge> > edges;

  for (edge_iterator iter(*this); iter.ok(); ++iter){
    edges.push_back(iter());
  }

  return edges;
}

RealVect& polygon::getNormal() {
  return m_normal;
}

const RealVect& polygon::getNormal() const {
  return m_normal;
}


RealVect& polygon::getCentroid() {
  return m_centroid;
}

const RealVect& polygon::getCentroid() const {
  return m_centroid;
}

Real& polygon::getArea() {
  return m_area;
}

const Real& polygon::getArea() const {
  return m_area;
}

RealVect& polygon::getBoundingBoxLo() {
  return m_lo;
}

const RealVect& polygon::getBoundingBoxLo() const {
  return m_lo;
}


RealVect& polygon::getBoundingBoxHi() {
  return m_hi;
}

const RealVect& polygon::getBoundingBoxHi() const {
  return m_hi;
}


Real polygon::signedDistance(const RealVect a_x0) {
#define bug_check 0
  Real retval = 1.234567E89;

  std::vector<std::shared_ptr<vertex> > vertices = this->getVertices();

#if bug_check // Debug, return shortest distance to vertex
  CH_assert(vertices.size() > 0);
  Real min = 1.E99;
  for (int i = 0; i < vertices.size(); i++){
    const Real d = (a_x0 - vertices[i]->getPosition()).vectorLength();
    min = (d < min) ? d : min;
  }

  return min;
#endif

  // Compute projection of x0 on the polygon plane
  const RealVect x1 = vertices[0]->getPosition();
  const Real ncomp  = PolyGeom::dot(a_x0-x1, m_normal);
  const RealVect xp = a_x0 - ncomp*m_normal;


  // Use angle rule to check if projected point lies inside the polygon
  Real anglesum = 0.0;
  const int n = vertices.size();
  for(int i = 0; i < n; i++){
    const RealVect p1 = vertices[i]->getPosition() - xp;
    const RealVect p2 = vertices[(i+1)%n]->getPosition() - xp;

    const Real m1 = p1.vectorLength();
    const Real m2 = p2.vectorLength();

    if(m1*m2 <= EPSILON) {// Projected point hits a vertex, return early. 
      anglesum = 0.;
      break;
    }
    else {
      const Real cos_theta = PolyGeom::dot(p1, p2)/(m1*m2);
      anglesum += acos(cos_theta);
    }
  }

  // Projected point is inside if angles sum to 2*pi
  bool inside = false;
  if(Abs(Abs(anglesum) - TWOPI) < EPSILON){
    inside = true;
  }

  // If projection is inside, shortest distance is the normal component of the point
  if(inside){
#if bug_check
    CH_assert(Abs(ncomp) <= min);
#endif
    retval = ncomp;
  }
  else{ // The projected point lies outside the triangle. Check distance to edges/vertices
    const std::vector<std::shared_ptr<edge> > edges = this->getEdges();
    for (int i = 0; i < edges.size(); i++){
      const Real cur_dist = edges[i]->signedDistance(a_x0);
      if(Abs(cur_dist) < Abs(retval)){
	retval = cur_dist;
      }
    }
  }

  return retval;
}
