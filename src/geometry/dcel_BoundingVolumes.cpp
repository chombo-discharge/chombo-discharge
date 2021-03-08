/*!
  @file   dcel_BoundingVolumes.cpp
  @brief  Implementation of dcel_BoundingSphere.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_BoundingVolumes.H"

using namespace dcel;

PointLocation BoundingVolume::getPointLocation(const RealVect& a_p) const{
  const bool inside  = this->isPointInside(a_p);
  const bool outside = this->isPointOutside(a_p);

  PointLocation ret;

  if(inside){
    ret = PointLocation::Inside;
  }
  else if(outside){
    ret = PointLocation::Outside;
  }
  else{
    ret = PointLocation::Boundary;
  }

  return ret;
}

bool BoundingVolume::isPointInside(const RealVect& a_x0) const{
  return this->getDistanceToPoint(a_x0) < 0.0;
}

bool BoundingVolume::isPointOutside(const RealVect& a_x0) const{
  return this->getDistanceToPoint(a_x0) > 0.0;
}

BoundingSphere::BoundingSphere(){
  m_radius = 0.0;
  m_center = RealVect::Zero;
}

BoundingSphere::BoundingSphere(const BoundingSphere& a_other){
  m_radius = a_other.m_radius;
  m_center = a_other.m_center;
}

BoundingSphere::BoundingSphere(const std::vector<RealVect>& a_points, const Algorithm& a_algorithm){
  this->define(a_points, a_algorithm);
}

BoundingSphere::~BoundingSphere(){

}

void BoundingSphere::define(const std::vector<RealVect>& a_points, const Algorithm& a_algorithm){
  switch(a_algorithm) {
  case Algorithm::Ritter:
    this->buildRitter(a_points);
    break;
  default:
    std::cerr << "BoundingSphere::define - unknown algorithm requested\n";
  }
}

Real BoundingSphere::getDistanceToPoint(const RealVect& a_x0) const {
  return (a_x0 - m_center).vectorLength() - m_radius;
}

Real& BoundingSphere::getRadius() {
  return m_radius;
}

const Real& BoundingSphere::getRadius() const {
  return m_radius;
}

RealVect& BoundingSphere::getCenter() {
  return m_center;
}

const RealVect& BoundingSphere::getCenter() const {
  return m_center;
}

void BoundingSphere::buildRitter(const std::vector<RealVect>& a_points){
  m_radius = 0.0;
  m_center = RealVect::Zero;

  constexpr Real half = 0.5;

  // INITIAL PASS
  std::vector<RealVect> min_coord(SpaceDim, a_points[0]); // [0] = Minimum x, [1] = Minimum y, [2] = Minimum z
  std::vector<RealVect> max_coord(SpaceDim, a_points[0]);
  
  for (int i = 1; i < a_points.size(); i++){
    for (int dir = 0; dir < SpaceDim; dir++){
      RealVect& min = min_coord[dir];
      RealVect& max = max_coord[dir];
      
      if(a_points[i][dir] < min[dir]){
	min = a_points[i];
      }
      if(a_points[i][dir] > max[dir]){
	max = a_points[i];
      }
    }
  }

  Real dist = -1;
  RealVect p1,p2;
  for (int dir = 0; dir < SpaceDim; dir++){
    const Real len = (max_coord[dir]-min_coord[dir]).vectorLength();
    if(len > dist ){
      dist = len;
      p1 = min_coord[dir];
      p2 = max_coord[dir];
    }
  }

  m_center = half*(p1+p2);
  m_radius = half*(p2-p1).vectorLength();


  // SECOND PASS
  for (int i = 0; i < a_points.size(); i++){
    const Real dist = (a_points[i]-m_center).vectorLength() - m_radius; 
    if(dist > 0){ // Point lies outside
      const RealVect v  = a_points[i] - m_center;
      const RealVect p1 = a_points[i];
      const RealVect p2 = m_center - m_radius*v/v.vectorLength();

      m_center = half*(p2+p1);
      m_radius = half*(p2-p1).vectorLength();
    }
  }

  // Ritter algorithm is very coarse and does not give an exact result anyways. Grow the dimension for safety. 
  m_radius *= (1.0 + 1.E-3);
}

AABB::AABB(){
  m_loCorner = RealVect::Zero;
  m_hiCorner = RealVect::Zero;
}

AABB::AABB(const AABB& a_other){
  m_loCorner = a_other.m_loCorner;
  m_hiCorner = a_other.m_hiCorner;
}

AABB::AABB(const std::vector<AABB>& a_others) {
  m_loCorner = a_others.front().getLowCorner();
  m_hiCorner = a_others.front().getHighCorner();

  for (const auto& other : a_others){
    const RealVect& otherLo = other.getLowCorner();
    const RealVect& otherHi = other.getHighCorner();
    
    for (int dir = 0; dir < SpaceDim; dir++){
      m_loCorner[dir] = std::min(m_loCorner[dir], otherLo[dir]);
      m_hiCorner[dir] = std::min(m_hiCorner[dir], otherHi[dir]);
    }
  }
}

AABB::AABB(const std::vector<RealVect>& a_points){
  this->define(a_points);
}

AABB::~AABB(){

}

void AABB::define(const std::vector<RealVect>& a_points){
  m_loCorner = a_points.front();
  m_hiCorner = a_points.front();

  for (const auto& p : a_points){
    for (int dir = 0; dir < SpaceDim; dir++){
      m_loCorner[dir] = std::min(m_loCorner[dir], p[dir]);
      m_hiCorner[dir] = std::max(m_hiCorner[dir], p[dir]);
    }
  }
}
    
Real AABB::getDistanceToPoint(const RealVect& a_point) const {
  const RealVect delta = RealVect(D_DECL(Max(m_loCorner[0] - a_point[0], a_point[0] - m_hiCorner[0]),
					 Max(m_loCorner[1] - a_point[1], a_point[1] - m_hiCorner[1]),
					 Max(m_loCorner[2] - a_point[2], a_point[2] - m_hiCorner[2])));

  const Real retval =  Min(0.0, delta[delta.maxDir(false)]) + max(RealVect::Zero, delta).vectorLength(); // This is negative inside. 
  
  return retval;
}

RealVect& AABB::getLowCorner(){
  return m_loCorner;
}

const RealVect& AABB::getLowCorner() const{
  return m_hiCorner;
}

RealVect& AABB::getHighCorner(){
  return m_loCorner;
}

const RealVect& AABB::getHighCorner() const{
  return m_hiCorner;
}
