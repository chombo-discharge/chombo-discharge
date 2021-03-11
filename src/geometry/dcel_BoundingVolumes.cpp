/*!
  @file   dcel_BoundingVolumes.cpp
  @brief  Implementation of dcel_BoundingSphere.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_BoundingVolumes.H"

#include <iostream>

using namespace dcel;

BoundingVolume::BoundingVolume(){

}

BoundingVolume::~BoundingVolume(){

}

auto BoundingVolume::getPointLocation(const Vec3& a_p) const noexcept -> PointLocation {
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

auto BoundingVolume::isPointInside(const Vec3& a_x0) const -> bool {
  return this->getDistanceToPoint(a_x0) < 0.0;
}

auto BoundingVolume::isPointOutside(const Vec3& a_x0) const -> bool {
  return this->getDistanceToPoint(a_x0) > 0.0;
}

BoundingSphere::BoundingSphere(){
  m_radius = 0.0;
  m_center = Vec3::Zero;
}

BoundingSphere::BoundingSphere(const BoundingSphere& a_other){
  m_radius = a_other.m_radius;
  m_center = a_other.m_center;
}

BoundingSphere::BoundingSphere(const std::vector<Vec3 >& a_points, const Algorithm& a_algorithm){
  this->define(a_points, a_algorithm);
}

BoundingSphere::~BoundingSphere(){

}

auto BoundingSphere::define(const std::vector<Vec3 >& a_points, const Algorithm& a_algorithm) noexcept -> void {
  switch(a_algorithm) {
  case Algorithm::Ritter:
    this->buildRitter(a_points);
    break;
  default:
    std::cerr << "BoundingSphere::define - unknown algorithm requested\n";
  }
}

auto BoundingSphere::intersects(const BoundingSphere& a_other) const noexcept -> bool {
  const Vec3 deltaV = m_center - a_other.getCenter();
  const double sumR       = m_radius + a_other.getRadius();

  return deltaV.dot(deltaV) < sumR*sumR;
}

auto BoundingSphere::getDistanceToPoint(const Vec3& a_x0) const noexcept -> double {
  return (a_x0 - m_center).length() - m_radius;
}

auto  BoundingSphere::getRadius() noexcept -> double& {
  return m_radius;
}

auto BoundingSphere::getRadius() const noexcept -> const double& {
  return m_radius;
}

auto BoundingSphere::getCenter() noexcept -> Vec3& {
  return m_center;
}

auto  BoundingSphere::getCenter() const noexcept -> const Vec3& {
  return m_center;
}

auto BoundingSphere::buildRitter(const std::vector<Vec3>& a_points) noexcept -> void {
  m_radius = 0.0;
  m_center = Vec3::Zero;

  constexpr double half = 0.5;

  constexpr int DIM = 3;

  // INITIAL PASS
  std::vector<Vec3 > min_coord(DIM, a_points[0]); // [0] = Minimum x, [1] = Minimum y, [2] = Minimum z
  std::vector<Vec3 > max_coord(DIM, a_points[0]);
  
  for (int i = 1; i < a_points.size(); i++){
    for (int dir = 0; dir < DIM; dir++){
      Vec3& min = min_coord[dir];
      Vec3& max = max_coord[dir];
      
      if(a_points[i][dir] < min[dir]){
	min = a_points[i];
      }
      if(a_points[i][dir] > max[dir]){
	max = a_points[i];
      }
    }
  }

  double dist = -1;
  Vec3 p1,p2;
  for (int dir = 0; dir < DIM; dir++){
    const double len = (max_coord[dir]-min_coord[dir]).length();
    if(len > dist ){
      dist = len;
      p1 = min_coord[dir];
      p2 = max_coord[dir];
    }
  }

  //  m_center = half*(p1+p2);
  m_center = (p1+p2)*half;
  m_radius = half*(p2-p1).length();


  // SECOND PASS
  for (int i = 0; i < a_points.size(); i++){
    const double dist = (a_points[i]-m_center).length() - m_radius; 
    if(dist > 0){ // Point lies outside
      const Vec3 v  = a_points[i] - m_center;
      const Vec3 p1 = a_points[i];
      const Vec3 p2 = m_center - m_radius*v/v.length();

      m_center = half*(p2+p1);
      m_radius = half*(p2-p1).length();
    }
  }

  // Ritter algorithm is very coarse and does not give an exact result anyways. Grow the dimension for safety. 
  m_radius *= (1.0 + 1.E-3);
}

AABB::AABB(){
  m_loCorner = Vec3::Zero;
  m_hiCorner = Vec3::Zero;
}

AABB::AABB(const AABB& a_other){
  m_loCorner = a_other.m_loCorner;
  m_hiCorner = a_other.m_hiCorner;
}

AABB::AABB(const std::vector<AABB>& a_others) {
  constexpr int DIM = 3;
  
  m_loCorner = a_others.front().getLowCorner();
  m_hiCorner = a_others.front().getHighCorner();

  for (const auto& other : a_others){
    const Vec3& otherLo = other.getLowCorner();
    const Vec3& otherHi = other.getHighCorner();
    
    for (int dir = 0; dir < DIM; dir++){
      m_loCorner[dir] = std::min(m_loCorner[dir], otherLo[dir]);
      m_hiCorner[dir] = std::min(m_hiCorner[dir], otherHi[dir]);
    }
  }
}

AABB::AABB(const std::vector<Vec3 >& a_points){
  this->define(a_points);
}

AABB::~AABB(){

}

auto AABB::define(const std::vector<Vec3 >& a_points) noexcept -> void{
  constexpr int DIM = 3;
  m_loCorner = a_points.front();
  m_hiCorner = a_points.front();

  for (const auto& p : a_points){
    for (int dir = 0; dir < 3; dir++){
      m_loCorner[dir] = std::min(m_loCorner[dir], p[dir]);
      m_hiCorner[dir] = std::max(m_hiCorner[dir], p[dir]);
    }
  }
}

auto AABB::intersects(const AABB& a_other) const noexcept -> bool {
  const Vec3& otherLo = a_other.getLowCorner();
  const Vec3& otherHi = a_other.getHighCorner();


  return (m_loCorner[0] < otherHi[0] && m_hiCorner[0] > otherLo[0]) &&
         (m_loCorner[1] < otherHi[1] && m_hiCorner[1] > otherLo[1]) &&
         (m_loCorner[2] < otherHi[2] && m_hiCorner[2] > otherLo[2]);
}
    
auto AABB::getDistanceToPoint(const Vec3& a_point) const noexcept -> double {
  const Vec3 delta = Vec3(std::max(m_loCorner[0] - a_point[0], a_point[0] - m_hiCorner[0]),
					  std::max(m_loCorner[1] - a_point[1], a_point[1] - m_hiCorner[1]),
					  std::max(m_loCorner[2] - a_point[2], a_point[2] - m_hiCorner[2]));

  const double retval =  std::min(0.0, delta[delta.maxDir(false)]) + std::max(Vec3::Zero, delta).length(); // This is negative inside. 
  
  return retval;
}

auto AABB::getLowCorner() noexcept -> Vec3& {
  return m_loCorner;
}

auto AABB::getLowCorner() const noexcept -> const Vec3& {
  return m_loCorner;
}

auto AABB::getHighCorner() noexcept -> Vec3& {
  return m_hiCorner;
}

auto AABB::getHighCorner() const noexcept -> const Vec3& {
  return m_hiCorner;
}
