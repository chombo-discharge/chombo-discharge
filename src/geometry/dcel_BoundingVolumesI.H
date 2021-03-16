/*!
  @file   dcel_BoundingVolumesI.H
  @brief  Implementation of dcel_BoundingVolumes.H
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_BOUNDINGVOLUMESI_H_
#define _DCEL_BOUNDINGVOLUMESI_H_

#include "dcel_BoundingVolumes.H"

#include <iostream>

namespace dcel {

  template <class T>
  BoundingVolumeT<T>::BoundingVolumeT(){

  }

  template <class T>
  BoundingVolumeT<T>::~BoundingVolumeT(){

  }

  template <class T>
  inline
  decltype(auto) BoundingVolumeT<T>::getPointLocation(const Vec3& a_p) const noexcept {
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

  template <class T>
  inline
  decltype(auto) BoundingVolumeT<T>::isPointInside(const Vec3& a_x0) const {
    constexpr T zero = 0.0;
    
    return this->getDistance(a_x0) < zero;
  }

  template <class T>
  inline
  decltype(auto) BoundingVolumeT<T>::isPointOutside(const Vec3& a_x0) const {
    constexpr T zero = 0.0;
    
    return this->getDistance(a_x0) > zero;
  }

  template <class T>
  inline
  BoundingSphereT<T>::BoundingSphereT(){
    m_radius = 0.0;
    m_center = Vec3::zero();
  }

  template <class T>
  BoundingSphereT<T>::BoundingSphereT(const BoundingSphereT& a_other){
    m_radius = a_other.m_radius;
    m_center = a_other.m_center;
  }

  template <class T>
  BoundingSphereT<T>::BoundingSphereT(const std::vector<Vec3>& a_points, const BoundingVolumeAlgorithm& a_algorithm){
    this->define(a_points, a_algorithm);
  }

  template <class T>
  BoundingSphereT<T>::~BoundingSphereT(){

  }

  template <class T>
  inline
  decltype(auto) BoundingSphereT<T>::define(const std::vector<Vec3>& a_points, const BoundingVolumeAlgorithm& a_algorithm) noexcept {
    switch(a_algorithm) {
    case BoundingVolumeAlgorithm::Ritter:
      this->buildRitter(a_points);
      break;
    default:
      std::cerr << "BoundingSphereT::define - unsupported algorithm requested\n";
    }
  }

  template <class T>
  inline
  decltype(auto) BoundingSphereT<T>::intersects(const BoundingSphereT& a_other) const noexcept {
    const Vec3 deltaV = m_center - a_other.getCenter();
    const T sumR       = m_radius + a_other.getRadius();

    return deltaV.dot(deltaV) < sumR*sumR;
  }

  template <class T>
  inline
  decltype(auto) BoundingSphereT<T>::getRadius() noexcept {
    return (m_radius);
  }

  template <class T>
  inline
  decltype(auto) BoundingSphereT<T>::getRadius() const noexcept {
    return (m_radius);
  }

  template <class T>
  inline
  decltype(auto) BoundingSphereT<T>::getCenter() noexcept {
    return (m_center);
  }

  template <class T>
  inline
  decltype(auto) BoundingSphereT<T>::getCenter() const noexcept {
    return (m_center);
  }

  template <class T>
  inline
  decltype(auto) BoundingSphereT<T>::buildRitter(const std::vector<Vec3>& a_points) noexcept {
    m_radius = 0.0;
    m_center = Vec3::zero();

    constexpr T half = 0.5;

    constexpr int DIM = 3;

    // INITIAL PASS
    std::vector<Vec3> min_coord(DIM, a_points[0]); // [0] = Minimum x, [1] = Minimum y, [2] = Minimum z
    std::vector<Vec3> max_coord(DIM, a_points[0]);
  
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

    T dist = -1;
    Vec3 p1,p2;
    for (int dir = 0; dir < DIM; dir++){
      const T len = (max_coord[dir]-min_coord[dir]).length();
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
      const T dist = (a_points[i]-m_center).length() - m_radius; 
      if(dist > 0){ // Point lies outside
	const Vec3 v  = a_points[i] - m_center;
	const Vec3 p1 = a_points[i];
	const Vec3 p2 = m_center - m_radius*v/v.length();

	m_center = half*(p2+p1);
	m_radius = half*(p2-p1).length();
      }
    }

    // Ritter algorithm is very coarse and does not give an exact result anyways. Grow the dimension for safety. 
    m_radius *= (1.0 + 1E-2);
  }

  template <class T>
  inline
  T BoundingSphereT<T>::getDistance(const Vec3& a_x0) const noexcept {
    constexpr T zero = 0.0;

    return std::max(zero, (a_x0-m_center).length() - m_radius);
  }

  template <class T>
  AABBT<T>::AABBT(){
    m_loCorner = Vec3::zero();
    m_hiCorner = Vec3::zero();
  }

  template <class T>
  AABBT<T>::AABBT(const AABBT<T>& a_other){
    m_loCorner = a_other.m_loCorner;
    m_hiCorner = a_other.m_hiCorner;
  }

  template <class T>
  AABBT<T>::AABBT(const std::vector<AABBT<T> >& a_others) {
    constexpr int DIM = 3;
  
    m_loCorner = a_others.front().getLowCorner();
    m_hiCorner = a_others.front().getHighCorner();

    for (const auto& other : a_others){
      const Vec3& otherLo = other.getLowCorner();
      const Vec3& otherHi = other.getHighCorner();
    
      for (int dir = 0; dir < DIM; dir++){
	m_loCorner[dir] = std::min(m_loCorner[dir], otherLo[dir]);
	m_hiCorner[dir] = std::max(m_hiCorner[dir], otherHi[dir]);
      }
    }
  }

  template <class T>
  AABBT<T>::AABBT(const std::vector<Vec3>& a_points){
    this->define(a_points);
  }

  template <class T>
  AABBT<T>::~AABBT(){

  }

  template <class T>
  inline
  decltype(auto) AABBT<T>::define(const std::vector<Vec3>& a_points) noexcept {
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

  template <class T>
  inline
  decltype(auto) AABBT<T>::intersects(const AABBT& a_other) const noexcept {
    const Vec3& otherLo = a_other.getLowCorner();
    const Vec3& otherHi = a_other.getHighCorner();

    return (m_loCorner[0] < otherHi[0] && m_hiCorner[0] > otherLo[0]) &&
           (m_loCorner[1] < otherHi[1] && m_hiCorner[1] > otherLo[1]) &&
           (m_loCorner[2] < otherHi[2] && m_hiCorner[2] > otherLo[2]);
  }

  template <class T>
  inline
  decltype(auto) AABBT<T>::getLowCorner() noexcept {
    return (m_loCorner);
  }

  template <class T>
  inline
  decltype(auto) AABBT<T>::getLowCorner() const noexcept {
    return (m_loCorner);
  }

  template <class T>
  inline
  decltype(auto) AABBT<T>::getHighCorner() noexcept {
    return (m_hiCorner);
  }

  template <class T>
  inline
  decltype(auto) AABBT<T>::getHighCorner() const noexcept {
    return (m_hiCorner);
  }

  template <class T>
  inline
  T AABBT<T>::getDistance(const Vec3& a_point) const noexcept {
    constexpr T zero = 0.0;
    
    const Vec3 delta = Vec3(std::max(m_loCorner[0] - a_point[0], a_point[0] - m_hiCorner[0]),
			    std::max(m_loCorner[1] - a_point[1], a_point[1] - m_hiCorner[1]),
			    std::max(m_loCorner[2] - a_point[2], a_point[2] - m_hiCorner[2]));

    const T retval =  std::min(zero, delta[delta.maxDir(false)]) + max(Vec3::zero(), delta).length(); 
  
    return std::max(zero, retval);
  }

}

#endif
