/*!
  @file   ritter_sphere.cpp
  @brief  Implementation of ritter_sphere.H
  @author Robert Marskar
  @date   Ap. 2018
*/

#include "ritter_sphere.H"

ritter_sphere::ritter_sphere(){
  m_radius = 0.0;
  m_center = RealVect::Zero;
}

ritter_sphere::ritter_sphere(const Vector<RealVect>& a_points){
  this->define(a_points);
}

ritter_sphere::~ritter_sphere(){

}

void ritter_sphere::define(const Vector<RealVect>& a_points){
  m_radius = 0.0;
  m_center = RealVect::Zero;


  // INITIAL PASS
  Vector<RealVect> min_coord(SpaceDim, a_points[0]);
  Vector<RealVect> max_coord(SpaceDim, a_points[0]);
  for (int i = 1; i < a_points.size(); i++){
    for (int dir = 0; dir < SpaceDim; dir++){
      if(a_points[i][dir] < min_coord[dir][dir]){
	min_coord[dir] = a_points[i];
      }
      if(a_points[i][dir] > max_coord[dir][dir]){
	max_coord[dir] = a_points[i];
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
      p2 = min_coord[dir];
    }
  }

  m_center = 0.5*(p1+p2);
  m_radius = 0.5*(p2-p1).vectorLength();


  // SECOND PASS
  for (int i = 0; i < a_points.size(); i++){
    const Real dist = (a_points[i]-m_center).vectorLength() - m_radius; // This could be sped up by using a squared distance test
    if(dist > 0){ // Point lies outside
      const RealVect p1 = a_points[i];
      const RealVect p2 = m_center - m_radius*(a_points[i]-m_center)/(a_points[i]-m_center).vectorLength();

      m_center = 0.5*(p2+p1);
      m_radius = 0.5*(p2-p1).vectorLength();
    }
  }
}

bool ritter_sphere::inside(const RealVect a_x0) const{
  return (a_x0 - m_center).vectorLength() <= m_radius;
}

Real ritter_sphere::get_radius() const {
  return m_radius;
}

RealVect ritter_sphere::get_center() const {
  return m_center;
}
