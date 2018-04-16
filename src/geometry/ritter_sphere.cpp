/*!
  @file   ritter_sphere.cpp
  @brief  Implementation of ritter_sphere.H
  @author Robert Marskar
  @date   Ap. 2018
*/

#include "ritter_sphere.H"

#include <ParmParse.H>

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
  Vector<RealVect> min_coord(SpaceDim, a_points[0]); // [0] = Minimum x, [1] = Minimum y, [2] = Minimum z
  Vector<RealVect> max_coord(SpaceDim, a_points[0]);
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

  m_center = 0.5*(p1+p2);
  m_radius = 0.5*(p2-p1).vectorLength();


  // SECOND PASS
  for (int i = 0; i < a_points.size(); i++){
    const Real dist = (a_points[i]-m_center).vectorLength() - m_radius; 
    if(dist > 0){ // Point lies outside
      const RealVect v  = a_points[i] - m_center;
      const RealVect p1 = a_points[i];
      const RealVect p2 = m_center - m_radius*v/v.vectorLength();

      m_center = 0.5*(p2+p1);
      m_radius = 0.5*(p2-p1).vectorLength();
    }
  }

  // Safety
  m_radius *= (1.0 + 1.E-10);

#if 1 // Debug
  for (int i = 0; i < a_points.size(); i++){
    const Real dist = (a_points[i] - m_center).vectorLength() - m_radius;
    if(dist > 0.0){
      pout() << "point = " << a_points[i] << endl;
      MayDay::Abort("ritter_sphere::define - point lies outside sphere!");
    }
  }
#endif
}

bool ritter_sphere::inside(const RealVect a_x0) const{
  return (a_x0 - m_center).vectorLength() <= m_radius;
}

Real ritter_sphere::get_radius() const {
  return m_radius;
}

Real ritter_sphere::dist(const RealVect a_x0) const{
  return (a_x0 - m_center).vectorLength() - m_radius;
}

RealVect ritter_sphere::get_center() const {
  return m_center;
}
