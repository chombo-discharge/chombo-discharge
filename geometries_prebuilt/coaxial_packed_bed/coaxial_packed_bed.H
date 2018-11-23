/*!
  @file   coaxial_packed_bed.H
  @brief  Coaxial cable geometry
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _COAXIAL_PACKED_BED_
#define _COAXIAL_PACKED_BED_

#include "computational_geometry.H"

/*!
  @brief This is a class that defines a coaxial cable geometry with packed bed spheres
*/
class coaxial_packed_bed : public computational_geometry {
public:
  coaxial_packed_bed(); 
  ~coaxial_packed_bed();
protected:
  void reseed(); 
  Real random();
  void build_conductors();
  void build_packed_bed();

  Real m_rad_outer;
  Real m_rad_inner;
  Real m_rad_sphere;
  Real m_min_dist;
  int  m_num_spheres;
  
};

#endif
