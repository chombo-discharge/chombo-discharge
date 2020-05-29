/*!
  @file   point_mass.cpp
  @brief  Implementation of point_mass.H
  @author Robert Marskar
  @date   May 2020
*/

#include "point_mass.H"

point_mass::point_mass(){

}

point_mass::point_mass(const RealVect a_pos, const Real a_mass){
  m_pos  = a_pos;
  m_mass = a_mass;
}

point_mass::~point_mass(){

}
