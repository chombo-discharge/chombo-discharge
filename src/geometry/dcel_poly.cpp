/*!
  @file   dcel_poly.cpp
  @brief  Implementation of dcel_poly.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_vert.H"
#include "dcel_edge.H"
#include "dcel_poly.H"

dcel_poly::dcel_poly(){
  m_edge = NULL;
}

dcel_poly::~dcel_poly(){

}

const dcel_edge* dcel_poly::get_edge() const{
  return m_edge;
}

RealVect dcel_poly::get_normal() const {
  return m_normal;
}

Real dcel_poly::signed_distance(const RealVect a_x0) const{
  return 0.0;
}
