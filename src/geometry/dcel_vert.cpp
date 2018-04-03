/*!
  @file   dcel_vert.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_vert.H"
#include "dcel_edge.H"
#include "dcel_poly.H"

dcel_vert::dcel_vert(){
  m_pos    = RealVect::Zero;
  m_normal = RealVect::Zero;
  m_edge   = NULL;
}

dcel_vert::~dcel_vert(){

}

const dcel_edge* dcel_vert::get_edge() const{
  return m_edge;
}

RealVect dcel_vert::get_pos() const {
  return m_pos;
}

RealVect dcel_vert::get_normal() const {
  return m_normal;
}
  
Real dcel_vert::signed_distance(const RealVect a_x0) const {
  return 0.0;
}
