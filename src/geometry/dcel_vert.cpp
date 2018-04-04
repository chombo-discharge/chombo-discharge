/*!
  @file   dcel_vert.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_vert.H"
#include "dcel_edge.H"
#include "dcel_poly.H"
#include "dcel_iterator.H"

dcel_vert::dcel_vert(){
  m_pos    = RealVect::Zero;
  m_normal = RealVect::Zero;
  m_edge   = NULL;
}

dcel_vert::~dcel_vert(){

}

void dcel_vert::define(const RealVect a_pos, const dcel_edge* const a_edge, const RealVect a_normal){
  m_pos    = a_pos;
  m_normal = a_normal;
  m_edge   = a_edge;
}

void dcel_vert::set_normal(const RealVect a_normal){
  m_normal = a_normal;
}

Vector<const dcel_poly*> dcel_vert::get_polygons() const{
  Vector<const dcel_poly*> polygons;
  for (edge_iterator iter(this); iter.ok(); ++iter){
    polygons.push_back(iter()->get_poly());
  }

  return polygons;
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
