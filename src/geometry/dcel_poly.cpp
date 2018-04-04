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
  m_normal = RealVect::Zero;
  m_edge   = NULL;
}

dcel_poly::~dcel_poly(){

}

void dcel_poly::define(const RealVect a_normal, const dcel_edge* const a_edge){
  m_normal = a_normal;
  m_edge   = a_edge;
}

const dcel_edge* dcel_poly::get_edge() const{
  return m_edge;
}

RealVect dcel_poly::get_normal() const {
  return m_normal;
}
