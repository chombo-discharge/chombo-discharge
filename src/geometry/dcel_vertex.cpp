/*!
  @file   vertex.cpp
  @brief  Implementation of vertex.H
  @author Robert Marskar
  @date   March 2021
*/

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_polygon.H"
#include "dcel_iterator.H"

using namespace dcel;

vertex::vertex(){
  m_pos    = RealVect::Zero;
  m_normal = RealVect::Zero;

  m_polycache.resize(0);
}

vertex::vertex(const RealVect& a_pos){
  m_pos    = a_pos;
  m_normal = RealVect::Zero;
}

vertex::~vertex(){

}

void vertex::define(const RealVect& a_pos, const std::shared_ptr<edge>& a_edge, const RealVect a_normal){
  this->set_pos(a_pos);
  this->set_edge(a_edge);
  this->set_normal(a_normal);
}

void vertex::set_pos(const RealVect& a_pos){
  m_pos = a_pos;
}

void vertex::set_edge(const std::shared_ptr<edge>& a_edge){
  m_edge = a_edge;
}

void vertex::set_normal(const RealVect& a_normal){
  m_normal = a_normal;
}

void vertex::add_polygon(const std::shared_ptr<polygon>& a_poly){
  m_polycache.push_back(a_poly);
}

std::vector<std::shared_ptr<polygon> > vertex::get_polygons() {
  std::vector<std::shared_ptr<polygon> > polygons;
  for (edge_iterator iter(*this); iter.ok(); ++iter){
    polygons.push_back(iter()->get_poly());
  }

  return polygons;
}

const std::vector<std::shared_ptr<polygon> >& vertex::get_polycache() const{
  return m_polycache;
}

std::vector<std::shared_ptr<polygon> >& vertex::get_polycache(){
  return m_polycache;
}

const std::shared_ptr<edge>& vertex::get_edge() const{
  return m_edge;
}

std::shared_ptr<edge>& vertex::get_edge(){
  return m_edge;
}

const RealVect& vertex::get_pos() const {
  return m_pos;
}

const RealVect& vertex::get_normal() const {
  return m_normal;
}

Real vertex::signed_distance(const RealVect a_x0) const {
  return 0.0;
}
