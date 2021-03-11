/*!
  @file   dcel_vertexI.H
  @brief  Implementaion of dcel_vertex.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_VERTEXI_
#define _DCEL_VERTEXI_

#include "dcel_vertex.H"
#include "dcel_edge.H"
#include "dcel_face.H"
#include "dcel_iterator.H"

namespace dcel {

  inline
  vertex::vertex(){
    m_pos    = RealVect::Zero;
    m_normal = RealVect::Zero;

    m_faceCache.resize(0);
  }

  inline
  vertex::vertex(const RealVect& a_pos){
    m_pos    = a_pos;
    m_normal = RealVect::Zero;
  }

  inline
  vertex::vertex(const RealVect& a_pos, const RealVect& a_normal){
    m_pos    = a_pos;
    m_normal = a_normal;
  }

  inline
  vertex::vertex(const vertex& a_otherVertex){
    this->define(a_otherVertex.getPosition(),
		 a_otherVertex.getEdge(),
		 a_otherVertex.getNormal());
  }

  inline
  vertex::~vertex(){

  }

  inline
  void vertex::define(const RealVect& a_pos, const std::shared_ptr<edge>& a_edge, const RealVect a_normal) noexcept {
    m_pos = a_pos;
    m_edge = a_edge;
    m_normal = a_normal;
  }

  inline
  void vertex::setPosition(const RealVect& a_pos) noexcept {
    m_pos = a_pos;
  }

  inline
  void vertex::setEdge(const std::shared_ptr<edge>& a_edge) noexcept {
    m_edge = a_edge;
  }

  inline
  void vertex::setNormal(const RealVect& a_normal) noexcept {
    m_normal = a_normal;
  }

  inline
  RealVect& vertex::getPosition() noexcept {
    return m_pos;
  }

  inline
  const RealVect& vertex::getPosition() const noexcept {
    return m_pos;
  }

  inline
  std::shared_ptr<edge>& vertex::getEdge() noexcept {
    return m_edge;
  }

  inline
  const std::shared_ptr<edge>& vertex::getEdge() const noexcept {
    return m_edge;
  }

  inline
  RealVect& vertex::getNormal() noexcept {
    return m_normal;
  }

  inline
  const RealVect& vertex::getNormal() const noexcept {
    return m_normal;
  }
  
}
#endif
