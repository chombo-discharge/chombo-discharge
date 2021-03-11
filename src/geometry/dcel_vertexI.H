/*!
  @file   dcel_vertexI.H
  @brief  Implementaion of dcel_vertex.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_VERTEXI_
#define _DCEL_VERTEXI_

#include "dcel_vertex.H"

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
