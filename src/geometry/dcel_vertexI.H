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
