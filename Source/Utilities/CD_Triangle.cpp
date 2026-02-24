/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Triangle.cpp
  @brief  Implementation of CD_Triangle.H
  @author Robert Marskar
*/

// Our includes
#include <CD_Triangle.H>
#include <CD_NamespaceHeader.H>

Triangle::Triangle() noexcept
{
  this->m_metaData = {0.0, 0.0, 0.0};
}

Triangle::Triangle(const std::array<Vec3, 3>& a_vertexPositions) noexcept
{
  this->setVertexPositions(a_vertexPositions);
  this->m_metaData = {0.0, 0.0, 0.0};
}

Triangle::~Triangle() noexcept {

}

void
Triangle::setVertexData(const std::array<Real, 3>& a_vertexData) noexcept {
  m_metaData = a_vertexData;
}

bool
Triangle::isInside(const Vec3& a_point) const noexcept {
  return false;
}

Real
Triangle::interpolate(const Vec3& a_point) const noexcept {
  return 0.0;
}

#include <CD_NamespaceFooter.H>
