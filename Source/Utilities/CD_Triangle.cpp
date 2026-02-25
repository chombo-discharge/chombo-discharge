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
  this->setVertexData({0.0, 0.0, 0.0});
  this->computeArea();
  this->computeNormal();
}

Triangle::~Triangle() noexcept {

}

void
Triangle::setVertexData(const std::array<Real, 3>& a_vertexData) noexcept
{
  m_metaData = a_vertexData;
}

void
Triangle::computeArea() noexcept
{
  const Vec3 x = m_vertexPositions[1] - m_vertexPositions[0];
  const Vec3 y = m_vertexPositions[2] - m_vertexPositions[0];

  this->m_area = 0.5 * (x.cross(y)).length();
}

bool
Triangle::isInside(const Vec3& a_point) const noexcept
{
  auto sgn = [](const Real x) -> int {
    return (x >= 0.0) ? 1 : -1;
  };

  const Vec3 v21 = m_vertexPositions[1] - m_vertexPositions[0];
  const Vec3 v32 = m_vertexPositions[2] - m_vertexPositions[1];
  const Vec3 v13 = m_vertexPositions[0] - m_vertexPositions[2];

  const Vec3 p1 = a_point - m_vertexPositions[0];
  const Vec3 p2 = a_point - m_vertexPositions[1];
  const Vec3 p3 = a_point - m_vertexPositions[2];

  const Real s0 = sgn(dot(v21.cross(m_triangleNormal), p1));
  const Real s1 = sgn(dot(v32.cross(m_triangleNormal), p2));
  const Real s2 = sgn(dot(v13.cross(m_triangleNormal), p3));

  return (s0 + s1 + s2 >= 2.0);
}

Real
Triangle::interpolate(const Vec3& a_point) const noexcept
{
  return 0.0;
}



#include <CD_NamespaceFooter.H>
