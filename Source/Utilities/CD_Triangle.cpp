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
}

Triangle::~Triangle() noexcept
{}

void
Triangle::setVertexPositions(const std::array<Vec3, 3>& a_vertexPositions) noexcept
{
  EBGeometry::Triangle<Real, std::array<Real, 3>>::setVertexPositions(a_vertexPositions);
}

void
Triangle::setVertexData(const std::array<Real, 3>& a_vertexData) noexcept
{
  m_metaData = a_vertexData;
}

void
Triangle::computeArea() noexcept
{
  m_area = computeTriangleArea(m_vertexPositions[0], m_vertexPositions[1], m_vertexPositions[2]);
}

EBGeometry::Vec3T<Real>
Triangle::projectToTrianglePlane(const Vec3& a_point) const noexcept
{
#warning "Remove later"
  std::cout << "meta = " << m_metaData[0] << "\t" << m_metaData[1] << "\t" << m_metaData[2] << "\n";

  return a_point - dot(a_point - m_vertexPositions[0], m_triangleNormal) * m_triangleNormal;
}

bool
Triangle::isInside(const Vec3& a_point) const noexcept
{
  auto sgn = [](const Real x) -> int {
    return (x >= 0.0) ? 1 : -1;
  };

  const Vec3 A = m_vertexPositions[0] - a_point;
  const Vec3 B = m_vertexPositions[1] - a_point;
  const Vec3 C = m_vertexPositions[2] - a_point;

  const Vec3 u = B.cross(C);
  const Vec3 v = C.cross(A);
  const Vec3 w = A.cross(B);

  return !((dot(u, v) < 0.0) || (dot(u, w) < 0.0));
}

Real
Triangle::computeTriangleArea(const Vec3 a_x1, const Vec3 a_x2, const Vec3 a_x3) const noexcept
{
  const Vec3 a = a_x2 - a_x1;
  const Vec3 b = a_x3 - a_x1;

  return 0.5 * (a.cross(b)).length();
}

Real
Triangle::interpolate(const Vec3& a_point) const noexcept
{
  if (m_area <= Real(0)) {
    return Real(0);
  }

  const Real u = computeTriangleArea(m_vertexPositions[1], m_vertexPositions[2], a_point) / m_area;
  const Real v = computeTriangleArea(m_vertexPositions[2], m_vertexPositions[0], a_point) / m_area;
  const Real w = computeTriangleArea(m_vertexPositions[0], m_vertexPositions[1], a_point) / m_area;

  return u * std::get<0>(m_metaData) + v * std::get<1>(m_metaData) + w * std::get<2>(m_metaData);
}

#include <CD_NamespaceFooter.H>
