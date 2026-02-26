/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Triangle.H
  @brief  Header file for a simple Triangle class.
  @author Robert Marskar
*/

#ifndef CD_Triangle_H
#define CD_Triangle_H

// Std includes
#include <array>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <EBGeometry.hpp>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that represents a single triangle together with associated metadata on the vertices (one float per vertex)
*/
class Triangle : public EBGeometry::Triangle<Real, std::array<Real, 3>>
{
public:
  /*!
    @brief Alias for always-3D vector type
  */
  using Vec3 = EBGeometry::Vec3T<Real>;

  /*!
    @brief Default constructor.
    @details Leaves vertex positions at their EBGeometry defaults and initializes vertex metadata
    to zero. The triangle area is zero-initialized. Vertex positions must be set before calling
    computeArea(), isInside(), or interpolate().
  */
  Triangle() noexcept;

  /*!
    @brief Full constructor.
    @details Sets vertex positions, computes the triangle normal and area, and zero-initializes
    vertex metadata.
    @param[in] a_vertexPositions Vertex positions, in CCW order when viewed from the normal direction.
  */
  Triangle(const std::array<Vec3, 3>& a_vertexPositions) noexcept;

  /*!
    @brief Destructor. 
  */
  virtual ~Triangle() noexcept;

  /*!
    @brief Compute and store the triangle area in m_area.
    @details Uses the cross-product formula: area = 0.5 * ||(v1-v0) x (v2-v0)||. This is called
    automatically by setVertexPositions(), so manual invocation is only needed if vertex positions
    were changed through the base class API directly.
  */
  void
  computeArea() noexcept;

  /*!
    @brief Set vertex positions and recompute the triangle area.
    @details This shadows EBGeometry::Triangle::setVertexPositions to ensure that m_area
    is kept in sync whenever vertex positions change.
    @param[in] a_vertexPositions Vertex positions.
  */
  void
  setVertexPositions(const std::array<Vec3, 3>& a_vertexPositions) noexcept;

  /*!
    @brief Set meta-data on the vertices
    @param[in] a_vertexData Data on the vertices. Must be in the same order as the vertex positions.
  */
  void
  setVertexData(const std::array<Real, 3>& a_vertexData) noexcept;

  /*!
    @brief Check if a point lies inside the triangle.
    @details Tests whether the projection of a_point onto the triangle plane falls within the
    triangle boundary. The test is based on the sign of the triple products of each edge direction
    with the triangle normal, so it is valid for arbitrary 3D points (not just coplanar ones).
    Points exactly on an edge are considered inside.
    @param[in] a_point Query point in 3D space.
    @return True if the projection of a_point onto the triangle plane is inside the triangle,
    false otherwise.
  */
  bool
  isInside(const Vec3& a_point) const noexcept;

  /*!
    @brief Interpolate vertex data at a point using barycentric coordinates.
    @details Computes barycentric weights from area ratios and returns the weighted sum of the
    three vertex data values. The result is exact when a_point lies on the triangle plane; for
    points off the plane the weights will not sum to one and the result is an approximation.
    No bounds check is performed: points outside the triangle produce barycentric coordinates
    outside [0,1] and the function extrapolates. Returns zero for degenerate triangles (m_area == 0).
    @param[in] a_point Query point in 3D space.
    @return Interpolated vertex data value at a_point.
  */
  Real
  interpolate(const Vec3& a_point) const noexcept;

protected:
  /*!
    @brief Triangle area.
  */
  Real m_area = Real(0);

  /*!
    @brief Compute the area of a triangle defined by three vertices.
    @param[in] a_x1 First vertex.
    @param[in] a_x2 Second vertex.
    @param[in] a_x3 Third vertex.
    @return Triangle area: 0.5 * ||(a_x2 - a_x1) x (a_x3 - a_x1)||. Returns zero for degenerate
    (collinear or coincident) vertices.
  */
  inline Real
  computeTriangleArea(const Vec3 a_x1, const Vec3 a_x2, const Vec3 a_x3) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
