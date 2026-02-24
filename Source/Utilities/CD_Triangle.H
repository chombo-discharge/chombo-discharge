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
    @brief Default constructor
  */
  Triangle() noexcept;

  /*!
    @brief Full constructor
  */
  Triangle(const std::array<Vec3, 3>& a_vertexPositions) noexcept;

  /*!
    @brief Destructor. 
  */
  virtual ~Triangle() noexcept;

  /*!
    @brief Set meta-data on the vertices
    @param[in] a_vertexData Data on the vertices. Must be in the same order as the vertex positions.
  */
  void
  setVertexData(const std::array<Real, 3>& a_vertexData) noexcept;

  /*!
    @brief Check if a point lies on the triangle
    @details Input point
  */
  bool
  isInside(const Vec3& a_point) const noexcept;

  /*!
    @brief Use barycentric interpolation within the triangle.
    @details Input point    
  */
  Real
  interpolate(const Vec3& a_point) const noexcept;

protected:
  /*!
    @brief Triangle area.
  */
  Real m_area;
};

#include <CD_NamespaceFooter.H>

#endif
