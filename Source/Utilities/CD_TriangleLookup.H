/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TriangleLookup.H
  @brief  Headers for required functionality for triangle lookup
  @author Robert Marskar
*/

#ifndef CD_TriangleLookup_H
#define CD_TriangleLookup_H

// Our includes
#include <EBGeometry.hpp>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that represents a single triangle together with associated metadata on the vertices.
*/
class Tri : public EBGeometry::Triangle<Real, std::array<Real, 3>>
{
public:
  /*!
    @brief Check if a point lies on the triangle
  */
  bool
  isInside(const Real x, const Real y, const Real z) const noexcept;

  /*!
    @brief Use barycentric interpolation within the triangle. 
  */
  Real
  interpolate(const Real x, const Real y, const Real z) const noexcept;

protected:
  /*!
    @brief Triangle area
  */
  Real m_area;
};

/*!
  @brief Class that represents a set of triangles.
  @details This class uses the bounding volume hierarchies supplied by EBGeometry.
*/
class TriangleCollection
{
public:
  /*!
    @brief Default constructor
  */
  TriangleCollection() noexcept;

  /*!
    @brief Full constructor. Takes the input triangles and stores them in a BVH-tree.
    @param[in] a_triangles Input triangles. These must already have their meta-data filled.
  */
  TriangleCollection(std::vector<std::shared_ptr<Tri>>& a_triangles) noexcept;

  /*!
    @brief Full constructor. Read the input file and stores the triangle in a BVH-tree.
    @note The input file must either be a VTK file or a PLY file. 
    @param[in] a_file File name
    @param[in] a_vertexData Identifier for the vertex data (must be a single property/float per vertex9
  */
  TriangleCollection(const std::string a_filename, const std::string a_vertexData);

  /*!
    @brief Destructor (does nothing). 
  */
  virtual ~TriangleCollection() noexcept;

  /*!
    @brief Full define function. Takes the input triangles and stores them in a BVH-tree.
    @param[in] a_triangles Input triangles. These must already have their meta-data filled.
  */
  virtual void
  define(std::vector<std::shared_ptr<Tri>>& a_triangles) noexcept;

protected:
  /*!
    @brief Triangle collection -- stored in a BVH-tree for faster lookup
  */
  std::shared_ptr<EBGeometry::BVH::LinearBVH<Real, Tri, EBGeometry::BoundingVolumes::AABBT<Real>, 4>> m_bvh;
};

#if 0
  /*!
    @brief Read triangles from a VTK file or a PLY file, and put them in a bounding volume hierarchy.
    @param[in] a_filename File name for the 
  */
  std::shared_ptr<LinearBVH<Real, Tri, EBGeometry::AABBT<Real>, 4>>
  readTriangles(const std::string a_filename, const std::string a_vertexData);
#endif

#include <CD_NamespaceFooter.H>

#endif
