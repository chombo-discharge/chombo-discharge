/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TriangleCollection.H
  @brief  File containing the declaration of TriangleCollection
  @author Robert Marskar
*/

#ifndef CD_TriangleCollection_H
#define CD_TriangleCollection_H

// Std includes
#include <vector>
#include <memory>

// Our includes
#include <CD_Triangle.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class that represents a set of triangles.
  @details This class uses the bounding volume hierarchies supplied by EBGeometry.
*/
class TriangleCollection
{
public:
  /*!
    @brief Vector alias. Inherited from Triangle class
  */
  using Vec3 = Triangle::Vec3;

  /*!
    @brief Branching factor for BVH tree
  */
  static constexpr int K = 4;

  /*!
    @brief Bounding volume type
  */
  using BV = EBGeometry::BoundingVolumes::AABBT<Real>;

  /*!
    @brief Default constructor
  */
  TriangleCollection() noexcept;

  /*!
    @brief Full constructor. Takes the input triangles and stores them in a BVH-tree.
    @param[in] a_triangles Input triangles. These must already have their meta-data filled.
  */
  TriangleCollection(std::vector<std::shared_ptr<Triangle>>& a_triangles) noexcept;

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
  define(std::vector<std::shared_ptr<Triangle>>& a_triangles) noexcept;

  /*!
    @brief Hierarchical search function that fetches the closest triangles
    @param[in] a_input Point
  */
  virtual std::vector<std::shared_ptr<const Triangle>>
  getClosestTriangles(const Vec3& a_point) const noexcept;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief Triangle collection -- stored in a BVH-tree for faster lookup
  */
  std::shared_ptr<EBGeometry::BVH::LinearBVH<Real, Triangle, BV, K>> m_bvh;
};

#include <CD_NamespaceFooter.H>

#endif
