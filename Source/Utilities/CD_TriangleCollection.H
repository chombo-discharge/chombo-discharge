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
    @brief 3D vector type, re-exported from Triangle for convenience.
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
    @brief Default constructor.
    @details Leaves the object in an undefined state (m_isDefined = false). define() must be called
    before getClosestTriangles() can be used.
  */
  TriangleCollection() noexcept;

  /*!
    @brief Full constructor. Takes the input triangles and stores them in a BVH-tree.
    @details Calls define(). The object is ready for use immediately after construction.
    @param[in] a_triangles Input triangles. Must be non-empty and must already have their vertex
    positions and meta-data set.
  */
  TriangleCollection(const std::vector<std::shared_ptr<Triangle>>& a_triangles) noexcept;

  /*!
    @brief Destructor (does nothing). 
  */
  virtual ~TriangleCollection() noexcept;

  /*!
    @brief Build the BVH tree from a set of triangles.
    @details Packs each triangle into an axis-aligned bounding box, constructs a K-ary BVH using
    top-down partitioning, and flattens it into a linear layout for cache-friendly traversal. Sets
    m_isDefined to true on success. This is the same operation performed by the full constructor;
    use it to (re-)initialize an object that was default-constructed or needs to be rebuilt.
    @param[in] a_triangles Input triangles. Must be non-empty and must already have their vertex
    positions and meta-data set.
  */
  virtual void
  define(const std::vector<std::shared_ptr<Triangle>>& a_triangles) noexcept;

  /*!
    @brief Hierarchical BVH search that returns candidate closest triangles to a query point.
    @details Traverses the BVH closest-node-first, pruning any node whose bounding volume distance
    exceeds the shortest triangle distance found so far. At each leaf the pruning threshold is
    tightened as each closer triangle is encountered, so the returned list contains every triangle
    that was the shortest-distance-so-far at the moment it was evaluated. As a result the list may
    include triangles that are not among the globally closest once the full traversal is complete;
    callers should use the first entry (index 0) as the definitive closest triangle.
    Requires define() to have been called first.
    @param[in] a_point Query point in 3D space.
    @return Candidate triangles paired with their unsigned distance to a_point, sorted ascending by
    distance. The first entry is the closest triangle found.
  */
  virtual std::vector<std::pair<std::shared_ptr<const Triangle>, Real>>
  getClosestTriangles(const Vec3& a_point) const noexcept;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined = false;

  /*!
    @brief Triangle collection -- stored in a BVH-tree for faster lookup
  */
  std::shared_ptr<EBGeometry::BVH::LinearBVH<Real, Triangle, BV, K>> m_bvh;
};

#include <CD_NamespaceFooter.H>

#endif
