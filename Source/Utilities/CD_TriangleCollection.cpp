/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TriangleCollection.cpp
  @brief  Implementation of CD_TriangleCollection.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_TriangleCollection.H>
#include <CD_NamespaceHeader.H>

TriangleCollection::TriangleCollection() noexcept
{
  m_isDefined = false;
}

TriangleCollection::TriangleCollection(std::vector<std::shared_ptr<Triangle>>& a_triangles) noexcept
{
  this->define(a_triangles);
}

TriangleCollection::~TriangleCollection() noexcept
{}

void
TriangleCollection::define(std::vector<std::shared_ptr<Triangle>>& a_triangles) noexcept
{
  CH_assert(a_triangles.size() > 0);

  // Pack triangles and their bounding volumes into a list.
  using TriAndBVList = std::vector<std::pair<std::shared_ptr<const Triangle>, BV>>;

  TriAndBVList triList;

  for (const auto& tri : a_triangles) {
    const auto& v = tri->getVertexPositions();

    std::vector<Vec3> vertices{v[0], v[1], v[2]};

    triList.emplace_back(std::make_pair(tri, BV(vertices)));
  }

  // Create the root BVH node.
  auto bvh = std::make_shared<EBGeometry::BVH::NodeT<Real, Triangle, BV, K>>(triList);

  // Sort the BVH tree.
  bvh->topDownSortAndPartition();

  // Flatten the tree onto a linear representation.
  m_bvh = bvh->flattenTree();

  m_isDefined = true;
}

std::vector<std::shared_ptr<const Triangle>>
TriangleCollection::getClosestTriangles(const Vec3& a_point) const noexcept
{
  CH_assert(m_isDefined);

  using BVHMeta    = Real;
  using TriAndDist = std::pair<std::shared_ptr<const Triangle>, Real>;
  using Node       = EBGeometry::BVH::LinearNodeT<Real, Triangle, BV, K>;

  std::vector<TriAndDist> candidates;

  BVHMeta shortestDistanceSoFar = std::numeric_limits<Real>::max();

  // Visitation pattern. Go into the node if the point is inside or the distance to the BV is shorter than the shortest distance
  // that we've found so far.
  EBGeometry::BVH::Visiter<Node, Real> visiter = [&shortestDistanceSoFar](const Node&    a_node,
                                                                          const BVHMeta& a_bvDist) noexcept -> bool {
    return a_bvDist <= 0.0 || a_bvDist <= shortestDistanceSoFar;
  };

  // Sort the candidate triangles based on their distance.
  std::sort(candidates.begin(), candidates.end(), [](const TriAndDist& a, const TriAndDist& b) {
    return a.second < b.second;
  });

  // Return the triangles without the distance attached.
  std::vector<std::shared_ptr<const Triangle>> ret;
  for (const auto& candidate : candidates) {
    ret.emplace_back(candidate.first);
  }

  return ret;
}

#include <CD_NamespaceFooter.H>
