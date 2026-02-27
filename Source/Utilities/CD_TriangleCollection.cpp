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
#include <CD_DataParser.H>
#include <CD_TriangleCollection.H>
#include <CD_NamespaceHeader.H>

TriangleCollection::TriangleCollection() noexcept
{
  m_isDefined = false;
}

TriangleCollection::TriangleCollection(const std::vector<std::shared_ptr<Triangle>>& a_triangles) noexcept
{
  this->define(a_triangles);
}

TriangleCollection::TriangleCollection(const std::string& a_filename, const std::string& a_vertexDataIdentifier)
{
  this->define(DataParser::readTriangles(a_filename, a_vertexDataIdentifier));
}

TriangleCollection::~TriangleCollection() noexcept
{}

void
TriangleCollection::define(const std::vector<std::shared_ptr<Triangle>>& a_triangles) noexcept
{
  CH_assert(a_triangles.size() > 0);
  if (a_triangles.empty()) {
    return;
  }

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

std::vector<std::pair<std::shared_ptr<const Triangle>, Real>>
TriangleCollection::getClosestTriangles(const Vec3& a_point) const noexcept
{
  CH_assert(m_isDefined);

  using BVHMeta    = Real;
  using TriAndDist = std::pair<std::shared_ptr<const Triangle>, Real>;
  using Node       = EBGeometry::BVH::LinearNodeT<Real, Triangle, BV, K>;

  BVHMeta shortestDistanceSoFar = std::numeric_limits<Real>::max();

  std::vector<TriAndDist> candidates;

  // Visitation pattern. Go into the node if the point is inside or the distance to the BV is shorter than the shortest distance
  // that we've found so far.
  EBGeometry::BVH::Visiter<Node, Real> visiter = [&shortestDistanceSoFar](const Node&    a_node,
                                                                          const BVHMeta& a_bvDist) noexcept -> bool {
    return a_bvDist <= 0.0 || a_bvDist <= shortestDistanceSoFar;
  };

  // Sorter for BVH nodes. When visiting an internal node, we sort its children based on the distance to the respective bounding volumes. In
  // the BVH traversel, we then visit the closest nodes first.
  EBGeometry::BVH::Sorter<Node, Real, K> sorter =
    [](std::array<std::pair<std::shared_ptr<const Node>, Real>, K>& a_leaves) noexcept -> void {
    std::sort(a_leaves.begin(),
              a_leaves.end(),
              [](const std::pair<std::shared_ptr<const Node>, Real>& n1,
                 const std::pair<std::shared_ptr<const Node>, Real>& n2) -> bool {
                return n1.second > n2.second;
              });
  };

  // Meta-data updater for the BVH nodes. This enters into the visitor pattern where we attach the distance to each
  // BVH node. This is important when we want to check if we should actually go into the node.
  EBGeometry::BVH::MetaUpdater<Node, BVHMeta> metaUpdater = [&a_point](const Node& a_node) noexcept -> BVHMeta {
    return a_node.getDistanceToBoundingVolume(a_point);
  };

  // Update rule for the BVH leaf nodes. This is called at each leaf node, and provides some logic as to what to do
  // when we finally get to the bottom of the tree. Here, we compute the distance to the triangles and if the distance
  // is shorter than the smallest distance we've found so far, we append those triangles to the list of triangles that
  // will be returned to the user.
  EBGeometry::BVH::Updater<Triangle> updater =
    [&shortestDistanceSoFar, &a_point, &candidates](
      const std::vector<std::shared_ptr<const Triangle>>& a_triangles) noexcept -> void {
    for (const auto& f : a_triangles) {
      const Real distToTri = std::abs(f->signedDistance(a_point));

      if (distToTri <= shortestDistanceSoFar) {
        candidates.emplace_back(f, distToTri);
        shortestDistanceSoFar = distToTri;
      }
    }
  };

  // Traverse the BVH tree with the above rules. This builds the candidate list.
  m_bvh->traverse(updater, visiter, sorter, metaUpdater);

  // Sort the candidate triangles based on their distance.
  std::sort(candidates.begin(), candidates.end(), [](const TriAndDist& a, const TriAndDist& b) {
    return a.second < b.second;
  });

  return candidates;
}

#include <CD_NamespaceFooter.H>
