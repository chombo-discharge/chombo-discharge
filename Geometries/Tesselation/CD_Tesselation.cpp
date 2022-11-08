/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Tesselation.cpp
  @brief  Implementation of CD_Tesselation.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// EBGeometry include
#include <EBGeometry.hpp>

// Our includes
#include <CD_Tesselation.H>
#include <CD_SignedDistanceBVH.H>
#include <CD_SignedDistanceDCEL.H>
#include <CD_NamespaceHeader.H>

using precision = float;

using Face   = EBGeometry::DCEL::FaceT<precision>;
using Mesh   = EBGeometry::DCEL::MeshT<precision>;
using AABB   = EBGeometry::BoundingVolumes::AABBT<precision>;
using Sphere = EBGeometry::BoundingVolumes::BoundingSphereT<precision>;

using BV = AABB;

Tesselation::Tesselation()
{

  std::string filename;
  std::string partitioner;

  // Binary tree

  Real zCoord = 0.0;

  ParmParse pp("Tesselation");

  pp.get("mesh_file", filename);
  pp.get("z_coord", zCoord);

  // Read the PLY file and put it in a DCEL mesh.
  auto m = EBGeometry::Parser::readIntoDCEL<precision>(filename);

  // Build the regular BVH tree. Use a quad-tree.
  constexpr int K    = 4;
  auto          root = std::make_shared<EBGeometry::BVH::NodeT<precision, Face, BV, K>>(m->getFaces());
  root->topDownSortAndPartitionPrimitives(EBGeometry::DCEL::defaultBVConstructor<precision, BV>,
                                          EBGeometry::DCEL::defaultPartitioner<precision, BV, K>,
                                          EBGeometry::DCEL::defaultStopFunction<precision, BV, K>);

  // Flatten the tree.
  auto linearNode = root->flattenTree();

  // Create our electrode.
  auto bif = RefCountedPtr<SignedDistanceBVH<precision, BV, K>>(
    new SignedDistanceBVH<precision, BV, K>(linearNode, false, zCoord));

  m_electrodes.push_back(Electrode(bif, true));
}

Tesselation::~Tesselation() {}

#include <CD_NamespaceFooter.H>
