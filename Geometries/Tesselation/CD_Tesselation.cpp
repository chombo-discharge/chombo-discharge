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

// Our includes
#include <CD_Tesselation.H>
#include <CD_BoundingVolumes.H>
#include <CD_BVH.H>
#include <CD_DcelMesh.H>
#include <CD_DcelParser.H>
#include <CD_DcelSdf.H>
#include <CD_DcelBVH.H>
#include <CD_BvhSdf.H>
#include <CD_NamespaceHeader.H>

using namespace Dcel;

using precision = float;

using Face   = FaceT<precision>;
using Mesh   = MeshT<precision>;
using AABB   = BoundingVolumes::AABBT<precision>;
using Sphere = BoundingVolumes::BoundingSphereT<precision>;

using BV = AABB;

Tesselation::Tesselation(){

  std::string filename;
  std::string partitioner;

  // Binary tree

  Real zCoord;

  ParmParse pp("Tesselation");

  pp.get("mesh_file",   filename);
  pp.get("partitioner", partitioner);
  pp.get("z_coord",     zCoord);

  // Build the dcel_mesh and the BVH
  auto m = std::make_shared<Mesh>();
  Parser::PLY<precision>::readASCII(*m, filename);
  m->reconcile(Dcel::MeshT<precision>::VertexNormalWeight::Angle);

  constexpr int K = 2;  
  
  auto root = std::make_shared<NodeT<precision, Face, BV, K> >(m->getFaces());

  // Build the BVH
  if(partitioner == "default"){
    root->topDownSortAndPartitionPrimitives(defaultStopFunction<precision, BV, K>,
					    spatialSplitPartitioner<precision, K>,
					    defaultBVConstructor<precision, BV>);
  }
  else if(partitioner == "binary"){
    root->topDownSortAndPartitionPrimitives(defaultStopFunction<precision, BV, K>,
					    spatialSplitBinaryPartitioner<precision, K>,
					    defaultBVConstructor<precision, BV>);
  }
  else
    MayDay::Error("Tesselation::Tesselation() -- badx partitioner requested");

  auto bif = RefCountedPtr<BvhSdf<precision, BV, K> > (new BvhSdf<precision, BV, K>(root,false,zCoord));

  m_electrodes.push_back(Electrode(bif, true));
  
}

Tesselation::~Tesselation(){
  
}

#include <CD_NamespaceFooter.H>
