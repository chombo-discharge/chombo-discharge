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

using namespace Dcel;

using precision = float;

using Face   = EBGeometry::Dcel::FaceT<precision>;
using Mesh   = EBGeometry::Dcel::MeshT<precision>;
using AABB   = EBGeometry::BoundingVolumes::AABBT<precision>;
using Sphere = EBGeometry::BoundingVolumes::BoundingSphereT<precision>;

using BV = AABB;

Tesselation::Tesselation(){

  std::string filename;
  std::string partitioner;

  // Binary tree

  Real zCoord = 0.0;
  bool useBVH = true;

  ParmParse pp("Tesselation");

  pp.get("mesh_file",   filename);
  pp.get("partitioner", partitioner);
  pp.get("z_coord",     zCoord);
  pp.get("use_bvh",     useBVH);

  // Build the dcel_mesh and the BVH
  auto m = std::make_shared<Mesh>();
  EBGeometry::Dcel::Parser::PLY<precision>::readASCII(*m, filename);
  //  m->reconcile(EBGeometry::Dcel::MeshT<precision>::VertexNormalWeight::Angle);

  constexpr int K = 2;  
  
  auto root = std::make_shared<EBGeometry::BVH::NodeT<precision, Face, BV, K> >(m->getFaces());

  // Build the BVH
  if(partitioner == "default"){
    root->topDownSortAndPartitionPrimitives(EBGeometry::Dcel::defaultStopFunction<precision, BV, K>,
					    EBGeometry::Dcel::spatialSplitPartitioner<precision, K>,
					    EBGeometry::Dcel::defaultBVConstructor<precision, BV>);
  }
  else if(partitioner == "binary"){
    root->topDownSortAndPartitionPrimitives(EBGeometry::Dcel::defaultStopFunction<precision, BV, K>,
					    EBGeometry::Dcel::spatialSplitBinaryPartitioner<precision, K>,
					    EBGeometry::Dcel::defaultBVConstructor<precision, BV>);
  }
  else
    MayDay::Error("Tesselation::Tesselation() -- badx partitioner requested");

  auto bif = RefCountedPtr<SignedDistanceBVH<precision, BV, K> > (new SignedDistanceBVH<precision, BV, K>(root,false,zCoord));

  m_electrodes.push_back(Electrode(bif, true));
  
}

Tesselation::~Tesselation(){
  
}

#include <CD_NamespaceFooter.H>
