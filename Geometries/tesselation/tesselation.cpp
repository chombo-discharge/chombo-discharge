/*!
  @file   tesselation.cpp
  @brief  Implementation of tesselation.H
  @author Robert Marskar
  @date   March 2021
*/

#include "tesselation.H"

#include <ParmParse.H>

#include <CD_BoundingVolumes.H>
#include <CD_BVH.H>
#include "dcel_mesh.H"
#include "dcel_parser.H"
#include "dcel_if.H"
#include "dcel_BVH.H"
#include "bvh_if.H"

#include "CD_NamespaceHeader.H"

using namespace dcel;

using precision = float;

using face   = faceT<precision>;
using mesh   = meshT<precision>;
using AABB   = BoundingVolumes::AABBT<precision>;
using Sphere = BoundingVolumes::BoundingSphereT<precision>;

using BV = AABB;

tesselation::tesselation(){

  std::string filename;
  std::string partitioner;

  ParmParse pp("tesselation");

  pp.get("mesh_file",   filename);
  pp.get("partitioner", partitioner);

  // Build the dcel_mesh and the BVH
  auto m = std::make_shared<mesh>();
  parser::PLY<precision>::readASCII(*m, filename);
  m->reconcile(VertexNormalWeight::Angle);

  auto root = std::make_shared<NodeT<precision, face, BV> >(m->getFaces());

  if(partitioner == "default"){
    root->topDownSortAndPartitionPrimitives(defaultStopFunction<precision, BV>,
					    defaultPartitionFunction<precision>,
					    defaultBVConstructor<precision, BV>);
  }
  else if(partitioner == "overlap"){
    root->topDownSortAndPartitionPrimitives(defaultStopFunction<precision, BV>,
					    partitionMinimumOverlap<precision, BV>,
					    defaultBVConstructor<precision, BV>);
  }
  else if(partitioner == "sah"){
    root->topDownSortAndPartitionPrimitives(defaultStopFunction<precision, BV>,
					    partitionSAH<precision, BV>,
					    defaultBVConstructor<precision, BV>);
  }
  else{
    MayDay::Abort("tesselation::tesselation() -- unknown partitioner requested");
  }


  auto bif = RefCountedPtr<bvh_if<precision, BV> > (new bvh_if<precision, BV>(root,false));

  m_electrodes.push_back(electrode(bif, true));
  
}

tesselation::~tesselation(){
  
}
#include "CD_NamespaceFooter.H"
