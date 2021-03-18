/*!
  @file   porsche.cpp
  @brief  Implementation of porsche.H
  @author Robert Marskar
  @date   Jan. 2019
*/

#include "porsche.H"

#include <ParmParse.H>

#include "dcel_BoundingVolumes.H"
#include "dcel_mesh.H"
#include "dcel_parser.H"
#include "dcel_if.H"
#include "bvh_if.H"

using namespace dcel;

using precision = float;

using face   = faceT<precision>;
using mesh   = meshT<precision>;
using AABB   = AABBT<precision>;
using Sphere = BoundingSphereT<precision>;

using BV = AABB;

porsche::porsche(){

  std::string filename;
  std::string partitioner;

  ParmParse pp("porsche");

  pp.get("mesh_file",   filename);
  pp.get("partitioner", partitioner);

  // Build the dcel_mesh and the BVH
  auto m = std::make_shared<mesh>();
  parser::PLY<precision>::readASCII(*m, filename);
  m->sanityCheck();
  m->reconcile();

  auto root = std::make_shared<NodeT<precision, face, BV> >(m->getFaces());

  if(partitioner == "default"){
    root->topDownSortAndPartitionPrimitives(bvh_if<precision, BV>::defaultStopFunction,
					    bvh_if<precision, BV>::defaultPartitionFunction,
					    bvh_if<precision, BV>::defaultBVConstructor);
  }
  else if(partitioner == "overlap"){
    root->topDownSortAndPartitionPrimitives(bvh_if<precision, BV>::defaultStopFunction,
					    bvh_if<precision, BV>::partitionMinimumOverlap,
					    bvh_if<precision, BV>::defaultBVConstructor);
  }
  else if(partitioner == "sah"){
    root->topDownSortAndPartitionPrimitives(bvh_if<precision, BV>::defaultStopFunction,
					    bvh_if<precision, BV>::partitionSAH,
					    bvh_if<precision, BV>::defaultBVConstructor);
  }
  else{
    MayDay::Abort("porsche::porsche() -- unknown partitioner requested");
  }


  auto bif = RefCountedPtr<bvh_if<precision, BV> > (new bvh_if<precision, BV>(root,false));

  m_electrodes.push_back(electrode(bif, true));
  
}

porsche::~porsche(){
  
}
