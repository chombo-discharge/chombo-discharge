/*!
  @file   porsche.cpp
  @brief  Implementation of porsche.H
  @author Robert Marskar
  @date   Jan. 2019
*/

#include "porsche.H"

#include <ParmParse.H>

#include "dcel_mesh.H"
#include "dcel_parser.H"
#include "dcel_if.H"


porsche::porsche(){

  std::string filename;
  int tree_depth;
  int max_elements;

  ParmParse pp("porsche");

  pp.get("mesh_file",    filename);
  pp.get("tree_depth",   tree_depth);
  pp.get("max_elements", max_elements);

  // Build the mesh
  std::shared_ptr<dcel::mesh> mesh = std::shared_ptr<dcel::mesh> (new dcel::mesh());
  dcel::parser::PLY::readASCII(*mesh, filename);
  mesh->sanityCheck();
  mesh->reconcile(dcel::VertexNormalWeight::Angle);

  // Build tree
  mesh->computeBoundingSphere();
  mesh->buildKdTree(tree_depth, max_elements);

  // Set algorithms
  //mesh->setAlgorithm(dcel::mesh::SearchAlgorithm::KdTree);
  mesh->setSearchAlgorithm(dcel::SearchAlgorithm::Direct);
  mesh->setInsideOutsideAlgorithm(dcel::InsideOutsideAlgorithm::CrossingNumber);

  // Create the if object
  bool flipNormal = false;
  RefCountedPtr<dcel_if> bif = RefCountedPtr<dcel_if>(new dcel_if(mesh, false));

  m_electrodes.push_back(electrode(bif, true));
  
}

porsche::~porsche(){
  
}
