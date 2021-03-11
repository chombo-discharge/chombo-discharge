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

  using prec = float;

  std::string filename;
  int tree_depth;
  int max_elements;

  ParmParse pp("porsche");

  pp.get("mesh_file",    filename);
  pp.get("tree_depth",   tree_depth);
  pp.get("max_elements", max_elements);

  // Build the mesh
  std::shared_ptr<dcel::mesh<prec> > m = std::shared_ptr<dcel::mesh<prec>> (new dcel::mesh<prec>());
  dcel::parser::PLY<prec>::readASCII(*m, filename);
  m->sanityCheck();
  m->reconcile(dcel::VertexNormalWeight::Angle);

  // Build tree
  //  mesh->computeBoundingSphere();
  //  mesh->buildKdTree(tree_depth, max_elements);

  // Set algorithms
  //  mesh->setSearchAlgorithm(dcel::SearchAlgorithm::KdTree);
  m->setSearchAlgorithm(dcel::SearchAlgorithm::Direct2);
  m->setInsideOutsideAlgorithm(dcel::InsideOutsideAlgorithm::CrossingNumber);

  // Create the if object
  bool flipNormal = false;
  RefCountedPtr<dcel_if<prec> > bif = RefCountedPtr<dcel_if<prec> >(new dcel_if<prec>(m,true));

  m_electrodes.push_back(electrode(bif, true));
  
}

porsche::~porsche(){
  
}
