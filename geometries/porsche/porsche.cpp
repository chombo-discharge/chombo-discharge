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

// Grouping in space. In principle we could use Morton ordering, even. 
template <class T>
auto comparator = [](const dcel::faceT<T>& f1, const dcel::faceT<T>& f2, const int a_dir){
  return f1.getCentroid()[a_dir] < f2.getCentroid()[a_dir];
};

using namespace dcel;

porsche::porsche(){

  using prec = float;
  using mesh = meshT<prec>;
  using AABB = AABBT<prec>;

  std::string filename;
  int tree_depth;
  int max_elements;

  ParmParse pp("porsche");

  pp.get("mesh_file",    filename);
  pp.get("tree_depth",   tree_depth);
  pp.get("max_elements", max_elements);

  // Build the mesh
  std::shared_ptr<mesh> m = std::shared_ptr<mesh> (new mesh());
  parser::PLY<prec>::readASCII(*m, filename);
  m->sanityCheck();
  m->reconcile(VertexNormalWeight::Angle);


  // Set algorithms
  //  mesh->setSearchAlgorithm(dcel::SearchAlgorithm::KdTree);
  m->setSearchAlgorithm(SearchAlgorithm::Direct2);
  m->setInsideOutsideAlgorithm(InsideOutsideAlgorithm::CrossingNumber);


  m->buildBVH(comparator<prec>);

  // Create the if object
  bool flipNormal = false;

  RefCountedPtr<dcel_if<prec, AABB> > bif = RefCountedPtr<dcel_if<prec, AABB> > (new dcel_if<prec, AABB>(m,true));

  m_electrodes.push_back(electrode(bif, true));
  
}

porsche::~porsche(){
  
}
