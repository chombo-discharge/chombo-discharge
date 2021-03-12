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

using namespace dcel;

using prec = double;
using face = faceT<prec>;
using mesh = meshT<prec>;
using AABB = AABBT<prec>;


#if 0
// In case you want to pass in a different lambda for sorting
SortingFunctor<prec> sortFunc = [](const face& f1, const face& f2, const int a_dir){
  return f1.getCentroid()[a_dir] > f2.getCentroid()[a_dir];
};

// In case you want to pass in a different cost function. 
CostFunctor<prec> costFunc = [](const std::vector<std::shared_ptr<face> >& sortedFaces, const int a_dir) {
  for (const auto& f : sortedFaces){
    std::cout << f->getCentroid()[0] << std::endl;
  }
    
  return std::pair<double, int>(0., 0);
};
#endif

porsche::porsche(){

  std::string filename;

  ParmParse pp("porsche");

  pp.get("mesh_file", filename);

  // Build the mesh and set default parameters. 
  auto m = std::make_shared<mesh>();

  // Read mesh from file and reconcile. 
  parser::PLY<prec>::readASCII(*m, filename);
  m->sanityCheck();
  m->reconcile(VertexNormalWeight::Angle);
  m->setSearchAlgorithm(SearchAlgorithm::Direct2);
  m->setInsideOutsideAlgorithm(InsideOutsideAlgorithm::CrossingNumber);


  // Creat the object and build the BVH. 
  RefCountedPtr<dcel_if<prec, AABB> > bif = RefCountedPtr<dcel_if<prec, AABB> > (new dcel_if<prec, AABB>(m,true));

  bif->buildBVH();

  m_electrodes.push_back(electrode(bif, true));
  
}

porsche::~porsche(){
  
}
