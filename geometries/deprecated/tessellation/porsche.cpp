/*!
  @file   tessellation.cpp
  @brief  Implementation of tessellation.H
  @author Robert Marskar
  @date   Jan. 2019
*/

#include "tessellation.H"

#include <string>
#include <iostream>
#include <fstream>

#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

#include "ply_reader.H"
#include "dcel_if.H"


tessellation::tessellation(){
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  std::string filename;
  int tree_depth;
  int max_elements;

  ParmParse pp("tessellation");

  pp.get("mesh_file", filename);
  pp.get("tree_depth", tree_depth);
  pp.get("max_elements", max_elements);

  // Build the mesh
  RefCountedPtr<dcel_mesh> mesh = RefCountedPtr<dcel_mesh> (new dcel_mesh());
  ply_reader::read_ascii(*mesh, filename);
  mesh->reconcile_polygons(true, false);
  mesh->build_tree(tree_depth, max_elements);

  // Create the if object
  RefCountedPtr<dcel_if> object = RefCountedPtr<dcel_if>(new dcel_if(mesh, true));

  m_electrodes.resize(1);
  m_electrodes[0].define(object, true);
  
}

tessellation::~tessellation(){
  
}
