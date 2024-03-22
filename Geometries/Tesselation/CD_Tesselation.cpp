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
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

using T = float;

Tesselation::Tesselation()
{

  std::string filename;
  Real        zCoord;
  bool        flipInside;

  // Get input options.
  ParmParse pp("Tesselation");

  pp.get("mesh_file", filename);
  pp.get("z_coord", zCoord);
  pp.get("flip_inside", flipInside);

  // Read the PLY file and put it in a linearized BVH hierarchy.
  auto implicitFunction = EBGeometry::Parser::readIntoLinearBVH<T>(filename);

  // Put our level-set into Chombo datastructures.
  RefCountedPtr<BaseIF> baseIF = RefCountedPtr<BaseIF>(new EBGeometryIF<T>(implicitFunction, flipInside, zCoord));

  m_electrodes.push_back(Electrode(baseIF, true));
}

Tesselation::~Tesselation()
{}

#include <CD_NamespaceFooter.H>
