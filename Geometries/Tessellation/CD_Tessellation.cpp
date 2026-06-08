/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Tessellation.cpp
  @brief  Implementation of CD_Tessellation.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// EBGeometry include
#include <EBGeometry.hpp>

// Our includes
#include <CD_Tessellation.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

/// @cond DOXYGEN_SKIP
using T = float;
/// @endcond

Tessellation::Tessellation()
{

  std::string filename;
  Real        zCoord;
  bool        flipInside;

  // Get input options.
  ParmParse pp("Tessellation");

  pp.get("mesh_file", filename);
  pp.get("z_coord", zCoord);
  pp.get("flip_inside", flipInside);

  // Read the PLY file and put it in a linearized BVH hierarchy.
  auto implicitFunction = EBGeometry::Parser::readIntoLinearBVH<T>(filename);

  // Put our level-set into Chombo datastructures.
  RefCountedPtr<BaseIF> baseIF = RefCountedPtr<BaseIF>(new EBGeometryIF<T>(implicitFunction, flipInside, zCoord));

  m_electrodes.push_back(Electrode(baseIF, true));
}

Tessellation::~Tessellation()
{}

#include <CD_NamespaceFooter.H>
