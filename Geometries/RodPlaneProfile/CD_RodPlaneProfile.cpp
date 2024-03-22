/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RodPlaneProfile.cpp
  @brief  Implementation of CD_RodPlaneProfile.H
  @author Robert Marskar
*/

// Std includes
#include <string>
#include <iostream>
#include <fstream>

// Chombo includes
#include <ParmParse.H>
#include <BaseIF.H>
#include <SphereIF.H>

// Our includes
#include <CD_RodPlaneProfile.H>
#include <CD_ProfilePlaneIF.H>
#include <CD_RodIF.H>
#include <CD_SphereSdf.H>
#include <CD_BoxSdf.H>
#include <CD_RoundedBoxIF.H>
#include <CD_NamespaceHeader.H>

RodPlaneProfile::RodPlaneProfile()
{
  if (SpaceDim == 3)
    MayDay::Abort("RodPlaneProfile::RodPlaneProfile - this is currently for 2D only");
  std::string  str;
  Vector<Real> vec(SpaceDim);

  bool        rod_live, has_rod, has_plane;
  Real        eps, rod_rad, xshift, curv, dist, yshift, width;
  int         numl, numr;
  RealVect    center1, center2, point, normal;
  std::string profile;

  ParmParse pp("RodPlaneProfile");

  pp.get("turn_on_rod", has_rod);
  pp.get("turn_on_plane", has_plane);
  pp.get("rod_live", rod_live);
  pp.get("rod_radius", rod_rad);
  pp.get("profile_num_left", numl);
  pp.get("profile_num_right", numr);
  pp.get("profile_curv", curv);
  pp.get("profile_dist", dist);
  pp.get("profile_xshift", xshift);
  pp.get("profile_yshift", yshift);
  pp.get("plane_width", width);
  pp.get("plane_eps", eps);
  pp.get("profile", str);
  pp.getarr("rod_center1", vec, 0, SpaceDim);
  center1 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("rod_center2", vec, 0, SpaceDim);
  center2 = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("plane_point", vec, 0, SpaceDim);
  point = RealVect(D_DECL(vec[0], vec[1], vec[2]));
  pp.getarr("plane_normal", vec, 0, SpaceDim);
  normal = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  // Set the shape.
  if (str == "circle") {
    m_profile = profile::circle;
  }
  else if (str == "square") {
    m_profile = profile::square;
  }
  else {
    MayDay::Abort("RodPlaneProfile::RodPlaneProfile - unknown profile requested");
  }

  // Build geometries
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  if (has_rod) {
    m_electrodes.resize(1);
    RefCountedPtr<BaseIF> rod = RefCountedPtr<BaseIF>(new RodIF(center1, center2, rod_rad, false));
    m_electrodes[0].define(rod, rod_live);
  }
  if (has_plane) {
    m_dielectrics.resize(1);
    BaseIF*               func  = this->getBaseIF();
    RefCountedPtr<BaseIF> plane = RefCountedPtr<BaseIF>(
      new ProfilePlaneIF(point, width, func, numl, numr, dist, xshift, yshift, curv, false));
    m_dielectrics[0].define(plane, eps);
  }

  this->setGasPermittivity(1.0);
}

RodPlaneProfile::~RodPlaneProfile()
{}

BaseIF*
RodPlaneProfile::getBaseIF()
{
  BaseIF* ret = nullptr;

  switch (m_profile) {
  case profile::circle:
    ret = this->getBaseIFCircle();
    break;
  case profile::square:
    ret = this->getBaseIFSquare();
    break;
  default:
    MayDay::Abort("RodPlaneProfile::getBaseIF - logic bust, unknown profile requested");
    break;
  }

  return ret;
}

BaseIF*
RodPlaneProfile::getBaseIFCircle()
{
  ParmParse pp("RodPlaneProfile");

  Vector<Real> vec(SpaceDim);
  RealVect     point;
  Real         rad;

  pp.get("circle_radius", rad);
  pp.getarr("plane_point", vec, 0, SpaceDim);

  point = RealVect(D_DECL(vec[0], vec[1], vec[2]));

  return new SphereSdf(point, rad, true);
}

BaseIF*
RodPlaneProfile::getBaseIFSquare()
{
  ParmParse pp("RodPlaneProfile");

  Real         width;
  Real         depth;
  Real         curv;
  Vector<Real> vec(SpaceDim);
  RealVect     point;

  pp.get("profile_curv", curv);
  pp.get("square_width", width);
  pp.get("square_depth", depth);
  pp.getarr("plane_point", vec, 0, SpaceDim);
  point = RealVect(D_DECL(vec[0], vec[1], vec[2]));

#if CH_SPACEDIM == 2
  const RealVect lo = point - 0.5 * width * BASISREALV(0) - depth * BASISREALV(1);
  const RealVect hi = point + 0.5 * width * BASISREALV(0) + (depth + curv) * BASISREALV(1);
#elif CH_SPACEDIM == 3
  const RealVect lo = point - 0.5 * width * (BASISREALV(0) + BASISREALV(1)) - depth * BASISREALV(2);
  const RealVect hi = point + 0.5 * width * (BASISREALV(0) + BASISREALV(1));
#endif

  return new RoundedBoxIF(lo, hi, curv, true);
}
#include "CD_NamespaceFooter.H"
