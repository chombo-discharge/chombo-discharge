/* chombo-discharge
 * Copyright © SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @brief CD_DiskProfiledPlane.cpp
  @brief Implementation of CD_DiskProfiledPlane.H
  @author Robert Marskar
*/

// Std includes
#include <string>
#include <iostream>

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_EBGeometryIF.H>
#include <CD_RoundedCylinderIF.H>
#include <CD_DiskProfiledPlane.H>
#include <CD_NamespaceHeader.H>

DiskProfiledPlane::DiskProfiledPlane() noexcept
{
  CH_TIME("DiskProfiledPlane::DiskProfiledPlane");

  bool useElectrode  = false;
  bool useDielectric = false;

  ParmParse pp("DiskProfiledPlane");

  pp.get("use_electrode", useElectrode);
  pp.get("use_dielectric", useDielectric);

  if (useElectrode) {
    this->defineElectrode();
  }
  else {
    m_electrodes.resize(0);
  }

  if (useDielectric) {
    this->defineDielectric();
  }
  else {
    m_dielectrics.resize(0);
  }
}

void
DiskProfiledPlane::defineElectrode() noexcept
{
  CH_TIME("DiskProfiledPlane::defineElectrode");

  using T       = Real;
  using Vec3    = EBGeometry::Vec3T<T>;
  using ImpFunc = EBGeometry::ImplicitFunction<T>;

  bool isLive         = false;
  Real wheelThickness = 0.0;
  Real wheelRadius    = 0.0;
  Real wheelCurvature = 0.0;

  Vector<Real> v;

  ParmParse pp("DiskProfiledPlane");

  // Get input parameters
  pp.get("wheel_live", isLive);
  pp.get("wheel_thickness", wheelThickness);
  pp.get("wheel_radius", wheelRadius);
  pp.get("wheel_curvature", wheelCurvature);
  pp.getarr("wheel_center", v, 0, SpaceDim);

  // Create the disk electrode using EBGeometry.
  RealVect lo = RealVect(D_DECL(v[0], v[1], v[2]));
  RealVect hi = RealVect(D_DECL(v[0], v[1], v[2]));

  lo[0] -= 0.5 * wheelThickness;
  hi[0] += 0.5 * wheelThickness;

  // std::shared_ptr<ImpFunc> wheel = std::make_shared<EBGeometry::CylinderSDF<T>>(lo, hi, wheelRadius);

  // if (wheelCurvature > 0.0) {
  //   wheel = EBGeometry::Mollify<T>(wheel, wheelCurvature);
  // }

  // Turn EBGeometry into a Chombo function
  auto wheelIF = RefCountedPtr<BaseIF>(new RoundedCylinderIF(lo, hi, wheelRadius, wheelCurvature, false));

  m_electrodes.push_back(Electrode(wheelIF, isLive));
}

void
DiskProfiledPlane::defineDielectric() noexcept
{
  CH_TIME("DiskProfiledPlane::defineDielectric");
}

#include <CD_NamespaceFooter.H>
