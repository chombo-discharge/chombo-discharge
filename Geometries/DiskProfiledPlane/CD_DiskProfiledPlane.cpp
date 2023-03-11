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

using namespace EBGeometry;

using Vec3    = EBGeometry::Vec3T<Real>;
using ImpFunc = EBGeometry::ImplicitFunction<Real>;

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

  bool isLive         = false;
  Real wheelThickness = 0.0;
  Real wheelRadius    = 0.0;
  Real wheelCurvature = 0.0;

  Vector<Real> v(3, 0.0);

  ParmParse pp("DiskProfiledPlane");

  // Get input parameters
  pp.get("wheel_live", isLive);
  pp.get("wheel_thickness", wheelThickness);
  pp.get("wheel_radius", wheelRadius);
  pp.get("wheel_curvature", wheelCurvature);
  pp.getarr("wheel_center", v, 0, SpaceDim);

  // Create disk electrode using EBGeometry. This is an elongation of the union of a torus and a cylinder. Created
  // in the xy-plane and then put into place later.
  auto torus    = std::make_shared<TorusSDF<Real>>(Vec3::zero(), wheelRadius, wheelCurvature);
  auto cylinder = std::make_shared<CylinderSDF<Real>>(Vec3::unit(2) * wheelCurvature,
                                                      -Vec3::unit(2) * wheelCurvature,
                                                      wheelRadius);

  // Make union and rotate and put into place.
  auto wheel = Union<Real>(torus, cylinder);
  if (wheelThickness > 0.0) {
    wheel = Elongate<Real>(wheel, 0.5 * wheelThickness * Vec3::unit(2));
  }
  wheel = Rotate<Real>(wheel, 90.0, 1);
  wheel = Translate<Real>(wheel, Vec3(v[0], v[1], v[2]));

  // Turn EBGeometry into Chombo
  auto wheelIF = RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(wheel, true));

  m_electrodes.push_back(Electrode(wheelIF, isLive));
}

void
DiskProfiledPlane::defineDielectric() noexcept
{
  CH_TIME("DiskProfiledPlane::defineDielectric");

  Vector<Real> boxDimensions(3, std::numeric_limits<Real>::max());
  Vector<Real> boxTranslation(3, 0.0);
  Vector<Real> squareDimensions(3, std::numeric_limits<Real>::max());
  Real         boxCurvature = 0.0;
  Real         permittivity = 1.0;

  // Get input parameters.
  ParmParse pp("DiskProfiledPlane");

  pp.get("box_permittivity", permittivity);
  pp.get("box_curvature", boxCurvature);
  pp.getarr("box_dimensions", boxDimensions, 0, SpaceDim);
  pp.getarr("box_translate", boxTranslation, 0, SpaceDim);
  pp.getarr("square_dimensions", squareDimensions, 0, SpaceDim);

  boxCurvature = std::max(boxCurvature, std::numeric_limits<Real>::min());

  // Basic rounded box.
  std::shared_ptr<ImpFunc> roundedBox = std::make_shared<RoundedBoxSDF<Real>>(Vec3(boxDimensions[0] - 2 * boxCurvature,
                                                                                   boxDimensions[1] - 2 * boxCurvature,
                                                                                   boxDimensions[2] - 2 * boxCurvature),
                                                                              boxCurvature);

  roundedBox = Translate<Real>(roundedBox, -0.5 * boxDimensions[2] * Vec3::unit(2));

  // Basic square profile shape.
  const Vec3               squareDims(squareDimensions[0], squareDimensions[1], squareDimensions[2]);
  std::shared_ptr<ImpFunc> square = std::make_shared<BoxSDF<Real>>(-0.5 * squareDims, 0.5 * squareDims);

  // Subtract the square
  roundedBox = Difference<Real>(roundedBox, square);

  // Translate into place.
  roundedBox = Translate<Real>(roundedBox, Vec3(boxTranslation[0], boxTranslation[1], boxTranslation[2]));

  // Hack for 2D.
  roundedBox = Elongate<Real>(roundedBox, std::numeric_limits<Real>::max() * Vec3::unit(1));

  // Turn EBGeometry into Chombo.
  auto profiledBoxIF = RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(roundedBox, true));

  m_dielectrics.push_back(Dielectric(profiledBoxIF, permittivity));
}

#include <CD_NamespaceFooter.H>
