/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MechanicalShaft.cpp
  @brief  Implementation of CD_MechanicalShaft.H
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
#include <CD_MechanicalShaft.H>
#include <CD_PerlinSphereSdf.H>
#include <CD_BoxSdf.H>
#include <CD_CylinderSdf.H>
#include <CD_ProfileCylinderIF.H>
#include <CD_PolygonRodIF.H>
#include <CD_HollowCylinderIF.H>
#include <CD_EBGeometryIF.H>
#include <CD_RoundedCylinderIF.H>
#include <CD_NamespaceHeader.H>

using Vec3    = EBGeometry::Vec3T<Real>;
using ImpFunc = EBGeometry::ImplicitFunction<Real>;

MechanicalShaft::MechanicalShaft() noexcept
{
  CH_TIME("MechanicalShaft::MechanicalShaft");

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  ParmParse pp("MechanicalShaft");

  Real eps0          = 1.0;
  bool useElectrode  = false;
  bool useDielectric = false;

  pp.get("eps0", eps0);
  pp.get("use_electrode", useElectrode);
  pp.get("use_dielectric", useDielectric);

  // Define shit.
  if (useElectrode) {
    this->defineElectrode();
  }
  if (useDielectric) {
    this->defineDielectric();
  }

  this->setGasPermittivity(eps0);
}

MechanicalShaft::~MechanicalShaft() noexcept { CH_TIME("MechanicalShaft::~MechanicalShaft"); }

void
MechanicalShaft::defineElectrode() noexcept
{
  CH_TIME("MechanicalShaft::defineElectrode");

  ParmParse pp("MechanicalShaft");

  Vector<Real> v;
  std::string  orientation = "+z";
  Vec3         translate   = Vec3::zero();

  bool live        = false;
  Real length      = -1.0;
  Real innerRadius = -1.0;
  Real outerRadius = -1.0;
  Real innerCurv   = -1.0;
  Real outerCurv   = -1.0;

  pp.get("electrode.live", live);
  pp.get("electrode.inner.radius", innerRadius);
  pp.get("electrode.outer.radius", outerRadius);
  pp.get("electrode.inner.curvature", innerCurv);
  pp.get("electrode.outer.curvature", outerCurv);
  pp.get("electrode.length", length);
  pp.get("electrode.orientation", orientation);
  pp.getarr("electrode.translate", v, 0, SpaceDim);

  for (int dir = 0; dir < SpaceDim; dir++) {
    translate[dir] = v[dir];
  }

  CH_assert(length > 0.0);
  CH_assert(innerRadius > 0.0);
  CH_assert(outerRadius > 0.0);
  CH_assert(innerCurv >= 0.0);
  CH_assert(outerCurv >= 0.0);
  CH_assert(outerRadius > innerRadius);

  std::shared_ptr<ImpFunc> outerCylinder;
  std::shared_ptr<ImpFunc> innerCylinder;
  std::shared_ptr<ImpFunc> torus;
  std::shared_ptr<ImpFunc> hollowCylinder;

  const Vec3 zHi(0.0, 0.0, 0.5 * outerCurv);
  const Vec3 zLo(0.0, 0.0, -0.5 * outerCurv);

  outerCylinder  = std::make_shared<EBGeometry::CylinderSDF<Real>>(zLo, zHi, outerRadius);
  innerCylinder  = std::make_shared<EBGeometry::CylinderSDF<Real>>(zLo, zHi, innerRadius);
  torus          = std::make_shared<EBGeometry::TorusSDF<Real>>(Vec3::zero(), outerRadius, 0.5 * outerCurv);
  hollowCylinder = EBGeometry::Union<Real>(outerCylinder, torus);
  hollowCylinder = EBGeometry::SmoothDifference<Real>(hollowCylinder, innerCylinder, innerCurv);

  if (length - outerCurv > 0.0) {
    hollowCylinder = EBGeometry::Elongate<Real>(hollowCylinder, Vec3(0.0, 0.0, 0.5 * (length - outerCurv)));
  }

  // Rotate into place and translate
  if (orientation == "+x") {
    hollowCylinder = EBGeometry::Rotate<Real>(hollowCylinder, 90.0, 1);
  }
  else if (orientation == "-x") {
    hollowCylinder = EBGeometry::Rotate<Real>(hollowCylinder, -90.0, 1);
  }
  else if (orientation == "+y") {
    hollowCylinder = EBGeometry::Rotate<Real>(hollowCylinder, -90.0, 0);
  }
  else if (orientation == "-y") {
    hollowCylinder = EBGeometry::Rotate<Real>(hollowCylinder, 90.0, 0);
  }
  else if (orientation == "+z") {
    hollowCylinder = EBGeometry::Rotate<Real>(hollowCylinder, 180.0, 0);
  }
  else if (orientation == "-z") {
  }
  else {
    MayDay::Error("MechanicalShaft::defineElectrode - logic bust in orientation. Must be -x, +x, -y, ....");
  }

  // Translate to user-specified position.
  hollowCylinder = EBGeometry::Translate<Real>(hollowCylinder, translate);

  // Stretch along z-axis if this is a 2D geometry.
  if (SpaceDim == 2) {
    hollowCylinder = EBGeometry::Elongate<Real>(hollowCylinder, std::numeric_limits<Real>::max() * Vec3::unit(2));
  }

  // Turn into Chombo-speak.
  RefCountedPtr<BaseIF> elec = RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(hollowCylinder, true));

  m_electrodes.resize(1);
  m_electrodes[0].define(elec, live);
}

void
MechanicalShaft::defineDielectric() noexcept
{
  CH_TIME("MechanicalShaft::defineDielectric");

  ParmParse pp("MechanicalShaft");

  std::string shape       = "";
  std::string orientation = "+z";

  Real eps       = -1.0;
  Vec3 translate = Vec3::zero();

  Vector<Real> v;

  pp.get("dielectric.shape", shape);
  pp.get("dielectric.permittivity", eps);
  pp.get("dielectric.orientation", orientation);
  pp.getarr("dielectric.translate", v, 0, SpaceDim);

  for (int dir = 0; dir < SpaceDim; dir++) {
    translate[dir] = v[dir];
  }

  std::shared_ptr<ImpFunc> shaft;

  if (shape == "cylinder") {
    shaft = this->getSimpleCylinder();
  }
  else if (shape == "polygon") {
    shaft = this->getPolygon();
  }
  else if (shape == "circular_profiles") {
    shaft = this->getCircularProfiles();
  }
  else {
    MayDay::Abort("MechanicalShaft::defineDielectric - unknown argument 'dielectric.shape'");
  }

  CH_assert(eps >= 1.0);
  CH_assert(shaft != nullptr);

  // Rotate into place and translate
  if (orientation == "+x") {
    shaft = EBGeometry::Rotate<Real>(shaft, 90.0, 1);
  }
  else if (orientation == "-x") {
    shaft = EBGeometry::Rotate<Real>(shaft, -90.0, 1);
  }
  else if (orientation == "+y") {
    shaft = EBGeometry::Rotate<Real>(shaft, -90.0, 0);
  }
  else if (orientation == "-y") {
    shaft = EBGeometry::Rotate<Real>(shaft, 90.0, 0);
  }
  else if (orientation == "+z") {
    shaft = EBGeometry::Rotate<Real>(shaft, 180.0, 0);
  }
  else if (orientation == "-z") {
  }
  else {
    MayDay::Error("MechanicalShaft::defineElectrode - logic bust in orientation. Must be -x, +x, -y, ....");
  }

  // Translate to user-specified position.
  shaft = EBGeometry::Translate<Real>(shaft, translate);

  // Stretch along z-axis if this is a 2D geometry.
  if (SpaceDim == 2) {
    shaft = EBGeometry::Elongate<Real>(shaft, std::numeric_limits<Real>::max() * Vec3::unit(2));
  }

  RefCountedPtr<BaseIF> diel = RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(shaft, true));

  m_dielectrics.resize(1);
  m_dielectrics[0].define(diel, eps);
}

std::shared_ptr<ImpFunc>
MechanicalShaft::getSimpleCylinder() const noexcept
{
  CH_TIME("MechanicalShaft::getSimpleCylinder");

  ParmParse pp("MechanicalShaft");

  Real radius = -1.0;

  pp.get("dielectric.cylinder.radius", radius);

  CH_assert(radius > 0.0);

  return std::make_shared<EBGeometry::InfiniteCylinderSDF<Real>>(Vec3::zero(), radius, 2);
}

std::shared_ptr<ImpFunc>
MechanicalShaft::getPolygon() const noexcept
{
  ParmParse pp("MechanicalShaft");

  int  numSides = -1;
  Real radius   = -1.0;
  Real curv     = -1.0;

  pp.get("dielectric.polygon.num_sides", numSides);
  pp.get("dielectric.polygon.radius", radius);
  pp.get("dielectric.polygon.curvature", curv);

  CH_assert(numSides > 2);
  CH_assert(radius > 0.0);
  CH_assert(curv > 0.0);
  CH_assert(curv < radius);

  const Real dTheta = 2.0 * M_PI / numSides;
  const Real alpha  = M_PI / numSides;
  const Real beta   = 0.5 * M_PI - alpha;
  const Real b      = sin(beta);
  const Real r      = curv;
  const Real R      = radius - r + r / b;

  std::vector<std::shared_ptr<ImpFunc>> planes;
  for (int iside = 0; iside < numSides; iside++) {
    const Real theta  = iside * dTheta;
    const Vec3 normal = Vec3(cos(theta + 0.5 * dTheta), sin(theta + 0.5 * dTheta), 0.0);
    const Vec3 point  = R * Vec3(cos(theta), sin(theta), 0.0);

    planes.emplace_back(std::make_shared<EBGeometry::PlaneSDF<Real>>(point, -normal));
  }

  return EBGeometry::SmoothUnion(planes, curv);
}

std::shared_ptr<EBGeometry::ImplicitFunction<Real>>
MechanicalShaft::getCircularProfiles() const noexcept
{
  CH_TIME("MechanicalShaft::getCircularProfiles");

  ParmParse pp("MechanicalShaft");

  Real cylinderRadius     = -1.0;
  Real profileMajorRadius = -1.0;
  Real profileMinorRadius = -1.0;
  Real profileTranslate   = -1.0;
  Real profilePeriod      = -1.0;
  Real profileSmooth      = -1.0;

  int profileRepLo = -1;
  int profileRepHi = -1;

  pp.get("dielectric.profile.circular.cylinder_radius", cylinderRadius);
  pp.get("dielectric.profile.circular.profile_major_radius", profileMajorRadius);
  pp.get("dielectric.profile.circular.profile_minor_radius", profileMinorRadius);
  pp.get("dielectric.profile.circular.profile_translate", profileTranslate);
  pp.get("dielectric.profile.circular.profile_period", profilePeriod);
  pp.get("dielectric.profile.circular.profile_repeat_lo", profileRepLo);
  pp.get("dielectric.profile.circular.profile_repeat_hi", profileRepHi);
  pp.get("dielectric.profile.circular.profile_smooth", profileSmooth);

  CH_assert(cylinderRadius > 0.0);
  CH_assert(profileMajorRadius > 0.0);
  CH_assert(profileMinorRadius > 0.0);
  CH_assert(profilePeriod > 0.0);
  CH_assert(profileRepLo > 0.0);
  CH_assert(profileRepHi > 0.0);
  CH_assert(profileSmooth >= 0.0);

  std::shared_ptr<ImpFunc> cylinder;
  std::shared_ptr<ImpFunc> profile;

  const Vec3 center = Vec3::zero();
  const Vec3 transl = profileTranslate * Vec3::unit(2);
  const Vec3 period = profilePeriod * Vec3::unit(2);
  const Vec3 repLo  = profileRepLo * Vec3::unit(2);
  const Vec3 repHi  = profileRepHi * Vec3::unit(2);

  cylinder = std::make_shared<EBGeometry::InfiniteCylinderSDF<Real>>(center, cylinderRadius, 2);
  profile  = std::make_shared<EBGeometry::TorusSDF<Real>>(center, profileMajorRadius, profileMinorRadius);
  profile  = EBGeometry::Translate<Real>(profile, transl);
  profile  = EBGeometry::FiniteRepetition<Real>(profile, period, repLo, repHi);

  return EBGeometry::SmoothDifference(cylinder, profile, profileSmooth);
}

#include <CD_NamespaceFooter.H>
