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

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  if (useElectrode) {
    this->defineElectrode();
  }
  if (useDielectric) {
    this->defineDielectric();
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
  Real stemRadius     = 0.0;
  Real wheelSmooth    = 0.0;

  Vector<Real> v(3, 0.0);

  ParmParse pp("DiskProfiledPlane");

  // Get input parameters
  pp.get("wheel_live", isLive);
  pp.get("wheel_extra_thickness", wheelThickness);
  pp.get("wheel_radius", wheelRadius);
  pp.get("wheel_curvature", wheelCurvature);
  pp.get("wheel_stem_radius", stemRadius);
  pp.get("wheel_smooth", wheelSmooth);
  pp.getarr("wheel_center", v, 0, SpaceDim);

  // Create disk electrode using EBGeometry. This is an elongation of the union of a torus and a cylinder. Created
  // in the xy-plane and then put into place later.

  auto torus    = std::make_shared<TorusSDF<Real>>(Vec3::zero(), wheelRadius, wheelCurvature);
  auto cylinder = std::make_shared<CylinderSDF<Real>>(Vec3::unit(2) * wheelCurvature,
                                                      -Vec3::unit(2) * wheelCurvature,
                                                      wheelRadius);

  // Make smooth union and put into place.
  auto wheel = Union<Real>(torus, cylinder);
  if (stemRadius > 0.0) {
    auto holder = std::make_shared<CapsuleSDF<Real>>(Vec3::zero(), 1.E10 * Vec3::unit(1), stemRadius);
    wheel       = SmoothUnion<Real>(wheel, holder, wheelSmooth);
  }
  if (wheelThickness > 0.0) {
    wheel = Elongate<Real>(wheel, 0.5 * wheelThickness * Vec3::unit(2));
  }
  wheel = Rotate<Real>(wheel, 90.0, 1);
  wheel = Translate<Real>(wheel, Vec3(v[0], v[1], v[2]));

  if (SpaceDim == 2) {
    wheel = Elongate<Real>(wheel, std::numeric_limits<Real>::max() * Vec3::unit(2));
  }

  // Turn EBGeometry into Chombo
  m_electrodes.push_back(Electrode(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(wheel, true)), isLive));
}

void
DiskProfiledPlane::defineDielectric() noexcept
{
  CH_TIME("DiskProfiledPlane::defineDielectric");

  std::string str;

  Real boxCurvature   = 0.0;
  Real permittivity   = 1.0;
  Real sphereRadius   = 0.0;
  Real cylinderRadius = 0.0;

  Vector<Real> boxDimensions(3, std::numeric_limits<Real>::max());
  Vector<Real> boxTranslation(3, 0.0);
  Vector<Real> profileTranslate(3, 0.0);
  Vector<Real> profileTranslate2(3, 0.0);
  Vector<Real> profileRepetitionLo(3, 0.0);
  Vector<Real> profileRepetitionHi(3, 0.0);
  Vector<Real> profilePeriod(3, 0.0);
  Vector<Real> squareDimensions(3, std::numeric_limits<Real>::max());
  Vector<Real> squareDimensions2(3, std::numeric_limits<Real>::max());

  // Get input parameters.
  ParmParse pp("DiskProfiledPlane");

  pp.get("profile_type", str);
  pp.get("box_permittivity", permittivity);
  pp.get("box_curvature", boxCurvature);
  pp.get("sphere_radius", sphereRadius);
  pp.get("cylinder_radius", cylinderRadius);

  pp.getarr("box_dimensions", boxDimensions, 0, 3);
  pp.getarr("box_translate", boxTranslation, 0, 3);
  pp.getarr("profile_translate", profileTranslate, 0, 3);
  pp.getarr("profile_translate2", profileTranslate2, 0, 3);
  pp.getarr("profile_repetition_lo", profileRepetitionLo, 0, 3);
  pp.getarr("profile_repetition_hi", profileRepetitionHi, 0, 3);
  pp.getarr("profile_period", profilePeriod, 0, 3);
  pp.getarr("square_dimensions", squareDimensions, 0, 3);
  pp.getarr("square_dimensions2", squareDimensions2, 0, 3);

  if (boxCurvature <= 0.0) {
    MayDay::Error("DiskProfiledPlane::defineDielectric - must have 'box_curvature' > 0.0");
  }

  const Vec3 boxDims(boxDimensions[0] - 2 * boxCurvature,
                     boxDimensions[1] - 2 * boxCurvature,
                     boxDimensions[2] - 2 * boxCurvature);
  const Vec3 profileTra(profileTranslate[0], profileTranslate[1], profileTranslate[2]);
  const Vec3 profileTra2(profileTranslate2[0], profileTranslate2[1], profileTranslate2[2]);
  const Vec3 profileRepLo(profileRepetitionLo[0], profileRepetitionLo[1], profileRepetitionLo[2]);
  const Vec3 profileRepHi(profileRepetitionHi[0], profileRepetitionHi[1], profileRepetitionHi[2]);
  const Vec3 profilePer(profilePeriod[0], profilePeriod[1], profilePeriod[2]);
  const Vec3 squareDim(squareDimensions[0] - 2 * boxCurvature,
                       squareDimensions[1] - 2 * boxCurvature,
                       squareDimensions[2] - 2 * boxCurvature);
  const Vec3 squareDim2(squareDimensions2[0] - 2 * boxCurvature,
                       squareDimensions2[1] - 2 * boxCurvature,
                       squareDimensions2[2] - 2 * boxCurvature);


  // Basic rounded box.
  std::shared_ptr<ImpFunc> roundedBox = std::make_shared<RoundedBoxSDF<Real>>(boxDims, boxCurvature);
  roundedBox                          = Translate<Real>(roundedBox, -0.5 * boxDimensions[1] * Vec3::unit(1));

  // Determine the requested profile type.
  std::shared_ptr<ImpFunc> profile;
  std::shared_ptr<ImpFunc> profile2; 
  if (str == "square") {
    profile = std::make_shared<RoundedBoxSDF<Real>>(squareDim, boxCurvature);
  }
  else if (str == "sphere") {
    profile = std::make_shared<SphereSDF<Real>>(Vec3::zero(), sphereRadius);
  }
  else if (str == "cylinder_x") {
    profile = std::make_shared<SphereSDF<Real>>(Vec3::zero(), cylinderRadius);
    profile = Elongate<Real>(profile, std::numeric_limits<Real>::max() * Vec3::unit(0));
  }
  else if (str == "cylinder_y") {
    profile = std::make_shared<SphereSDF<Real>>(Vec3::zero(), cylinderRadius);
    profile = Elongate<Real>(profile, std::numeric_limits<Real>::max() * Vec3::unit(1));
  }
  else if (str == "cylinder_z") {
    profile = std::make_shared<SphereSDF<Real>>(Vec3::zero(), cylinderRadius);
    profile = Elongate<Real>(profile, std::numeric_limits<Real>::max() * Vec3::unit(2));
  }
  else if (str == "l_shape") {
    profile = std::make_shared<BoxSDF<Real>>(-0.5*squareDim, 0.5*squareDim);
    profile2 = std::make_shared<BoxSDF<Real>>(-0.5*squareDim2, 0.5*squareDim2);
	// profile2 = std::make_shared<Union<Real>>({profile, profile2});
	profile2 = EBGeometry::Union<Real>(profile, profile2);
  }
  else if (str == "none") {
  }
  else {
    MayDay::Abort("DiskProfiledPlane::defineDielectric - unsupported profile type requested");
  }

  // Translate and repeat the profiles. Then do a smooth union
  // if (str != "none") {
    // profile    = FiniteRepetition<Real>(profile, profilePer, profileRepLo, profileRepHi);
	// profile    = Translate<Real>(profile, profileTra);
    // roundedBox = SmoothDifference<Real>(roundedBox, profile, boxCurvature);
  // }
  if (str == "l_shape"){
      profile2    = FiniteRepetition<Real>(profile2, profilePer, profileRepLo, profileRepHi);
	  profile2    = Translate<Real>(profile2, profileTra2);
      roundedBox = Difference<Real>(roundedBox, profile2);
  }
  // Translate box into place.
  roundedBox = Translate<Real>(roundedBox, Vec3(boxTranslation[0], boxTranslation[1], boxTranslation[2]));

  if (SpaceDim == 2) {
    roundedBox = Elongate<Real>(roundedBox, std::numeric_limits<Real>::max() * Vec3::unit(2));
  }

  // Finally, create our dielectric.
  m_dielectrics.push_back(Dielectric(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(roundedBox, true)), permittivity));
}

#include <CD_NamespaceFooter.H>
