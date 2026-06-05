/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SpherePlane.cpp
  @brief  High-voltage SpherePlane geometry
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_SpherePlane.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

SpherePlane::SpherePlane() noexcept
{
  CH_TIME("SpherePlane::SpherePlane");

  using ImpFunc = EBGeometry::ImplicitFunction<Real>;
  using Vec3    = EBGeometry::Vec3T<Real>;

  bool useSphere;
  bool useDisk;
  bool useSphereHolder;
  bool useDiskHolder;

  bool liveSphere;
  bool liveDisk;

  Real sphereCenter;
  Real sphereRadius;
  Real sphereHolderRadius;
  Real sphereHolderLength;
  Real sphereHolderSmooth;

  Real diskRadius;
  Real diskHolderRadius;
  Real diskHolderLength;
  Real diskThickness;
  Real diskCurvature;
  Real diskCenter;
  Real diskHolderSmooth;

  ParmParse pp("SpherePlane");

  pp.get("use_sphere", useSphere);
  pp.get("use_disk", useDisk);
  pp.get("use_sphere_holder", useSphereHolder);
  pp.get("use_disk_holder", useDiskHolder);

  pp.get("live_sphere", liveSphere);
  pp.get("live_disk", liveDisk);

  pp.get("sphere_radius", sphereRadius);
  pp.get("sphere_center", sphereCenter);
  pp.get("sphere_holder_radius", sphereHolderRadius);
  pp.get("sphere_holder_length", sphereHolderLength);
  pp.get("sphere_holder_smooth", sphereHolderSmooth);

  pp.get("disk_radius", diskRadius);
  pp.get("disk_center", diskCenter);
  pp.get("disk_holder_radius", diskHolderRadius);
  pp.get("disk_holder_length", diskHolderLength);
  pp.get("disk_holder_smooth", diskHolderSmooth);
  pp.get("disk_thickness", diskThickness);
  pp.get("disk_curvature", diskCurvature);

  if (useSphere) {
    std::shared_ptr<ImpFunc> sphere;
    std::shared_ptr<ImpFunc> holder;
    std::shared_ptr<ImpFunc> arrangement;

    sphere = std::make_shared<EBGeometry::SphereSDF<Real>>(Vec3::zero(), sphereRadius);
    holder = std::make_shared<EBGeometry::CylinderSDF<Real>>(Vec3::zero(),
                                                             sphereHolderLength * Vec3::unit(1),
                                                             sphereHolderRadius);

    arrangement = useSphereHolder ? EBGeometry::SmoothUnion<Real>(sphere, holder, sphereHolderSmooth) : sphere;
    arrangement = EBGeometry::Translate<Real>(arrangement, sphereCenter * Vec3::unit(1));

    m_electrodes.push_back(Electrode(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(arrangement, true)), liveSphere));
  }

  if (useDisk) {
    std::shared_ptr<ImpFunc> disk;
    std::shared_ptr<ImpFunc> holder;
    std::shared_ptr<ImpFunc> arrangement;

    disk   = std::make_shared<EBGeometry::RoundedCylinderSDF<Real>>(0.5 * (diskRadius + diskCurvature),
                                                                  diskCurvature,
                                                                  diskThickness);
    holder = std::make_shared<EBGeometry::CylinderSDF<Real>>(Vec3::zero(),
                                                             -diskHolderLength * Vec3::unit(1),
                                                             diskHolderRadius);

    arrangement = useDiskHolder ? EBGeometry::SmoothUnion<Real>(disk, holder, diskHolderSmooth) : disk;
    arrangement = EBGeometry::Translate<Real>(arrangement, -(diskCenter + 0.5 * diskThickness) * Vec3::unit(1));

    m_electrodes.push_back(Electrode(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(arrangement, true)), liveDisk));
  }
}

SpherePlane::~SpherePlane() noexcept
{
  CH_TIME("SpherePlane::~SpherePlane");
}

#include <CD_NamespaceFooter.H>
