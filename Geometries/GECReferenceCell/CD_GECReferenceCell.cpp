/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_GECReferenceCell.cpp
  @brief  Implementation of CD_GECReferenceCell.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBGeometryIF.H>
#include <CD_GECReferenceCell.H>
#include <CD_NamespaceHeader.H>

using Vec3    = EBGeometry::Vec3T<Real>;
using ImpFunc = std::shared_ptr<EBGeometry::ImplicitFunction<Real>>;

GECReferenceCell::GECReferenceCell()
{
  using namespace EBGeometry;

  const Real H   = 0.1016;
  const Real g   = 0.0254;
  const Real h   = 0.5 * (H - g);
  const Real r1  = 0.0508;
  const Real r2  = 0.0523;
  const Real r3  = 0.0538;
  const Real r4  = 0.1016;
  const Real c   = (r2 - r1) / 2.0;
  const Real eps = 2.1;

  // Construct the lower (non-powered) electrode parts.
  const Vec3 zero = Vec3::zero();
  const Vec3 yhat = Vec3::unit(1);

  ImpFunc innerElectrodeLo, insulatorLo, outerElectrodeLo;
  ImpFunc innerElectrodeHi, insulatorHi, outerElectrodeHi;
  ImpFunc vessel;

  // Lower electrode and insulator objects
  innerElectrodeLo = std::make_shared<RoundedCylinderSDF<Real>>(0.5 * r1 + c, c, 1.0);
  innerElectrodeLo = EBGeometry::Translate<Real>(innerElectrodeLo, -0.5 * yhat - 0.5 * g * yhat);
  insulatorLo      = std::make_shared<CylinderSDF<Real>>(-yhat, -0.5 * g * yhat, 0.5 * (r2 + r3));
  outerElectrodeLo = std::make_shared<TorusSDF<Real>>(zero, r3 - c, c);
  outerElectrodeLo = EBGeometry::Rotate<Real>(outerElectrodeLo, 90, 0);
  outerElectrodeLo = EBGeometry::Elongate<Real>(outerElectrodeLo, 0.5 * yhat);
  outerElectrodeLo = EBGeometry::Translate<Real>(outerElectrodeLo, -0.5 * yhat - 0.5 * (g + 2 * c) * yhat);

  // Upper objects -- simple reflections of the lower objects.
  innerElectrodeHi = EBGeometry::Reflect<Real>(innerElectrodeLo, 1);
  insulatorHi      = EBGeometry::Reflect<Real>(insulatorLo, 1);
  outerElectrodeHi = EBGeometry::Reflect<Real>(outerElectrodeLo, 1);

  // Cylindrical vessel.
  vessel = std::make_shared<CylinderSDF<Real>>(-0.5 * H * yhat, 0.5 * H * yhat, r4);

  // Create electrodes and dielectrics with associated BCs for chombo-discharge
  m_electrodes.push_back(Electrode(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(innerElectrodeLo, true)), true));
  m_dielectrics.push_back(Dielectric(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(insulatorLo, true)), eps));
  m_electrodes.push_back(Electrode(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(outerElectrodeLo, true)), false));

  m_electrodes.push_back(Electrode(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(innerElectrodeHi, true)), false));
  m_dielectrics.push_back(Dielectric(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(insulatorHi, true)), eps));
  m_electrodes.push_back(Electrode(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(outerElectrodeHi, true)), false));

  m_electrodes.push_back(Electrode(RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(vessel, false)), false));
}

GECReferenceCell::~GECReferenceCell()
{}

#include <CD_NamespaceFooter.H>
