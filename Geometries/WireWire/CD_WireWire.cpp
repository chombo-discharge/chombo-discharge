/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_WireWire.cpp
  @brief  Implementation of CD_WireWire.H
  @author Robert Marskar
*/

// Our includes
#include <CD_WireWire.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

using Vec3    = EBGeometry::Vec3T<Real>;
using ImpFunc = EBGeometry::ImplicitFunction<Real>;

WireWire::WireWire()
{
  CH_TIME("WireWire::WireWire");

  m_electrodes.resize(0);
  m_dielectrics.resize(0);

  ParmParse pp("WireWire");
  ParmParse pp1("WireWire.first");
  ParmParse pp2("WireWire.second");

  bool addWire1 = false;
  bool addWire2 = false;

  Real insulationPermittivity = -1.0;
  Real smoothLen              = 0.0;

  pp.get("insulation_permittivity", insulationPermittivity);
  pp.get("smooth_length", smoothLen);
  pp1.get("on", addWire1);
  pp2.get("on", addWire2);

  std::vector<std::pair<std::shared_ptr<ImpFunc>, Real>> electrodes;
  std::vector<std::shared_ptr<ImpFunc>>                  dielectrics;

  if (addWire1) {
    electrodes.emplace_back(this->addElectrode(pp1));
    dielectrics.emplace_back(this->addDielectric(pp1));
  }
  if (addWire2) {
    electrodes.emplace_back(this->addElectrode(pp2));
    dielectrics.emplace_back(this->addDielectric(pp2));
  }

  // Create electrodes
  for (int i = 0; i < electrodes.size(); i++) {
    const auto baseif = RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(electrodes[i].first, true));
    m_electrodes.push_back(Electrode(baseif, true, electrodes[i].second));
  }

  // Create the insulation object
  std::shared_ptr<ImpFunc> insulation;
  if (dielectrics.size() > 0 && insulationPermittivity > 0.0) {
    if (smoothLen > 0.0) {
      insulation = EBGeometry::SmoothUnion<Real>(dielectrics, smoothLen);
    }
    else {
      insulation = EBGeometry::Union<Real>(dielectrics);
    }

    const auto baseif = RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(insulation, true));

    m_dielectrics.push_back(Dielectric(baseif, insulationPermittivity));
  }
}

WireWire::~WireWire()
{
  CH_TIME("WireWire::~WireWire");
}

std::pair<std::shared_ptr<EBGeometry::ImplicitFunction<Real>>, Real>
WireWire::addElectrode(ParmParse& a_pp)
{
  CH_TIME("WireWire::addElectrode");

  Real         electrodeRadius;
  Real         potential;
  Vector<Real> v;

  a_pp.get("electrode_radius", electrodeRadius);
  a_pp.get("potential", potential);
  a_pp.getarr("center", v, 0, 2);

  std::pair<std::shared_ptr<EBGeometry::ImplicitFunction<Real>>, Real> ret;

  if (electrodeRadius > 0.0) {
    const auto cylinder = std::make_shared<EBGeometry::InfiniteCylinderSDF<Real>>(Vec3(v[0], v[1], 0.0),
                                                                                  electrodeRadius,
                                                                                  2);

    ret = std::make_pair(cylinder, potential);
  }

  return ret;
}

std::shared_ptr<EBGeometry::ImplicitFunction<Real>>
WireWire::addDielectric(ParmParse& a_pp)
{
  CH_TIME("WireWire::addDielectric");

  std::shared_ptr<EBGeometry::ImplicitFunction<Real>> ret = nullptr;

  ParmParse pp("WireWire");

  Real         electrodeRadius;
  Real         dielectricRadius;
  Vector<Real> v;

  a_pp.get("electrode_radius", electrodeRadius);
  a_pp.get("dielectric_radius", dielectricRadius);
  a_pp.getarr("center", v, 0, 2);

  if (dielectricRadius > electrodeRadius) {
    ret = std::make_shared<EBGeometry::InfiniteCylinderSDF<Real>>(Vec3(v[0], v[1], 0.0), dielectricRadius, 2);
  }

  return ret;
}

#include <CD_NamespaceFooter.H>
