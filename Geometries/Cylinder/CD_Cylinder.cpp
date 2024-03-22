/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Cylinder.cpp
  @brief  Implementation of CD_Cylinder.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_EBGeometryIF.H>
#include <CD_Cylinder.H>
#include <CD_NamespaceHeader.H>

Cylinder::Cylinder()
{
  bool turnOff      = false;
  bool fluidOutside = true;
  bool isElectrode  = true;
  bool isLive       = true;
  Real zCoord       = 0.0;
  Real radius       = 0.25;
  Real perm         = 1.0;

  Vector<Real> c1(3, 0.0);
  Vector<Real> c2(3, 0.0);

  std::string material;

  // Fetch input parameters
  ParmParse pp("Cylinder");
  pp.get("turn_off", turnOff);
  pp.get("fluid_outside", fluidOutside);
  pp.get("radius", radius);
  pp.get("z_coord", zCoord);
  pp.get("live", isLive);
  pp.get("permittivity", perm);
  pp.get("material", material);
  pp.getarr("center1", c1, 0, 3);
  pp.getarr("center2", c2, 0, 3);

  if (material == "electrode") {
    isElectrode = true;
  }
  else if (material == "dielectric") {
    isElectrode = false;
  }
  else {
    MayDay::Error("Cylinder::Cylinder() -- unrecognized material type");
  }

  if (!turnOff) {

    // Create the 3D cylinder
    const EBGeometry::Vec3T<Real> v1(c1[0], c1[1], c1[2]);
    const EBGeometry::Vec3T<Real> v2(c2[0], c2[1], c2[2]);

    auto cylinder = std::make_shared<EBGeometry::CylinderSDF<Real>>(v1, v2, radius);

    RefCountedPtr<BaseIF> baseIF = RefCountedPtr<BaseIF>(new EBGeometryIF<Real>(cylinder, fluidOutside, zCoord));

    if (isElectrode) {
      m_electrodes.push_back(Electrode(baseIF, isLive));
    }
    else {
      m_dielectrics.push_back(Dielectric(baseIF, perm));
    }
  }
}

Cylinder::~Cylinder()
{}

#include <CD_NamespaceFooter.H>
