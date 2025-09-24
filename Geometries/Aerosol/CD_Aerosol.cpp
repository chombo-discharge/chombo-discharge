/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Aerosol.cpp
  @brief  Implementation of CD_Aerosol.H
  @author Robert Marskar
*/

// Std includes
#include <string>
#include <iostream>
#include <fstream>

// Chombo includes
#include <ParmParse.H>
#include <BaseIF.H>

// Our includes
#include <CD_Aerosol.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

Aerosol::Aerosol()
{

  Real gasPermittivity   = 1.0;
  Real solidPermittivity = 1.0;
  bool invert;
  int  numSpheres;

  Real         radius;
  Vector<Real> v(SpaceDim);

  // Get a bunch of parameters from the input script
  ParmParse pp("Aerosol");

  pp.get("invert", invert);
  pp.get("eps0", gasPermittivity);
  pp.get("permittivity", solidPermittivity);
  pp.get("num_spheres", numSpheres);

  this->setGasPermittivity(gasPermittivity);

  if (gasPermittivity <= 0.0) {
    MayDay::Error("Aerosol::Aerosol() -- eps0 cannot be <= 0.0");
  }
  if (solidPermittivity <= 0.0) {
    MayDay::Error("Aerosol::Aerosol() -- permittivity cannot be <= 0.0");
  }

  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  // Construct tht CSG union of spheres
  if (numSpheres > 0) {
    std::vector<std::shared_ptr<EBGeometry::ImplicitFunction<Real>>> spheres;

    for (int i = 0; i < numSpheres; i++) {
      const std::string str = "Aerosol.sphere" + std::to_string(1 + i);

      // Get center and radius.
      ParmParse pp2(str);
      pp2.get("radius", radius);
      pp2.getarr("center", v, 0, SpaceDim);

      auto c = EBGeometry::Vec3T<Real>::zero();
      for (int dir = 0; dir < SpaceDim; dir++) {
        c[dir] = v[dir];
      }

      spheres.emplace_back(std::make_shared<EBGeometry::SphereSDF<Real>>(c, radius));
    }

    const auto unionChombo = RefCountedPtr<BaseIF>(new EBGeometryIF<>(EBGeometry::Union<Real>(spheres), !invert, 0.0));

    m_dielectrics.push_back(Dielectric(unionChombo, solidPermittivity));
  }
}

Aerosol::~Aerosol()
{}

#include <CD_NamespaceFooter.H>
