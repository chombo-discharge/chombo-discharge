/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RoughSphere.cpp
  @brief  Implementation of CD_RoughSphere.H
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
#include <CD_PerlinSphereSdf.H>
#include <CD_RoughSphere.H>
#include <CD_NamespaceHeader.H>

RoughSphere::RoughSphere()
{
  this->setGasPermittivity(1.0);

  ParmParse pp("RoughSphere");

  bool        useSphere;
  std::string whichMaterial;

  pp.get("on", useSphere);
  pp.get("material", whichMaterial);

  if (useSphere) {
    bool         reseed, live;
    RealVect     f, c;
    Real         r, persist, amp, eps;
    int          octaves;
    Vector<Real> v;

    pp.get("radius", r);
    pp.get("live", live);
    pp.get("eps", eps);
    pp.get("noise_reseed", reseed);
    pp.get("noise_amplitude", amp);
    pp.get("noise_octaves", octaves);
    pp.get("noise_persistence", persist);

    pp.getarr("noise_frequency", v, 0, SpaceDim);
    f = RealVect(D_DECL(v[0], v[1], v[2]));
    pp.getarr("center", v, 0, SpaceDim);
    c = RealVect(D_DECL(v[0], v[1], v[2]));

    RefCountedPtr<BaseIF> sph = RefCountedPtr<BaseIF>(
      new PerlinSphereSdf(r, c, false, amp, f, persist, octaves, reseed));

    if (whichMaterial == "electrode")
      m_electrodes.push_back(Electrode(sph, live));
    else if (whichMaterial == "dielectric")
      m_dielectrics.push_back(Dielectric(sph, eps));
    else
      MayDay::Abort("RoughSphere::RoughSphere - unknown material requested");
  }
}

RoughSphere::~RoughSphere()
{}

#include <CD_NamespaceFooter.H>
