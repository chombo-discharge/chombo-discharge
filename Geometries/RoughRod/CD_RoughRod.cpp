/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RoughRod.cpp
  @brief  Implementation of CD_RoughRod.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_RoughRod.H>
#include <CD_PerlinRodSdf.H>
#include <CD_NamespaceHeader.H>

RoughRod::RoughRod()
{
  this->setGasPermittivity(1.0);

  ParmParse pp("RoughRod");

  bool useRod;

  pp.get("on", useRod);

  if (useRod) {
    bool         live, reseed;
    RealVect     f, e1, e2;
    Real         r, persist, amp;
    int          octaves;
    Vector<Real> v;

    pp.get("radius", r);
    pp.get("live", live);
    pp.get("noise_reseed", reseed);
    pp.get("noise_amplitude", amp);
    pp.get("noise_octaves", octaves);
    pp.get("noise_persistence", persist);

    pp.getarr("noise_frequency", v, 0, SpaceDim);
    f = RealVect(D_DECL(v[0], v[1], v[2]));
    pp.getarr("endpoint1", v, 0, SpaceDim);
    e1 = RealVect(D_DECL(v[0], v[1], v[2]));
    pp.getarr("endpoint2", v, 0, SpaceDim);
    e2 = RealVect(D_DECL(v[0], v[1], v[2]));

    RefCountedPtr<BaseIF> rod = RefCountedPtr<BaseIF>(
      new PerlinRodSdf(r, e1, e2, false, amp, f, persist, octaves, reseed));

    m_electrodes.push_back(Electrode(rod, live));
  }
}

RoughRod::~RoughRod()
{}

#include <CD_NamespaceFooter.H>
