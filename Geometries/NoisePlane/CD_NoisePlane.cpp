/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NoisePlane.cpp
  @brief  Implementation of CD_NoisePlane.H
  @author Robert Marskar
*/

// Std includes
#include <string>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_NoisePlane.H>
#include <CD_BoundedNoisePlane.H>
#include <CD_NamespaceHeader.H>

NoisePlane::NoisePlane()
{

  this->setGasPermittivity(1.0);

  ParmParse pp("NoisePlane");

  bool usePlane;
  pp.get("use_plane", usePlane);

  if (usePlane) {

    std::string orientation;
    std::string material;

    int normal;
    int noiseOctaves;

    RealVect point;
    RealVect clampLo;
    RealVect clampHi;
    RealVect noiseFreq;

    Real eps;
    Real clampK;
    Real noiseAmp;
    Real noisePersist;

    bool live;
    bool reseed;

    auto get = [&](const std::string id) -> RealVect {
      Vector<Real> v;
      pp.getarr(id.c_str(), v, 0, SpaceDim);

      return RealVect(D_DECL(v[0], v[1], v[2]));
    };

    pp.get("live", live);
    pp.get("eps", eps);
    pp.get("orientation", orientation);
    pp.get("clamp_dx", clampK);
    pp.get("noise_amplitude", noiseAmp);
    pp.get("noise_persistence", noisePersist);
    pp.get("noise_reseed", reseed);
    pp.get("noise_octaves", noiseOctaves);
    pp.get("noise_reseed", reseed);
    pp.get("material", material);

    point     = get("point");
    clampLo   = get("clamp_lo");
    clampHi   = get("clamp_hi");
    noiseFreq = get("noise_frequency");

    RefCountedPtr<BaseIF> plane = RefCountedPtr<BaseIF>(new BoundedNoisePlane(orientation,
                                                                              point,
                                                                              clampLo,
                                                                              clampHi,
                                                                              1.0 / clampK,
                                                                              noiseAmp,
                                                                              noiseFreq,
                                                                              noisePersist,
                                                                              noiseOctaves,
                                                                              reseed));

    if (material == "electrode") {
      m_electrodes.push_back(Electrode(plane, live));
    }
    else if (material == "dielectric") {
      m_dielectrics.push_back(Dielectric(plane, eps));
    }
    else {
      MayDay::Error("NoisePlane::NoisePlane -- unknown mateiral requested");
    }
  }
}

NoisePlane::~NoisePlane() {}

#include <CD_NamespaceFooter.H>
