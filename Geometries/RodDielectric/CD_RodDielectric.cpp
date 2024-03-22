/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RodDielectric.cpp
  @brief  Implementation of CD_RodDielectric.H
  @author Robert Marskar
*/

// Chombo includes
#include <PlaneIF.H>
#include <ParmParse.H>

// Our includes
#include <CD_RodDielectric.H>
#include <CD_RodIF.H>
#include <CD_PerlinSlabSdf.H>
#include <CD_SphereSdf.H>
#include <CD_WedgeIF.H>
#include <CD_RoundedBoxIF.H>
#include <CD_NamespaceHeader.H>

RodDielectric::RodDielectric()
{
  m_dielectrics.resize(0);
  m_electrodes.resize(0);

  this->setGasPermittivity(1.0);

  ParmParse ppRod("RodDielectric.electrode");
  ParmParse ppIns("RodDielectric.dielectric");

  bool useRod;
  bool useIns;

  ppRod.get("on", useRod);
  ppIns.get("on", useIns);

  if (useRod)
    this->defineElectrode();
  if (useIns)
    this->defineInsulator();
}

RodDielectric::~RodDielectric()
{}

void
RodDielectric::defineElectrode()
{
  Vector<Real> v(SpaceDim);
  RealVect     e1, e2;
  Real         r;
  bool         live;

  ParmParse pp("RodDielectric.electrode");

  pp.get("radius", r);
  pp.get("live", live);
  pp.getarr("endpoint1", v, 0, SpaceDim);
  e1 = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("endpoint2", v, 0, SpaceDim);
  e2 = RealVect(D_DECL(v[0], v[1], v[2]));

  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF>(new RodIF(e1, e2, r, false));

  m_electrodes.push_back(Electrode(bif, live));
}

void
RodDielectric::defineInsulator()
{
  ParmParse pp("RodDielectric.dielectric");

  std::string str;
  Real        eps;

  pp.get("shape", str);
  pp.get("permittivity", eps);

  RefCountedPtr<BaseIF> bif;
  if (str == "plane") {
    bif = this->getPlane();
  }
  else if (str == "box") {
    bif = this->getBox();
  }
  else if (str == "perlin_box") {
    bif = this->getPerlinBox();
  }
  else if (str == "sphere") {
    bif = this->getSphere();
  }
  else {
    MayDay::Abort("RodDielectric.:defineInsulator - unsupported shape");
  }

  m_dielectrics.push_back(Dielectric(bif, eps));
}

RefCountedPtr<BaseIF>
RodDielectric::getBox()
{
  ParmParse pp("RodDielectric.box");

  Vector<Real> v(SpaceDim);
  RealVect     lo, hi;
  Real         curv;

  pp.get("curvature", curv);
  pp.getarr("lo_corner", v, 0, SpaceDim);
  lo = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("hi_corner", v, 0, SpaceDim);
  hi = RealVect(D_DECL(v[0], v[1], v[2]));

  return RefCountedPtr<BaseIF>(new RoundedBoxIF(lo, hi, curv, false));
}

RefCountedPtr<BaseIF>
RodDielectric::getPlane()
{
  ParmParse pp("RodDielectric.plane");

  Vector<Real> v(SpaceDim);
  RealVect     p, n;

  pp.getarr("point", v, 0, SpaceDim);
  p = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("normal", v, 0, SpaceDim);
  n = RealVect(D_DECL(v[0], v[1], v[2]));

  return RefCountedPtr<BaseIF>(new PlaneIF(n, p, true));
}

RefCountedPtr<BaseIF>
RodDielectric::getSphere()
{
  ParmParse pp("RodDielectric.sphere");

  Vector<Real> v(SpaceDim);
  RealVect     p;
  Real         r;

  pp.get("radius", r);
  pp.getarr("center", v, 0, SpaceDim);
  p = RealVect(D_DECL(v[0], v[1], v[2]));

  return RefCountedPtr<BaseIF>(new SphereSdf(p, r, false));
}

RefCountedPtr<BaseIF>
RodDielectric::getPerlinBox()
{
  ParmParse pp("RodDielectric.perlin_box");

  Vector<Real> v(SpaceDim);
  RealVect     p, n, xyz, freq;
  Real         amp, persist, curv;
  int          octaves;
  bool         reseed;

  pp.get("curvature", curv);
  pp.get("noise_octaves", octaves);
  pp.get("noise_amp", amp);
  pp.get("noise_reseed", reseed);
  pp.get("noise_persist", persist);

  pp.getarr("point", v, 0, SpaceDim);
  p = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("normal", v, 0, SpaceDim);
  n = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("dimensions", v, 0, SpaceDim);
  xyz = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("noise_freq", v, 0, SpaceDim);
  freq = RealVect(D_DECL(v[0], v[1], v[2]));

  return RefCountedPtr<BaseIF>(new PerlinSlabSdf(p, n, xyz, freq, octaves, amp, persist, curv, reseed, false));
}

#include <CD_NamespaceFooter.H>
