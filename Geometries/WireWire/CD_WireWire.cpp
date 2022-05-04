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
#include <CD_SphereSdf.H>
#include <CD_CylinderSdf.H>
#include <CD_NamespaceHeader.H>

WireWire::WireWire()
{
  ParmParse pp1("WireWire.first");
  ParmParse pp2("WireWire.second");

  bool useFirst;
  bool useSecond;

  pp1.get("on", useFirst);
  pp2.get("on", useSecond);

  if (useFirst)
    this->addWire(pp1);
  if (useSecond)
    this->addWire(pp2);
}

WireWire::~WireWire() {}

void
WireWire::addWire(ParmParse& a_pp)
{
  Real     r;
  RealVect e1, e2;
  bool     live;

  Vector<Real> v;

  a_pp.get("radius", r);
  a_pp.get("live", live);

  a_pp.getarr("endpoint1", v, 0, SpaceDim);
  e1 = RealVect(D_DECL(v[0], v[1], v[2]));
  a_pp.getarr("endpoint2", v, 0, SpaceDim);
  e2 = RealVect(D_DECL(v[0], v[1], v[2]));

  Real pot = live ? 0.5 : -0.5;

#if CH_SPACEDIM == 2
  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF>(new SphereSdf(e1, r, false));
#elif CH_SPACEDIM == 3
  RefCountedPtr<BaseIF> bif = RefCountedPtr<BaseIF>(new CylinderSdf(e1, e2, r, false));
#endif

  m_electrodes.push_back(Electrode(bif, true, pot));
}

#include <CD_NamespaceFooter.H>
