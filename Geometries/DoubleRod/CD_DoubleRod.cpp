/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DoubleRod.cpp
  @brief  Implementation of CD_DoubleRod.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_DoubleRod.H>
#include <CD_RodIF.H>
#include <CD_NamespaceHeader.H>

DoubleRod::DoubleRod()
{
  this->setGasPermittivity(1.0);

  ParmParse pp1("DoubleRod.rod1");
  ParmParse pp2("DoubleRod.rod2");

  bool use_rod1;
  bool use_rod2;

  pp1.get("on", use_rod1);
  pp2.get("on", use_rod2);

  Vector<Real> v(SpaceDim);
  Real         radius;
  RealVect     e1, e2;
  bool         live;

  if (use_rod1) {
    pp1.get("radius", radius);
    pp1.get("live", live);
    pp1.getarr("endpoint1", v, 0, SpaceDim);
    e1 = RealVect(D_DECL(v[0], v[1], v[2]));
    pp1.getarr("endpoint2", v, 0, SpaceDim);
    e2 = RealVect(D_DECL(v[0], v[1], v[2]));

    RefCountedPtr<BaseIF> rod1 = RefCountedPtr<BaseIF>(new RodIF(e1, e2, radius, false));

    m_electrodes.push_back(Electrode(rod1, live));
  }

  if (use_rod2) {
    pp2.get("radius", radius);
    pp2.get("live", live);
    pp2.getarr("endpoint1", v, 0, SpaceDim);
    e1 = RealVect(D_DECL(v[0], v[1], v[2]));
    pp2.getarr("endpoint2", v, 0, SpaceDim);
    e2 = RealVect(D_DECL(v[0], v[1], v[2]));

    RefCountedPtr<BaseIF> rod1 = RefCountedPtr<BaseIF>(new RodIF(e1, e2, radius, false));

    m_electrodes.push_back(Electrode(rod1, live));
  }
}

DoubleRod::~DoubleRod()
{}

#include <CD_NamespaceFooter.H>
