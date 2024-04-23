/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CoaxialCable.cpp
  @brief  Implementation of CD_CoaxialCable.H
  @author Robert Marskar
*/

// Std includes
#include <string>
#include <iostream>
#include <fstream>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CoaxialCable.H>
#include <CD_CylinderSdf.H>
#include <CD_SphereSdf.H>
#include <CD_NamespaceHeader.H>

CoaxialCable::CoaxialCable()
{

  ParmParse pp("CoaxialCable");
  ParmParse ppOuter("CoaxialCable.outer");
  ParmParse ppInner("CoaxialCable.inner");
  ParmParse ppMiddle("CoaxialCable.dielectric");

  // Get centers
  Vector<Real> v(SpaceDim);
  RealVect     e1, e2;
  pp.getarr("endpoint1", v, 0, SpaceDim);
  e1 = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("endpoint2", v, 0, SpaceDim);
  e2 = RealVect(D_DECL(v[0], v[1], v[2]));

  // Check which ones to use.
  bool outer, inner, middle;
  ppOuter.get("on", outer);
  ppInner.get("on", inner);
  ppMiddle.get("on", middle);

  // Reset
  m_electrodes.resize(0);
  m_dielectrics.resize(0);
  this->setGasPermittivity(1.0);

  if (outer) { // Add outer electrode
    Real rad;
    bool live;

    ppOuter.get("radius", rad);
    ppOuter.get("live", live);

#if CH_SPACEDIM == 2
    RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF>(new SphereSdf(e1, rad, true));
#elif CH_SPACEDIM == 3
    RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF>(new CylinderSdf(e1, e2, rad, true));
#endif
    m_electrodes.push_back(Electrode(baseif, live));
  }

  if (inner) { // Add inner electrode
    Real rad;
    bool live;

    ppInner.get("radius", rad);
    ppInner.get("live", live);

#if CH_SPACEDIM == 2
    RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF>(new SphereSdf(e1, rad, false));
#elif CH_SPACEDIM == 3
    RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF>(new CylinderSdf(e1, e2, rad, false));
#endif
    m_electrodes.push_back(Electrode(baseif, live));
  }

  if (middle) { // Add dielectric
    Real rad, eps;

    ppMiddle.get("radius", rad);
    ppMiddle.get("eps", eps);

#if CH_SPACEDIM == 2
    RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF>(new SphereSdf(e1, rad, false));
#elif CH_SPACEDIM == 3
    RefCountedPtr<BaseIF> baseif = RefCountedPtr<BaseIF>(new CylinderSdf(e1, e2, rad, false));
#endif
    m_dielectrics.push_back(Dielectric(baseif, eps));
  }
}

CoaxialCable::~CoaxialCable()
{}

#include <CD_NamespaceFooter.H>
