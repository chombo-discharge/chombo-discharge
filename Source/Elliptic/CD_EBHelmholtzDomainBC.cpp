/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDomainBC.cpp
  @brief  Implementation of CD_EBHelmholtzDomainBC.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int EBHelmholtzDomainBC::m_comp;
constexpr int EBHelmholtzDomainBC::m_nComp;

EBHelmholtzDomainBC::EBHelmholtzDomainBC()
{
  CH_TIME("EBHelmholtzDomainBC::EBHelmholtzDomainBC()");
}

EBHelmholtzDomainBC::~EBHelmholtzDomainBC()
{
  CH_TIME("EBHelmholtzDomainBC::~EBHelmholtzDomainBC()");
}

void
EBHelmholtzDomainBC::define(const Location::Cell a_dataLocation,
                            const EBLevelGrid&   a_eblg,
                            const RealVect&      a_probLo,
                            const Real           a_dx)
{
  CH_TIME("EBHelmholtzDomainBC::define(Location::Cell, EBLevelGrid, RefCountedPtr<LD<EBFluxFAB> >, RealVect, Real)");

  CH_assert(a_dx > 0.0);

  m_dataLocation = a_dataLocation;
  m_eblg         = a_eblg;
  m_probLo       = a_probLo;
  m_dx           = a_dx;
}

void
EBHelmholtzDomainBC::multiplyByBcoef(BaseFab<Real>&       a_flux,
                                     const BaseFab<Real>& a_bco,
                                     const int            a_dir,
                                     const Side::LoHiSide a_side) const
{
  CH_TIME("EBHelmholtzDomainBC::multiplyByBcoef");

  CH_assert(a_flux.nComp() == 1);
  CH_assert(a_bco.nComp() == 1);

  // For shifting the box over. a_flux is cell-centered and bco is face-centered so we need to get the indices right.
  IntVect shift;
  if (a_side == Side::Lo) {
    shift = IntVect::Zero;
  }
  else {
    shift = BASISV(a_dir);
  }

  // Kernel -- this just multiplies.
  auto kernel = [&](const IntVect& iv) {
    a_flux(iv, m_comp) *= a_bco(iv, m_comp);
  };

  // Execute kernel.
  BoxLoops::loop(a_flux.box(), kernel);
}

#include <CD_NamespaceFooter.H>
