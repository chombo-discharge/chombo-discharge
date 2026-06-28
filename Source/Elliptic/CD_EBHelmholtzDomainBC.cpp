/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
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

Real
EBHelmholtzDomainBC::getDiagWeight(const int /*a_dir*/, const Side::LoHiSide /*a_side*/) const
{
  return 1.0;
}

void
EBHelmholtzDomainBC::multiplyByBcoef(BaseFab<Real>&       a_flux,
                                     const BaseFab<Real>& a_bco,
                                     const int            a_dir,
                                     const Side::LoHiSide a_side)
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
  //
  // NOTE on box centering: a_flux is cell-centered (one cell per boundary face) while a_bco is
  // face-centered, so on the Hi side cell 'iv' actually maps to boundary face iv+BASISV(a_dir) (the
  // cell's hi face) -- i.e. the two are offset by one there, which is what 'shift' above represents.
  // Even so, we deliberately use the UNSHIFTED a_bco(iv) and 'shift' is left unused. The reason is the
  // operator path: EBHelmholtzOp::applyDomainFlux multiplies by this bco here and then divides it
  // straight back out with the SAME unshifted index (its ghost-cell trick needs the raw dphi/dn), so
  // the factor cancels and the true face b-coefficient is reintroduced later by the finite-volume
  // Laplacian stencil. Using a_bco(iv+shift) here would break that cancellation and corrupt the
  // operator for spatially-varying coefficients. (The only non-cancelling consumer, fillDomainFlux,
  // therefore scales its Hi-side boundary flux by bco(iv) rather than the boundary-face bco(iv+shift);
  // this is exact for constant coefficients and a known approximation for variable ones.)
  auto kernel = [&](const IntVect& iv) {
    a_flux(iv, m_comp) *= a_bco(iv, m_comp);
  };

  // Execute kernel.
  BoxLoops::loop<D_DECL(1, 1, 1)>(a_flux.box(), kernel);
}

#include <CD_NamespaceFooter.H>
