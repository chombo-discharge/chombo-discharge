/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzNaturalJumpBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzNaturalJumpBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzNaturalJumpBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzNaturalJumpBC::PetscHelmholtzNaturalJumpBC() noexcept
{
  CH_TIME("PetscHelmholtzNaturalJumpBC::PetscHelmholtzNaturalJumpBC");
}

PetscHelmholtzNaturalJumpBC::~PetscHelmholtzNaturalJumpBC() noexcept
{
  CH_TIME("PetscHelmholtzNaturalJumpBC::~PetscHelmholtzNaturalJumpBC");
  if (m_verbose) {
    pout() << "PetscHelmholtzNaturalJumpBC::~PetscHelmholtzNaturalJumpBC" << endl;
  }
}

void
PetscHelmholtzNaturalJumpBC::getInterfaceFluxStencil(Real&            a_dataWeight,
                                                     PetscStencil&    a_stencil,
                                                     const int        a_phase,
                                                     const int        a_level,
                                                     const DataIndex& a_din,
                                                     const IntVect&   a_iv,
                                                     const Real       a_bCoefPhase0,
                                                     const Real       a_bCoefPhase1) const noexcept
{
  CH_TIME("PetscHelmholtzNaturalJumpBC::getInterfaceFluxStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzNaturalJumpBC::getInterfaceFluxStencil" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(a_phase == 0 || a_phase == 1);

  // The one-sided normal derivatives D_k = w_k * phi_B + S_k on both sides.
  Real         w[2];
  PetscStencil S[2];

  for (int k = 0; k < 2; k++) {
    this->getNormalDerivativeStencil(w[k], S[k], k, a_level, a_din, VolIndex(a_iv, 0));
  }

  const Real b[2] = {a_bCoefPhase0, a_bCoefPhase1};
  const Real den  = b[0] * w[0] + b[1] * w[1];

  // Checked in every build: phi_B is not determined without it.
  if (den == 0.0) {
    MayDay::Abort("PetscHelmholtzNaturalJumpBC::getInterfaceFluxStencil -- b_0 * w_0 + b_1 * w_1 is zero");
  }

  // b_0 D_0 + b_1 D_1 = sigma gives phi_B = (sigma - b_0 S_0 - b_1 S_1)/den, and then, with o the other side,
  //
  //   b_p D_p = b_p w_p / den * sigma + b_p b_o (w_o S_p - w_p S_o) / den.
  const int p = a_phase;
  const int o = 1 - a_phase;

  PetscStencil stencilThis  = S[p];
  PetscStencil stencilOther = S[o];

  stencilThis *= b[p] * b[o] * w[o] / den;
  stencilOther *= -b[p] * b[o] * w[p] / den;

  a_stencil = stencilThis;
  a_stencil += stencilOther;
  a_stencil.compress();

  a_dataWeight = b[p] * w[p] / den;
}

Real
PetscHelmholtzNaturalJumpBC::getJumpData(const int          a_level,
                                         const DataIndex&   a_din,
                                         const IntVect&     a_iv,
                                         const EBAMRIVData& a_sigma) const noexcept
{
  CH_TIME("PetscHelmholtzNaturalJumpBC::getJumpData");
  if (m_verbose) {
    pout() << "PetscHelmholtzNaturalJumpBC::getJumpData" << endl;
  }

  // The surface charge density lives on the cut cells of the first phase.
  return (*a_sigma[a_level])[a_din](VolIndex(a_iv, 0), 0);
}

void
PetscHelmholtzNaturalJumpBC::updateJumpData(const Vec& /*a_phi*/, const MFAMRIVData& /*a_bCoefIrreg*/) noexcept
{
  CH_TIME("PetscHelmholtzNaturalJumpBC::updateJumpData");
  if (m_verbose) {
    pout() << "PetscHelmholtzNaturalJumpBC::updateJumpData" << endl;
  }
}

#include <CD_NamespaceFooter.H>

#endif
