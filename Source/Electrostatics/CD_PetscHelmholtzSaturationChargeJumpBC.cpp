/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzSaturationChargeJumpBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzSaturationChargeJumpBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzSaturationChargeJumpBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzSaturationChargeJumpBC::PetscHelmholtzSaturationChargeJumpBC(const phase::which_phase a_phase) noexcept
  : m_phase(a_phase)
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::PetscHelmholtzSaturationChargeJumpBC");
}

PetscHelmholtzSaturationChargeJumpBC::~PetscHelmholtzSaturationChargeJumpBC() noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::~PetscHelmholtzSaturationChargeJumpBC");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::~PetscHelmholtzSaturationChargeJumpBC" << endl;
  }
}

void
PetscHelmholtzSaturationChargeJumpBC::clear() noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::clear");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::clear" << endl;
  }

  m_sigma.clear();

  PetscHelmholtzJumpBC::clear();
}

void
PetscHelmholtzSaturationChargeJumpBC::getInterfaceFluxStencil(Real& /*a_dataWeight*/,
                                                              PetscStencil& /*a_stencil*/,
                                                              const int /*a_phase*/,
                                                              const int /*a_level*/,
                                                              const DataIndex& /*a_din*/,
                                                              const IntVect& /*a_iv*/,
                                                              const Real /*a_bCoefPhase0*/,
                                                              const Real /*a_bCoefPhase1*/) const noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::getInterfaceFluxStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::getInterfaceFluxStencil" << endl;
  }

  MayDay::Abort("PetscHelmholtzSaturationChargeJumpBC::getInterfaceFluxStencil -- not implemented yet");
}

Real
PetscHelmholtzSaturationChargeJumpBC::getJumpData(const int /*a_level*/,
                                                  const DataIndex& /*a_din*/,
                                                  const IntVect& /*a_iv*/,
                                                  const EBAMRIVData& /*a_sigma*/) const noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::getJumpData");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::getJumpData" << endl;
  }

  return 0.0;
}

void
PetscHelmholtzSaturationChargeJumpBC::updateJumpData(const Vec& /*a_phi*/, const MFAMRIVData& /*a_bCoefIrreg*/) noexcept
{
  CH_TIME("PetscHelmholtzSaturationChargeJumpBC::updateJumpData");
  if (m_verbose) {
    pout() << "PetscHelmholtzSaturationChargeJumpBC::updateJumpData" << endl;
  }

  MayDay::Abort("PetscHelmholtzSaturationChargeJumpBC::updateJumpData -- not implemented yet");
}

#include <CD_NamespaceFooter.H>

#endif
