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
PetscHelmholtzNaturalJumpBC::getInterfaceFluxStencil(Real& /*a_dataWeight*/,
                                                     PetscStencil& /*a_stencil*/,
                                                     const int /*a_phase*/,
                                                     const int /*a_level*/,
                                                     const DataIndex& /*a_din*/,
                                                     const IntVect& /*a_iv*/,
                                                     const Real /*a_bCoefPhase0*/,
                                                     const Real /*a_bCoefPhase1*/) const noexcept
{
  CH_TIME("PetscHelmholtzNaturalJumpBC::getInterfaceFluxStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzNaturalJumpBC::getInterfaceFluxStencil" << endl;
  }

  MayDay::Abort("PetscHelmholtzNaturalJumpBC::getInterfaceFluxStencil -- not implemented yet");
}

Real
PetscHelmholtzNaturalJumpBC::getJumpData(const int /*a_level*/,
                                         const DataIndex& /*a_din*/,
                                         const IntVect& /*a_iv*/,
                                         const EBAMRIVData& /*a_sigma*/) const noexcept
{
  CH_TIME("PetscHelmholtzNaturalJumpBC::getJumpData");
  if (m_verbose) {
    pout() << "PetscHelmholtzNaturalJumpBC::getJumpData" << endl;
  }

  MayDay::Abort("PetscHelmholtzNaturalJumpBC::getJumpData -- not implemented yet");

  return 0.0;
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
