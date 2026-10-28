/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzJumpBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzJumpBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_PetscHelmholtzJumpBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzJumpBC::PetscHelmholtzJumpBC() noexcept
{
  CH_TIME("PetscHelmholtzJumpBC::PetscHelmholtzJumpBC");
}

PetscHelmholtzJumpBC::~PetscHelmholtzJumpBC() noexcept
{
  CH_TIME("PetscHelmholtzJumpBC::~PetscHelmholtzJumpBC");
  if (m_verbose) {
    pout() << "PetscHelmholtzJumpBC::~PetscHelmholtzJumpBC" << endl;
  }
}

#include <CD_NamespaceFooter.H>

#endif
