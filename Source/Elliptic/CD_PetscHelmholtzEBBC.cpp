/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzEBBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzEBBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_PetscHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzEBBC::PetscHelmholtzEBBC() noexcept
{
  CH_TIME("PetscHelmholtzEBBC::PetscHelmholtzEBBC");
}

PetscHelmholtzEBBC::~PetscHelmholtzEBBC() noexcept
{
  CH_TIME("PetscHelmholtzEBBC::~PetscHelmholtzEBBC");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBBC::~PetscHelmholtzEBBC" << endl;
  }
}

#include <CD_NamespaceFooter.H>

#endif
