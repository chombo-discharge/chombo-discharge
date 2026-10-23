/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzElectrostaticEBBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzElectrostaticEBBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzElectrostaticEBBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzElectrostaticEBBC::PetscHelmholtzElectrostaticEBBC(const ElectrostaticEbBc& a_electrostaticBCs) noexcept
  : m_electrostaticBCs(a_electrostaticBCs)
{
  CH_TIME("PetscHelmholtzElectrostaticEBBC::PetscHelmholtzElectrostaticEBBC");
}

PetscHelmholtzElectrostaticEBBC::~PetscHelmholtzElectrostaticEBBC() noexcept
{
  CH_TIME("PetscHelmholtzElectrostaticEBBC::~PetscHelmholtzElectrostaticEBBC");
  if (m_verbose) {
    pout() << "PetscHelmholtzElectrostaticEBBC::~PetscHelmholtzElectrostaticEBBC" << endl;
  }
}

Real
PetscHelmholtzElectrostaticEBBC::getBoundaryData(const int /*a_phase*/,
                                                 const int /*a_level*/,
                                                 const DataIndex& /*a_din*/,
                                                 const VolIndex& /*a_vof*/) const noexcept
{
  CH_TIME("PetscHelmholtzElectrostaticEBBC::getBoundaryData");
  if (m_verbose) {
    pout() << "PetscHelmholtzElectrostaticEBBC::getBoundaryData" << endl;
  }

  MayDay::Abort("PetscHelmholtzElectrostaticEBBC::getBoundaryData -- not implemented yet");

  return 0.0;
}

Real
PetscHelmholtzElectrostaticEBBC::getElectrodePotential(const RealVect& /*a_pos*/) const noexcept
{
  CH_TIME("PetscHelmholtzElectrostaticEBBC::getElectrodePotential");
  if (m_verbose) {
    pout() << "PetscHelmholtzElectrostaticEBBC::getElectrodePotential" << endl;
  }

  MayDay::Abort("PetscHelmholtzElectrostaticEBBC::getElectrodePotential -- not implemented yet");

  return 0.0;
}

#include <CD_NamespaceFooter.H>

#endif
