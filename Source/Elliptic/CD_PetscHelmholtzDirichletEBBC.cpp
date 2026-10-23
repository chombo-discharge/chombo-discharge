/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzDirichletEBBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzDirichletEBBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzDirichletEBBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzDirichletEBBC::PetscHelmholtzDirichletEBBC() noexcept
{
  CH_TIME("PetscHelmholtzDirichletEBBC::PetscHelmholtzDirichletEBBC");

  m_useConstant   = false;
  m_useFunction   = false;
  m_constantValue = 0.0;
}

PetscHelmholtzDirichletEBBC::~PetscHelmholtzDirichletEBBC() noexcept
{
  CH_TIME("PetscHelmholtzDirichletEBBC::~PetscHelmholtzDirichletEBBC");
  if (m_verbose) {
    pout() << "PetscHelmholtzDirichletEBBC::~PetscHelmholtzDirichletEBBC" << endl;
  }
}

void
PetscHelmholtzDirichletEBBC::setValue(const Real a_value) noexcept
{
  CH_TIME("PetscHelmholtzDirichletEBBC::setValue(Real)");
  if (m_verbose) {
    pout() << "PetscHelmholtzDirichletEBBC::setValue(Real)" << endl;
  }

  m_useConstant   = true;
  m_useFunction   = false;
  m_constantValue = a_value;
}

void
PetscHelmholtzDirichletEBBC::setValue(const std::function<Real(const RealVect& a_pos)>& a_value) noexcept
{
  CH_TIME("PetscHelmholtzDirichletEBBC::setValue(std::function)");
  if (m_verbose) {
    pout() << "PetscHelmholtzDirichletEBBC::setValue(std::function)" << endl;
  }

  m_useConstant   = false;
  m_useFunction   = true;
  m_functionValue = a_value;
}

void
PetscHelmholtzDirichletEBBC::getEBFluxStencil(Real& /*a_boundaryWeight*/,
                                              PetscStencil& /*a_stencil*/,
                                              const int /*a_phase*/,
                                              const int /*a_level*/,
                                              const DataIndex& /*a_din*/,
                                              const VolIndex& /*a_vof*/) const noexcept
{
  CH_TIME("PetscHelmholtzDirichletEBBC::getEBFluxStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzDirichletEBBC::getEBFluxStencil" << endl;
  }

  MayDay::Abort("PetscHelmholtzDirichletEBBC::getEBFluxStencil -- not implemented yet");
}

Real
PetscHelmholtzDirichletEBBC::getBoundaryData(const int /*a_phase*/,
                                             const int /*a_level*/,
                                             const DataIndex& /*a_din*/,
                                             const VolIndex& /*a_vof*/) const noexcept
{
  CH_TIME("PetscHelmholtzDirichletEBBC::getBoundaryData");
  if (m_verbose) {
    pout() << "PetscHelmholtzDirichletEBBC::getBoundaryData" << endl;
  }

  MayDay::Abort("PetscHelmholtzDirichletEBBC::getBoundaryData -- not implemented yet");

  return 0.0;
}

#include <CD_NamespaceFooter.H>

#endif
