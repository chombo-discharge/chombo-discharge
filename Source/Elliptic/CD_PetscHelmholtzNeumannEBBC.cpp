/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzNeumannEBBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzNeumannEBBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzNeumannEBBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzNeumannEBBC::PetscHelmholtzNeumannEBBC() noexcept
{
  CH_TIME("PetscHelmholtzNeumannEBBC::PetscHelmholtzNeumannEBBC");

  m_useConstant    = false;
  m_useFunction    = false;
  m_constantDphiDn = 0.0;
}

PetscHelmholtzNeumannEBBC::~PetscHelmholtzNeumannEBBC() noexcept
{
  CH_TIME("PetscHelmholtzNeumannEBBC::~PetscHelmholtzNeumannEBBC");
  if (m_verbose) {
    pout() << "PetscHelmholtzNeumannEBBC::~PetscHelmholtzNeumannEBBC" << endl;
  }
}

void
PetscHelmholtzNeumannEBBC::setDphiDn(const Real a_DphiDn) noexcept
{
  CH_TIME("PetscHelmholtzNeumannEBBC::setDphiDn(Real)");
  if (m_verbose) {
    pout() << "PetscHelmholtzNeumannEBBC::setDphiDn(Real)" << endl;
  }

  m_useConstant    = true;
  m_useFunction    = false;
  m_constantDphiDn = a_DphiDn;
}

void
PetscHelmholtzNeumannEBBC::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn) noexcept
{
  CH_TIME("PetscHelmholtzNeumannEBBC::setDphiDn(std::function)");
  if (m_verbose) {
    pout() << "PetscHelmholtzNeumannEBBC::setDphiDn(std::function)" << endl;
  }

  m_useConstant    = false;
  m_useFunction    = true;
  m_functionDphiDn = a_DphiDn;
}

void
PetscHelmholtzNeumannEBBC::getEBFluxStencil(Real& /*a_boundaryWeight*/,
                                            PetscStencil& /*a_stencil*/,
                                            const int /*a_phase*/,
                                            const int /*a_level*/,
                                            const DataIndex& /*a_din*/,
                                            const VolIndex& /*a_vof*/) const noexcept
{
  CH_TIME("PetscHelmholtzNeumannEBBC::getEBFluxStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzNeumannEBBC::getEBFluxStencil" << endl;
  }

  MayDay::Abort("PetscHelmholtzNeumannEBBC::getEBFluxStencil -- not implemented yet");
}

Real
PetscHelmholtzNeumannEBBC::getBoundaryData(const int /*a_phase*/,
                                           const int /*a_level*/,
                                           const DataIndex& /*a_din*/,
                                           const VolIndex& /*a_vof*/) const noexcept
{
  CH_TIME("PetscHelmholtzNeumannEBBC::getBoundaryData");
  if (m_verbose) {
    pout() << "PetscHelmholtzNeumannEBBC::getBoundaryData" << endl;
  }

  MayDay::Abort("PetscHelmholtzNeumannEBBC::getBoundaryData -- not implemented yet");

  return 0.0;
}

#include <CD_NamespaceFooter.H>

#endif
