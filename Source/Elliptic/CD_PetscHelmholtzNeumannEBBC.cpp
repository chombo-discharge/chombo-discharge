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
PetscHelmholtzNeumannEBBC::getEBFluxStencil(Real&         a_boundaryWeight,
                                            PetscStencil& a_stencil,
                                            const int /*a_phase*/,
                                            const int /*a_level*/,
                                            const DataIndex& /*a_din*/,
                                            const VolIndex& /*a_vof*/) const noexcept
{
  CH_TIME("PetscHelmholtzNeumannEBBC::getEBFluxStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzNeumannEBBC::getEBFluxStencil" << endl;
  }

  // The flux is the given normal derivative g.
  a_boundaryWeight = 1.0;
  a_stencil.clear();
}

Real
PetscHelmholtzNeumannEBBC::getBoundaryData(const int        a_phase,
                                           const int        a_level,
                                           const DataIndex& a_din,
                                           const VolIndex&  a_vof) const noexcept
{
  CH_TIME("PetscHelmholtzNeumannEBBC::getBoundaryData");
  if (m_verbose) {
    pout() << "PetscHelmholtzNeumannEBBC::getBoundaryData" << endl;
  }

  CH_assert(m_isDefined);

  Real DphiDn = 0.0;

  if (m_useConstant) {
    DphiDn = m_constantDphiDn;
  }
  else if (m_useFunction) {
    DphiDn = m_functionDphiDn(this->getBoundaryPosition(a_phase, a_level, a_din, a_vof));
  }
  else {
    // Checked in every build: without a value the right-hand side would silently use dphi/dn = 0.
    MayDay::Abort("PetscHelmholtzNeumannEBBC::getBoundaryData -- the normal derivative has not been set");
  }

  return DphiDn;
}

#include <CD_NamespaceFooter.H>

#endif
