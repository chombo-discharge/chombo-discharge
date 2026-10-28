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
PetscHelmholtzDirichletEBBC::getEBFluxStencil(Real&            a_boundaryWeight,
                                              PetscStencil&    a_stencil,
                                              const int        a_phase,
                                              const int        a_level,
                                              const DataIndex& a_din,
                                              const VolIndex&  a_vof) const noexcept
{
  CH_TIME("PetscHelmholtzDirichletEBBC::getEBFluxStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzDirichletEBBC::getEBFluxStencil" << endl;
  }

  CH_assert(m_isDefined);

  // The flux is the one-sided normal derivative itself, with g = phi_B.
  this->getNormalDerivativeStencil(a_boundaryWeight, a_stencil, a_phase, a_level, a_din, a_vof);
}

Real
PetscHelmholtzDirichletEBBC::getBoundaryData(const int        a_phase,
                                             const int        a_level,
                                             const DataIndex& a_din,
                                             const VolIndex&  a_vof) const noexcept
{
  CH_TIME("PetscHelmholtzDirichletEBBC::getBoundaryData");
  if (m_verbose) {
    pout() << "PetscHelmholtzDirichletEBBC::getBoundaryData" << endl;
  }

  CH_assert(m_isDefined);

  Real value = 0.0;

  if (m_useConstant) {
    value = m_constantValue;
  }
  else if (m_useFunction) {
    value = m_functionValue(this->getBoundaryPosition(a_phase, a_level, a_din, a_vof));
  }
  else {
    // Checked in every build: without a value the right-hand side would silently use phi_B = 0.
    MayDay::Abort("PetscHelmholtzDirichletEBBC::getBoundaryData -- the value has not been set");
  }

  return value;
}

#include <CD_NamespaceFooter.H>

#endif
