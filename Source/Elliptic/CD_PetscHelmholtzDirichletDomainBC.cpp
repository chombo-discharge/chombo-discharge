/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzDirichletDomainBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzDirichletDomainBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzDirichletDomainBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzDirichletDomainBC::PetscHelmholtzDirichletDomainBC() : m_useConstant(false), m_useFunction(false)
{
  CH_TIME("PetscHelmholtzDirichletDomainBC::PetscHelmholtzDirichletDomainBC()");
}

PetscHelmholtzDirichletDomainBC::PetscHelmholtzDirichletDomainBC(const Real a_value)
{
  CH_TIME("PetscHelmholtzDirichletDomainBC::PetscHelmholtzDirichletDomainBC(Real)");

  this->setValue(a_value);
}

PetscHelmholtzDirichletDomainBC::PetscHelmholtzDirichletDomainBC(
  const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("PetscHelmholtzDirichletDomainBC::PetscHelmholtzDirichletDomainBC(std::function)");

  this->setValue(a_value);
}

PetscHelmholtzDirichletDomainBC::~PetscHelmholtzDirichletDomainBC()
{
  CH_TIME("PetscHelmholtzDirichletDomainBC::~PetscHelmholtzDirichletDomainBC");
}

void
PetscHelmholtzDirichletDomainBC::setValue(const Real a_value)
{
  CH_TIME("PetscHelmholtzDirichletDomainBC::setValue(Real)");

  m_useConstant   = true;
  m_useFunction   = false;
  m_constantValue = a_value;
}

void
PetscHelmholtzDirichletDomainBC::setValue(const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("PetscHelmholtzDirichletDomainBC::setValue(std::function)");

  m_useConstant   = false;
  m_useFunction   = true;
  m_functionValue = a_value;
}

std::pair<VoFStencil, VoFStencil>
PetscHelmholtzDirichletDomainBC::getFaceFluxStencil(const VolIndex&      a_vof,
                                                    const int            a_dir,
                                                    const Side::LoHiSide a_side,
                                                    const int            a_level,
                                                    const DataIndex& /*a_din*/) const
{
  CH_assert(m_isDefined);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  // The flux is s * (phi_B - phi)/(dx/2), with s = -1 on the low side and +1 on the high side.
  VoFStencil stencil;
  stencil.add(a_vof, -sign(a_side) * 2.0 / m_dx[a_level]);

  return std::make_pair(stencil, VoFStencil());
}

Real
PetscHelmholtzDirichletDomainBC::getFaceFluxConstant(const VolIndex&      a_vof,
                                                     const int            a_dir,
                                                     const Side::LoHiSide a_side,
                                                     const int            a_level,
                                                     const DataIndex& /*a_din*/) const
{
  CH_assert(m_isDefined);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  Real value = 0.0;

  if (m_useConstant) {
    value = m_constantValue;
  }
  else if (m_useFunction) {
    value = m_functionValue(this->getBoundaryPosition(a_vof.gridIndex(), a_dir, a_side, a_level));
  }
  else {
    MayDay::Abort("PetscHelmholtzDirichletDomainBC::getFaceFluxConstant -- value not set");
  }

  return sign(a_side) * 2.0 * value / m_dx[a_level];
}

#include <CD_NamespaceFooter.H>

#endif
