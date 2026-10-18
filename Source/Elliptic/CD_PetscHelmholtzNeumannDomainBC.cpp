/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzNeumannDomainBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzNeumannDomainBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzNeumannDomainBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzNeumannDomainBC::PetscHelmholtzNeumannDomainBC() : m_useConstant(false), m_useFunction(false)
{
  CH_TIME("PetscHelmholtzNeumannDomainBC::PetscHelmholtzNeumannDomainBC()");
}

PetscHelmholtzNeumannDomainBC::PetscHelmholtzNeumannDomainBC(const Real a_DphiDn)
{
  CH_TIME("PetscHelmholtzNeumannDomainBC::PetscHelmholtzNeumannDomainBC(Real)");

  this->setDphiDn(a_DphiDn);
}

PetscHelmholtzNeumannDomainBC::PetscHelmholtzNeumannDomainBC(const std::function<Real(const RealVect& a_pos)>& a_DphiDn)
{
  CH_TIME("PetscHelmholtzNeumannDomainBC::PetscHelmholtzNeumannDomainBC(std::function)");

  this->setDphiDn(a_DphiDn);
}

PetscHelmholtzNeumannDomainBC::~PetscHelmholtzNeumannDomainBC()
{
  CH_TIME("PetscHelmholtzNeumannDomainBC::~PetscHelmholtzNeumannDomainBC");
}

void
PetscHelmholtzNeumannDomainBC::setDphiDn(const Real a_DphiDn)
{
  CH_TIME("PetscHelmholtzNeumannDomainBC::setDphiDn(Real)");

  m_useConstant    = true;
  m_useFunction    = false;
  m_constantDphiDn = a_DphiDn;
}

void
PetscHelmholtzNeumannDomainBC::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn)
{
  CH_TIME("PetscHelmholtzNeumannDomainBC::setDphiDn(std::function)");

  m_useConstant    = false;
  m_useFunction    = true;
  m_functionDphiDn = a_DphiDn;
}

std::pair<VoFStencil, VoFStencil>
PetscHelmholtzNeumannDomainBC::getFaceFluxStencil(const VolIndex& /*a_vof*/,
                                                  const int a_dir,
                                                  const Side::LoHiSide /*a_side*/,
                                                  const DataIndex& /*a_din*/) const
{
  CH_assert(m_isDefined);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  return std::make_pair(VoFStencil(), VoFStencil());
}

Real
PetscHelmholtzNeumannDomainBC::getFaceFluxConstant(const VolIndex&      a_vof,
                                                   const int            a_dir,
                                                   const Side::LoHiSide a_side,
                                                   const DataIndex& /*a_din*/) const
{
  CH_assert(m_isDefined);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  Real DphiDn = 0.0;

  if (m_useConstant) {
    DphiDn = m_constantDphiDn;
  }
  else if (m_useFunction) {
    DphiDn = m_functionDphiDn(this->getBoundaryPosition(a_vof.gridIndex(), a_dir, a_side));
  }
  else {
    MayDay::Abort("PetscHelmholtzNeumannDomainBC::getFaceFluxConstant -- dphi/dn not set");
  }

  // A positive dphi/dn is a flux into the domain, as in EBHelmholtzNeumannDomainBC.
  return -sign(a_side) * DphiDn;
}

#include <CD_NamespaceFooter.H>

#endif
