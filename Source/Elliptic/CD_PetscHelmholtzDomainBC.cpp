/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzDomainBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzDomainBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzDomainBC::PetscHelmholtzDomainBC() : m_isDefined(false), m_dx(-1.0)
{
  CH_TIME("PetscHelmholtzDomainBC::PetscHelmholtzDomainBC");
}

PetscHelmholtzDomainBC::~PetscHelmholtzDomainBC()
{
  CH_TIME("PetscHelmholtzDomainBC::~PetscHelmholtzDomainBC");
}

void
PetscHelmholtzDomainBC::define(const Location::Cell a_dataLocation,
                               const EBLevelGrid&   a_eblg,
                               const RealVect&      a_probLo,
                               const Real           a_dx)
{
  CH_TIME("PetscHelmholtzDomainBC::define");

  CH_assert(a_dx > 0.0);

  if (a_dataLocation != Location::Cell::Center) {
    MayDay::Abort("PetscHelmholtzDomainBC::define -- only cell-centered data is supported");
  }

  m_dataLocation = a_dataLocation;
  m_eblg         = a_eblg;
  m_probLo       = a_probLo;
  m_dx           = a_dx;
  m_isDefined    = true;
}

RealVect
PetscHelmholtzDomainBC::getBoundaryPosition(const IntVect&       a_iv,
                                            const int            a_dir,
                                            const Side::LoHiSide a_side) const noexcept
{
  CH_assert(m_isDefined);

  return m_probLo + (0.5 * RealVect::Unit + RealVect(a_iv)) * m_dx +
         0.5 * m_dx * RealVect(BASISV(a_dir)) * sign(a_side);
}

#include <CD_NamespaceFooter.H>

#endif
