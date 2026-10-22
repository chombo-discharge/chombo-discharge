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

// Our includes
#include <CD_PetscHelmholtzDomainBC.H>
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzDomainBC::PetscHelmholtzDomainBC() : m_isDefined(false)
{
  CH_TIME("PetscHelmholtzDomainBC::PetscHelmholtzDomainBC");
}

PetscHelmholtzDomainBC::~PetscHelmholtzDomainBC()
{
  CH_TIME("PetscHelmholtzDomainBC::~PetscHelmholtzDomainBC");
}

void
PetscHelmholtzDomainBC::define(const RefCountedPtr<PetscGrid>& a_petscGrid,
                               const RealVect&                 a_probLo,
                               const Vector<Real>&             a_dx)
{
  CH_TIME("PetscHelmholtzDomainBC::define");

  CH_assert(!a_petscGrid.isNull());
  CH_assert(a_petscGrid->isDefined());
  CH_assert(a_dx.size() >= a_petscGrid->getMFLevelGrids().size());

  m_petscGrid = a_petscGrid;
  m_probLo    = a_probLo;
  m_dx        = a_dx;
  m_isDefined = true;
}

RealVect
PetscHelmholtzDomainBC::getBoundaryPosition(const int            a_phase,
                                            const int            a_level,
                                            const DataIndex&     a_din,
                                            const VolIndex&      a_vof,
                                            const int            a_dir,
                                            const Side::LoHiSide a_side) const noexcept
{
  CH_assert(m_isDefined);
  CH_assert(a_level >= 0 && a_level < static_cast<int>(m_dx.size()));

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg    = m_petscGrid->getMFLevelGrids();
  const EBLevelGrid&                        eblg    = mflg[a_level]->getEBLevelGrid(a_phase);
  const EBISLayout&                         ebisl   = eblg.getEBISL();
  const EBISBox&                            ebisbox = ebisl[a_din];

  const Vector<FaceIndex> faces = ebisbox.getFaces(a_vof, a_dir, a_side);

  CH_assert(faces.size() == 1);

  return m_probLo + Location::position(Location::Face::Centroid, faces[0], ebisbox, m_dx[a_level]);
}

#include <CD_NamespaceFooter.H>

#endif
