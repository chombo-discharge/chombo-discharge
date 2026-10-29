/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzRobinDomainBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzRobinDomainBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzRobinDomainBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzRobinDomainBC::PetscHelmholtzRobinDomainBC() : m_useConstant(false), m_useFunction(false)
{
  CH_TIME("PetscHelmholtzRobinDomainBC::PetscHelmholtzRobinDomainBC");
}

PetscHelmholtzRobinDomainBC::~PetscHelmholtzRobinDomainBC()
{
  CH_TIME("PetscHelmholtzRobinDomainBC::~PetscHelmholtzRobinDomainBC");
}

void
PetscHelmholtzRobinDomainBC::setCoefficients(const Real a_A, const Real a_B, const Real a_C)
{
  CH_TIME("PetscHelmholtzRobinDomainBC::setCoefficients(Real, Real, Real)");

  m_useConstant = true;
  m_useFunction = false;
  m_constantA   = a_A;
  m_constantB   = a_B;
  m_constantC   = a_C;
}

void
PetscHelmholtzRobinDomainBC::setCoefficients(const std::function<Real(const RealVect& a_pos)>& a_A,
                                             const std::function<Real(const RealVect& a_pos)>& a_B,
                                             const std::function<Real(const RealVect& a_pos)>& a_C)
{
  CH_TIME("PetscHelmholtzRobinDomainBC::setCoefficients(3x std::function)");

  m_useConstant = false;
  m_useFunction = true;
  m_functionA   = a_A;
  m_functionB   = a_B;
  m_functionC   = a_C;
}

void
PetscHelmholtzRobinDomainBC::getCoefficients(Real&                a_A,
                                             Real&                a_B,
                                             Real&                a_C,
                                             const int            a_phase,
                                             const int            a_level,
                                             const DataIndex&     a_din,
                                             const VolIndex&      a_vof,
                                             const int            a_dir,
                                             const Side::LoHiSide a_side) const
{
  if (m_useConstant) {
    a_A = m_constantA;
    a_B = m_constantB;
    a_C = m_constantC;
  }
  else if (m_useFunction) {
    const RealVect pos = this->getBoundaryPosition(a_phase, a_level, a_din, a_vof, a_dir, a_side);

    a_A = m_functionA(pos);
    a_B = m_functionB(pos);
    a_C = m_functionC(pos);
  }
  else {
    MayDay::Abort("PetscHelmholtzRobinDomainBC::getCoefficients -- coefficients not set");
  }
}

void
PetscHelmholtzRobinDomainBC::getFaceFluxStencil(Real&                a_boundaryWeight,
                                                VoFStencil&          a_stencilCoar,
                                                VoFStencil&          a_stencil,
                                                VoFStencil&          a_stencilFine,
                                                const int            a_phase,
                                                const int            a_level,
                                                const DataIndex&     a_din,
                                                const VolIndex&      a_vof,
                                                const int            a_dir,
                                                const Side::LoHiSide a_side) const
{
  CH_assert(m_isDefined);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  Real A;
  Real B;
  Real C;

  this->getCoefficients(A, B, C, a_phase, a_level, a_din, a_vof, a_dir, a_side);

  const Vector<RefCountedPtr<MFLevelGrid>>& mflg    = m_petscGrid->getMFLevelGrids();
  const EBLevelGrid&                        eblg    = mflg[a_level]->getEBLevelGrid(a_phase);
  const EBISLayout&                         ebisl   = eblg.getEBISL();
  const EBISBox&                            ebisbox = ebisl[a_din];
  const BaseFab<PetscAMRCell>&              cells   = (*m_petscGrid->getAMRToPetsc()[a_level])[a_din];

  const Real    factor = sign(a_side) * A / B;
  const IntVect ivNear = a_vof.gridIndex() - sign(a_side) * BASISV(a_dir);

  // The interior neighbor must be a cell on this level: inside the domain, not covered by the EB, and not covered by a
  // finer grid.
  const bool hasNear = eblg.getDomain().contains(ivNear) && !ebisbox.isCovered(ivNear) &&
                       !cells(ivNear).isCoveredByFinerGrid();

  a_stencilCoar.clear();
  a_stencil.clear();
  a_stencilFine.clear();

  // phi on the face is extrapolated linearly from the cell and its interior neighbor. Without the neighbor, the
  // extrapolation needs a least-squares reconstruction over the valid cells on this level and the next finer level.
  if (hasNear) {
    a_stencil.add(a_vof, 1.5 * factor);
    a_stencil.add(VolIndex(ivNear, 0), -0.5 * factor);
  }
  else {
    MayDay::Abort("PetscHelmholtzRobinDomainBC::getFaceFluxStencil -- interior neighbor is not a cell on this level, "
                  "and the least-squares reconstruction is not implemented");
  }

  a_boundaryWeight = -sign(a_side) / B;
}

Real
PetscHelmholtzRobinDomainBC::getBoundaryData(const int            a_phase,
                                             const int            a_level,
                                             const DataIndex&     a_din,
                                             const VolIndex&      a_vof,
                                             const int            a_dir,
                                             const Side::LoHiSide a_side) const
{
  CH_assert(m_isDefined);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  Real A;
  Real B;
  Real C;

  this->getCoefficients(A, B, C, a_phase, a_level, a_din, a_vof, a_dir, a_side);

  return C;
}

#include <CD_NamespaceFooter.H>

#endif
