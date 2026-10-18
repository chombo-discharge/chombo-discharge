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
                                             const IntVect&       a_iv,
                                             const int            a_dir,
                                             const Side::LoHiSide a_side,
                                             const int            a_level) const
{
  if (m_useConstant) {
    a_A = m_constantA;
    a_B = m_constantB;
    a_C = m_constantC;
  }
  else if (m_useFunction) {
    const RealVect pos = this->getBoundaryPosition(a_iv, a_dir, a_side, a_level);

    a_A = m_functionA(pos);
    a_B = m_functionB(pos);
    a_C = m_functionC(pos);
  }
  else {
    MayDay::Abort("PetscHelmholtzRobinDomainBC::getCoefficients -- coefficients not set");
  }
}

std::pair<VoFStencil, VoFStencil>
PetscHelmholtzRobinDomainBC::getFaceFluxStencil(const VolIndex&      a_vof,
                                                const int            a_dir,
                                                const Side::LoHiSide a_side,
                                                const int            a_level,
                                                const DataIndex&     a_din) const
{
  CH_assert(m_isDefined);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  Real A;
  Real B;
  Real C;

  this->getCoefficients(A, B, C, a_vof.gridIndex(), a_dir, a_side, a_level);

  const Real    factor  = sign(a_side) * A / B;
  const IntVect ivNear  = a_vof.gridIndex() - sign(a_side) * BASISV(a_dir);
  const bool    hasNear = m_eblg[a_level].getDomain().contains(ivNear) &&
                       !m_eblg[a_level].getEBISL()[a_din].isCovered(ivNear);

  // phi on the face is extrapolated linearly from the cell and its interior neighbor, or set to the cell value if
  // there is no neighbor.
  VoFStencil stencil;

  if (hasNear) {
    stencil.add(a_vof, 1.5 * factor);
    stencil.add(VolIndex(ivNear, 0), -0.5 * factor);
  }
  else {
    stencil.add(a_vof, factor);
  }

  return std::make_pair(stencil, VoFStencil());
}

Real
PetscHelmholtzRobinDomainBC::getFaceFluxConstant(const VolIndex&      a_vof,
                                                 const int            a_dir,
                                                 const Side::LoHiSide a_side,
                                                 const int            a_level,
                                                 const DataIndex& /*a_din*/) const
{
  CH_assert(m_isDefined);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  Real A;
  Real B;
  Real C;

  this->getCoefficients(A, B, C, a_vof.gridIndex(), a_dir, a_side, a_level);

  return C / B;
}

#include <CD_NamespaceFooter.H>

#endif
