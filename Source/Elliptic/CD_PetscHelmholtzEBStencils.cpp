/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzEBStencils.cpp
 * @brief  Implementation of CD_PetscHelmholtzEBStencils.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>
#include <ParmParse.H>

// Our includes
#include <CD_PetscHelmholtzEBStencils.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzEBStencils::PetscHelmholtzEBStencils() noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::PetscHelmholtzEBStencils");

  m_isDefined = false;
  m_verbose   = false;
  m_order     = -1;
  m_weight    = -1;

  ParmParse pp("PetscHelmholtzEBStencils");
  pp.query("verbose", m_verbose);
}

PetscHelmholtzEBStencils::~PetscHelmholtzEBStencils() noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::~PetscHelmholtzEBStencils");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::~PetscHelmholtzEBStencils" << endl;
  }
}

void
PetscHelmholtzEBStencils::setOrder(const int a_order) noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::setOrder");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::setOrder" << endl;
  }

  m_order = a_order;
}

void
PetscHelmholtzEBStencils::setWeight(const int a_weight) noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::setWeight");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::setWeight" << endl;
  }

  m_weight = a_weight;
}

void
PetscHelmholtzEBStencils::clear() noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::clear");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::clear" << endl;
  }

  for (int iphase = 0; iphase < 2; iphase++) {
    m_derivStencils[iphase].resize(0);
    m_derivWeights[iphase].resize(0);
  }

  m_petscGrid = RefCountedPtr<PetscGrid>();
  m_probLo    = RealVect::Zero;
  m_isDefined = false;

  m_dx.resize(0);
}

void
PetscHelmholtzEBStencils::define(const RefCountedPtr<PetscGrid>& /*a_petscGrid*/,
                                 const RealVect& /*a_probLo*/,
                                 const Vector<Real>& /*a_dx*/) noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::define");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::define" << endl;
  }

  this->clear();

  MayDay::Abort("PetscHelmholtzEBStencils::define -- not implemented yet");
}

void
PetscHelmholtzEBStencils::getNormalDerivativeStencil(Real& /*a_boundaryWeight*/,
                                                     PetscStencil& /*a_stencil*/,
                                                     const int /*a_phase*/,
                                                     const int /*a_level*/,
                                                     const DataIndex& /*a_din*/,
                                                     const VolIndex& /*a_vof*/) const noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::getNormalDerivativeStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::getNormalDerivativeStencil" << endl;
  }

  MayDay::Abort("PetscHelmholtzEBStencils::getNormalDerivativeStencil -- not implemented yet");
}

void
PetscHelmholtzEBStencils::computeNormalDerivativeStencil(Real& /*a_boundaryWeight*/,
                                                         PetscStencil& /*a_stencil*/,
                                                         const int /*a_phase*/,
                                                         const int /*a_level*/,
                                                         const DataIndex& /*a_din*/,
                                                         const VolIndex& /*a_vof*/) const noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::computeNormalDerivativeStencil");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::computeNormalDerivativeStencil" << endl;
  }

  MayDay::Abort("PetscHelmholtzEBStencils::computeNormalDerivativeStencil -- not implemented yet");
}

RealVect
PetscHelmholtzEBStencils::getBoundaryPosition(const int /*a_phase*/,
                                              const int /*a_level*/,
                                              const DataIndex& /*a_din*/,
                                              const VolIndex& /*a_vof*/) const noexcept
{
  CH_TIME("PetscHelmholtzEBStencils::getBoundaryPosition");
  if (m_verbose) {
    pout() << "PetscHelmholtzEBStencils::getBoundaryPosition" << endl;
  }

  MayDay::Abort("PetscHelmholtzEBStencils::getBoundaryPosition -- not implemented yet");

  return RealVect::Zero;
}

#include <CD_NamespaceFooter.H>

#endif
