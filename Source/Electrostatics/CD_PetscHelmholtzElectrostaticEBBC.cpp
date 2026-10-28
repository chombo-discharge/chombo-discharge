/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PetscHelmholtzElectrostaticEBBC.cpp
 * @brief  Implementation of CD_PetscHelmholtzElectrostaticEBBC.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Std includes
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_PetscHelmholtzElectrostaticEBBC.H>
#include <CD_NamespaceHeader.H>

PetscHelmholtzElectrostaticEBBC::PetscHelmholtzElectrostaticEBBC(const ElectrostaticEbBc& a_electrostaticBCs) noexcept
  : m_electrostaticBCs(a_electrostaticBCs)
{
  CH_TIME("PetscHelmholtzElectrostaticEBBC::PetscHelmholtzElectrostaticEBBC");
}

PetscHelmholtzElectrostaticEBBC::~PetscHelmholtzElectrostaticEBBC() noexcept
{
  CH_TIME("PetscHelmholtzElectrostaticEBBC::~PetscHelmholtzElectrostaticEBBC");
  if (m_verbose) {
    pout() << "PetscHelmholtzElectrostaticEBBC::~PetscHelmholtzElectrostaticEBBC" << endl;
  }
}

Real
PetscHelmholtzElectrostaticEBBC::getBoundaryData(const int        a_phase,
                                                 const int        a_level,
                                                 const DataIndex& a_din,
                                                 const VolIndex&  a_vof) const noexcept
{
  CH_TIME("PetscHelmholtzElectrostaticEBBC::getBoundaryData");
  if (m_verbose) {
    pout() << "PetscHelmholtzElectrostaticEBBC::getBoundaryData" << endl;
  }

  CH_assert(m_isDefined);

  return this->getElectrodePotential(this->getBoundaryPosition(a_phase, a_level, a_din, a_vof));
}

Real
PetscHelmholtzElectrostaticEBBC::getElectrodePotential(const RealVect& a_pos) const noexcept
{
  CH_TIME("PetscHelmholtzElectrostaticEBBC::getElectrodePotential");
  if (m_verbose) {
    pout() << "PetscHelmholtzElectrostaticEBBC::getElectrodePotential" << endl;
  }

  const std::vector<std::pair<Electrode, ElectrostaticEbBc::BcFunction>>& electrodeBCs = m_electrostaticBCs.getBcs();

  // Checked in every build: an electrode face without electrodes has no potential.
  if (electrodeBCs.empty()) {
    MayDay::Abort("PetscHelmholtzElectrostaticEBBC::getElectrodePotential -- there are no electrodes");
  }

  // The closest electrode is the one with the smallest absolute level-set value.
  int  closestElectrode = 0;
  Real minDist          = std::numeric_limits<Real>::infinity();

  for (int i = 0; i < static_cast<int>(electrodeBCs.size()); i++) {
    const Real curDist = std::abs(electrodeBCs[i].first.getImplicitFunction()->value(a_pos));

    if (curDist < minDist) {
      closestElectrode = i;
      minDist          = curDist;
    }
  }

  // The potentials capture the time by reference and ignore their time argument, as in MFHelmholtzElectrostaticEBBC.
  constexpr Real dummyDt = 0.0;

  return electrodeBCs[closestElectrode].second(a_pos, dummyDt);
}

#include <CD_NamespaceFooter.H>

#endif
