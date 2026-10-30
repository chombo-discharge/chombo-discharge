/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleMassRedistribution.cpp
 * @brief  Implementation of CD_ParticleMassRedistribution.H
 * @author Robert Marskar
 */

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_ParticleMassRedistribution.H>
#include <CD_NamespaceHeader.H>

ParticleMassRedistribution::ParticleMassRedistribution(
  const Vector<RefCountedPtr<EBFluxRedistribution>>& a_fluxRedistribution) noexcept
{
  CH_TIME("ParticleMassRedistribution::ParticleMassRedistribution");

  m_fluxRedistribution = a_fluxRedistribution;
}

ParticleMassRedistribution::~ParticleMassRedistribution() noexcept
{
  CH_TIME("ParticleMassRedistribution::~ParticleMassRedistribution");
}

void
ParticleMassRedistribution::redistribute(EBAMRCellData& a_phi, const EBAMRIVData& a_deltaM) const noexcept
{
  CH_TIME("ParticleMassRedistribution::redistribute");

  const int finestLevel = a_phi.size() - 1;

  CH_assert(a_deltaM.size() == a_phi.size());
  CH_assert(m_fluxRedistribution.size() > finestLevel);

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    const Real     scale     = 1.0;
    const Interval variables = Interval(0, 0);
    const bool     hasCoar   = lvl > 0;
    const bool     hasFine   = lvl < finestLevel;

    CH_assert(a_phi[lvl]->nComp() == 1);
    CH_assert(a_deltaM[lvl]->nComp() == 1);

    if (hasCoar) {
      m_fluxRedistribution[lvl]->redistributeCoar(*a_phi[lvl - 1], *a_deltaM[lvl], scale, variables);
    }

    m_fluxRedistribution[lvl]->redistributeLevel(*a_phi[lvl], *a_deltaM[lvl], scale, variables);

    if (hasFine) {
      m_fluxRedistribution[lvl]->redistributeFine(*a_phi[lvl + 1], *a_deltaM[lvl], scale, variables);
    }
  }
}

#include <CD_NamespaceFooter.H>
