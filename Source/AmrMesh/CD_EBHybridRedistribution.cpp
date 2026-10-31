/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBHybridRedistribution.cpp
 * @brief  Implementation of CD_EBHybridRedistribution.H
 * @author Robert Marskar
 */

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHybridRedistribution.H>
#include <CD_EBFluxRedistribution.H>
#include <CD_DataOps.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

EBHybridRedistribution::EBHybridRedistribution(const RefCountedPtr<AmrMesh>& a_amr,
                                               const std::string&            a_realm,
                                               const phase::which_phase      a_phase,
                                               const bool                    a_blendConservation,
                                               const bool                    a_redistribute) noexcept
{
  CH_TIME("EBHybridRedistribution::EBHybridRedistribution");

  CH_assert(!a_amr.isNull());

  m_amr               = a_amr;
  m_realm             = a_realm;
  m_phase             = a_phase;
  m_blendConservation = a_blendConservation;
  m_redistribute      = a_redistribute;

  m_amr->allocate(m_nonConservativeDivergence, m_realm, m_phase, 1);
  m_amr->allocate(m_massDifference, m_realm, m_phase, 1);

  DataOps::setValue(m_nonConservativeDivergence, 0.0);
  DataOps::setValue(m_massDifference, 0.0);
}

EBHybridRedistribution::~EBHybridRedistribution() noexcept
{
  CH_TIME("EBHybridRedistribution::~EBHybridRedistribution");
}

void
EBHybridRedistribution::redistribute(EBAMRCellData& a_divG,
                                     const EBAMRCellData& /*a_phi*/,
                                     const Real /*a_dt*/) noexcept
{
  CH_TIME("EBHybridRedistribution::redistribute");

  CH_assert(a_divG[0]->nComp() == 1);

  const int finestLevel = m_amr->getFinestLevel();

  // The non-conservative divergence D_nc = sum(kappa*D)/sum(kappa) over the neighbourhood of each cut cell. Without
  // blending it is zero.
  if (m_blendConservation) {
    m_amr->nonConservativeDivergence(m_nonConservativeDivergence, a_divG, m_realm, m_phase);
  }
  else {
    DataOps::setValue(m_nonConservativeDivergence, 0.0);
  }

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    this->hybridDivergence(*a_divG[lvl], *m_massDifference[lvl], *m_nonConservativeDivergence[lvl], lvl);
  }

  // Redistribute the missed mass to the neighbouring cells on this level and across refinement boundaries.
  if (m_redistribute) {
    Vector<RefCountedPtr<EBFluxRedistribution>>& redistOps = m_amr->getRedistributionOp(m_realm, m_phase);

    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      const Real     scale     = 1.0;
      const Interval variables = Interval(0, 0);
      const bool     hasCoar   = lvl > 0;
      const bool     hasFine   = lvl < finestLevel;

      if (hasCoar) {
        redistOps[lvl]->redistributeCoar(*a_divG[lvl - 1], *m_massDifference[lvl], scale, variables);
      }

      redistOps[lvl]->redistributeLevel(*a_divG[lvl], *m_massDifference[lvl], scale, variables);

      if (hasFine) {
        redistOps[lvl]->redistributeFine(*a_divG[lvl + 1], *m_massDifference[lvl], scale, variables);
      }
    }
  }
}

void
EBHybridRedistribution::hybridDivergence(LevelData<EBCellFAB>&             a_hybridDivergence,
                                         LevelData<BaseIVFAB<Real>>&       a_massDifference,
                                         const LevelData<BaseIVFAB<Real>>& a_nonConservativeDivergence,
                                         const int                         a_level) const noexcept
{
  CH_TIME("EBHybridRedistribution::hybridDivergence");

  CH_assert(a_hybridDivergence.nComp() == 1);
  CH_assert(a_massDifference.nComp() == 1);
  CH_assert(a_nonConservativeDivergence.nComp() == 1);

  // The hybrid divergence is kappa*D + (1-kappa)*D_nc, where a_hybridDivergence holds kappa*D on input.
  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[a_level];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[a_level];
  const DataIterator&      dit   = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    EBCellFAB&             divH   = a_hybridDivergence[din];
    BaseIVFAB<Real>&       deltaM = a_massDifference[din];
    const BaseIVFAB<Real>& divNC  = a_nonConservativeDivergence[din];

    const EBISBox& ebisbox = ebisl[din];

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[din];

    auto kernel = [&](const VolIndex& vof) -> void {
      const Real kappa = ebisbox.volFrac(vof);
      const Real dc    = divH(vof, 0);
      const Real dnc   = divNC(vof, 0);

      // dc is kappa*D, so the missed mass is (1-kappa)*(dc - kappa*dnc).
      divH(vof, 0)   = dc + (1 - kappa) * dnc;
      deltaM(vof, 0) = (1 - kappa) * (dc - kappa * dnc);
    };

    BoxLoops::loop(vofit, kernel);
  }
}

#include <CD_NamespaceFooter.H>
