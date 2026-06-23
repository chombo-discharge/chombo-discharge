/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoSolverImplem.H
  @brief  Implementation of CD_ItoSolver.H
  @author Robert Marskar
  @date   April 2020
*/

#ifndef CD_ITOSOLVERIMPLEM_H
#define CD_ITOSOLVERIMPLEM_H

// Std includes
#include <chrono>

// Chombo includes
#include <EBAlias.H>
#include <PolyGeom.H>

// Our includes
#include <CD_ItoSolver.H>
#include <CD_DataOps.H>
#include <CD_Random.H>
#include <CD_EBParticleMeshSoA.H>
#include <CD_NamespaceHeader.H>

inline RealVect
ItoSolver::randomGaussian() const
{
  // TLDR: We draw a random number from a Gaussian distribution for each coordinate, and truncate the distribution at m_normalDistributionTruncation.

  auto sign = [](const Real& a) -> Real {
    return (a > 0) - (a < 0);
  };

  RealVect r = RealVect::Zero;
  for (int i = 0; i < SpaceDim; i++) {
    r[i] = Random::getNormal01();

    r[i] = sign(r[i]) * std::min(std::abs(r[i]), m_normalDistributionTruncation);
  }

  return r;
}

template <typename Gather>
void
ItoSolver::depositGatheredNGP(LevelData<EBCellFAB>&                    a_output,
                              const ParticleContainerSoA<ItoParticle>& a_particles,
                              const int                                a_level,
                              Gather                                   a_gather) const noexcept
{
  CH_TIME("ItoSolver::depositGatheredNGP");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositGatheredNGP" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  const ProblemDomain&     domain = m_amr->getDomains()[a_level];
  const DisjointBoxLayout& dbl    = m_amr->getGrids(a_particles.getRealm())[a_level];
  const DataIterator&      dit    = dbl.dataIterator();
  const EBISLayout&        ebisl  = m_amr->getEBISLayout(a_particles.getRealm(), m_phase)[a_level];
  const Real               dx     = m_amr->getDx()[a_level];
  const RealVect           probLo = m_amr->getProbLo();

  CH_assert(a_output.disjointBoxLayout() == dbl);

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box      cellBox = dbl[din];
    const EBISBox& ebisbox = ebisl[din];

    EBParticleMeshSoA particleMesh(domain, cellBox, ebisbox, dx * RealVect::Unit, probLo);

    EBCellFAB&                      output = a_output[din];
    const ParticleSoA<ItoParticle>& leaf   = a_particles[a_level][din];

    // The per-patch deposit INCREMENTS, so start from a clean slate.
    output.setVal(0.0);

    particleMesh
      .depositGathered(output, 0, leaf, DepositionType::NGP, 1.0, true, [&](const std::size_t a_i, Real* a_out) {
        a_out[0] = a_gather(leaf, a_i);
      });
  }
}

template <typename Gather>
void
ItoSolver::depositGathered(EBAMRCellData&                           a_phi,
                           const ParticleContainerSoA<ItoParticle>& a_particles,
                           const DepositionType                     a_deposition,
                           const CoarseFineDeposition               a_coarseFineDeposition,
                           Gather                                   a_gather) const
{
  CH_TIME("ItoSolver::depositGathered");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositGathered" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isOrganizedByCell());

  // Deposit the gathered quantity onto the mesh (resets a_phi internally), then run the shared
  // redistribution + average-down + ghost-interpolation tail.
  m_amr->depositGathered(a_phi,
                         m_realm,
                         m_phase,
                         a_particles,
                         a_deposition,
                         a_coarseFineDeposition,
                         m_forceIrregDepositionNGP,
                         a_gather);

  this->finalizeDeposit(a_phi);
}

#include <CD_NamespaceFooter.H>

#endif
