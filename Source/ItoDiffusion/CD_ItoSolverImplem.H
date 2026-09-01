/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ItoSolverImplem.H
 * @brief  Implementation of CD_ItoSolver.H
 * @author Robert Marskar
 * @date   April 2020
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
#include <CD_MirrorDeposition.H>
#include <CD_ParallelOps.H>
#include <CD_ParticleOps.H>
#include <CD_Random.H>
#include <CD_EBParticleMesh.H>
#include <CD_NamespaceHeader.H>

inline RealVect
ItoSolver::randomGaussian() const
{
  // TLDR: We draw a random number from a Gaussian distribution for each coordinate, and truncate the distribution at
  // m_normalDistributionTruncation.

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

template <typename P, typename Traits>
void
ItoSolver::depositWeight(EBAMRCellData&                      a_phi,
                         const ParticleContainer<P, Traits>& a_particles,
                         const DepositionType                a_deposition,
                         const CoarseFineDeposition          a_coarseFineDeposition) const
{
  CH_TIME("ItoSolver::depositWeight");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositWeight" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(a_phi.getRealm() == m_realm);
  CH_assert(a_particles.getRealm() == m_realm);
  CH_assert(!a_particles.isOrganizedByCell());

  // Deposit the weight column onto the mesh (resets a_phi internally), then redistribute. Redistribution belongs
  // here because it decides what the valid-cell values are; coarsening and ghost filling are the caller's, see
  // coarsenAndFillGhosts().
  m_amr
    ->depositWeight(a_phi, m_realm, m_phase, a_particles, a_deposition, a_coarseFineDeposition, m_irregularDeposition);

  if (m_irregularDeposition == IrregularDeposition::Mirror) {
    this->mirrorPass(a_phi,
                     a_particles,
                     a_deposition,
                     a_coarseFineDeposition,
                     [](const ParticleSoA<P, Traits>& a_leaf, const std::size_t a_index) -> Real {
                       return a_leaf.weight(a_index);
                     });
  }

  this->redistributeAMR(a_phi);
}

template <typename Gather>
void
ItoSolver::depositGatheredNGP(LevelData<EBCellFAB>&                 a_output,
                              const ParticleContainer<ItoParticle>& a_particles,
                              const int                             a_level,
                              Gather                                a_gather) const noexcept
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

    EBParticleMesh particleMesh(domain, cellBox, ebisbox, dx * RealVect::Unit, probLo);

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
ItoSolver::depositGathered(EBAMRCellData&                        a_phi,
                           const ParticleContainer<ItoParticle>& a_particles,
                           const DepositionType                  a_deposition,
                           const CoarseFineDeposition            a_coarseFineDeposition,
                           Gather                                a_gather) const
{
  CH_TIME("ItoSolver::depositGathered");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositGathered" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isOrganizedByCell());

  // Deposit the gathered quantity onto the mesh (resets a_phi internally), then redistribute. Coarsening and ghost
  // filling are the caller's, see coarsenAndFillGhosts().
  m_amr->depositGathered(a_phi,
                         m_realm,
                         m_phase,
                         a_particles,
                         a_deposition,
                         a_coarseFineDeposition,
                         m_irregularDeposition,
                         a_gather);

  if (m_irregularDeposition == IrregularDeposition::Mirror) {
    this->mirrorPass(a_phi, a_particles, a_deposition, a_coarseFineDeposition, a_gather);
  }

  // A no-op unless the selector is one of the two redistributing ones, which Mirror is not -- the enum is what makes
  // that combination unrepresentable rather than merely discouraged.
  this->redistributeAMR(a_phi);
}

template <typename P, typename Traits, typename Strength>
void
ItoSolver::mirrorPass(EBAMRCellData&                      a_phi,
                      const ParticleContainer<P, Traits>& a_particles,
                      const DepositionType                a_deposition,
                      const CoarseFineDeposition          a_coarseFineDeposition,
                      Strength                            a_strength) const
{
  CH_TIME("ItoSolver::mirrorPass");
  if (m_verbosity > 5) {
    pout() << m_name + "::mirrorPass" << endl;
  }

  const Vector<int>& hasCutCells = m_amr->getMirrorHasCutCells(m_realm, m_phase);

  // The early-out has to be a GLOBAL decision, not a per-rank one: remap() below is collective, so a rank that
  // skipped it while another entered it deadlocks rather than misbehaves. getMirrorHasCutCells is all-reduced at
  // regrid precisely so this test is safe to make independently on every rank.
  bool anyCutCells = false;
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    anyCutCells = anyCutCells || (hasCutCells[lvl] != 0);
  }

  if (!anyCutCells) {
    return;
  }

  const EBAMRCellData& surfaceData = m_amr->getMirrorSurfaceData(m_realm, m_phase);
  const RealVect       probLo      = m_amr->getProbLo();
  const Real           maxJacobian = m_amr->getMirrorMaxJacobian();
  const Real           minDenom    = m_amr->getMirrorMinDenominator();

  // Cleared at the START of the pass, so that a pass which aborts partway cannot leak its images into the next
  // deposit.
  m_mirrorImages.clearParticles();

  // Diagnostics. Counted unconditionally because it is a handful of increments, but reduced and printed only when
  // asked for: the reduction is a collective and this pass runs on every deposit of every species, so doing it
  // always would put a global synchronization in the hot path. m_verbosity is parsed from the same input on every
  // rank, so gating a collective on it is safe.
  //
  // These matter more than they look. A refused image still deposits, with J = 1, which is the FLAT mirror -- and a
  // flat mirror on a curved surface is biased, not neutral. So a run whose images are nearly all refused has quietly
  // degraded to a worse model rather than failed, and without these counters it looks identical to a healthy one.
  long long numImages       = 0;
  long long numNoData       = 0;
  long long numInside       = 0;
  long long numSmallDenom   = 0;
  long long numLargeJac     = 0;
  long long numNonPositiveJ = 0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const Real               dx  = m_amr->getDx()[lvl];

    const int nbox = dit.size();

    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const EBCellFAB&     surf    = (*surfaceData[lvl])[din];
      const BaseFab<Real>& surfReg = surf.getSingleValuedFAB();

      const ParticleSoA<P, Traits>& leaf   = a_particles[lvl][din];
      ParticleSoA<NoPayload>&       images = m_mirrorImages[lvl][din];

      const std::size_t numParticles = leaf.size();

      for (std::size_t i = 0; i < numParticles; i++) {
        const RealVect pos = leaf.position(i);
        const IntVect  iv  = ParticleOps::getParticleCellIndex(pos, probLo, dx);

        // Outside the band the surface data was never written, and a particle there has nothing to reflect across.
        if (surfReg(iv, MirrorDeposition::compStatus) < 0.5) {
          numNoData++;

          continue;
        }

        // reflect() takes a contiguous array. A BaseFab stores the component index LAST, so &surfReg(iv, 0) would
        // walk into the next cell rather than the next component.
        Real surfaceComponents[MirrorDeposition::numComp];
        for (int c = 0; c < MirrorDeposition::numComp; c++) {
          surfaceComponents[c] = surfReg(iv, c);
        }

        RealVect                  image;
        Real                      jacobian = 1.0;
        Real                      d        = 0.0;
        MirrorDeposition::Refusal refusal  = MirrorDeposition::Refusal::None;

        // A refused image is still deposited, with J = 1. Dropping it would remove the correction in exactly the
        // cells where the correction is largest.
        MirrorDeposition::reflect(pos, surfaceComponents, maxJacobian, minDenom, image, jacobian, d, refusal);

        switch (refusal) {
        case MirrorDeposition::Refusal::SmallDenominator: {
          numSmallDenom++;

          break;
        }
        case MirrorDeposition::Refusal::LargeJacobian: {
          numLargeJac++;

          break;
        }
        case MirrorDeposition::Refusal::NonPositiveJacobian: {
          numNonPositiveJ++;

          break;
        }
        case MirrorDeposition::Refusal::None: {
          break;
        }
        }

        // A particle at or inside the surface has no image to add: the even extension is defined on the fluid side.
        if (d <= 0.0) {
          numInside++;

          continue;
        }

        numImages++;

        images.append(image, a_strength(leaf, i) * jacobian);
      }
    }
  }

  // The ORDINARY remap. A level-preserving one was evaluated and refuted in both directions -- pinning an image to
  // the coarse level puts it outside the coarse-fine band where it would simply be discarded.
  m_mirrorImages.remap();

  // Native, because the images are the cut-cell treatment; asking for another one here would apply it twice.
  m_amr->depositWeight(m_mirrorScratch,
                       m_realm,
                       m_phase,
                       m_mirrorImages,
                       a_deposition,
                       a_coarseFineDeposition,
                       IrregularDeposition::Native);

  // No covered-cell reset. Every image is inside the solid by construction -- that is what reflection means -- and
  // the density is only defined where kappa > 0. The deposit's contract is that the field is correct where it is
  // read.
  DataOps::incr(a_phi, m_mirrorScratch, 1.0);

  if (m_verbosity > 5) {
    numImages       = ParallelOps::sum(numImages);
    numNoData       = ParallelOps::sum(numNoData);
    numInside       = ParallelOps::sum(numInside);
    numSmallDenom   = ParallelOps::sum(numSmallDenom);
    numLargeJac     = ParallelOps::sum(numLargeJac);
    numNonPositiveJ = ParallelOps::sum(numNonPositiveJ);

    pout() << m_name + "::mirrorPass - " << numImages << " images; skipped " << numNoData << " outside the band, "
           << numInside << " at d <= 0; refused (J = 1, the FLAT mirror) " << numSmallDenom << " small denominator, "
           << numLargeJac << " large Jacobian, " << numNonPositiveJ << " non-positive Jacobian" << endl;
  }
}

#include <CD_NamespaceFooter.H>

#endif
