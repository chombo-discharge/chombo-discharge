/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoSolverImplem.H
  @brief  Implementation of CD_ItoSolver.H
  @author Robert Marskar
  @date   April 2020
*/

#ifndef CD_ItoSolverImplem_H
#define CD_ItoSolverImplem_H

// Std includes
#include <chrono>

// Chombo includes
#include <EBAlias.H>
#include <PolyGeom.H>

// Our includes
#include <CD_ItoSolver.H>
#include <CD_DataOps.H>
#include <CD_Random.H>
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

template <ItoSolver::WhichContainer C>
void
ItoSolver::addParticles(ListBox<ItoParticle>& a_inputParticles,
                        const int             a_lvl,
                        const DataIndex       a_dit,
                        const bool            a_destructive)
{
  CH_TIME("ItoSolver::addParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::addParticles" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(C);

  ListBox<ItoParticle>& my_particles = particles[a_lvl][a_dit];
  if (a_destructive) {
    my_particles.addItemsDestructive(a_inputParticles.listItems());
  }
  else {
    my_particles.addItems(a_inputParticles.listItems());
  }
}

template <class P, const Real& (P::*particleScalarField)() const>
void
ItoSolver::depositParticles(EBAMRCellData& a_phi, ParticleContainer<P>& a_particles) const
{
  CH_TIME("ItoSolver::depositParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles" << endl;
  }

  this->depositParticles<P, particleScalarField>(a_phi, a_particles, m_deposition, m_coarseFineDeposition);
}

template <class P, const Real& (P ::*particleScalarField)() const>
void
ItoSolver::depositParticlesNGP(LevelData<EBCellFAB>&       a_output,
                               const ParticleContainer<P>& a_particles,
                               const int                   a_level) const noexcept
{
  CH_TIME("ItoSolver::depositParticles(NGP)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles(NGP)" << endl;
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

    EBCellFAB&     output    = a_output[din];
    const List<P>& particles = a_particles[a_level][din].listItems();

    particleMesh.deposit<P, particleScalarField>(particles, output, DepositionType::NGP, true);
  }
}

template <class P, Real (P ::*particleScalarField)() const>
void
ItoSolver::depositParticlesNGP(LevelData<EBCellFAB>&       a_output,
                               const ParticleContainer<P>& a_particles,
                               const int                   a_level) const noexcept
{
  CH_TIME("ItoSolver::depositParticles(NGP)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles(NGP)" << endl;
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

    EBCellFAB&     output    = a_output[din];
    const List<P>& particles = a_particles[a_level][din].listItems();

    particleMesh.deposit<P, particleScalarField>(particles, output, DepositionType::NGP, true);
  }
}

template <class P, const Real& (P::*particleScalarField)() const>
void
ItoSolver::depositParticles(EBAMRCellData&             a_phi,
                            ParticleContainer<P>&      a_particles,
                            const DepositionType       a_deposition,
                            const CoarseFineDeposition a_coarseFineDeposition) const
{
  CH_TIME("ItoSolver::depositParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isOrganizedByCell());

  // TLDR: First, deposit onto the mesh as usual (as if the EB wasn't there). If the user asks for it, he can redistribute mass in order to
  //       conserve total mass (if that is important). But the corresponding scheme will be O(1) accurate.
  this->depositKappaConservative<P, particleScalarField>(a_phi, a_particles, a_deposition, a_coarseFineDeposition);

  // Redistribution magic.
  this->redistributeAMR(a_phi);

  // Average down and interpolate
  m_amr->conservativeAverage(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);
}

template <class P, Real (P::*particleScalarField)() const>
void
ItoSolver::depositParticles(EBAMRCellData&             a_phi,
                            ParticleContainer<P>&      a_particles,
                            const DepositionType       a_deposition,
                            const CoarseFineDeposition a_coarseFineDeposition) const
{
  CH_TIME("ItoSolver::depositParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isOrganizedByCell());

  // TLDR: First, deposit onto the mesh as usual (as if the EB wasn't there). If the user asks for it, he can redistribute mass in order to
  //       conserve total mass (if that is important). But the corresponding scheme will be O(1) accurate.
  this->depositKappaConservative<P, particleScalarField>(a_phi, a_particles, a_deposition, a_coarseFineDeposition);

  // Redistribution magic.
#if 0
  this->redistributeAMR(a_phi);

  // Average down and interpolate
  m_amr->conservativeAverage(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);
#endif
}

template <class P, const Real& (P::*particleScalarField)() const>
void
ItoSolver::depositKappaConservative(EBAMRCellData&             a_phi,
                                    ParticleContainer<P>&      a_particles,
                                    const DepositionType       a_deposition,
                                    const CoarseFineDeposition a_coarseFineDeposition) const
{
  CH_TIME("ItoSolver::depositKappaConservative");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositKappaConservative" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);

#warning "This needs to go"
  {
    DataOps::setValue(a_phi, 0.0);
  }

  // Now do the deposition. Recall that when we deposit with "halos", we need to fetch the subset of coarse-level particles that surround the
  // refinement boundary.
  switch (a_coarseFineDeposition) {
  case CoarseFineDeposition::Interp: {
    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    DepositionType::TSC,
                                                    CoarseFineDeposition::Interp,
                                                    m_forceIrregDepositionNGP);
    break;
  }
  case CoarseFineDeposition::Halo: {

#warning "Need to try to deposit the halo particles using HALF the grid width."
    // Copy the mask particles.
    const AMRMask& mask = m_amr->getMask(s_outer_particle_halo, m_outerHaloBuffer, m_realm);
    a_particles.copyMaskParticles(mask);

    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::Halo,
                                                    m_forceIrregDepositionNGP);

    a_particles.clearMaskParticles();

    break;
  }
  case CoarseFineDeposition::HaloNGP: {

    // Transfer mask particles
    const AMRMask& mask = m_amr->getMask(s_outer_particle_halo, m_outerHaloBuffer, m_realm);
    a_particles.transferMaskParticles(mask);

    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::HaloNGP,
                                                    m_forceIrregDepositionNGP);

    a_particles.transferParticles(a_particles.getMaskParticles());

    break;
  }
  case CoarseFineDeposition::Buffer: {
#warning "Some of this code should be moved to EBAMRParticleMesh"
    // Transfer mask particles
    const AMRMask& mask = m_amr->getMask(s_inner_particle_halo, m_innerHaloBuffer, m_realm);

    a_particles.transferMaskParticles(mask);

    for (int lvl = m_amr->getFinestLevel(); lvl > 0; lvl--) {
      ParticleData<P>& coarParticles = a_particles[lvl - 1];
      ParticleData<P>& maskParticles = *a_particles.getMaskParticles()[lvl];

#warning "Could there be a bug here related to if there aren't enough MPI ranks...?"
      List<P>& coarOutcasts = coarParticles.outcast();

      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        coarOutcasts.catenate(maskParticles[dit()].listItems());
      }

      coarParticles.gatherOutcast();
      coarParticles.remapOutcast();
    }

    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::Interp,
                                                    m_forceIrregDepositionNGP);

    a_particles.remap();

#warning "Debug code here to test the halo"
#if 0
    DataOps::setValue(a_phi, 0.0);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const auto& amrMask = m_amr->getMask(s_inner_particle_halo, m_innerHaloBuffer, m_realm);
        const auto& mask    = (*amrMask[lvl])[dit()];

        FArrayBox& phi = (*a_phi[lvl])[dit()].getFArrayBox();

        auto kernel = [&](const IntVect& iv) -> void {
          if (mask(iv, 0)) {
            phi(iv) = 1.0;
          }
        };

        BoxLoops::loop(dbl[dit()], kernel);
      }
    }
#endif

    break;
  }
  default: {
    MayDay::Error("ItoSolverImplem.H in function ItoSolver::depositKappaConservative -- logic bust!");

    break;
  }
  }
}
template <class P, Real (P::*particleScalarField)() const>
void
ItoSolver::depositKappaConservative(EBAMRCellData&             a_phi,
                                    ParticleContainer<P>&      a_particles,
                                    const DepositionType       a_deposition,
                                    const CoarseFineDeposition a_coarseFineDeposition) const
{
  CH_TIME("ItoSolver::depositKappaConservative");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositKappaConservative" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);

  // Now do the deposition. Recall that when we deposit with "halos", we need to fetch the subset of coarse-level particles that surround the
  // refinement boundary.
  switch (a_coarseFineDeposition) {
  case CoarseFineDeposition::Interp: {
    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::Interp,
                                                    m_forceIrregDepositionNGP);
    break;
  }
  case CoarseFineDeposition::Halo: {

    // Copy the mask particles.
    const AMRMask& mask = m_amr->getMask(s_outer_particle_halo, m_outerHaloBuffer, m_realm);
    a_particles.copyMaskParticles(mask);

    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::Halo,
                                                    m_forceIrregDepositionNGP);

    a_particles.clearMaskParticles();

    break;
  }
  case CoarseFineDeposition::HaloNGP: {

    // Transfer mask particles
    const AMRMask& mask = m_amr->getMask(s_outer_particle_halo, m_outerHaloBuffer, m_realm);
    a_particles.transferMaskParticles(mask);

    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::HaloNGP,
                                                    m_forceIrregDepositionNGP);

    a_particles.transferParticles(a_particles.getMaskParticles());

    break;
  }
  case CoarseFineDeposition::Buffer: {
    MayDay::Warning("ItoSolverImplem.H -- copy code from other function once we're done");

    break;
  }
  default: {
    MayDay::Error("ItoSolverImplem.H in function ItoSolver::depositKappaConservative -- logic bust!");

    break;
  }
  }
}

#include <CD_NamespaceFooter.H>

#endif
