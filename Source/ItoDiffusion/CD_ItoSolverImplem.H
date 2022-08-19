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

inline Real
ItoSolver::sign(const Real& a) const
{
  return (a > 0) - (a < 0);
}

inline RealVect
ItoSolver::randomGaussian() const
{
  // TLDR: We draw a random number from a Gaussian distribution for each coordinate, and truncate the distribution at m_normalDistributionTruncation.

  RealVect r = RealVect::Zero;
  for (int i = 0; i < SpaceDim; i++) {
    r[i] = Random::getNormal01();
    r[i] = sign(r[i]) * std::min(std::abs(r[i]), m_normalDistributionTruncation);
  }

  return r;
}

inline RealVect
ItoSolver::randomPosition(const RealVect a_cellPos,
                          const RealVect a_lo,
                          const RealVect a_hi,
                          const RealVect a_bndryCentroid,
                          const RealVect a_bndryNormal,
                          const Real     a_dx,
                          const Real     a_kappa) const
{

  RealVect pos;

  if (a_kappa < 1.0) { // Rejection sampling.
    pos = this->randomPosition(a_lo, a_hi, a_bndryCentroid, a_bndryNormal);
  }
  else { // Regular cell. Get a position.
    pos = this->randomPosition(a_lo, a_hi);
  }

  // Convert from unit cell coordinates to physical coordinates.
  pos = a_cellPos + pos * a_dx;

  return pos;
}

inline RealVect
ItoSolver::randomPosition(const RealVect a_lo,
                          const RealVect a_hi,
                          const RealVect a_bndryCentroid,
                          const RealVect a_bndryNormal) const
{

  // Draw a position.
  RealVect pos   = this->randomPosition(a_lo, a_hi);
  bool     valid = PolyGeom::dot(pos - a_bndryCentroid, a_bndryNormal) >= 0.0;

  // Rejectino sampling.
  while (!valid) {
    pos   = this->randomPosition(a_lo, a_hi);
    valid = PolyGeom::dot(pos - a_bndryCentroid, a_bndryNormal) >= 0.0;
  }

  return pos;
}

inline RealVect
ItoSolver::randomPosition(const RealVect a_lo, const RealVect a_hi) const
{
  RealVect pos;

  for (int dir = 0; dir < SpaceDim; dir++) {
    pos[dir] = a_lo[dir] + 0.5 * (1.0 + Random::getUniformReal11()) * (a_hi[dir] - a_lo[dir]);
  }

  return pos;
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
ItoSolver::depositParticles(EBAMRCellData&              a_phi,
                            const ParticleContainer<P>& a_particles,
                            const DepositionType        a_deposition) const
{
  CH_TIME("ItoSolver::depositParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isCellSorted());

  // TLDR: First, deposit onto the mesh as usual (as if the EB wasn't there). If the user asks for it, he can redistribute mass in order to
  //       conserve total mass (if that is important). But the corresponding scheme will be O(1) accurate.

  this->depositKappaConservative<
    P,
    particleScalarField>(a_phi, a_particles, a_deposition); // a_phi contains only weights, i.e. not divided by kappa
  this->redistributeAMR(a_phi);                             // Redistribution magic, if you want it.

  // Average down and interpolate
  m_amr->averageDown(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);
}

template <class P, Real (P::*particleScalarField)() const>
void
ItoSolver::depositParticles(EBAMRCellData&              a_phi,
                            const ParticleContainer<P>& a_particles,
                            const DepositionType        a_deposition) const
{
  CH_TIME("ItoSolver::depositParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isCellSorted());

  // TLDR: First, deposit onto the mesh as usual (as if the EB wasn't there). If the user asks for it, he can redistribute mass in order to
  //       conserve total mass (if that is important). But the corresponding scheme will be O(1) accurate.

  this->depositKappaConservative<
    P,
    particleScalarField>(a_phi, a_particles, a_deposition); // a_phi contains only weights, i.e. not divided by kappa
  this->redistributeAMR(a_phi);                             // Redistribution magic, if you want it.

  // Average down and interpolate
  m_amr->averageDown(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);
}

template <class P, const Real& (P::*particleScalarField)() const>
void
ItoSolver::depositKappaConservative(EBAMRCellData&              a_phi,
                                    const ParticleContainer<P>& a_particles,
                                    const DepositionType        a_deposition) const
{
  CH_TIME("ItoSolver::depositKappaConservative");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositKappaConservative" << endl;
  }

  CH_assert(m_haloBuffer == 0 || m_pvrBuffer == 0);
  CH_assert(a_phi[0]->nComp() == 1);

  if (m_haloBuffer > 0 && m_pvrBuffer > 0) {
    MayDay::Error(
      "ItoSolver::depositKappaConservative - conflicting inputs, halo/pvr buffers both > 0 but one of them must be 0");
  }

  // Figure out how to handle the coarse-fine deposition.
  CoarseFineDeposition coarseFineDeposition;
  if (m_pvrBuffer > 0) {
    coarseFineDeposition = CoarseFineDeposition::Interp;
  }
  else if (m_haloBuffer > 0) {
    if (m_forceHaloNGP) {
      coarseFineDeposition = CoarseFineDeposition::HaloNGP;
    }
    else {
      coarseFineDeposition = CoarseFineDeposition::Halo;
    }
  }
  else {
    MayDay::Error("ItoSolver::depositKappaConservative -- logic bust");
  }

  // Now do the deposition. Recall that when we deposit with "halos", we need to fetch the subset of coarse-level particles that surround the
  // refinement boundary.
  switch (coarseFineDeposition) {
  case CoarseFineDeposition::Interp:
    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::Interp,
                                                    m_forceIrregDepositionNGP);
    break;
  case CoarseFineDeposition::Halo: {
    const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);
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

    MayDay::Error("CD_ItoSolverImplem.H -- HaloNGP not supported");

    break;
  }
  default:
    MayDay::Error("ItoSolverImplem.H in function ItoSolver::depositKappaConservative -- logic bust!");
    break;
  }
}
template <class P, Real (P::*particleScalarField)() const>
void
ItoSolver::depositKappaConservative(EBAMRCellData&              a_phi,
                                    const ParticleContainer<P>& a_particles,
                                    const DepositionType        a_deposition) const
{
  CH_TIME("ItoSolver::depositKappaConservative");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositKappaConservative" << endl;
  }

  CH_assert(m_haloBuffer == 0 || m_pvrBuffer == 0);
  CH_assert(a_phi[0]->nComp() == 1);

  if (m_haloBuffer > 0 && m_pvrBuffer > 0) {
    MayDay::Error(
      "ItoSolver::depositKappaConservative - conflicting inputs, halo/pvr buffers both > 0 but one of them must be 0");
  }

  // Figure out how to handle the coarse-fine deposition.
  CoarseFineDeposition coarseFineDeposition;
  if (m_pvrBuffer > 0) {
    coarseFineDeposition = CoarseFineDeposition::Interp;
  }
  else if (m_haloBuffer > 0) {
    if (m_forceHaloNGP) {
      coarseFineDeposition = CoarseFineDeposition::HaloNGP;
    }
    else {
      coarseFineDeposition = CoarseFineDeposition::Halo;
    }
  }
  else {
    MayDay::Error("ItoSolver::depositKappaConservative -- logic bust");
  }

  // Now do the deposition. Recall that when we deposit with "halos", we need to fetch the subset of coarse-level particles that surround the
  // refinement boundary.
  switch (coarseFineDeposition) {
  case CoarseFineDeposition::Interp:
    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::Interp,
                                                    m_forceIrregDepositionNGP);
    break;
  case CoarseFineDeposition::Halo: {
    const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);
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

    MayDay::Error("CD_ItoSolverImplem.H -- HaloNGP not supported");

    break;
  }
  default:
    MayDay::Error("ItoSolverImplem.H in function ItoSolver::depositKappaConservative -- logic bust!");
    break;
  }
}

#include <CD_NamespaceFooter.H>

#endif
