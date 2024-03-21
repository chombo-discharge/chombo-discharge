/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_McPhotoImplem.H
  @brief  Implementation of CD_McPhoto.H
  @author Robert Marskar
*/

#ifndef CD_McPhotoImplem_H
#define CD_McPhotoImplem_H

// Our includes
#include <CD_McPhoto.H>
#include <CD_NamespaceHeader.H>

template <class P, const Real& (P::*particleScalarField)() const>
void
McPhoto::depositPhotons(EBAMRCellData&        a_phi,
                        ParticleContainer<P>& a_photons,
                        const DepositionType& a_deposition) const noexcept
{
  CH_TIME("McPhoto::depositPhotons(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositPhotons(ParticleContainer)" << endl;
  }

  // a_phi contains only weights, i.e. not divided by kappa
  this->depositKappaConservative<P, particleScalarField>(a_phi, a_photons, a_deposition, m_coarseFineDeposition);

  // Compute m_depositionNC = sum(kappa*Wc)/sum(kappa)
  this->depositNonConservative(m_depositionNC, a_phi);

  // Compute hybrid deposition, including mass differnce
  this->depositHybrid(a_phi, m_massDiff, m_depositionNC);

  // Redistribute
  if (m_blendConservation) {
    Vector<RefCountedPtr<EBFluxRedistribution>>& redistOps = m_amr->getRedistributionOp(m_realm, m_phase);
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real     scale     = 1.0;
      const Interval variables = Interval(0, 0);
      const bool     hasCoar   = lvl > 0;
      const bool     hasFine   = lvl < m_amr->getFinestLevel();

      if (hasCoar) {
        redistOps[lvl]->redistributeCoar(*a_phi[lvl - 1], *m_massDiff[lvl], scale, variables);
      }

      redistOps[lvl]->redistributeLevel(*a_phi[lvl], *m_massDiff[lvl], scale, variables);

      if (hasFine) {
        redistOps[lvl]->redistributeFine(*a_phi[lvl + 1], *m_massDiff[lvl], scale, variables);
      }
    }
  }

  // Average down and interpolate
  m_amr->conservativeAverage(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);
}

template <class P, const Real& (P::*particleScalarField)() const>
void
McPhoto::depositKappaConservative(EBAMRCellData&             a_phi,
                                  ParticleContainer<P>&      a_particles,
                                  const DepositionType       a_deposition,
                                  const CoarseFineDeposition a_coarseFineDeposition) const noexcept
{
  CH_TIME("McPhoto::depositKappaConservative");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositKappaConservative" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);

  switch (a_coarseFineDeposition) {
  case CoarseFineDeposition::Interp: {
    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::Interp,
                                                    false);

    break;
  }
  case CoarseFineDeposition::Halo: {

    // Copy particles living on the mask.
    const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);
    a_particles.copyMaskParticles(mask);

    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::Halo,
                                                    false);

    // Clear out the mask particles.
    a_particles.clearMaskParticles();

    break;
  }
  case CoarseFineDeposition::HaloNGP: {
    const AMRMask& mask = m_amr->getMask(s_particle_halo, m_haloBuffer, m_realm);

    // Transfer particles living on the mask.
    a_particles.transferMaskParticles(mask);

    m_amr->depositParticles<P, particleScalarField>(a_phi,
                                                    m_realm,
                                                    m_phase,
                                                    a_particles,
                                                    a_deposition,
                                                    CoarseFineDeposition::HaloNGP,
                                                    false);

    // Transfer them back.
    a_particles.transferParticles(a_particles.getMaskParticles());

    break;
  }
  default: {
    MayDay::Error("McPhoto::depositKappaConservative -- logic bust due to unsupported coarse-fine deposition");
  }
  }
}

#include <CD_NamespaceFooter.H>

#endif
