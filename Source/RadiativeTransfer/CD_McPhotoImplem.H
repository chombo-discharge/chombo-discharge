/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_McPhotoImplem.H
  @brief  Implementation of CD_McPhoto.H
  @author Robert Marskar
*/

#ifndef CD_MCPHOTOIMPLEM_H
#define CD_MCPHOTOIMPLEM_H

// Our includes
#include <CD_McPhoto.H>
#include <CD_NamespaceHeader.H>

template <class P, class Ret, Ret (P::*MemberFunc)() const>
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
  this->depositKappaConservative<P, Ret, MemberFunc>(a_phi, a_photons, a_deposition, m_coarseFineDeposition);

  // Blend with the non-conservative divergence, redistribute, average down, and interpolate ghosts.
  this->depositHybridDivergence(a_phi);
}

template <class P, class Ret, Ret (P::*MemberFunc)() const>
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

  m_amr->depositParticles<P, Ret, MemberFunc>(a_phi,
                                              m_realm,
                                              m_phase,
                                              a_particles,
                                              a_deposition,
                                              a_coarseFineDeposition,
                                              false);
}

#include <CD_NamespaceFooter.H>

#endif
