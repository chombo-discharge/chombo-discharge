/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_DemoParticle.H
  @brief  A demonstration particle type used to exercise ParticleSoA.
  @author Robert Marskar
  @details DemoParticle is a drifting computational super-particle carrying a
  mobility and a derived conductivity field. It serves as a concrete example of the
  recipe documented on ParticleTraits: a plain struct of byte-copyable members plus
  a ParticleTraits specialization that names the columns by member pointer.
*/

#ifndef CD_DEMOPARTICLE_H
#define CD_DEMOPARTICLE_H

// Std includes
#include <tuple>

// Chombo includes
#include <REAL.H>
#include <RealVect.H>

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Demonstration particle: position, weight, drift velocity, mobility, and a
         conductivity field filled by SoA kernels.
*/
struct DemoParticle
{
  RealVect pos          = RealVect::Zero; ///< Position (mandatory vector field)
  Real     weight       = 0.0;            ///< Number of physical particles (mandatory scalar; acts as "mass")
  RealVect velocity     = RealVect::Zero; ///< Drift velocity
  Real     mobility     = 0.0;            ///< Mobility coefficient
  Real     conductivity = 0.0;            ///< weight * mobility (filled by a kernel)
};

/**
  @brief ParticleTraits specialization describing the SoA layout of DemoParticle.
*/
template <>
struct ParticleTraits<DemoParticle>
{
  /** @brief One pointer-to-member per SoA column (single source of truth for the layout). */
  static constexpr auto columns = std::make_tuple(&DemoParticle::pos,
                                                  &DemoParticle::weight,
                                                  &DemoParticle::velocity,
                                                  &DemoParticle::mobility,
                                                  &DemoParticle::conductivity);

  /** @brief Pointer to the mandatory RealVect position member. */
  static constexpr auto positionPtr = &DemoParticle::pos;

  /** @brief Pointer to the mandatory Real weight member. */
  static constexpr auto weightPtr = &DemoParticle::weight;
};

#include <CD_NamespaceFooter.H>

#endif
