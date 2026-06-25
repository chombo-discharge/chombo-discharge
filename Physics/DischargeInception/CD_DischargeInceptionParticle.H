/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_DischargeInceptionParticle.H
  @brief  SoA payload for the discharge-inception tracer particles
  @author Robert Marskar
*/

#ifndef CD_DISCHARGEINCEPTIONPARTICLE_H
#define CD_DISCHARGEINCEPTIONPARTICLE_H

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief SoA payload for the discharge-inception integrators.
  @details Per-component scalar columns (one member per spatial dimension; the z-component exists only when
  CH_SPACEDIM == 3). vx/vy/vz hold the interpolated electron velocity (required by TracerParticleSolver).
  x0 is the position snapshot used to rewind particles after an integration; vk is the velocity snapshot
  from the previous Runge-Kutta stage; ga is the interpolated grad(alpha) field. alphaEff and dtStep cache
  the effective ionization coefficient and step size from the predictor stage. Weight (the inception
  integral K) and position are container-owned columns.
*/
struct DischargeInceptionParticle
{
  ParticleReal vx = 0.0; ///< Interpolated velocity, x-component.
  ParticleReal vy = 0.0; ///< Interpolated velocity, y-component.
#if CH_SPACEDIM == 3
  ParticleReal vz = 0.0; ///< Interpolated velocity, z-component.
#endif

  ParticleReal x0x = 0.0; ///< Position snapshot for rewind, x-component.
  ParticleReal x0y = 0.0; ///< Position snapshot for rewind, y-component.
#if CH_SPACEDIM == 3
  ParticleReal x0z = 0.0; ///< Position snapshot for rewind, z-component.
#endif

  ParticleReal vkx = 0.0; ///< Previous-stage velocity snapshot, x-component.
  ParticleReal vky = 0.0; ///< Previous-stage velocity snapshot, y-component.
#if CH_SPACEDIM == 3
  ParticleReal vkz = 0.0; ///< Previous-stage velocity snapshot, z-component.
#endif

  ParticleReal gax = 0.0; ///< Interpolated grad(alpha), x-component.
  ParticleReal gay = 0.0; ///< Interpolated grad(alpha), y-component.
#if CH_SPACEDIM == 3
  ParticleReal gaz = 0.0; ///< Interpolated grad(alpha), z-component.
#endif

  ParticleReal alphaEff = 0.0; ///< Cached effective ionization coefficient.
  ParticleReal dtStep   = 0.0; ///< Cached predictor step size.
};

/**
  @brief ParticleTraits specialization listing the payload columns of DischargeInceptionParticle.
*/
template <>
struct ParticleTraits<DischargeInceptionParticle>
{
  /// @brief Payload column member-pointers, listed one field group at a time (SpaceDim entries per group).
  static constexpr auto columns = std::make_tuple(
    D_DECL(&DischargeInceptionParticle::vx, &DischargeInceptionParticle::vy, &DischargeInceptionParticle::vz),
    D_DECL(&DischargeInceptionParticle::x0x, &DischargeInceptionParticle::x0y, &DischargeInceptionParticle::x0z),
    D_DECL(&DischargeInceptionParticle::vkx, &DischargeInceptionParticle::vky, &DischargeInceptionParticle::vkz),
    D_DECL(&DischargeInceptionParticle::gax, &DischargeInceptionParticle::gay, &DischargeInceptionParticle::gaz),
    &DischargeInceptionParticle::alphaEff,
    &DischargeInceptionParticle::dtStep);
};

#include <CD_NamespaceFooter.H>

#endif
