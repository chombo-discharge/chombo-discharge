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
  @details Per-component scalar columns (D_DECL). vx/vy/vz hold the interpolated electron velocity
  (required by TracerParticleSolver). x0 is the position snapshot used to rewind particles after an
  integration; vk is the velocity snapshot from the previous Runge-Kutta stage; ga is the interpolated
  grad(alpha) field. alphaEff and dtStep cache the effective ionization coefficient and step size from
  the predictor stage. Weight (the inception integral K) and position are container-owned columns.
*/
struct DischargeInceptionParticle
{
  ParticleReal D_DECL(vx = 0.0, vy = 0.0, vz = 0.0);    ///< Interpolated velocity (per component).
  ParticleReal D_DECL(x0x = 0.0, x0y = 0.0, x0z = 0.0); ///< Position snapshot for rewind (per component).
  ParticleReal D_DECL(vkx = 0.0, vky = 0.0, vkz = 0.0); ///< Previous-stage velocity snapshot (per component).
  ParticleReal D_DECL(gax = 0.0, gay = 0.0, gaz = 0.0); ///< Interpolated grad(alpha) (per component).
  ParticleReal alphaEff = 0.0;                          ///< Cached effective ionization coefficient.
  ParticleReal dtStep   = 0.0;                          ///< Cached predictor step size.
};

/**
  @brief ParticleTraits specialization listing the payload columns of DischargeInceptionParticle.
*/
template <>
struct ParticleTraits<DischargeInceptionParticle>
{
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
