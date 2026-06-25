/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_TracerParticle.H
  @brief  SoA payload for tracer particles (velocity + Runge-Kutta stage scratch)
  @author Robert Marskar
*/

#ifndef CD_TRACERPARTICLE_H
#define CD_TRACERPARTICLE_H

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief SoA payload for tracer particles: interpolated velocity plus Runge-Kutta stage scratch.
  @details Per-component scalar columns (one member per spatial dimension; the z-component exists only when
  CH_SPACEDIM == 3). vx/vy/vz hold the interpolated velocity (required by TracerParticleSolver). xk is the
  position snapshot x^k; k1/k2/k3 hold the Runge-Kutta stage velocities. Weight and position are
  container-owned columns, NOT payload members.
*/
struct TracerParticle
{
  ParticleReal vx = 0.0; ///< Interpolated velocity, x-component.
  ParticleReal vy = 0.0; ///< Interpolated velocity, y-component.
#if CH_SPACEDIM == 3
  ParticleReal vz = 0.0; ///< Interpolated velocity, z-component.
#endif

  ParticleReal xkx = 0.0; ///< Position snapshot x^k, x-component.
  ParticleReal xky = 0.0; ///< Position snapshot x^k, y-component.
#if CH_SPACEDIM == 3
  ParticleReal xkz = 0.0; ///< Position snapshot x^k, z-component.
#endif

  ParticleReal k1x = 0.0; ///< Runge-Kutta stage k1, x-component.
  ParticleReal k1y = 0.0; ///< Runge-Kutta stage k1, y-component.
#if CH_SPACEDIM == 3
  ParticleReal k1z = 0.0; ///< Runge-Kutta stage k1, z-component.
#endif

  ParticleReal k2x = 0.0; ///< Runge-Kutta stage k2, x-component.
  ParticleReal k2y = 0.0; ///< Runge-Kutta stage k2, y-component.
#if CH_SPACEDIM == 3
  ParticleReal k2z = 0.0; ///< Runge-Kutta stage k2, z-component.
#endif

  ParticleReal k3x = 0.0; ///< Runge-Kutta stage k3, x-component.
  ParticleReal k3y = 0.0; ///< Runge-Kutta stage k3, y-component.
#if CH_SPACEDIM == 3
  ParticleReal k3z = 0.0; ///< Runge-Kutta stage k3, z-component.
#endif
};

/**
  @brief ParticleTraits specialization listing the payload columns of TracerParticle.
*/
template <>
struct ParticleTraits<TracerParticle>
{
  /// @brief Payload column member-pointers, listed one field group at a time (SpaceDim entries per group).
  static constexpr auto columns = std::make_tuple(
    D_DECL(&TracerParticle::vx, &TracerParticle::vy, &TracerParticle::vz),
    D_DECL(&TracerParticle::xkx, &TracerParticle::xky, &TracerParticle::xkz),
    D_DECL(&TracerParticle::k1x, &TracerParticle::k1y, &TracerParticle::k1z),
    D_DECL(&TracerParticle::k2x, &TracerParticle::k2y, &TracerParticle::k2z),
    D_DECL(&TracerParticle::k3x, &TracerParticle::k3y, &TracerParticle::k3z));
};

#include <CD_NamespaceFooter.H>

#endif
