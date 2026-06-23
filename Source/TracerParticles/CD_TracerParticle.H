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
  @details Per-component scalar columns (D_DECL). vx/vy/vz hold the interpolated velocity (required by
  TracerParticleSolver). xk is the position snapshot x^k; k1/k2/k3 hold the Runge-Kutta stage
  velocities. Weight and position are container-owned columns, NOT payload members.
*/
struct TracerParticle
{
  ParticleReal D_DECL(vx = 0.0, vy = 0.0, vz = 0.0);    ///< Interpolated velocity (per component).
  ParticleReal D_DECL(xkx = 0.0, xky = 0.0, xkz = 0.0); ///< Position snapshot x^k (per component).
  ParticleReal D_DECL(k1x = 0.0, k1y = 0.0, k1z = 0.0); ///< Runge-Kutta stage k1 (per component).
  ParticleReal D_DECL(k2x = 0.0, k2y = 0.0, k2z = 0.0); ///< Runge-Kutta stage k2 (per component).
  ParticleReal D_DECL(k3x = 0.0, k3y = 0.0, k3z = 0.0); ///< Runge-Kutta stage k3 (per component).
};

/**
  @brief ParticleTraits specialization listing the payload columns of TracerParticle.
*/
template <>
struct ParticleTraits<TracerParticle>
{
  static constexpr auto columns = std::make_tuple(
    D_DECL(&TracerParticle::vx, &TracerParticle::vy, &TracerParticle::vz),
    D_DECL(&TracerParticle::xkx, &TracerParticle::xky, &TracerParticle::xkz),
    D_DECL(&TracerParticle::k1x, &TracerParticle::k1y, &TracerParticle::k1z),
    D_DECL(&TracerParticle::k2x, &TracerParticle::k2y, &TracerParticle::k2z),
    D_DECL(&TracerParticle::k3x, &TracerParticle::k3y, &TracerParticle::k3z));
};

#include <CD_NamespaceFooter.H>

#endif
