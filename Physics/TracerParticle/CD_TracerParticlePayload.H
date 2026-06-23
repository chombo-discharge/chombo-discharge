/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_TracerParticlePayload.H
  @brief  SoA payload for tracer particles (velocity + Runge-Kutta stage scratch)
  @author Robert Marskar
*/

#ifndef CD_TRACERPARTICLEPAYLOAD_H
#define CD_TRACERPARTICLEPAYLOAD_H

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief SoA payload for tracer particles: interpolated velocity plus Runge-Kutta stage scratch.
  @details Per-component scalar columns (D_DECL). vx/vy/vz hold the interpolated velocity (required by
  TracerParticleSolver). xk is the position snapshot x^k; k1/k2/k3 hold the Runge-Kutta stage
  velocities. Weight and position are container-owned columns, NOT payload members.
*/
struct TracerParticlePayload
{
  ParticleReal D_DECL(vx = 0.0, vy = 0.0, vz = 0.0);    ///< Interpolated velocity (per component).
  ParticleReal D_DECL(xkx = 0.0, xky = 0.0, xkz = 0.0); ///< Position snapshot x^k (per component).
  ParticleReal D_DECL(k1x = 0.0, k1y = 0.0, k1z = 0.0); ///< Runge-Kutta stage k1 (per component).
  ParticleReal D_DECL(k2x = 0.0, k2y = 0.0, k2z = 0.0); ///< Runge-Kutta stage k2 (per component).
  ParticleReal D_DECL(k3x = 0.0, k3y = 0.0, k3z = 0.0); ///< Runge-Kutta stage k3 (per component).
};

/**
  @brief ParticleTraits specialization listing the payload columns of TracerParticlePayload.
*/
template <>
struct ParticleTraits<TracerParticlePayload>
{
  static constexpr auto columns = std::make_tuple(
    D_DECL(&TracerParticlePayload::vx, &TracerParticlePayload::vy, &TracerParticlePayload::vz),
    D_DECL(&TracerParticlePayload::xkx, &TracerParticlePayload::xky, &TracerParticlePayload::xkz),
    D_DECL(&TracerParticlePayload::k1x, &TracerParticlePayload::k1y, &TracerParticlePayload::k1z),
    D_DECL(&TracerParticlePayload::k2x, &TracerParticlePayload::k2y, &TracerParticlePayload::k2z),
    D_DECL(&TracerParticlePayload::k3x, &TracerParticlePayload::k3y, &TracerParticlePayload::k3z));
};

#include <CD_NamespaceFooter.H>

#endif
