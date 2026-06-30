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
  CH_SPACEDIM == 3). v_x/v_y/v_z hold the interpolated velocity (required by TracerParticleSolver). xk is the
  position snapshot x^k; k1/k2/k3 hold the Runge-Kutta stage velocities. Weight and position are
  container-owned columns, NOT payload members.
*/
struct TracerParticle
{
  ParticleReal v_x = 0.0; ///< Interpolated velocity, x-component.
  ParticleReal v_y = 0.0; ///< Interpolated velocity, y-component.
#if CH_SPACEDIM == 3
  ParticleReal v_z = 0.0; ///< Interpolated velocity, z-component.
#endif

  double xk_x = 0.0; ///< Position snapshot x^k, x-component (double: position-like).
  double xk_y = 0.0; ///< Position snapshot x^k, y-component (double: position-like).
#if CH_SPACEDIM == 3
  double xk_z = 0.0; ///< Position snapshot x^k, z-component (double: position-like).
#endif

  ParticleReal k1_x = 0.0; ///< Runge-Kutta stage k1, x-component.
  ParticleReal k1_y = 0.0; ///< Runge-Kutta stage k1, y-component.
#if CH_SPACEDIM == 3
  ParticleReal k1_z = 0.0; ///< Runge-Kutta stage k1, z-component.
#endif

  ParticleReal k2_x = 0.0; ///< Runge-Kutta stage k2, x-component.
  ParticleReal k2_y = 0.0; ///< Runge-Kutta stage k2, y-component.
#if CH_SPACEDIM == 3
  ParticleReal k2_z = 0.0; ///< Runge-Kutta stage k2, z-component.
#endif

  ParticleReal k3_x = 0.0; ///< Runge-Kutta stage k3, x-component.
  ParticleReal k3_y = 0.0; ///< Runge-Kutta stage k3, y-component.
#if CH_SPACEDIM == 3
  ParticleReal k3_z = 0.0; ///< Runge-Kutta stage k3, z-component.
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
    D_DECL(&TracerParticle::v_x, &TracerParticle::v_y, &TracerParticle::v_z),
    D_DECL(&TracerParticle::xk_x, &TracerParticle::xk_y, &TracerParticle::xk_z),
    D_DECL(&TracerParticle::k1_x, &TracerParticle::k1_y, &TracerParticle::k1_z),
    D_DECL(&TracerParticle::k2_x, &TracerParticle::k2_y, &TracerParticle::k2_z),
    D_DECL(&TracerParticle::k3_x, &TracerParticle::k3_y, &TracerParticle::k3_z));

  /**
     @brief Payload columns written to H5Part visualization files: the interpolated velocity vector (id,
            position and weight are written automatically; the position snapshot and Runge-Kutta stage
            columns are internal and not plotted).
  */
  static constexpr auto h5PartColumns = std::make_tuple(
    H5Part::Vector{"v", D_DECL(&TracerParticle::v_x, &TracerParticle::v_y, &TracerParticle::v_z)});
};

#include <CD_NamespaceFooter.H>

#endif
