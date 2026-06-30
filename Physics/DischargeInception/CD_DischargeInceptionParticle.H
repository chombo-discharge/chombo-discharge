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
  CH_SPACEDIM == 3). v_x/v_y/v_z hold the interpolated electron velocity (required by TracerParticleSolver).
  x0 is the position snapshot used to rewind particles after an integration; vk is the velocity snapshot
  from the previous Runge-Kutta stage; ga is the interpolated grad(alpha) field. alphaEff and dtStep cache
  the effective ionization coefficient and step size from the predictor stage. Weight (the inception
  integral K) and position are container-owned columns.
*/
struct DischargeInceptionParticle
{
  ParticleReal v_x = 0.0; ///< Interpolated velocity, x-component.
  ParticleReal v_y = 0.0; ///< Interpolated velocity, y-component.
#if CH_SPACEDIM == 3
  ParticleReal v_z = 0.0; ///< Interpolated velocity, z-component.
#endif

  double x0_x = 0.0; ///< Position snapshot for rewind, x-component (double: position-like).
  double x0_y = 0.0; ///< Position snapshot for rewind, y-component (double: position-like).
#if CH_SPACEDIM == 3
  double x0_z = 0.0; ///< Position snapshot for rewind, z-component (double: position-like).
#endif

  ParticleReal vk_x = 0.0; ///< Previous-stage velocity snapshot, x-component.
  ParticleReal vk_y = 0.0; ///< Previous-stage velocity snapshot, y-component.
#if CH_SPACEDIM == 3
  ParticleReal vk_z = 0.0; ///< Previous-stage velocity snapshot, z-component.
#endif

  ParticleReal ga_x = 0.0; ///< Interpolated grad(alpha), x-component.
  ParticleReal ga_y = 0.0; ///< Interpolated grad(alpha), y-component.
#if CH_SPACEDIM == 3
  ParticleReal ga_z = 0.0; ///< Interpolated grad(alpha), z-component.
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
    D_DECL(&DischargeInceptionParticle::v_x, &DischargeInceptionParticle::v_y, &DischargeInceptionParticle::v_z),
    D_DECL(&DischargeInceptionParticle::x0_x, &DischargeInceptionParticle::x0_y, &DischargeInceptionParticle::x0_z),
    D_DECL(&DischargeInceptionParticle::vk_x, &DischargeInceptionParticle::vk_y, &DischargeInceptionParticle::vk_z),
    D_DECL(&DischargeInceptionParticle::ga_x, &DischargeInceptionParticle::ga_y, &DischargeInceptionParticle::ga_z),
    &DischargeInceptionParticle::alphaEff,
    &DischargeInceptionParticle::dtStep);

  /**
     @brief Payload columns written to H5Part visualization files: the interpolated velocity vector (id,
            position and weight are written automatically; the remaining columns are internal integrator
            state and not plotted).
  */
  static constexpr auto h5PartColumns = std::make_tuple(H5Part::Vector{
    "v",
    D_DECL(&DischargeInceptionParticle::v_x, &DischargeInceptionParticle::v_y, &DischargeInceptionParticle::v_z)});
};

#include <CD_NamespaceFooter.H>

#endif
