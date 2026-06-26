/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Photon.H
  @brief  SoA payload for Monte Carlo radiative-transfer photons
  @author Robert Marskar
*/

#ifndef CD_PHOTON_H
#define CD_PHOTON_H

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief SoA payload for Monte Carlo radiative-transfer photons.
  @details kappa is the mean absorption coefficient (inverse mean free path); velocity is the photon
  propagation velocity (one member per spatial dimension; the z-component exists only when
  CH_SPACEDIM == 3). Weight and position are container-owned columns, NOT payload members. All payload
  columns (kappa and velocity) are checkpointed; position and weight are always checkpointed by the
  container, so no h5Columns subset is declared.
*/
struct Photon
{
  ParticleReal kappa = 0.0; ///< Mean absorption coefficient (inverse mean free path).

  ParticleReal v_x = 0.0; ///< Velocity, x-component.
  ParticleReal v_y = 0.0; ///< Velocity, y-component.
#if CH_SPACEDIM == 3
  ParticleReal v_z = 0.0; ///< Velocity, z-component.
#endif
};

/**
  @brief ParticleTraits specialization listing the payload columns of Photon.
*/
template <>
struct ParticleTraits<Photon>
{
  /// @brief Payload column member-pointers: kappa followed by the per-component velocity.
  static constexpr auto columns = std::make_tuple(&Photon::kappa, D_DECL(&Photon::v_x, &Photon::v_y, &Photon::v_z));
};

#include <CD_NamespaceFooter.H>

#endif
