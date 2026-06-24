/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoParticle.H
  @brief  SoA payload for Ito diffusion particles (drifting Brownian walkers)
  @author Robert Marskar
*/

#ifndef CD_ITOPARTICLE_H
#define CD_ITOPARTICLE_H

// Std includes
#include <string>
#include <vector>

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief SoA payload for ItoSolver particles, i.e. drifting Brownian walkers.
  @details Per-component scalar columns (D_DECL). mobility/diffusion/energy are the particle transport
  coefficients and average energy; tmpReal is scratch scalar storage. oldPosition holds the previous
  particle position (double precision, because it subtracts from the double-precision container position);
  velocity is the interpolated drift velocity; tmpVect is scratch vector storage. Weight and position are
  container-owned columns, NOT payload members. Only the velocity columns are checkpointed (position and
  weight are always checkpointed by the container).
*/
struct ItoParticle
{
  /**
    @brief Naming convention for scalar fields (used by ItoKMC plotting).
  */
  static std::vector<std::string> s_realVariables;

  /**
    @brief Naming convention for vector fields (used by ItoKMC plotting).
  */
  static std::vector<std::string> s_vectVariables;

  ParticleReal mobility  = 0.0; ///< Mobility coefficient.
  ParticleReal diffusion = 0.0; ///< Diffusion coefficient.
  ParticleReal energy    = 0.0; ///< Average particle energy.
  ParticleReal tmpReal   = 0.0; ///< Scratch scalar storage.

  double D_DECL(oldPositionX = 0.0, oldPositionY = 0.0, oldPositionZ = 0.0); ///< Previous position (per component).
  ParticleReal D_DECL(velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0);    ///< Interpolated velocity (per component).
  ParticleReal D_DECL(tmpVectX = 0.0, tmpVectY = 0.0, tmpVectZ = 0.0); ///< Scratch vector storage (per component).
};

/**
  @brief ParticleTraits specialization listing the payload columns of ItoParticle.
  @details Only the velocity columns are written to HDF5 (h5Columns); position and weight are always
  checkpointed by the container.
*/
template <>
struct ParticleTraits<ItoParticle>
{
  static constexpr auto columns = std::make_tuple(
    &ItoParticle::mobility,
    &ItoParticle::diffusion,
    &ItoParticle::energy,
    &ItoParticle::tmpReal,
    D_DECL(&ItoParticle::oldPositionX, &ItoParticle::oldPositionY, &ItoParticle::oldPositionZ),
    D_DECL(&ItoParticle::velocityX, &ItoParticle::velocityY, &ItoParticle::velocityZ),
    D_DECL(&ItoParticle::tmpVectX, &ItoParticle::tmpVectY, &ItoParticle::tmpVectZ));

  static constexpr auto h5Columns = std::make_tuple(
    D_DECL(&ItoParticle::velocityX, &ItoParticle::velocityY, &ItoParticle::velocityZ));
};

#include <CD_NamespaceFooter.H>

#endif
