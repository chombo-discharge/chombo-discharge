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

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief SoA payload for ItoSolver particles, i.e. drifting Brownian walkers.
  @details mobility/diffusion/energy are the particle transport coefficients and average energy; scratch is
  scratch scalar storage. old_x/old_y/old_z holds the previous particle position (double precision, because it
  subtracts from the double-precision container position); velocity is the interpolated drift velocity;
  tmpVect is scratch vector storage. The vector quantities carry one member per spatial dimension (the
  z-component exists only when CH_SPACEDIM == 3). Weight and position are container-owned columns, NOT
  payload members. Only the velocity columns are checkpointed (position and weight are always checkpointed
  by the container).
*/
struct ItoParticle
{
  ParticleReal mobility  = 0.0; ///< Mobility coefficient.
  ParticleReal diffusion = 0.0; ///< Diffusion coefficient.
  ParticleReal energy    = 0.0; ///< Average particle energy.
  ParticleReal scratch   = 0.0; ///< Scratch scalar storage.

  double old_x = 0.0; ///< Previous position, x-component.
  double old_y = 0.0; ///< Previous position, y-component.
#if CH_SPACEDIM == 3
  double old_z = 0.0; ///< Previous position, z-component.
#endif

  ParticleReal vx = 0.0; ///< Interpolated velocity, x-component.
  ParticleReal vy = 0.0; ///< Interpolated velocity, y-component.
#if CH_SPACEDIM == 3
  ParticleReal vz = 0.0; ///< Interpolated velocity, z-component.
#endif

  ParticleReal scratch_x = 0.0; ///< Scratch vector storage, x-component.
  ParticleReal scratch_y = 0.0; ///< Scratch vector storage, y-component.
#if CH_SPACEDIM == 3
  ParticleReal scratch_z = 0.0; ///< Scratch vector storage, z-component.
#endif
};

/**
  @brief ParticleTraits specialization listing the payload columns of ItoParticle.
  @details Only the velocity columns are written to HDF5 (h5Columns); position and weight are always
  checkpointed by the container.
*/
template <>
struct ParticleTraits<ItoParticle>
{
  /**
     @brief Payload column member-pointers, listed one field group at a time (SpaceDim entries per vector group).
  */
  static constexpr auto columns = std::make_tuple(
    &ItoParticle::mobility,
    &ItoParticle::diffusion,
    &ItoParticle::energy,
    &ItoParticle::scratch,
    D_DECL(&ItoParticle::old_x, &ItoParticle::old_y, &ItoParticle::old_z),
    D_DECL(&ItoParticle::vx, &ItoParticle::vy, &ItoParticle::vz),
    D_DECL(&ItoParticle::scratch_x, &ItoParticle::scratch_y, &ItoParticle::scratch_z));

  /**
     @brief Subset of payload columns written to HDF5 checkpoints: the per-component velocity.
  */
  static constexpr auto h5Columns = std::make_tuple(D_DECL(&ItoParticle::vx, &ItoParticle::vy, &ItoParticle::vz));

  /**
     @brief Payload columns written to H5Part visualization files: the transport coefficients, the average
            energy and the velocity vector (id, position and weight are written automatically; the scratch
            and previous-position columns are internal and not plotted).
  */
  static constexpr auto h5PartColumns = std::make_tuple(
    H5Part::Scalar{"mobility", &ItoParticle::mobility},
    H5Part::Scalar{"diffusion", &ItoParticle::diffusion},
    H5Part::Scalar{"energy", &ItoParticle::energy},
    H5Part::Vector{"v", D_DECL(&ItoParticle::vx, &ItoParticle::vy, &ItoParticle::vz)});
};

#include <CD_NamespaceFooter.H>

#endif
