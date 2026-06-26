/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoKMCFieldParticle.H
  @brief  SoA payload for the transient E-dot-J energy computation in ItoKMCStepper
  @author Robert Marskar
*/

#ifndef CD_ITOKMCFIELDPARTICLE_H
#define CD_ITOKMCFIELDPARTICLE_H

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief SoA payload for the transient particles used by ItoKMCStepper::computeEdotJSource.
  @details The energy increase as a particle moves from position A to position B is q*(V(B)-V(A)). phiA and phiB
  hold the interpolated electrostatic potential at the two positions; x0 holds the "other"
  position (A or B), one member per spatial dimension (the z-component exists only when CH_SPACEDIM == 3),
  so the single position-based interpolation routine can be reused at both locations.
  Weight and position are container-owned columns, NOT payload members.
*/
struct ItoKMCFieldParticle
{
  ParticleReal phiA = 0.0; ///< Electrostatic potential at position A.
  ParticleReal phiB = 0.0; ///< Electrostatic potential at position B.

  double x0_x = 0.0; ///< Alternate (other) position, x-component (double: position-like).
  double x0_y = 0.0; ///< Alternate (other) position, y-component (double: position-like).
#if CH_SPACEDIM == 3
  double x0_z = 0.0; ///< Alternate (other) position, z-component (double: position-like).
#endif
};

/**
  @brief ParticleTraits specialization listing the payload columns of ItoKMCFieldParticle.
*/
template <>
struct ParticleTraits<ItoKMCFieldParticle>
{
  /// @brief Payload column member-pointers: phiA, phiB, then the per-component alternate position.
  static constexpr auto columns = std::make_tuple(
    &ItoKMCFieldParticle::phiA,
    &ItoKMCFieldParticle::phiB,
    D_DECL(&ItoKMCFieldParticle::x0_x, &ItoKMCFieldParticle::x0_y, &ItoKMCFieldParticle::x0_z));
};

#include <CD_NamespaceFooter.H>

#endif
