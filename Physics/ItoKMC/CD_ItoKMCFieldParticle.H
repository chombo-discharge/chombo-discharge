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
  @details The energy increase as a particle moves from position A to position B is q*(V(B)-V(A)). vA and vB
  hold the interpolated electrostatic potential at the two positions; alt holds the per-component "other"
  position (A or B) so the single position-based interpolation routine can be reused at both locations.
  Weight and position are container-owned columns, NOT payload members.
*/
struct ItoKMCFieldParticle
{
  ParticleReal vA = 0.0; ///< Electrostatic potential at position A.
  ParticleReal vB = 0.0; ///< Electrostatic potential at position B.

  ParticleReal D_DECL(altX = 0.0, altY = 0.0, altZ = 0.0); ///< Alternate (other) position (per component).
};

/**
  @brief ParticleTraits specialization listing the payload columns of ItoKMCFieldParticle.
*/
template <>
struct ParticleTraits<ItoKMCFieldParticle>
{
  static constexpr auto columns = std::make_tuple(
    &ItoKMCFieldParticle::vA,
    &ItoKMCFieldParticle::vB,
    D_DECL(&ItoKMCFieldParticle::altX, &ItoKMCFieldParticle::altY, &ItoKMCFieldParticle::altZ));
};

#include <CD_NamespaceFooter.H>

#endif
