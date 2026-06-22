/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_DemoParticle.H
  @brief  A demonstration payload type used to exercise ParticleSoA.
  @author Robert Marskar
  @details DemoPayload is the EXTRA-fields struct of a drifting computational
  super-particle: a drift velocity, a mobility, and a derived conductivity. Position,
  weight, and the id/rank metadata are NOT declared here -- ParticleSoA owns those
  columns. This is a concrete example of the recipe documented on ParticleTraits: a
  plain struct of trivially-copyable members plus a ParticleTraits specialization that
  names the payload columns by member pointer. Note the velocity is stored PER-COMPONENT
  (vx/vy/vz scalar columns, guarded by D_DECL) rather than as a RealVect column: vector
  fields always go per-component so the field-update kernels vectorize -- promote to a
  RealVect in downstream code if you need one.
*/

#ifndef CD_DEMOPARTICLE_H
#define CD_DEMOPARTICLE_H

// Std includes
#include <tuple>

// Chombo includes
#include <REAL.H>
#include <RealVect.H>

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Demonstration payload: drift velocity, mobility, and a conductivity field filled
         by SoA kernels. (Position + weight + id/rank are container-owned, not declared here.)
*/
struct DemoPayload
{
  /** @brief Drift velocity, stored per-component (vx/vy/vz scalar columns). */
  ParticleReal D_DECL(vx = 0.0, vy = 0.0, vz = 0.0);
  ParticleReal mobility     = 0.0; ///< Mobility coefficient
  ParticleReal conductivity = 0.0; ///< weight * mobility (filled by a kernel)
};

/**
  @brief ParticleTraits specialization describing the payload columns of DemoPayload.
*/
template <>
struct ParticleTraits<DemoPayload>
{
  /** @brief One pointer-to-member per payload column (single source of truth for the layout). */
  static constexpr auto columns = std::make_tuple(D_DECL(&DemoPayload::vx, &DemoPayload::vy, &DemoPayload::vz),
                                                  &DemoPayload::mobility,
                                                  &DemoPayload::conductivity);

  /** @brief Persist velocity + mobility to HDF5; conductivity is recomputed on restart. */
  static constexpr auto h5Columns = std::make_tuple(D_DECL(&DemoPayload::vx, &DemoPayload::vy, &DemoPayload::vz),
                                                    &DemoPayload::mobility);
};

#include <CD_NamespaceFooter.H>

#endif
