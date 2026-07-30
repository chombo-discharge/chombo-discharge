/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleMergers.H
  @brief  Umbrella header pulling in every ParticleManagement super-particle merge algorithm.
  @author Robert Marskar
*/

#ifndef CD_PARTICLEMERGERS_H
#define CD_PARTICLEMERGERS_H

// Our includes
#include <CD_ParticleManagement.H>
#include <CD_NearestNeighborParticleMerge.H>

/**
  @brief Convenience header that makes ALL ParticleManagement merge algorithms available in one
  include.
  @details Include THIS header only when a translation unit needs all merge families; a consumer
  that uses only the per-cell mergers should include CD_ParticleManagement.H directly to avoid the
  heavy transitive includes caused by some of the other merge classes.
*/

#endif
