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
  @details The merge algorithms live in two header pairs split by dependency weight (see the
  ParticleManagement namespace docs in CD_ParticleManagement.H): the lightweight per-cell mergers in
  CD_ParticleManagement.H, and the AMR-granularity distributed nearest-neighbor merge in
  CD_NearestNeighborParticleMerge.H (which additionally pulls in EBGeometry, AmrMesh, and MPI).
  Include THIS header only when a translation unit genuinely needs both families; a consumer that
  uses only the per-cell mergers should include CD_ParticleManagement.H directly to avoid the heavy
  transitive includes.
*/

#endif
