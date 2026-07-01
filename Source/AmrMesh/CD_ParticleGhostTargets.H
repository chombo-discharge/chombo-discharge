/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleGhostTargets.H
  @brief  Declaration of ParticleGhostTargets, a per-box CSR lookup of particle ghost scatter targets
  @author Robert Marskar
*/

#ifndef CD_PARTICLEGHOSTTARGETS_H
#define CD_PARTICLEGHOSTTARGETS_H

// Std includes
#include <vector>

// Chombo includes
#include <BaseFab.H>
#include <LayoutData.H>
#include <RefCountedPtr.H>
#include <Vector.H>

// Our includes
#include <CD_LevelTiles.H>
#include <CD_NamespaceHeader.H>

/**
  @brief CSR per-box lookup from a source cell to the destination boxes whose ghosted region covers it.
  @details For each cell of the (source) box, m_start gives the offset into m_flat (-1 if the cell has no
  target) and m_count the number of targets; the targets themselves are packed contiguously in m_flat. A
  cell may map to several targets (overlapping ghosted boxes, e.g. diagonal coarse-fine cells), which the
  CSR stores without a fixed cap. Each target is a LevelTiles::BoxIDs pair = (grid index on the target
  level, rank that owns and will receive that box); the sending rank is intentionally not stored.
*/
struct ParticleGhostTargets
{
  BaseFab<int>                    m_start; ///< per cell: first index into m_flat, or -1 if no targets
  BaseFab<int>                    m_count; ///< per cell: number of targets
  std::vector<LevelTiles::BoxIDs> m_flat;  ///< packed targets, size == sum(m_count)
};

/**
  @brief AMR-wide particle ghost-target lookup: one per-box CSR table per level, indexed by SOURCE level.
*/
using AMRParticleGhostTargets = Vector<RefCountedPtr<LayoutData<ParticleGhostTargets>>>;

#include <CD_NamespaceFooter.H>

#endif
