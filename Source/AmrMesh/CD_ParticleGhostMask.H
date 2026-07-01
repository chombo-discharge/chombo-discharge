/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleGhostMask.H
  @brief  Declaration of ParticleGhostMask, a per-box CSR lookup of particle ghost scatter targets
  @author Robert Marskar
*/

#ifndef CD_PARTICLEGHOSTMASK_H
#define CD_PARTICLEGHOSTMASK_H

// Std includes
#include <vector>

// Chombo includes
#include <CH_assert.H>
#include <Box.H>
#include <BaseFab.H>
#include <BoxIterator.H>
#include <LayoutData.H>
#include <RefCountedPtr.H>
#include <Vector.H>

// Our includes
#include <CD_LevelTiles.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Per-box CSR lookup from a source cell to the destination boxes whose ghosted region covers it.
  @details Each cell of the box maps to zero or more scatter targets, a target being a
  LevelTiles::BoxIDs pair = (grid index on the destination level, rank that owns and will receive that
  box). Storage is compressed-sparse-row: for a cell, m_start holds the offset into the packed target
  array m_flat and m_count the number of targets. A cell may map to several targets (overlapping ghosted
  boxes, e.g. diagonal coarse-fine cells) with no fixed cap; the sending rank is not stored.

  Optionally (coarse-to-fine only) each target also carries an acceptance box: the destination box's
  ghost region measured in DESTINATION cells (grown by the ghost width, periodic-shifted into the source
  cell's index space). It lets the scatter prune whole-coarse-cell over-communication down to the exact
  sub-cell shell -- a candidate particle is scattered only if its destination cell lies in the acceptance
  box. This parallel array m_flatBox is populated only when allocateTargetBoxes()/addTarget(...,box) are
  used; same-level and fine-to-coarse masks leave it empty (no per-target box, no memory cost).

  The build API is a two-pass fill (define -> incrementCount per contribution -> allocate [->
  allocateTargetBoxes] -> addTarget per contribution -> finalize). The query API (numTargets/target/
  targetGridIndex/targetRank[/targetBox]) is inlined so downstream scatter code never touches the CSR
  internals.
*/
class ParticleGhostMask
{
public:
  /**
    @brief Default constructor. Leaves the table undefined; call define() before use.
  */
  ParticleGhostMask() = default;

  /**
    @brief Destructor.
  */
  ~ParticleGhostMask() = default;

  /**
    @brief Reset to an empty table over a box (every cell has zero targets).
    @param[in] a_box Box that the source cells live on.
  */
  inline void
  define(const Box& a_box) noexcept;

  /**
    @brief Build pass 1: register that cell a_iv gains one more target.
    @param[in] a_iv Source cell.
  */
  inline void
  incrementCount(const IntVect& a_iv) noexcept;

  /**
    @brief Build step between the passes: turn the per-cell counts into CSR offsets and size the packed
    target array. Must be called after all incrementCount() and before any addTarget().
  */
  inline void
  allocate() noexcept;

  /**
    @brief Optional build step: size the per-target acceptance-box array to match the packed target array.
    @details Call after allocate() and before any addTarget(...,box), only for masks that carry acceptance
    boxes (coarse-to-fine). Masks that skip it leave m_flatBox empty.
  */
  inline void
  allocateTargetBoxes() noexcept;

  /**
    @brief Build pass 2: append one target for cell a_iv.
    @details Uses m_start as the advancing write head; call exactly numTargets(a_iv) times per cell (the
    same cells/counts as pass 1). finalize() restores m_start afterwards.
    @param[in] a_iv     Source cell (must have been counted in pass 1).
    @param[in] a_target Destination (grid index, receiving rank).
  */
  inline void
  addTarget(const IntVect& a_iv, const LevelTiles::BoxIDs& a_target) noexcept;

  /**
    @brief Build pass 2 with an acceptance box: append one target and its acceptance box for cell a_iv.
    @details Requires allocateTargetBoxes() to have been called. Writes the acceptance box at the same
    packed slot as the target, then delegates to addTarget() to store the target and advance the head.
    @param[in] a_iv        Source cell (must have been counted in pass 1).
    @param[in] a_target    Destination (grid index, receiving rank).
    @param[in] a_targetBox Acceptance box in destination cells (particle's dest cell must lie in it).
  */
  inline void
  addTarget(const IntVect& a_iv, const LevelTiles::BoxIDs& a_target, const Box& a_targetBox) noexcept;

  /**
    @brief Finish the build: rewind m_start (advanced by addTarget) back to the CSR start offsets.
  */
  inline void
  finalize() noexcept;

  /**
    @brief Get the box that the source cells live on.
    @return m_box
  */
  inline const Box&
  box() const noexcept;

  /**
    @brief Get the number of scatter targets for a cell.
    @param[in] a_iv Source cell (must lie in box()).
    @return Number of targets (0 if the cell scatters nowhere in this direction).
  */
  inline int
  numTargets(const IntVect& a_iv) const noexcept;

  /**
    @brief Get the i-th scatter target for a cell.
    @param[in] a_iv Source cell (must lie in box()).
    @param[in] a_i  Target index in [0, numTargets(a_iv)).
    @return (grid index on the destination level, receiving rank).
  */
  inline const LevelTiles::BoxIDs&
  target(const IntVect& a_iv, const int a_i) const noexcept;

  /**
    @brief Get the destination grid index of the i-th target for a cell.
    @param[in] a_iv Source cell (must lie in box()).
    @param[in] a_i  Target index in [0, numTargets(a_iv)).
    @return Grid index on the destination level.
  */
  inline unsigned int
  targetGridIndex(const IntVect& a_iv, const int a_i) const noexcept;

  /**
    @brief Get the receiving rank of the i-th target for a cell.
    @param[in] a_iv Source cell (must lie in box()).
    @param[in] a_i  Target index in [0, numTargets(a_iv)).
    @return Rank that owns (and will receive) the destination box.
  */
  inline unsigned int
  targetRank(const IntVect& a_iv, const int a_i) const noexcept;

  /**
    @brief Get the acceptance box of the i-th target for a cell.
    @details Only valid when the mask was built with acceptance boxes (hasTargetBoxes()).
    @param[in] a_iv Source cell (must lie in box()).
    @param[in] a_i  Target index in [0, numTargets(a_iv)).
    @return Acceptance box in destination cells.
  */
  inline const Box&
  targetBox(const IntVect& a_iv, const int a_i) const noexcept;

  /**
    @brief Whether this mask carries per-target acceptance boxes.
    @return True if allocateTargetBoxes() was used (coarse-to-fine masks), false otherwise.
  */
  inline bool
  hasTargetBoxes() const noexcept;

protected:
  /**
    @brief Box that the source cells (and m_start/m_count) live on.
  */
  Box m_box;

  /**
    @brief Per cell: first index into m_flat, or -1 if the cell has no targets.
  */
  BaseFab<int> m_start;

  /**
    @brief Per cell: number of targets.
  */
  BaseFab<int> m_count;

  /**
    @brief Packed targets, size == sum(m_count).
  */
  std::vector<LevelTiles::BoxIDs> m_flat;

  /**
    @brief Packed per-target acceptance boxes, parallel to m_flat. Empty unless the mask carries them.
  */
  std::vector<Box> m_flatBox;
};

/**
  @brief AMR-wide particle ghost-target lookup: one per-box CSR table per level, indexed by SOURCE level.
*/
using AMRParticleGhostMask = Vector<RefCountedPtr<LayoutData<ParticleGhostMask>>>;

#include <CD_NamespaceFooter.H>

#include <CD_ParticleGhostMaskImplem.H>

#endif
