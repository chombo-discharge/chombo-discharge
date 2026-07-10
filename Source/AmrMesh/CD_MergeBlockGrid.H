/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MergeBlockGrid.H
  @brief  Declaration of MergeBlockGrid
  @author Robert Marskar
*/

#ifndef CD_MERGEBLOCKGRID_H
#define CD_MERGEBLOCKGRID_H

// Chombo includes
#include <IntVect.H>
#include <RealVect.H>
#include <Box.H>
#include <DisjointBoxLayout.H>
#include <DataIndex.H>
#include <LayoutData.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Class for storing, per patch, the static block-index-space geometry a nearest-neighbor
  particle-merge spatial-index walk needs for AABB pruning.
  @details This is a purely geometric structure (patch box -> block-index-space bounds plus a
  1-block halo), built once per regrid and never touched again until the next one -- the same
  lifecycle as LevelTiles, which this class is deliberately modeled on (weak-construct + define(),
  stored/aliased the same way). It answers a DIFFERENT question than LevelTiles does, though:
  LevelTiles maps a tile to its OWNING box/rank (a cross-patch ownership/routing question, used by
  ParticleContainer::findDestination()); this class answers a purely-local geometric question for
  whichever patch a nearest-neighbor search is currently walking ("what are my own block-index-space
  bounds, plus a 1-block halo, and what is a given block's physical AABB"). Kept as its own type
  rather than folded into LevelTiles to avoid overloading that class's established meaning.
  @note The block size is always the AMR hierarchy's own minimum block size (see
  AmrMesh::getMinBlockSize()) -- not an independent tunable -- since patches are already guaranteed
  to be unions of minBlockSize-aligned boxes (the same invariant LevelTiles itself relies on), so
  coarsening a patch's box by that size always divides it with no ragged remainder.
*/
class MergeBlockGrid
{
public:
  /**
    @brief Per-patch block-index-space geometry.
  */
  struct PatchBlocks
  {
    /**
      @brief This patch's own valid box, in cell-index space.
    */
    Box cellBox;

    /**
      @brief This patch's own box, in block-index space (coarsen(cellBox, blockSize)) -- no halo.
    */
    Box blockBox;

    /**
      @brief blockBox grown by one block in every direction.
      @details A nearest-neighbor walk over this patch's own gathered particles (native valids plus
      every ghost already shipped in, per the merge algorithm's own ghost-fill) never needs to look
      further than one block outside this patch's own block range: the registered ghost width is
      one CELL, and the minimum block size is always at least that wide. Given this feature's
      history of subtle geometric-assumption bugs, callers relying on this bound should still
      defensively assert a computed block index never falls outside this box rather than trusting
      it silently.
    */
    Box haloBlockBox;
  };

  /**
    @brief Weak constructor. Must subsequently call the define function.
  */
  MergeBlockGrid() noexcept;

  /**
    @brief Full constructor -- puts object in usable state.
    @param[in] a_dbl        Grids.
    @param[in] a_blockSize  Block size (see AmrMesh::getMinBlockSize()).
  */
  MergeBlockGrid(const DisjointBoxLayout& a_dbl, const int a_blockSize) noexcept;

  /**
    @brief Destructor (does nothing).
  */
  virtual ~MergeBlockGrid() noexcept;

  /**
    @brief Define function. Puts object in usable state.
    @param[in] a_dbl        Grids.
    @param[in] a_blockSize  Block size (see AmrMesh::getMinBlockSize()).
  */
  virtual void
  define(const DisjointBoxLayout& a_dbl, const int a_blockSize) noexcept;

  /**
    @brief Get the block size this object was defined with.
    @return Block size.
  */
  virtual int
  getBlockSize() const noexcept;

  /**
    @brief Get the block-index-space geometry for one of this rank's own patches.
    @param[in] a_din DataIndex of the patch.
    @return This patch's PatchBlocks.
  */
  virtual const PatchBlocks&
  getPatchBlocks(const DataIndex& a_din) const noexcept;

  /**
    @brief Component-wise floor-division block index for a cell.
    @details Pure arithmetic, translation-invariant -- no lookup needed, hence static. Duplicated
    (not shared) in CD_NearestNeighborParticleMergeImplem.H's own nnMergeBlockIndexOf(), since that
    header deliberately stays AmrMesh-agnostic and does not depend on this class -- see this
    class's own docs and that header's architectural-boundary note. Keep the two in sync if this
    formula ever changes.
    @param[in] a_cell      Cell index.
    @param[in] a_blockSize Block size in cells.
    @return floor(a_cell / a_blockSize), component-wise.
  */
  static IntVect
  blockIndexOf(const IntVect& a_cell, const int a_blockSize) noexcept;

  /**
    @brief Physical axis-aligned bounding box of a block.
    @param[in]  a_blockIndex Block index (see blockIndexOf()).
    @param[in]  a_blockSize  Block size in cells.
    @param[in]  a_probLo     Physical position of the domain's low corner.
    @param[in]  a_dx         Cell spacing (per direction).
    @param[out] a_lo         Block's physical low corner.
    @param[out] a_hi         Block's physical high corner.
  */
  static void
  blockAABB(const IntVect&  a_blockIndex,
            const int       a_blockSize,
            const RealVect& a_probLo,
            const RealVect& a_dx,
            RealVect&       a_lo,
            RealVect&       a_hi) noexcept;

protected:
  /**
    @brief Is defined or not.
  */
  bool m_isDefined;

  /**
    @brief Block size this object was defined with.
  */
  int m_blockSize;

  /**
    @brief Per-patch block-index-space geometry, one entry per this rank's own patch.
  */
  LayoutData<PatchBlocks> m_patchBlocks;
};

#include <CD_NamespaceFooter.H>

#endif
