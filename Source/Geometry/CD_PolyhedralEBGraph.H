/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PolyhedralEBGraph.H
 * @brief  Declaration of the cut-cell connectivity of one level of the polyhedral geometry.
 * @author Robert Marskar
 */

#ifndef CD_POLYHEDRALEBGRAPH_H
#define CD_POLYHEDRALEBGRAPH_H

// Chombo includes
#include <BaseFab.H>
#include <BaseIF.H>
#include <DisjointBoxLayout.H>
#include <IVSFAB.H>
#include <IntVectSet.H>
#include <LayoutData.H>
#include <LevelData.H>
#include <ProblemDomain.H>
#include <RealVect.H>

// Our includes
#include <CD_CutCellSurface.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief The cut cells of one level of the polyhedral geometry, and what they connect to.
 * @details Holds, over the level's cut tiles, the classification of every cell, the surface every cut cell is
 * reconstructed from, and for every face of every cell what lies across it: nothing, a cell of this level, a cell
 * of the finer level, a cell of the coarser level, or the domain boundary. The surface is the source of everything
 * else about a cut cell -- its polyhedron and moments are rebuilt from it on demand -- and it regrids: a coarse
 * surface's corners are fine nodes and its crossings are the children's, and children are reconstructed from the
 * implicit function. A cell whose face lies across a level boundary is not described by its own surface alone:
 * the face is described by the level across it, which is what the face state records.
 *
 * The graph covers the cut tiles only. A cell outside every tile is regular or covered, never cut, and which of
 * the two is the geometry's box classification rather than the graph's business.
 *
 * Built in two steps: define makes one level from its implicit function and its tiles, and knows about the
 * domain boundary and about faces that open onto no tile of its own level; link, given two adjacent levels,
 * marks the coarse cells the fine level carries and the coarse faces that open onto them.
 */
class PolyhedralEBGraph
{
public:
  /**
   * @brief Cell states, in the convention the polyhedral generator hands the index space.
   */
  static constexpr int s_covered = -1;

  /**
   * @brief Cell states, in the convention the polyhedral generator hands the index space.
   */
  static constexpr int s_cut = 0;

  /**
   * @brief Cell states, in the convention the polyhedral generator hands the index space.
   */
  static constexpr int s_regular = 1;

  /**
   * @brief Face states: the face is closed, no fluid crosses it.
   */
  static constexpr int s_faceClosed = 0;

  /**
   * @brief Face states: the face opens onto a cell of this level.
   */
  static constexpr int s_faceSameLevel = 1;

  /**
   * @brief Face states: the face opens onto cells the finer level carries, and is described by them.
   */
  static constexpr int s_faceFiner = 2;

  /**
   * @brief Face states: the face opens onto a cell no tile of this level carries, which the coarser level does.
   */
  static constexpr int s_faceCoarser = 3;

  /**
   * @brief Face states: the face lies in the domain boundary and is open.
   */
  static constexpr int s_faceBoundary = 4;

  /**
   * @brief Constructor. Leaves the graph undefined.
   */
  PolyhedralEBGraph();

  /**
   * @brief Destructor.
   */
  virtual ~PolyhedralEBGraph();

  /**
   * @brief Build the graph of one level over its cut tiles.
   * @details Load balances the tiles by cell count, reconstructs every cell of every tile and its ghost region
   * from the implicit function, classifies it as the polyhedral generator does (including the volume threshold
   * and the dust rule), keeps the surface of every cut cell, and sets the face states this level can know by
   * itself: closed or open from the apertures, open onto the domain boundary, and open onto a cell no tile of
   * this level carries. Collective.
   * @param[in] a_function        Implicit function of the phase.
   * @param[in] a_cutTiles        The level's cut tiles, identical on every rank.
   * @param[in] a_domain          Domain of the level.
   * @param[in] a_probLo          Lower-left corner of the domain.
   * @param[in] a_dx              Grid spacing of the level.
   * @param[in] a_numGhost        Ghost cells the cell states are kept for.
   * @param[in] a_volumeThreshold Volume fraction below which a cut cell is covered, and the dust rule's threshold.
   */
  virtual void
  define(const BaseIF&        a_function,
         const Vector<Box>&   a_cutTiles,
         const ProblemDomain& a_domain,
         const RealVect&      a_probLo,
         const Real           a_dx,
         const int            a_numGhost,
         const Real           a_volumeThreshold);

  /**
   * @brief Build the graph of a level over a given layout, from another graph of the same level.
   * @details The other graph's cell states, face states and surfaces are copied onto the new layout's valid cells, and
   * exchanged into its ghost cells; the ghost cells no tile of the new layout carries are classified from the implicit
   * function by their corners, as define classifies them. The new layout must cover the same cells as the other
   * graph's, since the surfaces of the cut cells come from it alone; a layout of another cell count stops the run. The
   * refined mask is not copied: it is the finer level's to set, so link the copy to its finer level afterwards as the
   * original was. Collective.
   * @param[in] a_source   Graph to copy from.
   * @param[in] a_grids    Layout to copy onto, covering the same cells as a_source's.
   * @param[in] a_function Implicit function of the phase, for the ghost cells no tile carries.
   */
  virtual void
  define(const PolyhedralEBGraph& a_source, const DisjointBoxLayout& a_grids, const BaseIF& a_function);

  /**
   * @brief Whether another graph over the same layout holds the same cells, faces, mask and surfaces, bit for bit.
   * @details Compares valid and ghost cells alike, on every rank; the answer is reduced so that every rank gets
   * the same one. Collective.
   * @param[in] a_other Graph over the same layout.
   * @return True if nothing differs.
   */
  bool
  equals(const PolyhedralEBGraph& a_other) const;

  /**
   * @brief Mark, on the coarse level, the cells the fine level carries and the faces that open onto them.
   * @details The fine tiles are coarsened onto the coarse layout the way Realm builds its valid-cell mask: a
   * marker on the coarsened fine layout is copied onto the coarse layout and its ghost cells. A coarse cell the
   * marker reaches is refined; an open face of an unrefined coarse cell onto a refined neighbour is described by
   * the fine level. The fine level's own faces onto the coarse level were set when it was defined. Collective.
   * @param[in,out] a_coarse Graph of the coarser level.
   * @param[in]     a_fine   Graph of the finer level, a factor-two refinement of it.
   */
  static void
  link(PolyhedralEBGraph& a_coarse, const PolyhedralEBGraph& a_fine);

  /**
   * @brief Whether define has run.
   * @return True after define.
   */
  bool
  isDefined() const noexcept;

  /**
   * @brief The level's domain.
   * @return The domain.
   */
  const ProblemDomain&
  getDomain() const noexcept;

  /**
   * @brief The level's grid spacing.
   * @return The spacing.
   */
  Real
  getDx() const noexcept;

  /**
   * @brief The layout over the cut tiles.
   * @return The layout.
   */
  const DisjointBoxLayout&
  getGrids() const noexcept;

  /**
   * @brief The cut cells of every box, its own and those of its ghost region that some tile carries.
   * @return One set per box.
   */
  const LayoutData<IntVectSet>&
  getCutCells() const noexcept;

  /**
   * @brief The state of every cell, s_covered, s_cut or s_regular, with ghost cells filled from the geometry. A signed
   * char per cell, which is what the three values need. A ghost cell no tile carries is classified from its corners
   * alone and has no surface even when cut; the face onto it is s_faceCoarser.
   * @return One component per cell.
   */
  const LevelData<BaseFab<signed char>>&
  getCellStates() const noexcept;

  /**
   * @brief The state of every face of every cell, one of the s_face constants; component 2 * dir + side.
   * @return 2 * SpaceDim components per cell.
   */
  const LevelData<BaseFab<signed char>>&
  getFaceStates() const noexcept;

  /**
   * @brief Whether the finer level carries a cell, one where it does, with two ghost cells; all zero until link has
   * run.
   * @return One component per cell.
   */
  const LevelData<BaseFab<signed char>>&
  getRefinedMask() const noexcept;

  /**
   * @brief The surface every cut cell is reconstructed from, in the box's own cells and in its ghost cells.
   * @return One surface per cut cell of every box's set.
   */
  const LevelData<IVSFAB<PolyhedralEB::CutCellSurface>>&
  getSurfaces() const noexcept;

protected:
  /**
   * @brief Whether define has run.
   */
  bool m_isDefined;

  /**
   * @brief Ghost cells the cell states are kept for.
   */
  int m_numGhost;

  /**
   * @brief Grid spacing of the level.
   */
  Real m_dx;

  /**
   * @brief Lower-left corner of the domain.
   */
  RealVect m_probLo;

  /**
   * @brief Domain of the level.
   */
  ProblemDomain m_domain;

  /**
   * @brief Layout over the cut tiles.
   */
  DisjointBoxLayout m_grids;

  /**
   * @brief Cut cells of every box.
   */
  LayoutData<IntVectSet> m_cutCells;

  /**
   * @brief State of every cell, with ghost cells.
   */
  LevelData<BaseFab<signed char>> m_cellStates;

  /**
   * @brief State of every face of every cell, component 2 * dir + side.
   */
  LevelData<BaseFab<signed char>> m_faceStates;

  /**
   * @brief One where the finer level carries the cell.
   */
  LevelData<BaseFab<signed char>> m_refined;

  /**
   * @brief Surface of every cut cell.
   */
  LevelData<IVSFAB<PolyhedralEB::CutCellSurface>> m_surfaces;

  /**
   * @brief Build the layout over the tiles, load balanced by cell count.
   * @param[in] a_cutTiles The tiles.
   */
  void
  defineGrids(const Vector<Box>& a_cutTiles);

  /**
   * @brief Allocate the per-cell containers over m_grids: the cut-cell sets, the cell states, the face states and the
   * refined mask. The surfaces wait for the sets.
   */
  void
  defineData();

  /**
   * @brief Classify the ghost cells no tile of this level carries from their corners, and take the ghost cut cells some
   * tile carries into each box's cut-cell set.
   * @details The cells another tile carries are expected to hold that tile's state already, by exchange.
   * @param[in] a_function Implicit function.
   * @param[in] a_carried  The marker from markCarried.
   */
  void
  defineGhostCells(const BaseIF& a_function, const LevelData<BaseFab<signed char>>& a_carried);

  /**
   * @brief Mark the cells some tile of this level carries: one in every valid cell and in every ghost cell another tile
   * covers, zero elsewhere.
   * @param[out] a_carried The marker, defined here over the layout with the cell states' ghost width.
   */
  void
  markCarried(LevelData<BaseFab<signed char>>& a_carried) const;

  /**
   * @brief Reconstruct and classify every cell of every box and its ghost region, keep the cut cells' surfaces with
   * ghost cells, and set the face states the apertures decide.
   * @param[in] a_function        Implicit function.
   * @param[in] a_volumeThreshold Volume-fraction threshold.
   * @param[in] a_carried         The marker from markCarried.
   */
  void
  defineCells(const BaseIF& a_function, const Real a_volumeThreshold, const LevelData<BaseFab<signed char>>& a_carried);

  /**
   * @brief Mark the open faces that lie in the domain boundary, and those onto cells no tile of this level carries.
   * @param[in] a_carried The marker from markCarried.
   */
  void
  defineOuterFaces(const LevelData<BaseFab<signed char>>& a_carried);
};

#include <CD_NamespaceFooter.H>

#endif
