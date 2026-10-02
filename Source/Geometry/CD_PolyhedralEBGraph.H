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
#include <CD_CutCellFaceOverrides.H>
#include <CD_CutCellSurface.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief The cut cells of one level of the polyhedral geometry, and what they connect to.
 * @details Holds, over the level's cut tiles, the classification of every cell, the surface every cut cell is
 * reconstructed from, and for every face of every cell what lies across it: nothing, a cell of this level, a cell
 * of the finer level, a cell of the coarser level, or the domain boundary. The surface is the source of everything
 * else about a cut cell -- its polyhedron and moments are rebuilt from it on demand -- and it regrids: a coarse
 * surface's corners are fine nodes and its crossings are the children's. A cell whose face lies across a level
 * boundary, or opens onto a covered cell, is not described by its own surface alone: PolyhedralGeometryShop
 * stitches it once, from the stored cells of the level across, writes the crossings it takes from there into the
 * surface and keeps the faces it changes as face overrides, so that the cell's body is rebuilt from its own surface
 * and its own overrides with no call to the implicit function.
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
   * @brief Why a cell was filled: it was not.
   */
  static constexpr int s_notFilled = 0;

  /**
   * @brief Why a cell was filled: its fluid falls into more than one piece, which one body cannot hold.
   */
  static constexpr int s_filledPieces = 1;

  /**
   * @brief Why a cell was filled: its one piece of fluid is the smaller part of it, a gap between pieces of solid.
   */
  static constexpr int s_filledGap = 2;

  /**
   * @brief Why a cell was filled: its fluid reaches no regular cell.
   */
  static constexpr int s_filledSealed = 3;

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
   * from the implicit function, classifies it as the polyhedral generator does, keeping every cell that holds
   * anything -- a body with nothing in it is a regular cell, and nothing else is discarded, since a cell's
   * interface is what closes the surface against its neighbours -- keeps the surface of every cell whose corners
   * disagree, cut or not, and of every cell beside one the finer level carries, and sets the face states this level can
   * know by itself: closed or open from the apertures, open onto the domain boundary, and open onto a cell no tile of
   * this level carries. Collective.
   * @param[in] a_function        Implicit function of the phase.
   * @param[in] a_cutTiles        The level's cut tiles, identical on every rank.
   * @param[in] a_domain          Domain of the level.
   * @param[in] a_probLo          Lower-left corner of the domain.
   * @param[in] a_dx              Grid spacing of the level.
   * @param[in] a_numGhost        Ghost cells the cell states are kept for.
   * @param[in] a_covered         Cells of this level the finer level will carry, coarsened onto this level. A cell
   * in them is described up there and is not filled down here.
   * @param[in] a_coarser         Graph of the coarser level, already defined, or null where there is none. A face onto
   * a cell it filled is closed here, and a cell with such a face is kept cut, as beside a cell this level fills.
   */
  virtual void
  define(const BaseIF&            a_function,
         const Vector<Box>&       a_cutTiles,
         const ProblemDomain&     a_domain,
         const RealVect&          a_probLo,
         const Real               a_dx,
         const int                a_numGhost,
         const Vector<Box>&       a_covered,
         const PolyhedralEBGraph* a_coarser);

  /**
   * @brief Build the graph of a level over a given layout, from another graph of the same level.
   * @details The other graph's cell states, face states and surfaces are copied onto the new layout's valid cells, and
   * exchanged into its ghost cells; the ghost cells no tile of the new layout carries are classified from the implicit
   * function by their corners, as define classifies them. The new layout must cover the same cells as the other
   * graph's, since the surfaces come from it alone; a layout of another cell count stops the run. The refined mask
   * is not copied: it is the finer level's to set, so link the copy to its finer level afterwards as the original
   * was. Nor are the face overrides, which are stitched afresh on the copy. Collective.
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
   * @brief Mark, on the coarse level, the faces that open onto cells the fine level carries.
   * @details The fine tiles are coarsened onto the coarse layout the way Realm builds its valid-cell mask: a
   * marker on the coarsened fine layout is copied onto the coarse layout and its ghost cells. An open face of an
   * unrefined coarse cell onto a cell the marker reaches is described by the fine level. The coarse level already
   * knows which of its cells are refined -- it was told when it was defined -- so the marker is checked against
   * that rather than written over it, and a disagreement between the coarse level's covered region and the fine
   * level's tiles stops the run. The fine level's own faces onto the coarse level were set when it was defined.
   * Collective.
   * @param[in,out] a_coarse Graph of the coarser level.
   * @param[in]     a_fine   Graph of the finer level, a factor-two refinement of it.
   */
  static void
  link(PolyhedralEBGraph& a_coarse, const PolyhedralEBGraph& a_fine);

  /**
   * @brief Spread anchoring through the open faces of this level until it stops.
   * @details Filling a cell closes the faces onto it, and fluid that reached the rest of the fluid only through such
   * faces is left sealed off: beside a filled gap thinner than a cell, the slivers at the gap's corners are joined
   * to one another and to nothing else. No stencil fits in fluid that holds no regular cell, so such fluid is
   * filled. A cut cell is anchored when an open face joins it to an anchored cell; the spread runs box by box with
   * an exchange in between until no box anchors anything new. A face on the domain boundary anchors nothing, and a
   * cell the finer level carries passes nothing through itself: what it holds is the finer level's. Collective.
   * @return Number of cells anchored, over all ranks.
   */
  long long
  spreadAnchors();

  /**
   * @brief Carry anchoring across a level boundary, both ways.
   * @details A coarse cell the fine level carries is anchored when any of its children is, and a cut ghost cell of
   * the fine level that the coarse level describes is anchored when the coarse cell over it is. Where the coarse
   * level holds nothing over such a ghost cell, the ghost cell is anchored. Collective.
   * @param[in,out] a_coarse Graph of the coarser level.
   * @param[in,out] a_fine   Graph of the finer level, a factor-two refinement of it.
   * @return Number of cells newly anchored, over all ranks.
   */
  static long long
  anchorAcrossLevels(PolyhedralEBGraph& a_coarse, PolyhedralEBGraph& a_fine);

  /**
   * @brief Anchor every cut ghost cell no tile carries, for a level with no coarser graph to answer for them.
   */
  void
  anchorUncarriedCells();

  /**
   * @brief Fill every cut cell anchoring did not reach, and let go of the anchoring.
   * @details The cell becomes covered with every face closed and leaves the cut cells; its surface stays, as a
   * filled cell's does. A sealed region has no open face onto anchored fluid, so filling it changes no face of a
   * cell that stays. Run once spreadAnchors and anchorAcrossLevels have settled on every level. Collective.
   * @return Number of cells filled, over all ranks.
   */
  long long
  fillSealedCells();

  /**
   * @brief Find the cells of a region that one body cannot describe, or that are mostly a gap of fluid.
   * @details A cell the surface enters as more than one sheet holds more than one piece of fluid or more than one
   * piece of solid. Two pieces of fluid -- a plate between them -- are more than one body can hold, and the cell is
   * filled, which is the only reading that stays single valued and keeps a barrier a barrier instead of letting
   * the fluid run through it. One piece of fluid with several pieces of solid is either a gap of fluid thinner than
   * the cell between two solid lobes, or pieces that merely sit apart in it, such as two solid tips at opposite
   * corners of a face. The two share their topology and their corners, and a pattern that repeats at every scale --
   * a surface lying in a node plane -- looks the same in every subcell as well, so it is the volume that decides:
   * the cell's body is built from its crossings, and the cell is filled when the fluid is the smaller part of it,
   * erring toward blocking a gap, and kept, as one body with several sheets of interface, when the fluid fills
   * most of it. Fluid pieces are counted the way the body joins them: across a face by its face decision, never
   * through the interior. Filling is a
   * classification and changes no value another level reads, so the children a coarse cell restricts against still
   * agree with it about where every edge is crossed -- which is what a snap of the nodes themselves cannot promise.
   *
   * A cell the finer level carries is left alone: the feature is resolved up there, and this level's reading of
   * it is not the one that gets used. Only the combinatorics are read, so no root is solved for here.
   *
   * Nothing is counted here. The region a caller asks about reaches past the cells it owns, so that a cell can
   * see whether its neighbours were filled, and those reaches overlap between neighbouring boxes; a count taken
   * here would report a cell near a box boundary once for every box that reaches it. The caller counts over the
   * cells it owns.
   * @param[in]  a_function   Implicit function, read at face centres and along the edges of multi-sheet cells.
   * @param[in]  a_nodeValues Node values of the region.
   * @param[in]  a_region     Cells to ask.
   * @param[in]  a_covered    Cells a finer level carries.
   * @param[in]  a_probLo     Lower-left corner of the domain.
   * @param[in]  a_dx         Grid spacing.
   * @param[out] a_unresolved For each cell of the region, why it must be filled -- s_filledPieces or s_filledGap -- or
   * s_notFilled.
   */
  static void
  findUnresolvedCells(const BaseIF&         a_function,
                      const BaseFab<Real>&  a_nodeValues,
                      const Box&            a_region,
                      const Vector<Box>&    a_covered,
                      const RealVect&       a_probLo,
                      const Real            a_dx,
                      BaseFab<signed char>& a_unresolved);

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
   * @brief Why each valid cell was filled: s_notFilled, s_filledPieces, s_filledGap or s_filledSealed. Kept by the
   * graph that defined the cells, not by a copy of it.
   * @return One component per cell, no ghost cells.
   */
  const LevelData<BaseFab<signed char>>&
  getFillReasons() const noexcept;

  /**
   * @brief The surface every cell of the surface sets is reconstructed from, in the box's own cells and in its ghost
   * cells. A stitched level's surfaces carry the crossings its refined edges take from the finer level.
   * @return One surface per cell of every box's surface set.
   */
  const LevelData<IVSFAB<PolyhedralEB::CutCellSurface>>&
  getSurfaces() const noexcept;

  /**
   * @brief The surfaces, for the shop to write the crossings it takes from the finer level into.
   * @return One surface per cell of every box's surface set.
   */
  LevelData<IVSFAB<PolyhedralEB::CutCellSurface>>&
  getSurfaces() noexcept;

  /**
   * @brief The cells of every box that hold a surface: every cut cell of its set, and every cell, own or ghost, whose
   * corners disagree though it was dropped as dust or filled, or that lies beside a cell the finer level carries.
   * @details A cell that is not cut keeps its surface because its crossings lie on edges other cells share, and the
   * coarser level stitches from them; one beside the finer level keeps it because stitching may make it cut.
   * @return One set per box.
   */
  const LayoutData<IntVectSet>&
  getSurfaceCells() const noexcept;

  /**
   * @brief The faces of every box's cut cells, own and ghost, that their surfaces do not describe.
   * @details Empty until PolyhedralGeometryShop has stitched the level.
   * @return One set of overrides per box.
   */
  const LayoutData<PolyhedralEB::CutCellFaceOverrides>&
  getFaceOverrides() const noexcept;

  /**
   * @brief The face overrides, for the shop to write.
   * @return One set of overrides per box.
   */
  LayoutData<PolyhedralEB::CutCellFaceOverrides>&
  getFaceOverrides() noexcept;

  /**
   * @brief Mark the cells some tile of this level carries: one in every valid cell and in every ghost cell another tile
   * covers, zero elsewhere. Collective.
   * @param[out] a_carried The marker, defined here over the layout with the cell states' ghost width.
   */
  void
  markCarried(LevelData<BaseFab<signed char>>& a_carried) const;

  /**
   * @brief Mark the cells of this level, valid and ghost, that lie under a cell the coarser level filled.
   * @details Such a cell is never carried here, since the coarser level fills no cell this level describes, so the
   * marker lands on ghost cells the coarser level answers for. A face of this level onto one of them opens onto
   * nothing and has to be closed from this side. Collective.
   * @param[out] a_filled  The marker, defined here over the layout with a_ghost ghost cells.
   * @param[in]  a_coarser Graph of the coarser level, a factor two coarser.
   * @param[in]  a_ghost   Ghost cells of the marker.
   */
  void
  markCoarserFilled(LevelData<BaseFab<signed char>>& a_filled,
                    const PolyhedralEBGraph&         a_coarser,
                    const int                        a_ghost) const;

  /**
   * @brief Make a regular or covered cell of a box, own or ghost, a cut cell of this level.
   * @details For a cell on a refinement boundary whose face the finer level describes otherwise than its corners do: a
   * regular cell whose face the finer level cuts or closes, or a covered one onto whose face the finer level opens. The
   * cell becomes cut and enters the box's cut set, and every face of it onto a cell the finer level carries opens onto
   * the finer level. Its surface is already kept, since every cell on a refinement boundary keeps one.
   * @param[in] a_dit  Box index.
   * @param[in] a_cell The cell, in the box or its ghost ring.
   */
  void
  promoteToCut(const DataIndex& a_dit, const IntVect& a_cell);

  /**
   * @brief Mark the cells that hold a surface.
   * @details The cell states cannot say which cells that are not cut keep a surface, so a copy of the graph carries
   * this mask across first. Ghost cells are marked from the tiles that own them. Collective.
   * @param[out] a_mask One where a cell holds a surface, zero elsewhere, over the graph's layout and its ghost cells.
   */
  void
  markSurfaces(LevelData<BaseFab<signed char>>& a_mask) const;

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
   * @brief Cells of this level the finer level carries, coarsened onto this level.
   * @details Known before any cell is built, since the tiles are the geometry's and not the graph's. A cell in
   * them is described by the finer level: this level classifies it from its corners but builds no surface or body
   * for it, and it is not one of this level's cut cells.
   */
  Vector<Box> m_covered;

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
   * @brief One where the finer level carries the cell, with at least two ghost cells.
   * @details Set from m_covered when the graph is defined, before any cell is built, so that nothing on this
   * level is built for a cell the finer level describes. link checks it against the finer level's tiles.
   */
  LevelData<BaseFab<signed char>> m_refined;

  /**
   * @brief Why each valid cell was filled, for inspection.
   */
  LevelData<BaseFab<signed char>> m_fillReasons;

  /**
   * @brief Anchoring, from beginAnchoring until fillSealedCells: one where a cell's fluid reaches a regular cell, and
   * in a second component one on the cut ghost cells no tile carries, which the coarser level answers for.
   */
  LevelData<BaseFab<signed char>> m_anchored;

  /**
   * @brief Cells of every box that hold a surface.
   */
  LayoutData<IntVectSet> m_surfaceCells;

  /**
   * @brief Surface of every cell of m_surfaceCells.
   */
  LevelData<IVSFAB<PolyhedralEB::CutCellSurface>> m_surfaces;

  /**
   * @brief Faces of every box's cut cells that their surfaces do not describe.
   */
  LayoutData<PolyhedralEB::CutCellFaceOverrides> m_faceOverrides;

  /**
   * @brief Put into every box's surface set the ghost cells whose owning tile keeps a surface for them.
   * @details Run once the valid cells' sets are complete and before the surfaces are exchanged. Collective.
   */
  void
  addGhostSurfaceCells();

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
   * @brief Record which cells anchor fluid before any of it is spread, for fillSealedCells.
   * @details Regular cells anchor. A cell the finer level carries anchors when the finer level's fluid in it does,
   * and a cut ghost cell no tile carries -- the coarser level describes it -- when the coarser level's fluid there
   * does; both start unanchored and are told by anchorAcrossLevels. A regular ghost cell no tile carries anchors, a
   * covered one does not.
   * @param[in] a_carried The marker from markCarried.
   */
  void
  beginAnchoring(const LevelData<BaseFab<signed char>>& a_carried);

  /**
   * @brief Set m_refined from m_covered, over every box and its ghost cells.
   */
  void
  markRefined();

  /**
   * @brief Reconstruct and classify every cell of every box and its ghost region, keep the cut cells' surfaces with
   * ghost cells, and set the face states the apertures decide.
   * @param[in] a_function      Implicit function.
   * @param[in] a_carried       The marker from markCarried.
   * @param[in] a_coarserFilled The marker from markCoarserFilled, with one ghost cell, or null where there is no
   * coarser level.
   */
  void
  defineCells(const BaseIF&                          a_function,
              const LevelData<BaseFab<signed char>>& a_carried,
              const LevelData<BaseFab<signed char>>* a_coarserFilled);

  /**
   * @brief Mark the open faces that lie in the domain boundary, and those onto cells no tile of this level carries.
   * @param[in] a_carried The marker from markCarried.
   */
  void
  defineOuterFaces(const LevelData<BaseFab<signed char>>& a_carried);
};

#include <CD_NamespaceFooter.H>

#endif
