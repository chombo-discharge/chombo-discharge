/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_AMRPolyhedralEBGraph.H
 * @brief  Declaration of the polyhedral graph as an AMR hierarchy of any layout holds it.
 * @author Robert Marskar
 */

#ifndef CD_AMRPOLYHEDRALEBGRAPH_H
#define CD_AMRPOLYHEDRALEBGRAPH_H

// Std includes
#include <vector>

// Chombo includes
#include <DisjointBoxLayout.H>
#include <LayoutData.H>
#include <ProblemDomain.H>
#include <RefCountedPtr.H>

// Our includes
#include <CD_CutCellDescription.H>
#include <CD_CutCellFaceOverrides.H>
#include <CD_CutCellSurface.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief The polyhedral graph as an AMR hierarchy of any layout holds it: every cell of every box, ghost cells
 * included, described by the geometry level that owns it.
 * @details The geometry's graph lives on its own levels and its own tiles. A consumer lives on its own hierarchy --
 * one DisjointBoxLayout per level, each box filled to a ghost width, less what the next finer level covers -- and a
 * cell of it may be owned by the geometry level at its own resolution, by a finer one, or by a coarser one. Each box
 * keeps the part of the graph that intersects it and nothing else: per geometry level, the pieces it received -- boxes
 * of cells one tile of that level owns, with their states -- and the surfaces and face overrides of the cells among
 * them that hold any. Nothing here is recomputed from the implicit function; PolyhedralGeometryShop::copyGraphs fills
 * it from the graphs alone.
 *
 * A consumer asks a destination cell for the cut cells that intersect it, getCutCells, each as a CutCellDescription
 * from which CutCellBody::define builds its polyhedron, and asks getState what a cell with none is.
 */
class AMRPolyhedralEBGraph
{
public:
  /**
   * @brief State of a destination cell, or of a geometry cell, that no geometry level owns here.
   */
  static constexpr int s_unowned = 2;

  /**
   * @brief Constructor. Holds nothing until define.
   */
  AMRPolyhedralEBGraph() noexcept;

  /**
   * @brief Destructor.
   */
  virtual ~AMRPolyhedralEBGraph() noexcept;

  /**
   * @brief Set the destination hierarchy, holding nothing of the graph yet.
   * @param[in] a_grids           Destination layout of every level, coarsest first.
   * @param[in] a_domains         Domain of every destination level.
   * @param[in] a_geometryLevels  Geometry level of every destination level.
   * @param[in] a_numSourceLevels Number of geometry levels.
   * @param[in] a_ghost           Ghost width every destination box is filled to.
   */
  void
  define(const std::vector<DisjointBoxLayout>& a_grids,
         const std::vector<ProblemDomain>&     a_domains,
         const std::vector<int>&               a_geometryLevels,
         const int                             a_numSourceLevels,
         const int                             a_ghost);

  /**
   * @brief Whether define has run.
   * @return True after define.
   */
  bool
  isDefined() const noexcept;

  /**
   * @brief Number of destination levels.
   * @return The count.
   */
  int
  getNumLevels() const noexcept;

  /**
   * @brief Number of geometry levels.
   * @return The count.
   */
  int
  getNumSourceLevels() const noexcept;

  /**
   * @brief Ghost width every destination box is filled to.
   * @return The width.
   */
  int
  getGhost() const noexcept;

  /**
   * @brief Layout of a destination level.
   * @param[in] a_level Destination level.
   * @return The layout.
   */
  const DisjointBoxLayout&
  getGrids(const int a_level) const noexcept;

  /**
   * @brief Domain of a destination level.
   * @param[in] a_level Destination level.
   * @return The domain.
   */
  const ProblemDomain&
  getDomain(const int a_level) const noexcept;

  /**
   * @brief Geometry level a destination level lives on.
   * @param[in] a_level Destination level.
   * @return The geometry level.
   */
  int
  getGeometryLevel(const int a_level) const noexcept;

  /**
   * @brief The cut cells of the geometry that intersect a cell of a destination box, each with everything its
   * polyhedron is built from.
   * @details A cell the geometry owns at the destination's resolution is its own cut cell; one owned by a finer level
   * is intersected by every cut cell of that level under it; one owned by a coarser level is intersected by the cut
   * ancestor, whose polyhedron the consumer subdivides down to the cell. Each comes with its face overrides. A cell
   * none of whose owners is cut gives none; getState says what it is.
   * @param[out] a_cells The cut cells, appended to an emptied list.
   * @param[in]  a_level Destination level.
   * @param[in]  a_dit   Destination box.
   * @param[in]  a_cell  The cell, on the destination level, in the box or its ghost cells.
   */
  void
  getCutCells(std::vector<PolyhedralEB::CutCellDescription>& a_cells,
              const int                                      a_level,
              const DataIndex&                               a_dit,
              const IntVect&                                 a_cell) const;

  /**
   * @brief What a cell of a destination box is: cut if any geometry cell owning part of it is, regular or covered if
   * all of them are, and s_unowned if no geometry level owns any of it here.
   * @param[in] a_level Destination level.
   * @param[in] a_dit   Destination box.
   * @param[in] a_cell  The cell, on the destination level.
   * @return PolyhedralEBGraph::s_cut, s_regular or s_covered, or s_unowned.
   */
  int
  getState(const int a_level, const DataIndex& a_dit, const IntVect& a_cell) const;

  /**
   * @brief The pieces a destination box holds of a geometry level: boxes of cells, each owned by one of its tiles.
   * @param[in] a_level       Destination level.
   * @param[in] a_dit         Destination box.
   * @param[in] a_sourceLevel Geometry level.
   * @return The pieces, in the geometry level's index space.
   */
  const std::vector<Box>&
  getPieces(const int a_level, const DataIndex& a_dit, const int a_sourceLevel) const noexcept;

  /**
   * @brief The state of a geometry cell, as a destination box holds it.
   * @param[in] a_level       Destination level.
   * @param[in] a_dit         Destination box.
   * @param[in] a_sourceLevel Geometry level of the cell.
   * @param[in] a_sourceCell  The cell, on that level.
   * @return Its state, or s_unowned if the box holds no piece with it.
   */
  int
  getSourceState(const int a_level, const DataIndex& a_dit, const int a_sourceLevel, const IntVect& a_sourceCell) const;

  /**
   * @brief The geometry cells of a level that a destination box holds a surface for, in increasing order.
   * @param[in] a_level       Destination level.
   * @param[in] a_dit         Destination box.
   * @param[in] a_sourceLevel Geometry level.
   * @return The cells.
   */
  const std::vector<IntVect>&
  getSurfaceCells(const int a_level, const DataIndex& a_dit, const int a_sourceLevel) const noexcept;

  /**
   * @brief The surface of a geometry cell, as a destination box holds it.
   * @param[in] a_level       Destination level.
   * @param[in] a_dit         Destination box.
   * @param[in] a_sourceLevel Geometry level of the cell.
   * @param[in] a_sourceCell  The cell, on that level.
   * @return The surface, or null if the box holds none for the cell.
   */
  const PolyhedralEB::CutCellSurface*
  getSourceSurface(const int        a_level,
                   const DataIndex& a_dit,
                   const int        a_sourceLevel,
                   const IntVect&   a_sourceCell) const;

  /**
   * @brief The face overrides of a geometry level's cells, as a destination box holds them.
   * @param[in] a_level       Destination level.
   * @param[in] a_dit         Destination box.
   * @param[in] a_sourceLevel Geometry level.
   * @return The overrides, keyed by cells of that level.
   */
  const PolyhedralEB::CutCellFaceOverrides&
  getFaceOverrides(const int a_level, const DataIndex& a_dit, const int a_sourceLevel) const noexcept;

  /**
   * @brief Add a piece to a destination box.
   * @param[in] a_level       Destination level.
   * @param[in] a_dit         Destination box.
   * @param[in] a_sourceLevel Geometry level of the piece.
   * @param[in] a_region      The piece, on that level.
   * @param[in] a_states      Its cells' states, in the order a BoxIterator meets them.
   */
  void
  addPiece(const int          a_level,
           const DataIndex&   a_dit,
           const int          a_sourceLevel,
           const Box&         a_region,
           const signed char* a_states);

  /**
   * @brief Add the surface of a geometry cell to a destination box. Surfaces may come in any order; finish puts them in
   * order before any is read.
   * @param[in] a_level       Destination level.
   * @param[in] a_dit         Destination box.
   * @param[in] a_sourceLevel Geometry level of the cell.
   * @param[in] a_sourceCell  The cell, on that level.
   * @param[in] a_surface     Its surface.
   */
  void
  addSurface(const int                           a_level,
             const DataIndex&                    a_dit,
             const int                           a_sourceLevel,
             const IntVect&                      a_sourceCell,
             const PolyhedralEB::CutCellSurface& a_surface);

  /**
   * @brief Put every box's surfaces in cell order and index its pieces, after the last addPiece and addSurface and
   * before anything is read.
   * @details The index is a coarse grid over each box's pieces, of bins as wide as its widest piece, each listing the
   * pieces that reach it: a cell is found among the few pieces of its bin rather than among all of a box's pieces,
   * which a coarse destination box gathering a finer level's tiles has many of.
   */
  void
  finish();

  /**
   * @brief The face overrides of a geometry level's cells in a destination box, to be filled.
   * @param[in] a_level       Destination level.
   * @param[in] a_dit         Destination box.
   * @param[in] a_sourceLevel Geometry level.
   * @return The overrides.
   */
  PolyhedralEB::CutCellFaceOverrides&
  getFaceOverrides(const int a_level, const DataIndex& a_dit, const int a_sourceLevel) noexcept;

protected:
  /**
   * @brief Whether define has run.
   */
  bool m_isDefined;

  /**
   * @brief Ghost width every destination box is filled to.
   */
  int m_ghost;

  /**
   * @brief Number of geometry levels.
   */
  int m_numSourceLevels;

  /**
   * @brief Destination layout of every level.
   */
  std::vector<DisjointBoxLayout> m_grids;

  /**
   * @brief Domain of every destination level.
   */
  std::vector<ProblemDomain> m_domains;

  /**
   * @brief Geometry level of every destination level.
   */
  std::vector<int> m_geometryLevels;

  /**
   * @brief The pieces of every box, per destination level and geometry level.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<std::vector<Box>>>>> m_pieces;

  /**
   * @brief Where each piece's states begin in m_states, and one past the last piece's end.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<std::vector<int>>>>> m_pieceOffsets;

  /**
   * @brief The pieces' states, one piece after another, each in the order a BoxIterator meets its cells.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<std::vector<signed char>>>>> m_states;

  /**
   * @brief The cells that hold a surface, in increasing order once finish has run.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<std::vector<IntVect>>>>> m_surfaceCells;

  /**
   * @brief Their surfaces, in the same order.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<std::vector<PolyhedralEB::CutCellSurface>>>>> m_surfaces;

  /**
   * @brief The face overrides of the cells whose surfaces are overridden.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<PolyhedralEB::CutCellFaceOverrides>>>> m_faceOverrides;

  /**
   * @brief Width of the bins each box's pieces are indexed by.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<int>>>> m_binSize;

  /**
   * @brief The bins of each box's index, on the geometry level coarsened by the bin width.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<Box>>>> m_binBox;

  /**
   * @brief Where each bin's pieces begin in m_binPieces, the bins in the order a BoxIterator meets them, and one past
   * the last bin's end.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<std::vector<int>>>>> m_binOffsets;

  /**
   * @brief The pieces that reach each bin, one bin after another.
   */
  std::vector<std::vector<RefCountedPtr<LayoutData<std::vector<int>>>>> m_binPieces;

  /**
   * @brief The pieces of a box that may meet a region: those of the bins the region reaches, each once.
   * @param[out] a_pieces      Indices into the box's pieces.
   * @param[in]  a_level       Destination level.
   * @param[in]  a_dit         Destination box.
   * @param[in]  a_sourceLevel Geometry level.
   * @param[in]  a_region      Region on that level.
   */
  void
  piecesNear(std::vector<int>& a_pieces,
             const int         a_level,
             const DataIndex&  a_dit,
             const int         a_sourceLevel,
             const Box&        a_region) const;

  /**
   * @brief Append the description of a cut geometry cell this box holds.
   * @param[out] a_cells       Descriptions to append to.
   * @param[in]  a_level       Destination level.
   * @param[in]  a_dit         Destination box.
   * @param[in]  a_sourceLevel Geometry level of the cell.
   * @param[in]  a_sourceCell  The cell, on that level.
   */
  void
  describe(std::vector<PolyhedralEB::CutCellDescription>& a_cells,
           const int                                      a_level,
           const DataIndex&                               a_dit,
           const int                                      a_sourceLevel,
           const IntVect&                                 a_sourceCell) const;
};

#include <CD_NamespaceFooter.H>

#endif
