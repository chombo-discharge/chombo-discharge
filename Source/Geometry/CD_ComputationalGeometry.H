/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ComputationalGeometry.H
 * @brief  Declaration of base class for defining geometries.
 * @author Robert Marskar
 */

#ifndef CD_COMPUTATIONALGEOMETRY_H
#define CD_COMPUTATIONALGEOMETRY_H

// Std includes
#include <string>

// Chombo includes
#include <BaseIF.H>
#include <GeometryService.H>
#include <MFIndexSpace.H>
#include <Box.H>
#include <RealVect.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_Dielectric.H>
#include <CD_Electrode.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_RealBox.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Abstract base class for geometries
 * @details This class encapsulates computational geometries in chombo-discharge. If you construct this object as-is,
 * you will get a blank geometry. To include EBs one must set the electrodes and dielectrics. This is not a pure
 * function, so you can set those objects directly from a ComputationalGeometry object. However, in almost all cases one
 * will want to derive from ComputationalGeometry and create a parametrized geometry (that is what
 * $DISCHARGE_HOME/Geometries is for!).
 */
class ComputationalGeometry
{
public:
  /**
   * @brief Constructor. Sets a blank geometry.
   */
  ComputationalGeometry();

  /**
   * @brief Destructor
   */
  virtual ~ComputationalGeometry();

  /**
   * @brief Get dielectrics
   * @return Dielectrics (m_dielectrics)
   */
  const Vector<Dielectric>&
  getDielectrics() const;

  /**
   * @brief Get electrodes
   * @return Electrodes (m_electrodes)
   */
  const Vector<Electrode>&
  getElectrodes() const;

  /**
   * @brief Get the background gas permittivity
   * @return Background gas permittivity
   */
  Real
  getGasPermittivity() const;

  /**
   * @brief Calls for ComputationalGeometry to use ScanShop rather than Chombo's default geometry generation tool.
   * @param[in] a_beginDomain Coarse domain where ScanShop begins the load balancing recursion process.
   */
  void
  useScanShop(const ProblemDomain& a_beginDomain);

  /**
   * @brief Calls for ComputationalGeometry to use Chombo's geometry generation tool
   */
  void
  useChomboShop();

  /**
   * @brief Generate the geometry with PolyhedralGeometryShop.
   * @details Same topology as the chombo-discharge generator, but the moments are exact
   * integrals of an explicitly reconstructed surface rather than quadratures of the implicit
   * function.
   * @param[in] a_beginDomain Level on which to initiate the load balancing sequence
   */
  void
  usePolyhedralShop(const ProblemDomain& a_beginDomain);

  /**
   * @brief Set dielectrics
   * @param[in] a_dielectrics Dielectris
   */
  void
  setDielectrics(const Vector<Dielectric>& a_dielectrics);

  /**
   * @brief Set electrodes
   * @param[in] a_electrodes Electrodes
   */
  void
  setElectrodes(const Vector<Electrode>& a_electrodes);

  /**
   * @brief Set the background permittivity
   * @param[in] a_eps0 Gas permittivity
   */
  void
  setGasPermittivity(const Real a_eps0);

  /**
   * @brief Get the multifluid index space.
   * @return m_multiFluidIndexSpace
   */
  const RefCountedPtr<MultiFluidIndexSpace>&
  getMfIndexSpace() const;

  /**
   * @brief Get the implicit function used to generate the gas-phase EBIS.
   * @return Gas-phase implicit function.
   */
  const RefCountedPtr<BaseIF>&
  getGasImplicitFunction() const;

  /**
   * @brief Get the implicit function used to generate the solid-phase EBIS.
   * @return Solid-phase implicit function.
   */
  const RefCountedPtr<BaseIF>&
  getSolidImplicitFunction() const;

  /**
   * @brief Get implicit function for the specified phase.
   * @param[in] a_phase Phase identifier.
   * @return Implicit function for the requested phase.
   */
  const RefCountedPtr<BaseIF>&
  getImplicitFunction(const phase::which_phase a_phase) const;

  /**
   * @brief Build the grids the index space is generated over.
   * @details Built before any index space exists, from the implicit functions alone. The start domain is built
   * whole, as is every coarser domain down to the coarsest that can still be coarsened by two; above it the
   * levels are factor-two refinements up to and including the stop domain. Levels are
   * indexed from the coarsest domain; getLevel finds the index of a domain. There is one box hierarchy for both
   * phases, built from both implicit functions together, and every box carries one classification per phase.
   * Without any implicit function every level is one regular box. The tile and super-tile sizes are the
   * class's own, ComputationalGeometry.min_block_size and ComputationalGeometry.max_block_size, read in the
   * constructor.
   * @param[in] a_startDomain  Domain the upward pass starts from.
   * @param[in] a_stopDomain   Finest domain built; a factor-two refinement of a_startDomain some number of times.
   * @param[in] a_probLo       Lower-left corner of the domain.
   * @param[in] a_startDx      Grid spacing on a_startDomain.
   * @param[in] a_refineAngle  Angle, in degrees, between neighbouring normals above which a box is split.
   * @param[in] a_maxGhostEB   Ghost cells a box is grown by when it is classified. At most half the tile size.
   */
  virtual void
  makeGrids(const ProblemDomain& a_startDomain,
            const ProblemDomain& a_stopDomain,
            const RealVect&      a_probLo,
            const Real           a_startDx,
            const Real           a_refineAngle,
            const int            a_maxGhostEB);

  /**
   * @brief Number of levels makeGrids built.
   * @return Number of levels, zero before makeGrids has run.
   */
  int
  getNumGridLevels() const noexcept;

  /**
   * @brief Level of a domain among the levels makeGrids built.
   * @param[in] a_domain Domain to look up.
   * @return Its level, 0 being the coarsest, or -1 if there is no such level.
   */
  int
  getLevel(const ProblemDomain& a_domain) const noexcept;

  /**
   * @brief Boxes of one classification for a phase on a level.
   * @param[in] a_phase Phase.
   * @param[in] a_level Level, 0 being the coarsest.
   * @param[in] a_type  Regular, covered or irregular.
   * @return The boxes.
   */
  Vector<Box>
  getBoxes(const phase::which_phase a_phase, const int a_level, const GeometryService::InOut a_type) const noexcept;

  /**
   * @brief Every box on a level, both phases' boxes being the same.
   * @param[in] a_level Level, 0 being the coarsest.
   * @return The boxes.
   */
  const Vector<Box>&
  getBoxes(const int a_level) const noexcept;

  /**
   * @brief The cut tiles on a level: the boxes of the level that carry cut cells or nest the level above, which is what
   * a simulation grid on the level is made of. Empty on the start level and below, which are whole.
   * @param[in] a_level Level, 0 being the coarsest.
   * @return The cut tiles.
   */
  const Vector<Box>&
  getCutTiles(const int a_level) const noexcept;

  /**
   * @brief The boxes that split in the upward pass on a level, with why: 1 a pair inside the box, 2 a pair reaching
   * into the ring outside it, 3 two surfaces facing each other across the band.
   * @param[in]  a_level   Level, 0 being the coarsest.
   * @param[out] a_reasons One reason per box, parallel to the returned boxes.
   * @return The boxes.
   */
  const Vector<Box>&
  getSplitBoxes(const int a_level, Vector<int>& a_reasons) const noexcept;

  /**
   * @brief Classification of every box on a level for a phase, parallel to getBoxes(level).
   * @param[in] a_phase Phase.
   * @param[in] a_level Level, 0 being the coarsest.
   * @return One classification per box.
   */
  const Vector<GeometryService::InOut>&
  getTypes(const phase::which_phase a_phase, const int a_level) const noexcept;

  /**
   * @brief Domain of a level.
   * @param[in] a_level Level, 0 being the coarsest.
   * @return The domain.
   */
  const ProblemDomain&
  getDomain(const int a_level) const noexcept;

  /**
   * @brief Grid spacing of a level.
   * @param[in] a_level Level, 0 being the coarsest.
   * @return The spacing.
   */
  Real
  getDx(const int a_level) const noexcept;

  /**
   * @brief Build geometries and the MFIndexSpace.
   * @param[in] a_finestDomain Finest domain
   * @param[in] a_probLo       Lower-left corner
   * @param[in] a_finestDx     Finest grid resolution
   * @param[in] a_nCellMax     Patch size
   * @param[in] a_maxGhostEB   Maximum number of EB ghosts that will be encountered.
   * @param[in] a_maxCoarsen   Max coarsenings to run. If = -1 then coarsen all the way down.
   * @details This will build the gas and solid phase domain. The input is the finest-level stuff and you can control
   * the division into boxes as well as the maximum number of coarsenings. The computed domain is (a_probLo, a_probLo +
   * a_box*a_dx)
   */
  virtual void
  buildGeometries(const ProblemDomain& a_finestDomain,
                  const RealVect&      a_probLo,
                  const Real           a_finestDx,
                  const int            a_nCellMax,
                  const int            a_maxGhostEB,
                  const int            a_maxCoarsen = -1);

protected:
  /**
   * @brief Threshold for Vof computation
   */
  constexpr static Real s_thresh = 1.E-15;

  /**
   * @brief Whether a cut cell whose body will not close stops the run.
   * @details Loud by choice while the polyhedral generator is being developed: a cell it cannot
   * close is an edge case worth seeing rather than one worth approximating around. No geometry
   * tested has produced one. Not exposed as an option, so that turning it off is a deliberate
   * edit rather than something an input file can do by accident.
   */
  constexpr static bool s_strictGeometry = true;

  /**
   * @brief Which generator to build the geometry with.
   * @details Each enumerator names the class that the geometry is generated with.
   */
  enum class Generator
  {
    GeometryShop,
    ScanShop,
    PolyhedralShop
  };

  /**
   * @brief Generator selected by the user.
   */
  Generator m_generator;

  /**
   * @brief Multifluid index spaces
   */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /**
   * @brief The gas-phase implicit function (i.e. outside electrodes and dielectrics).
   */
  RefCountedPtr<BaseIF> m_implicitFunctionGas;

  /**
   * @brief The solid-phase implicit function (i.e. the inside of the dielectrics).
   */
  RefCountedPtr<BaseIF> m_implicitFunctionSolid;

  /**
   * @brief List of dielectrics
   */
  Vector<Dielectric> m_dielectrics;

  /**
   * @brief List of electrodes
   */
  Vector<Electrode> m_electrodes;

  /**
   * @brief Grid level where we begin using ScanShop
   */
  ProblemDomain m_scanDomain;

  /**
   * @brief Lower-left corner of the domain.
   */
  RealVect m_probLo;

  /**
   * @brief Background permittivity
   */
  Real m_eps0;

  /**
   * @brief Angle, in degrees, between neighbouring normals above which an irregular box is split.
   */
  Real m_refineAngle;

  /**
   * @brief Maximum number of ghost cells that we will ever need
   */
  int m_maxGhostEB;

  /**
   * @brief Level of the start domain; every level below it is built whole.
   */
  int m_startLevel;

  /**
   * @brief Level of the stop domain, the finest level.
   */
  int m_stopLevel;

  /**
   * @brief Tile size, in cells. ComputationalGeometry.min_block_size, 8 when not given.
   */
  int m_minBlockSize;

  /**
   * @brief Super-tile size, in cells; the size every box makeGrids makes is split to.
   * ComputationalGeometry.max_block_size, 8 when not given.
   */
  int m_maxBlockSize;

  /**
   * @brief Domains of the levels makeGrids builds, coarsest first.
   */
  Vector<ProblemDomain> m_domains;

  /**
   * @brief Grid spacings of the levels makeGrids builds, coarsest first.
   */
  Vector<Real> m_dx;

  /**
   * @brief Cut-cell tiles per level, common to both phases.
   */
  Vector<Vector<Box>> m_cutTiles;

  /**
   * @brief Boxes per level, common to both phases.
   */
  Vector<Vector<Box>> m_boxes;

  /**
   * @brief Per level, the number of boxes that split for each SplitReason, and per level the number of tiles that lie
   * in an irregular box (tagged) versus elsewhere (nesting), for the report.
   */
  Vector<Vector<int>> m_splitCounts;

  /**
   * @brief Per level, the boxes that split in the upward pass, and why (SplitReason as an integer), parallel.
   */
  Vector<Vector<Box>> m_splitBoxes;

  /**
   * @brief Per level, the reason each box in m_splitBoxes split.
   */
  Vector<Vector<int>> m_splitReasons;

  /**
   * @brief Gas-phase classification of each box in m_boxes, parallel to it.
   */
  Vector<Vector<GeometryService::InOut>> m_gasTypes;

  /**
   * @brief Solid-phase classification of each box in m_boxes, parallel to it.
   */
  Vector<Vector<GeometryService::InOut>> m_solidTypes;

  /**
   * @brief Report the boxes and tiles of every level to pout.
   */
  void
  reportGrids() const;

  /**
   * @brief Build the composite implicit functions of the two phases from the electrodes and dielectrics.
   * @details The gas phase is the region outside every object; the solid phase is the region inside the
   * dielectrics and outside the electrodes, and is left null when there are no dielectrics.
   */
  void
  buildImplicitFunctions();

  /**
   * @brief The per-level classifications of one phase.
   * @param[in] a_phase Phase.
   * @return m_gasTypes or m_solidTypes.
   */
  Vector<Vector<GeometryService::InOut>>&
  types(const phase::which_phase a_phase) noexcept;

  /**
   * @brief The per-level classifications of one phase.
   * @param[in] a_phase Phase.
   * @return m_gasTypes or m_solidTypes.
   */
  const Vector<Vector<GeometryService::InOut>>&
  types(const phase::which_phase a_phase) const noexcept;

  /**
   * @brief Step 0: build and classify the start level.
   * @details Pseudocode:
   *
   *   m_boxes[m_startLevel] = domainSplit(m_domains[m_startLevel], m_maxBlockSize, block factor m_minBlockSize)
   *   classifyBoxes(m_boxes[m_startLevel], m_startLevel) -> m_gasTypes, m_solidTypes on m_startLevel
   */
  void
  buildStartLevel();

  /**
   * @brief Step 1: the upward pass, from the start level to the stop level.
   * @details Pseudocode:
   *
   *   for lvl = m_startLevel .. m_stopLevel - 1:
   *     for each box on lvl that is regular or covered in both phases:
   *       append refine(box, 2) on lvl + 1 with both tags                   -- never a hole beneath these
   *     for the boxes on lvl that are irregular in some phase:
   *       flags = splitFlags(...): exceedsCurvature in any phase where the box is irregular
   *       for each box with a flag: pieces = domainSplit(refine(box, 2), m_maxBlockSize, m_minBlockSize)
   *       classifyBoxes(all pieces, lvl + 1) in both phases, append on lvl + 1
   *       a box without a flag is a leaf: nothing is appended, the region above it is a hole on lvl + 1
   *
   *   The loop always runs to m_stopLevel: once no box splits, the curvature test has nothing to do, but the
   *   boxes that are regular or covered in both phases still refine whole onto every remaining level.
   *
   * The links from each box to its children are handed back for classifyTiles to descend by; they describe
   * the lists as the upward pass leaves them and are stale once decimateBoxes has rewritten a level.
   *
   * @param[out] a_firstChild  Per level and box, the index on the next level of the first child, or -1 for a leaf.
   * @param[out] a_numChildren Per level and box, the number of children: one for a whole refinement, the number
   *                           of pieces for a split, zero for a leaf.
   */
  void
  buildFinerLevels(Vector<Vector<int>>& a_firstChild, Vector<Vector<int>>& a_numChildren);

  /**
   * @brief Classify a list of boxes in both phases, the work shared between the ranks.
   * @details Rank r classifies every box whose index is r modulo the number of ranks; the answers are
   * combined with one all-reduce, so every rank returns the full lists in the input order.
   * @param[in]  a_boxes      Boxes to classify, identical on every rank.
   * @param[in]  a_level      Level the boxes live on.
   * @param[out] a_gasTypes   Gas-phase classification of each box.
   * @param[out] a_solidTypes Solid-phase classification of each box.
   */
  void
  classifyBoxes(const Vector<Box>&              a_boxes,
                const int                       a_level,
                Vector<GeometryService::InOut>& a_gasTypes,
                Vector<GeometryService::InOut>& a_solidTypes) const;

  /**
   * @brief Whether each box must split, the work shared between the ranks.
   * @details Same distribution as classifyBoxes. A box splits if exceedsCurvature holds in any phase in which
   * the box is irregular; a box irregular in neither phase never splits.
   * @param[in] a_boxes      Boxes to test, identical on every rank.
   * @param[in] a_level      Level the boxes live on.
   * @param[in] a_gasTypes   Gas-phase classification of each box.
   * @param[in] a_solidTypes Solid-phase classification of each box.
   * @return One flag per box, parallel to a_boxes.
   */
  Vector<int>
  splitFlags(const Vector<Box>&                    a_boxes,
             const int                             a_level,
             const Vector<GeometryService::InOut>& a_gasTypes,
             const Vector<GeometryService::InOut>& a_solidTypes) const;

  /**
   * @brief Classify a box of one phase on one level as regular, covered or irregular.
   * @details ScanShop::isRegular/isCovered on the implicit function directly; a phase without an implicit
   * function is regular everywhere. Pseudocode:
   *
   *   grown = grow(box, m_maxGhostEB) & domain
   *   halfDiagonal = 0.5 * dx * sqrt(SpaceDim)
   *   regular = every cell centre in grown has f(x) < -halfDiagonal
   *   covered = every cell centre in grown has f(x) >  halfDiagonal
   *   return Regular, Covered, or Irregular if neither
   *
   * @param[in] a_box   Box to classify.
   * @param[in] a_level Level the box lives on.
   * @param[in] a_phase Phase whose implicit function is asked.
   * @return The classification.
   */
  GeometryService::InOut
  classifyBox(const Box& a_box, const int a_level, const phase::which_phase a_phase) const;

  /**
   * @brief Why a box split, for the report: no split, a pair inside the box, a pair reaching into the ring one cell
   * outside it, or a pair whose normals face each other (two surfaces closer than the band).
   */
  enum class SplitReason
  {
    None,
    Interior,
    Ring,
    Medial
  };

  /**
   * @brief Whether the interface turns by more than m_refineAngle between neighbouring cut cells of a box.
   * @details No index space exists at this point; the cut cells are reconstructed from the implicit function
   * on the level, as the polyhedral shop reconstructs them, and the normal of each is that of its interface
   * polygon. That normal depends on the edge roots alone, so it is a property of the zero set: the implicit
   * function's gradient is not used, because inside a body built by CSG the function is not a distance and
   * its gradient carries the kinks and medial shells of the construction. A phase without an implicit function
   * never asks for a split. Pseudocode:
   *
   *   grown = grow(box, 1) & domain
   *   node values over grown, once per node; edge crossings bisected once and shared
   *   for each cell in grown: if it is cut and closes, n(cell) = polygon normal
   *   for each cut cell in the box and each cut neighbour in its 3^D block:
   *     angle > m_refineAngle: Interior if the neighbour is in the box, Ring if in the grown ring;
   *     normals facing (dot < 0): Medial
   *   the first pair inside the box decides; a ring pair decides only if the box has none of its own
   *
   * Across a sharp edge the angle never shrinks with dx, so edges refine to the stop level. That is intended.
   *
   * @param[in] a_box   Box to test.
   * @param[in] a_level Level the box lives on.
   * @param[in] a_phase Phase whose implicit function is asked.
   * @return True if the box must be split.
   */
  SplitReason
  exceedsCurvature(const Box& a_box, const int a_level, const phase::which_phase a_phase) const;

  /**
   * @brief Step 2: tile the boxes that are irregular in either phase into one properly nested set per level.
   * @details Pseudocode:
   *
   *   for lvl = m_startLevel + 1 .. m_stopLevel:
   *     tags[lvl - 1] = every box on lvl irregular in some phase, coarsened by two
   *                     -- TiledMeshRefine takes tags on the level below the one it tiles
   *   TiledMeshRefine(m_domains[m_startLevel], ratio 2, m_minBlockSize, m_maxBlockSize).regrid(tiles, tags)
   *   m_cutTiles[lvl] = tiles[lvl - m_startLevel] for lvl > m_startLevel; the start level and below are whole and
   *   are not tiled
   *
   * Every box is a super-tile already, so the tags are exact in tile units. The nesting buffer is one tile per
   * level, which is at least m_maxGhostEB when m_minBlockSize >= 2 * m_maxGhostEB. No collar is added: a
   * simulation ghost cell that reaches past the tiles into a hole is the index space's to fill, not carried here.
   */
  void
  makeTiles();

  /**
   * @brief Step 3: place every tile of every level in the boxes and classify it in both phases.
   * @details Pseudocode:
   *
   *   for lvl = m_startLevel + 1 .. m_stopLevel:
   *     for each tile t on lvl (every box is a super-tile or a whole refinement of one, so t lies in at most one):
   *       if t lies in box i:
   *         a_tileHosts[lvl][t] = i
   *         per phase: regular or covered box -> t inherits; irregular box -> classifyBox(t, lvl, phase),
   *                    since the box is conservative and t may be clear of the surface
   *       else (above a leaf, a hole):
   *         a_tileHosts[lvl][t] = -1
   *         per phase: classifyBox(t, lvl, phase) -- a carried tile is generated from the implicit function at
   *                    its own level, so that is the consistent answer
   *
   * Nothing is cut here; the hosts are what decimateBoxes reads. The box a tile lies in is found by descending
   * the child links from the start level, where the box containing a point is lattice arithmetic.
   *
   * @param[in]  a_firstChild     Child links from buildFinerLevels.
   * @param[in]  a_numChildren    Child counts from buildFinerLevels.
   * @param[out] a_gasTileTypes   Gas-phase classification of each tile, indexed [level][tile].
   * @param[out] a_solidTileTypes Solid-phase classification of each tile, indexed [level][tile].
   * @param[out] a_tileHosts      Index of the box each tile lies in, or -1, indexed [level][tile].
   */
  void
  classifyTiles(const Vector<Vector<int>>&              a_firstChild,
                const Vector<Vector<int>>&              a_numChildren,
                Vector<Vector<GeometryService::InOut>>& a_gasTileTypes,
                Vector<Vector<GeometryService::InOut>>& a_solidTileTypes,
                Vector<Vector<int>>&                    a_tileHosts) const;

  /**
   * @brief Index of the box of a whole level that contains a point.
   * @details A whole level is a domainSplit of its domain to m_maxBlockSize, so the box holding a point is
   * found from the point's position on that lattice; the table maps lattice cells to box indices.
   * @param[in] a_level Level, at or below the start level.
   * @return One box index per lattice cell, in lexicographic order of the lattice.
   */
  Vector<int>
  latticeTable(const int a_level) const;

  /**
   * @brief Look a point up in a lattice table.
   * @param[in] a_table Table from latticeTable.
   * @param[in] a_level Level the table was made for.
   * @param[in] a_point Cell on that level.
   * @return Index of the box containing the cell.
   */
  int
  latticeLookup(const Vector<int>& a_table, const int a_level, const IntVect& a_point) const;

  /**
   * @brief Key of the super-tile a cell lies in, for finding tiles.
   * @param[in] a_cell Cell on some level.
   * @return Packed super-tile coordinates, 21 bits per direction.
   */
  long long
  superTileKey(const IntVect& a_cell) const noexcept;

  /**
   * @brief Tiles of a level intersecting a box, found through their super-tile keys.
   * @details The tiles of a level are indexed by the key of the super-tile they lie in, so the tiles meeting a
   * box are found by looking up the keys of the super-tiles the box touches.
   * @param[in] a_level  Level, above the start level.
   * @param[in] a_order  Tile indices of the level sorted by super-tile key, from sortTilesByKey.
   * @param[in] a_box    Box on that level.
   * @return Indices into m_cutTiles[a_level] of the tiles that intersect the box.
   */
  Vector<int>
  tilesMeeting(const int a_level, const Vector<int>& a_order, const Box& a_box) const;

  /**
   * @brief Tile indices of a level sorted by super-tile key.
   * @param[in] a_level Level, above the start level.
   * @return The indices.
   */
  Vector<int>
  sortTilesByKey(const int a_level) const;

  /**
   * @brief Step 4: cut every hit box down to what the tiles left of it, and replace the lists on the tiled
   * levels with the result.
   * @details Pseudocode:
   *
   *   for lvl = m_startLevel + 1 .. m_stopLevel:
   *     new list = every tile on lvl with its two classifications
   *     for each box b on lvl:
   *       tiles = { t : a_tileHosts[lvl][t] == b }
   *       if tiles is empty:                         new list += b with its two classifications
   *       else if b is irregular in some phase:      nothing, the tiles cover it
   *       else:                                      new list += (b minus tiles) with b's two classifications
   *                            -- first implementation: TreeIntVectSet in tile coordinates, then createBoxes;
   *                            -- correct and octree-graded, not tight; the packer is replaceable here
   *     m_boxes, m_gasTypes, m_solidTypes on lvl = new list
   *
   * @param[in] a_gasTileTypes   Gas-phase classification of each tile, from classifyTiles.
   * @param[in] a_solidTileTypes Solid-phase classification of each tile, from classifyTiles.
   * @param[in] a_tileHosts      Host of each tile, from classifyTiles.
   */
  void
  decimateBoxes(const Vector<Vector<GeometryService::InOut>>& a_gasTileTypes,
                const Vector<Vector<GeometryService::InOut>>& a_solidTileTypes,
                const Vector<Vector<int>>&                    a_tileHosts);

  /**
   * @brief Step 5: build the levels coarser than the start level, and push irregularity down onto every whole
   * level.
   * @details The coarser levels are built as ScanShop builds them, whole and each box classified on its own.
   * On its own a box's classification agrees with the level above only because classifyBox is conservative for
   * a signed-distance function; the second half makes the agreement hold by construction instead. Pseudocode:
   *
   *   for lvl = m_startLevel - 1 down to 0:
   *     m_boxes[lvl] = domainSplit(m_domains[lvl], m_maxBlockSize)   -- no block factor: these levels are not
   *                                                                  -- tiled and the coarsest may be smaller
   *                                                                  -- than a tile
   *     classifyBoxes(m_boxes[lvl], lvl) in both phases
   *
   *   for lvl = m_stopLevel - 1 down to 0, per phase:
   *     for each box b on lvl + 1 irregular in that phase:
   *       the box on lvl containing coarsen(b, 2) is irregular in that phase
   *       -- exactly one box contains it: every box is a whole super-tile, a whole refinement of one, or a union
   *       -- of whole tiles, and coarsen(b, 2) is at most half a tile wide on the same lattice.
   *
   * Runs after decimation: the coarser levels depend on nothing above them, and the push-down needs the final
   * lists of every level above the one it marks.
   */
  void
  buildCoarserLevels();

  /**
   * @brief Set up the geometry generation tool for the gas phase.
   * @param[in,out] a_geoserver     Geometry service object which is later used for making build the EB information.
   * @param[in]    a_finestDomain  Finest domain which will be used
   * @param[in]    a_probLo        Lower-left corner of simulation domain.
   * @param[in]    a_finestDx      Finest resolution which will be used.
   */
  void
  buildGasGeometry(GeometryService*&    a_geoserver,
                   const ProblemDomain& a_finestDomain,
                   const RealVect&      a_probLo,
                   const Real           a_finestDx);

  /**
   * @brief Set up the geometry generation tool for the solid phase, i.e. the part inside the dielectrics.
   * @param[in,out] a_geoserver     Geometry service object which is later used for making build the EB information.
   * @param[in]    a_finestDomain  Finest domain which will be used
   * @param[in]    a_probLo        Lower-left corner of simulation domain.
   * @param[in]    a_finestDx      Finest resolution which will be used.
   */
  void
  buildSolidGeometry(GeometryService*&    a_geoserver,
                     const ProblemDomain& a_finestDomain,
                     const RealVect&      a_probLo,
                     const Real           a_finestDx);
};

#include <CD_NamespaceFooter.H>

#endif
