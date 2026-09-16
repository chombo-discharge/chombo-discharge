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
   * @brief Build the grids the index space is generated over, per phase and per level.
   * @details The grids are built before any index space exists, from the implicit functions alone. The start
   * domain is built whole; every finer level is a factor-two refinement of it, up to a_maxEbDepth levels, the
   * last of which has spacing a_finestDx; every coarser level, down to the coarsest domain that can still be
   * coarsened by two, is built whole as well, as ScanShop builds them. Levels are indexed from that coarsest
   * domain, and getLevel finds the index of a domain. Every box made here is at most a_maxBlockSize wide and
   * lies on the a_maxBlockSize lattice, so the cut-cell grids of both phases are unions of the same super-tiles.
   * The result, per phase and per level, is a box list with a classification for every box, read through
   * getGrids/getGridTypes; the cut-cell tiles common to both phases are read through getTiles. Without an
   * implicit function a phase is one regular box on every level.
   * @param[in] a_startDomain  Domain the upward pass starts from. Built whole in both phases.
   * @param[in] a_probLo       Lower-left corner of the domain.
   * @param[in] a_finestDx     Grid spacing on the finest level, a_maxEbDepth refinements above a_startDomain.
   * @param[in] a_refineAngle  Angle, in degrees, between neighbouring normals above which a box is split.
   * @param[in] a_maxEbDepth   Number of factor-two levels above a_startDomain that may be built.
   * @param[in] a_minBlockSize Tile size in cells: the finest granularity of the cut-cell grids, and the nesting
   *                           buffer between levels. Must be at least twice a_maxGhostEB.
   * @param[in] a_maxBlockSize Super-tile size in cells: the size every box made here is split to. A multiple of
   *                           a_minBlockSize.
   * @param[in] a_maxGhostEB   Ghost cells a box is grown by when it is classified.
   */
  virtual void
  makeGrids(const ProblemDomain& a_startDomain,
            const RealVect&      a_probLo,
            const Real           a_finestDx,
            const Real           a_refineAngle,
            const int            a_maxEbDepth,
            const int            a_minBlockSize,
            const int            a_maxBlockSize,
            const int            a_maxGhostEB);

  /**
   * @brief Number of levels makeGrids built, level 0 being the coarsest domain.
   * @return Number of levels; zero before makeGrids has run.
   */
  int
  getNumGridLevels() const noexcept;

  /**
   * @brief Level index of a domain among the levels makeGrids built.
   * @param[in] a_domain Domain to look up.
   * @return Its level, 0 being the coarsest domain, or -1 if makeGrids built no such level.
   */
  int
  getLevel(const ProblemDomain& a_domain) const noexcept;

  /**
   * @brief Final grids of one phase on one level, after decimation.
   * @param[in] a_phase Phase.
   * @param[in] a_level Level, 0 being the coarsest domain.
   * @return Boxes of the level, parallel to getGridTypes.
   */
  const Vector<Box>&
  getGrids(const phase::which_phase a_phase, const int a_level) const noexcept;

  /**
   * @brief Classification of each box in getGrids.
   * @param[in] a_phase Phase.
   * @param[in] a_level Level, 0 being the coarsest domain.
   * @return One classification per box, parallel to getGrids.
   */
  const Vector<GeometryService::InOut>&
  getGridTypes(const phase::which_phase a_phase, const int a_level) const noexcept;

  /**
   * @brief Boxes of one classification that the upward pass made for a phase on a level, before tiling.
   * @details These are the regular, covered and irregular boxes as the upward pass left them; the irregular
   * ones are what the tiles are made from, and what will be shown as a simulation grid for inspection.
   * @param[in] a_phase Phase.
   * @param[in] a_level Level, 0 being the coarsest domain.
   * @param[in] a_type  Classification to select.
   * @return The boxes with that classification.
   */
  Vector<Box>
  getBoxes(const phase::which_phase a_phase, const int a_level, const GeometryService::InOut a_type) const noexcept;

  /**
   * @brief The properly nested cut-cell tiles on one level, common to both phases.
   * @details This is the coverage the simulation regrids onto.
   * @param[in] a_level Level, 0 being the coarsest domain.
   * @return Tiles of the level.
   */
  const Vector<Box>&
  getTiles(const int a_level) const noexcept;

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
   * @brief Multifluid index spaces
   */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /**
   * @brief Background permittivity
   */
  Real m_eps0;

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
   * @brief Grid level where we begin using ScanShop
   */
  ProblemDomain m_scanDomain;

  /**
   * @brief Maximum number of ghost cells that we will ever need
   */
  int m_maxGhostEB;

  /**
   * @brief List of dielectrics
   */
  Vector<Dielectric> m_dielectrics;

  /**
   * @brief List of electrodes
   */
  Vector<Electrode> m_electrodes;

  /**
   * @brief The gas-phase implicit function (i.e. outside electrodes and dielectrics).
   */
  RefCountedPtr<BaseIF> m_implicitFunctionGas;

  /**
   * @brief The solid-phase implicit function (i.e. the inside of the dielectrics).
   */
  RefCountedPtr<BaseIF> m_implicitFunctionSolid;

  /**
   * @brief Build the composite implicit functions of the two phases from the electrodes and dielectrics.
   * @details The gas phase is the region outside every object; the solid phase is the region inside the
   * dielectrics and outside the electrodes, and is left null when there are no dielectrics.
   */
  void
  buildImplicitFunctions();

  /**
   * @brief Step 0: build and classify the start level of one phase.
   * @details Pseudocode:
   *
   *   boxes = domainSplit(m_gridDomains[m_startLevel], m_maxBlockSize)
   *   for each box (this rank's share): type = classifyBox(box, m_startLevel, phase)
   *   gather boxes and types onto every rank, in a fixed order
   *   m_levelBoxes[phase][m_startLevel], m_levelBoxTypes[phase][m_startLevel] = boxes, types
   *
   * @param[in] a_phase Phase to build.
   */
  void
  buildStartLevel(const phase::which_phase a_phase);

  /**
   * @brief Step 1: the upward pass of one phase, from the start level to the finest level.
   * @details Pseudocode:
   *
   *   for lvl = m_startLevel .. m_startLevel + m_maxEbDepth - 1:
   *     for each box on lvl (this rank's share):
   *       if type is Regular or Covered:
   *         push refine(box, 2) on lvl + 1 with the same type          -- never a hole beneath these
   *       else if exceedsCurvature(box, lvl, phase):
   *         for each piece of domainSplit(refine(box, 2), m_maxBlockSize):
   *           push piece on lvl + 1 with classifyBox(piece, lvl + 1, phase)
   *       else:
   *         leaf: nothing is pushed; the region above it is a hole on lvl + 1
   *     gather lvl + 1 onto every rank, in a fixed order
   *     stop if no rank split a box
   *
   * @param[in] a_phase Phase to build.
   */
  void
  buildFinerLevels(const phase::which_phase a_phase);

  /**
   * @brief Classify a box of one phase on one level as regular, covered or irregular.
   * @details Pseudocode, ScanShop::isRegular/isCovered on the implicit function directly:
   *
   *   grown = grow(box, m_maxGhostEB) & domain
   *   halfDiagonal = 0.5 * dx * sqrt(SpaceDim)
   *   regular = every cell centre in grown has f(x) < -halfDiagonal   -- with the signed-distance scan skip
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
   * @brief Whether the implicit function's normal turns by more than m_refineAngle inside a box.
   * @details Cut cells do not exist at this point and are not looked for; the test runs on the implicit
   * function in the band of cells near its zero set. Pseudocode:
   *
   *   grown = grow(box, 1) & domain
   *   for each cell centre x in grown with |f(x)| <= sqrt(SpaceDim) * dx      -- scan skip prunes the rest
   *     n(x) = central difference of f at x, normalised
   *   for each such cell and each of its neighbours in the band:
   *     if angle(n(cell), n(neighbour)) > m_refineAngle: return true      -- early exit
   *   return false
   *
   * Across a sharp edge the angle never shrinks with dx, so edges refine to m_maxEbDepth. That is intended.
   *
   * @param[in] a_box   Box to test.
   * @param[in] a_level Level the box lives on.
   * @param[in] a_phase Phase whose implicit function is asked.
   * @return True if the box must be split.
   */
  bool
  exceedsCurvature(const Box& a_box, const int a_level, const phase::which_phase a_phase) const;

  /**
   * @brief Step 2: tile the union of both phases' irregular boxes into one properly nested set per level.
   * @details Runs once, after both phases have finished the upward pass. Pseudocode:
   *
   *   for lvl = m_startLevel + 1 .. m_startLevel + m_maxEbDepth:
   *     tags[lvl - 1] = union over phases of the irregular boxes on lvl, coarsened by two
   *                     -- TiledMeshRefine takes tags on the level below the one it tiles
   *   TiledMeshRefine(m_gridDomains[m_startLevel], ratio 2, m_minBlockSize, m_maxBlockSize).regrid(tiles, tags)
   *   m_tiles[lvl] = tiles[lvl - m_startLevel] for lvl > m_startLevel; the start level and below are whole and
   *   are not tiled
   *
   * Every irregular box is a super-tile already, so the union is exact in tile units. The nesting buffer is one
   * tile per level, which is at least m_maxGhostEB when m_minBlockSize >= 2 * m_maxGhostEB. No collar is added
   * beyond the union: a simulation ghost cell that reaches past the tiles into a hole is filled by the index
   * space from the irregular leaf beneath it, not carried here.
   */
  void
  makeTiles();

  /**
   * @brief Step 3: place every tile of every level in this phase's boxes and classify it.
   * @details Pseudocode:
   *
   *   for lvl = m_startLevel + 1 .. m_startLevel + m_maxEbDepth:
   *     walk the BVH over m_levelBoxes[phase][lvl] against the BVH over m_tiles[lvl]
   *     for each tile t (every box is a super-tile or a whole refinement of one, so t meets at most one box):
   *       if t lies in box b:
   *         m_tileHost[phase][lvl][t] = b
   *         m_tileTypes[phase][lvl][t] = type of b                -- Irregular stays Irregular; Regular/Covered inherit
   *       else:                                                   -- above a leaf of this phase on lvl - 1
   *         m_tileHost[phase][lvl][t] = -1
   *         m_tileTypes[phase][lvl][t] = classifyBox(t, lvl, phase)
   *
   * Nothing is cut here; which tiles each regular or covered box hosts is what decimateBoxes reads.
   *
   * @param[in] a_phase Phase to classify the tiles for.
   */
  void
  classifyTiles(const phase::which_phase a_phase);

  /**
   * @brief Step 4: cut every hit regular or covered box down to what the tiles left of it, and assemble the
   * final grids of the phase.
   * @details Pseudocode:
   *
   *   for lvl = m_startLevel .. m_startLevel + m_maxEbDepth:
   *     grids = every tile on lvl with its type from m_tileTypes
   *     for each box b on lvl:
   *       tiles = { t : m_tileHost[phase][lvl][t] == b }
   *       if tiles is empty:   grids += b with its type
   *       else if b is Irregular: nothing                          -- the tiles already cover it
   *       else:                grids += (b minus tiles) with b's type
   *                             -- first implementation: TreeIntVectSet in tile coordinates, then createBoxes;
   *                             -- correct and octree-graded, not tight; the packer is replaceable here
   *     m_grids[phase][lvl], m_gridTypes[phase][lvl] = grids
   *
   * @param[in] a_phase Phase to assemble.
   */
  void
  decimateBoxes(const phase::which_phase a_phase);

  /**
   * @brief Step 5: build the levels coarser than the start level, and push irregularity down onto every whole
   * level.
   * @details The coarser levels are built as ScanShop builds them, whole and each box classified on its own.
   * On its own a box's classification agrees with the level above only because classifyBox is conservative for
   * a signed-distance function; the second half makes the agreement hold by construction instead. Pseudocode:
   *
   *   for lvl = m_startLevel - 1 down to 0, per phase:
   *     boxes = domainSplit(m_gridDomains[lvl], m_maxBlockSize)
   *     for each box (this rank's share): type = classifyBox(box, lvl, phase)
   *     gather onto every rank, in a fixed order
   *     m_levelBoxes/m_levelBoxTypes and m_grids/m_gridTypes on lvl = boxes, types    -- nothing to decimate
   *
   *   for lvl = finest - 1 down to 0, per phase:
   *     for each Irregular box b in m_grids[phase][lvl + 1]:
   *       the box of m_grids[phase][lvl] containing coarsen(b, 2) is Irregular
   *       -- exactly one box contains it: every box is a whole super-tile, a whole refinement of one, or a union
   *       -- of whole tiles, and coarsen(b, 2) is at most half a tile wide on the same lattice. The box it flips
   *       -- is a tile or a super-tile, never a decimated remainder, since above the start level a finer
   *       -- irregular tile always lies inside a tile of this level by nesting.
   *
   * Runs after decimation: the coarser levels depend on nothing above them, and the push-down needs the final
   * lists of every level above the one it marks.
   */
  void
  buildCoarserLevels();

  /**
   * @brief Domains of the levels makeGrids builds, level 0 being the coarsest domain.
   */
  Vector<ProblemDomain> m_gridDomains;

  /**
   * @brief Index into m_gridDomains of the start domain; every level below it is built whole.
   */
  int m_startLevel;

  /**
   * @brief Grid spacings of the levels makeGrids builds.
   */
  Vector<Real> m_gridDx;

  /**
   * @brief Lower-left corner of the domain, as makeGrids was given it.
   */
  RealVect m_gridProbLo;

  /**
   * @brief Number of levels above the start domain that makeGrids may build.
   */
  int m_maxEbDepth;

  /**
   * @brief Tile size, in cells.
   */
  int m_minBlockSize;

  /**
   * @brief Super-tile size, in cells; the size every box makeGrids makes is split to.
   */
  int m_maxBlockSize;

  /**
   * @brief Angle, in degrees, between neighbouring normals above which an irregular box is split.
   */
  Real m_refineAngle;

  /**
   * @brief Boxes the upward pass made, indexed [phase][level], before tiling and decimation.
   * @details Kept as made: the irregular ones are what the tiles come from, and the three classifications
   * are what getBoxes selects among.
   */
  Vector<Vector<Vector<Box>>> m_levelBoxes;

  /**
   * @brief Classification of each box in m_levelBoxes, parallel to it.
   */
  Vector<Vector<Vector<GeometryService::InOut>>> m_levelBoxTypes;

  /**
   * @brief The properly nested cut-cell tiles, indexed [level], common to both phases.
   */
  Vector<Vector<Box>> m_tiles;

  /**
   * @brief For each tile, indexed [phase][level][tile]: the index into m_levelBoxes of the box it lies in, or -1
   * when it lies above a leaf of that phase.
   */
  Vector<Vector<Vector<int>>> m_tileHost;

  /**
   * @brief Classification of each tile in each phase, indexed [phase][level][tile].
   */
  Vector<Vector<Vector<GeometryService::InOut>>> m_tileTypes;

  /**
   * @brief Final grids, indexed [phase][level]: the tiles, the decimated remainders, and the untouched boxes.
   */
  Vector<Vector<Vector<Box>>> m_grids;

  /**
   * @brief Classification of each box in m_grids, parallel to it.
   */
  Vector<Vector<Vector<GeometryService::InOut>>> m_gridTypes;

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
