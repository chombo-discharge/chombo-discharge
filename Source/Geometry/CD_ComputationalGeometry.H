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
   * @brief Build the grids the index space is generated over.
   * @details Built before any index space exists, from the implicit functions alone. The start domain is built
   * whole, as is every coarser domain down to the coarsest that can still be coarsened by two; above it the
   * levels are factor-two refinements up to and including the stop domain. Levels are
   * indexed from the coarsest domain; getLevel finds the index of a domain. Every box lies on the
   * a_maxBlockSize lattice, so the cut-cell grids of both phases are unions of the same super-tiles. Without an
   * implicit function a phase is one regular box on every level.
   * @param[in] a_startDomain  Domain the upward pass starts from.
   * @param[in] a_stopDomain   Finest domain built; a factor-two refinement of a_startDomain some number of times.
   * @param[in] a_probLo       Lower-left corner of the domain.
   * @param[in] a_startDx      Grid spacing on a_startDomain.
   * @param[in] a_refineAngle  Angle, in degrees, between neighbouring normals above which a box is split.
   * @param[in] a_minBlockSize Tile size in cells. At least twice a_maxGhostEB.
   * @param[in] a_maxBlockSize Super-tile size in cells. A multiple of a_minBlockSize.
   * @param[in] a_maxGhostEB   Ghost cells a box is grown by when it is classified.
   */
  virtual void
  makeGrids(const ProblemDomain& a_startDomain,
            const ProblemDomain& a_stopDomain,
            const RealVect&      a_probLo,
            const Real           a_startDx,
            const Real           a_refineAngle,
            const int            a_minBlockSize,
            const int            a_maxBlockSize,
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
  const Vector<Box>&
  getBoxes(const phase::which_phase a_phase, const int a_level, const GeometryService::InOut a_type) const noexcept;

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
   * @brief Tile size, in cells.
   */
  int m_minBlockSize;

  /**
   * @brief Super-tile size, in cells; the size every box makeGrids makes is split to.
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
  Vector<Vector<Box>> m_tiles;

  /**
   * @brief Regular boxes of the gas phase, per level.
   */
  Vector<Vector<Box>> m_gasRegularBoxes;

  /**
   * @brief Covered boxes of the gas phase, per level.
   */
  Vector<Vector<Box>> m_gasCoveredBoxes;

  /**
   * @brief Irregular boxes of the gas phase, per level.
   */
  Vector<Vector<Box>> m_gasIrregularBoxes;

  /**
   * @brief Regular boxes of the solid phase, per level.
   */
  Vector<Vector<Box>> m_solidRegularBoxes;

  /**
   * @brief Covered boxes of the solid phase, per level.
   */
  Vector<Vector<Box>> m_solidCoveredBoxes;

  /**
   * @brief Irregular boxes of the solid phase, per level.
   */
  Vector<Vector<Box>> m_solidIrregularBoxes;

  /**
   * @brief Build the composite implicit functions of the two phases from the electrodes and dielectrics.
   * @details The gas phase is the region outside every object; the solid phase is the region inside the
   * dielectrics and outside the electrodes, and is left null when there are no dielectrics.
   */
  void
  buildImplicitFunctions();

  /**
   * @brief The per-level box list of one phase and one classification.
   * @param[in] a_phase Phase.
   * @param[in] a_type  Regular, covered or irregular.
   * @return The list, indexed by level.
   */
  Vector<Vector<Box>>&
  boxes(const phase::which_phase a_phase, const GeometryService::InOut a_type) noexcept;

  /**
   * @brief The per-level box list of one phase and one classification.
   * @param[in] a_phase Phase.
   * @param[in] a_type  Regular, covered or irregular.
   * @return The list, indexed by level.
   */
  const Vector<Vector<Box>>&
  boxes(const phase::which_phase a_phase, const GeometryService::InOut a_type) const noexcept;

  /**
   * @brief Step 0: build and classify the start level of one phase.
   * @details Pseudocode:
   *
   *   for each box of domainSplit(m_domains[m_startLevel], m_maxBlockSize) (this rank's share):
   *     append it to the phase's regular, covered or irregular list on m_startLevel, per classifyBox
   *   gather the three lists onto every rank, in a fixed order
   *
   * @param[in] a_phase Phase to build.
   */
  void
  buildStartLevel(const phase::which_phase a_phase);

  /**
   * @brief Step 1: the upward pass of one phase, from the start level to the finest level.
   * @details Pseudocode:
   *
   *   for lvl = m_startLevel .. m_stopLevel - 1:
   *     for each regular or covered box on lvl (this rank's share):
   *       append refine(box, 2) to the same list on lvl + 1                -- never a hole beneath these
   *     for each irregular box on lvl (this rank's share):
   *       if exceedsCurvature(box, lvl, phase):
   *         for each piece of domainSplit(refine(box, 2), m_maxBlockSize):
   *           append it to the list on lvl + 1 that classifyBox(piece, lvl + 1, phase) names
   *       else:
   *         leaf: nothing is appended; the region above it is a hole on lvl + 1
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
   * Across a sharp edge the angle never shrinks with dx, so edges refine to the finest level. That is intended.
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
   *   for lvl = m_startLevel + 1 .. m_stopLevel:
   *     tags[lvl - 1] = union over phases of the irregular boxes on lvl, coarsened by two
   *                     -- TiledMeshRefine takes tags on the level below the one it tiles
   *   TiledMeshRefine(m_domains[m_startLevel], ratio 2, m_minBlockSize, m_maxBlockSize).regrid(tiles, tags)
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
   *   for lvl = m_startLevel + 1 .. m_stopLevel:
   *     walk the BVH over the phase's boxes on lvl against the BVH over m_tiles[lvl]
   *     for each tile t (every box is a super-tile or a whole refinement of one, so t meets at most one box):
   *       if t lies in an irregular box:            a_tileTypes[lvl][t] = Irregular,  a_tileHosts[lvl][t] = -1
   *       else if t lies in regular/covered box i:  a_tileTypes[lvl][t] = that type,  a_tileHosts[lvl][t] = i
   *       else (above a leaf of this phase):        a_tileTypes[lvl][t] = classifyBox(t, lvl, phase), host -1
   *
   * Nothing is cut here; the hosts are what decimateBoxes reads.
   *
   * @param[in]  a_phase     Phase to classify the tiles for.
   * @param[out] a_tileTypes Classification of each tile, indexed [level][tile].
   * @param[out] a_tileHosts Index of the regular or covered box each tile lies in, or -1, indexed [level][tile].
   */
  void
  classifyTiles(const phase::which_phase                a_phase,
                Vector<Vector<GeometryService::InOut>>& a_tileTypes,
                Vector<Vector<int>>&                    a_tileHosts) const;

  /**
   * @brief Step 4: cut every hit regular or covered box down to what the tiles left of it, and replace the
   * phase's lists on the tiled levels with the result.
   * @details Pseudocode:
   *
   *   for lvl = m_startLevel + 1 .. m_stopLevel:
   *     irregular on lvl = the tiles typed Irregular                       -- the tiles cover every irregular box
   *     for regular and for covered:
   *       new list = the tiles typed so
   *       for each box b in the old list:
   *         tiles = { t : a_tileHosts[lvl][t] == b }
   *         if tiles is empty: new list += b
   *         else:              new list += b minus tiles
   *                            -- first implementation: TreeIntVectSet in tile coordinates, then createBoxes;
   *                            -- correct and octree-graded, not tight; the packer is replaceable here
   *
   * @param[in] a_phase     Phase to assemble.
   * @param[in] a_tileTypes Classification of each tile, from classifyTiles.
   * @param[in] a_tileHosts Host of each tile, from classifyTiles.
   */
  void
  decimateBoxes(const phase::which_phase                      a_phase,
                const Vector<Vector<GeometryService::InOut>>& a_tileTypes,
                const Vector<Vector<int>>&                    a_tileHosts);

  /**
   * @brief Step 5: build the levels coarser than the start level, and push irregularity down onto every whole
   * level.
   * @details The coarser levels are built as ScanShop builds them, whole and each box classified on its own.
   * On its own a box's classification agrees with the level above only because classifyBox is conservative for
   * a signed-distance function; the second half makes the agreement hold by construction instead. Pseudocode:
   *
   *   for lvl = m_startLevel - 1 down to 0, per phase:
   *     for each box of domainSplit(m_domains[lvl], m_maxBlockSize) (this rank's share):
   *       append it to the phase's regular, covered or irregular list on lvl, per classifyBox
   *       -- no block factor: these levels are not tiled and the coarsest may be smaller than a tile
   *     gather onto every rank, in a fixed order
   *
   *   for lvl = m_stopLevel - 1 down to 0, per phase:
   *     for each irregular box b on lvl + 1:
   *       move the box on lvl containing coarsen(b, 2) from the regular or covered list to the irregular one
   *       -- exactly one box contains it: every box is a whole super-tile, a whole refinement of one, or a union
   *       -- of whole tiles, and coarsen(b, 2) is at most half a tile wide on the same lattice. The box it moves
   *       -- is a tile or a super-tile, never a decimated remainder, since above the start level a finer
   *       -- irregular tile always lies inside a tile of this level by nesting.
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
