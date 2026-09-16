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
#include <MFIndexSpace.H>
#include <Box.H>
#include <RealVect.H>
#include <ProblemDomain.H>
#include <IntVectSet.H>

// Our includes
#include <CD_Dielectric.H>
#include <CD_Electrode.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_PolyhedralEBUtils.H>
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
   * @brief Write the embedded boundary the polyhedral generator reconstructs, as STL.
   * @details One file per phase, named for the base name given here with the phase appended. It
   *          holds the interface of every cut cell that no finer level also carries, so it is the
   *          finest description of each part of the geometry rather than one surface per level.
   *          An empty name writes nothing.
   * @param[in] a_fileName Base name for the output, or empty to write nothing.
   */
  void
  setGeometrySurfaceFile(const std::string& a_fileName) noexcept;

  /**
   * @brief The generator that built one phase's index space.
   * @details Null when the geometry was not built, or was built by something other than the
   *          polyhedral generator. It is kept alive for as long as the index space it built, so a
   *          cell can be reconstructed from it long after generation.
   * @param[in] a_phase Phase to ask about.
   * @return The generator, or null.
   */
  const RefCountedPtr<GeometryService>&
  getGeometryGenerator(const phase::which_phase a_phase) const noexcept;

  /**
   * @brief Generate the geometry with PolyhedralGeometryShop.
   * @details Same topology as the chombo-discharge generator, but the moments are exact
   * integrals of an explicitly reconstructed surface rather than quadratures of the implicit
   * function.
   * @param[in] a_beginDomain Level on which to initiate the load balancing sequence
   * @param[in] a_refinement  How many times the generated cells refine the cells the surface is
   * reconstructed on. One reconstructs on the generated cells themselves.
   */
  void
  usePolyhedralShop(const ProblemDomain& a_beginDomain, const int a_refinement);

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

  /**
   * @brief Assemble the implicit functions for the two phases from the objects set on this class.
   * @details Called by buildGeometries before it makes the generators, and callable ahead of it
   * by anything that needs to read the geometry before the index spaces exist. Doing it twice is
   * harmless: the result depends only on the electrodes and dielectrics.
   */
  void
  buildImplicitFunctions();

  /**
   * @brief Say how much of each level the generator should carry.
   * @details Must be set before buildGeometries, since that is when the generator is made. The
   * boxes are the ones getCurvatureTags clustered its tags into, coarsest level first.
   * @param[in] a_regions        Boxes each level was swept over, which is how much of it to carry.
   * @param[in] a_coarsestDomain Domain of the first entry in a_regions.
   * @param[in] a_buffer         Cells to carry beyond what the boxes ask for, on every level.
   */
  void
  setCoverage(const Vector<Vector<Box>>& a_regions, const ProblemDomain& a_coarsestDomain, const int a_buffer);

  /**
   * @brief Get the boxes the geometry was generated over, coarsest level first.
   * @details Empty unless the pre-pass ran. Where the index space is carried only over part of a
   * level, these are the boxes it was carried over, so grids built anywhere else ask for embedded
   * boundary data that was never generated. Indexed as the AMR hierarchy is, coarsest first,
   * which is the opposite of how the index space indexes its own levels.
   * @return Boxes each level was generated over.
   */
  const Vector<Vector<Box>>&
  getCoverageRegions() const noexcept;

  /**
   * @brief Where the interface bends too sharply for the resolution it sits on.
   * @details A pre-pass over the implicit function alone. Each candidate cell's interface normal
   * is fitted to the positions of its edge crossings, and a cell whose normal turns by more than
   * a_angle against any of its neighbours' is tagged. No moments, no connectivity and no index
   * space are built, so the tags can be had before any geometry is committed and nothing has to
   * be generated and thrown away.
   *
   * Normals come from crossingNormal rather than from derivatives of the implicit function,
   * because root finding only asks where the function changes sign. Implicit functions that are
   * not differentiable at their own zero set are common enough that a derivative-based criterion
   * fails on geometries already in the test suite.
   *
   * Each level is swept only where the level above it tagged, with TiledMeshRefine turning tags
   * into the regions to descend into, so the work follows the surface rather than the volume and
   * proper nesting is left to the code that already does it.
   * @param[out] a_regions       Boxes each level was swept over, coarsest first. These are what
   * the grids of a run are built from, so a level carried only where the geometry needs it has to
   * be carried over exactly these and not over the tags with a margin guessed around them.
   * @param[in] a_coarsestDomain Domain of the coarsest level.
   * @param[in] a_refRatios      Refinement ratio between each level and the next. Must reach
   * a_maxDepth, which is the pre-pass's own cap and need not be how deep the run may refine.
   * @param[in] a_tileSize       Tile size for the clustering, in cells.
   * @param[in] a_maxBlockSize   Largest box the clustering may emit, in cells.
   * @param[in] a_probLo         Lower left corner of the domain.
   * @param[in] a_coarsestDx     Grid spacing on the coarsest level.
   * @param[in] a_angle          Turn, in degrees, above which a cell is tagged.
   * @param[in] a_growth         Cells by which each level's tags are grown before the next level
   * is swept. The refined region then reaches past the cells that asked for it, which is what
   * keeps the boundary between coarsened and freshly generated cells away from the features that
   * asked for the refinement in the first place.
   * @param[in] a_maxDepth       Hard cap on how many levels deep the pre-pass will go. Needed
   * rather than optional: where two surfaces meet in a crease the turn does not fall off under
   * refinement, because a crease is a crease at every resolution, so nothing but the cap stops
   * those cells refining forever. A feature thinner than a few cells also keeps asking for more,
   * until its two sides stop reaching each other through a neighbour stencil, which takes it to
   * around three or four cells thick.
   * @return Tags on each level, coarsest first.
   */
  Vector<IntVectSet>
  getCurvatureTags(Vector<Vector<Box>>& a_regions,
                   const ProblemDomain& a_coarsestDomain,
                   const Vector<int>&   a_refRatios,
                   const IntVect&       a_tileSize,
                   const IntVect&       a_maxBlockSize,
                   const RealVect&      a_probLo,
                   const Real&          a_coarsestDx,
                   const Real&          a_angle,
                   const int            a_growth,
                   const int            a_maxDepth) const;

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
   * @brief Refinement between the cells the surface is reconstructed on and the ones generated.
   */
  int m_geometryRefinement;

  /**
   * @brief Base name for writing the reconstructed surface, or empty for none.
   */
  std::string m_geometrySurfaceFile;

  /**
   * @brief Boxes each level was swept over, coarsest first, and so how much of it to carry.
   */
  Vector<Vector<Box>> m_coverageRegions;

  /**
   * @brief Domain the first entry of m_coverageRegions belongs to.
   */
  ProblemDomain m_coverageDomain;

  /**
   * @brief Cells carried beyond what the coverage boxes ask for, on every level.
   */
  int m_coverageBuffer;

  /**
   * @brief Whether each level of the index space is generated in turn rather than coarsened down
   * from the finest one.
   */
  bool m_generateEveryLevel;

  /**
   * @brief The generators, one per phase, kept for as long as the index space they built.
   * @details The index space borrows a raw pointer to its generator and keeps it: EBISLevel holds
   * a const GeometryService* that it never owns. Everything Chombo does with it today happens
   * while the index space is being defined, so a generator that lived only that long was enough.
   * Refining a cell that the index space did not carry has to ask the generator afterwards, which
   * means it has to outlive the call that built the geometry.
   */
  Vector<RefCountedPtr<GeometryService>> m_geoservers;

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
   * @brief Tag the cells in a region that the resolution does not do justice to.
   * @details Two things are looked for. A cell whose interface turns by more than a_angle against
   * a neighbour's is under-resolved because the surface bends inside it. A cell that no corner
   * reports as cut but that the interface passes through anyway is under-resolved because the
   * feature is thinner than the cell, and only sampling the edge midpoints and the cell centre
   * finds it. That sampling is confined to cells whose smallest corner value does not exceed the
   * largest change across any of their edges, which bounds how far the interface could reach from
   * the values themselves rather than from assuming a signed distance function.
   *
   * Node values and edge crossings are evaluated once each and shared between the cells that meet
   * on them, so a cell costs one implicit function evaluation plus a root find on each of its
   * edges that carries a sign change. Normals are wanted one cell beyond the region so that every
   * cell in it can be compared against all of its neighbours.
   * @param[in,out] a_tags             Tagged cells, added to.
   * @param[in]     a_region           Cells to tag over.
   * @param[in]     a_domain           Domain of the level the region sits on.
   * @param[in]     a_implicitFunction Implicit function to reconstruct from.
   * @param[in]     a_probLo           Lower left corner of the domain.
   * @param[in]     a_dx               Grid spacing on this level.
   * @param[in]     a_angle            Turn, in degrees, above which a cell is tagged.
   */
  void
  tagUnderResolvedCells(IntVectSet&                  a_tags,
                        const Box&                   a_region,
                        const ProblemDomain&         a_domain,
                        const RefCountedPtr<BaseIF>& a_implicitFunction,
                        const RealVect&              a_probLo,
                        const Real&                  a_dx,
                        const Real&                  a_angle) const;

  /**
   * @brief Position along a cell edge where the implicit function changes sign.
   * @details Bisection, addressed by the edge rather than by either of the cells sharing it, so
   * that every cell reaching the edge is handed the same interval and gets the same root back.
   * @param[in] a_implicitFunction Implicit function to find the root of.
   * @param[in] a_lowPoint         Position of the edge's low end.
   * @param[in] a_dir              Direction the edge runs along.
   * @param[in] a_lowValue         Implicit function at the edge's low end.
   * @param[in] a_dx               Length of the edge.
   * @return Root, as a fraction of the edge measured from its low end.
   */
  Real
  edgeRoot(const RefCountedPtr<BaseIF>& a_implicitFunction,
           const RealVect&              a_lowPoint,
           const int                    a_dir,
           const Real&                  a_lowValue,
           const Real&                  a_dx) const;

  /**
   * @brief Set up the geometry generation tool for the gas phase.
   * @param[out] a_geoserver     Generator, owned here and borrowed by the index space.
   * @param[in]  a_finestDomain  Finest domain which will be used
   * @param[in]  a_probLo        Lower-left corner of simulation domain.
   * @param[in]  a_finestDx      Finest resolution which will be used.
   */
  void
  buildGasGeometry(RefCountedPtr<GeometryService>& a_geoserver,
                   const ProblemDomain&            a_finestDomain,
                   const RealVect&                 a_probLo,
                   const Real                      a_finestDx);

  /**
   * @brief Set up the geometry generation tool for the solid phase, i.e. the part inside the dielectrics.
   * @param[out] a_geoserver     Generator, owned here and borrowed by the index space.
   * @param[in]  a_finestDomain  Finest domain which will be used
   * @param[in]  a_probLo        Lower-left corner of simulation domain.
   * @param[in]  a_finestDx      Finest resolution which will be used.
   */
  void
  buildSolidGeometry(RefCountedPtr<GeometryService>& a_geoserver,
                     const ProblemDomain&            a_finestDomain,
                     const RealVect&                 a_probLo,
                     const Real                      a_finestDx);
};

#include <CD_NamespaceFooter.H>

#endif
