/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PolyhedralGeometryShop.H
 * @brief  Declaration of a geometry generator whose moments are integrals of a stored surface.
 * @author Robert Marskar
 */

#ifndef CD_POLYHEDRALGEOMETRYSHOP_H
#define CD_POLYHEDRALGEOMETRYSHOP_H

// Std includes
#include <string>

// Chombo includes
#include <BaseFab.H>
#include <IntVectSet.H>

// Our includes
#include <CD_ScanShop.H>
#include <CD_CutCellBody.H>
#include <CD_CutCellSurface.H>
#include <CD_PolyhedralEBGraph.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_NamespaceHeader.H>

class ComputationalGeometry;

/**
 * @brief Geometry generation whose moments are exact integrals of an explicit surface.
 * @details Takes the same topology as GeometryShop -- the same classification, the same graph and the same covered-cell
 * fix-ups -- and replaces only the quadrature. Where GeometryShop integrates the arc the interface really
 * follows, this assembles the closed polyhedron whose faces are bounded by straight chords between the edge
 * crossings and integrates that.
 *
 * The difference is deliberate and costs accuracy on curved geometry, where the apertures drop
 * from third order to second. What it buys is that the moments are the moments of a body, so
 * they satisfy the relations Chombo's own coarsening applies to them, and a subdivision of that
 * body sums back to it exactly.
 *
 * Cells are declined rather than approximated when the assembled polygons do not close. That is
 * fatal by default: it has never fired on any geometry tested, and a cell that trips it is an
 * edge case worth seeing rather than absorbing.
 */
class PolyhedralGeometryShop : public ScanShop
{
public:
  /**
   * @brief Disallowed default constructor.
   */
  PolyhedralGeometryShop() = delete;

  /**
   * @brief Constructor.
   * @param[in] a_localGeom    The implicit function for the geometry
   * @param[in] a_verbosity    Verbosity
   * @param[in] a_dx           Resolution on the finest level
   * @param[in] a_probLo       Physical coordinates of the lower-left corner of the domain
   * @param[in] a_finestDomain Finest grid level
   * @param[in] a_scanLevel    Level on which to initiate the load balancing sequence
   * @param[in] a_ebGhost      Number of EB ghost cells
   * @param[in] a_thrshdVoF    Volume fraction below which a cut cell is dropped
   * @param[in] a_strict       Whether a declined cell aborts the run
   */
  PolyhedralGeometryShop(const BaseIF&        a_localGeom,
                         const int            a_verbosity,
                         const Real           a_dx,
                         const RealVect&      a_probLo,
                         const ProblemDomain& a_finestDomain,
                         const ProblemDomain& a_scanLevel,
                         const int            a_ebGhost,
                         const Real           a_thrshdVoF,
                         const bool           a_strict);

  /**
   * @brief Destructor.
   */
  ~PolyhedralGeometryShop() override;

  /**
   * @brief Disallow copy construction.
   */
  PolyhedralGeometryShop(const PolyhedralGeometryShop& a_other) = delete;

  /**
   * @brief Disallow copy assignment.
   */
  PolyhedralGeometryShop&
  operator=(const PolyhedralGeometryShop& a_other) = delete;

  /**
   * @brief Fill the graph and the nodes for one box.
   * @details Classifies the cells, reuses GeometryShop's fix-up for the regular cells bordering a covered one, and then
   * builds each remaining irregular cell's node from the polyhedron rather than from a quadrature of the
   * implicit function.
   * @param[in,out] a_regIrregCovered Regular/Covered/Irregular cells
   * @param[in,out] a_nodes           Nodes
   * @param[in]     a_validRegion     Grid region
   * @param[in]     a_ghostRegion     Grid region, including ghost cells
   * @param[in]     a_domain          Domain
   * @param[in]     a_probLo          Lower-left corner of computational domain
   * @param[in]     a_dx              Grid resolution
   * @param[in]     a_di              Grid index
   */
  void
  fillGraph(BaseFab<int>&        a_regIrregCovered,
            Vector<IrregNode>&   a_nodes,
            const Box&           a_validRegion,
            const Box&           a_ghostRegion,
            const ProblemDomain& a_domain,
            const RealVect&      a_probLo,
            const Real&          a_dx,
            const DataIndex&     a_di) const override;

  /**
   * @brief Say which grids the polyhedra are built over, and which phase this generator serves on them.
   * @param[in] a_compGeom Builder holding the boxes and their classifications per level.
   * @param[in] a_phase Phase whose classifications and implicit function this generator uses.
   */
  void
  setGrids(const ComputationalGeometry& a_compGeom, const phase::which_phase a_phase) noexcept;

  /**
   * @brief Build the graph of every tiled level: the cut cells, their surfaces, and what every face opens onto.
   * @details One PolyhedralEBGraph per level above the geometry's start level, over that level's cut tiles, each
   * defined from this phase's implicit function and then linked to the level above it so that the coarse side of
   * every level boundary knows which of its cells and faces the finer level describes. The levels at and below
   * the start level carry no graph. Collective. Requires setGrids.
   */
  void
  buildGraphs();

  /**
   * @brief The graph of a level, undefined at and below the start level and before buildGraphs.
   * @param[in] a_level Level, 0 being the coarsest of the geometry's levels.
   * @return The graph.
   */
  const PolyhedralEBGraph&
  getGraph(const int a_level) const noexcept;

  /**
   * @brief Build the polyhedra over the grids, check that they form a surface, and write it if asked.
   * @details Every cut cell of every irregular box is reconstructed at its own level, cells the next finer level
   * carries are skipped, and a face shared with such a cell is restricted to the finer chords and the
   * interface closed. The interface triangles of every level are then checked by sanityCheck, and written as
   * one STL per level and one composite when PolyhedralGeometryShop.write_stl is set. Collective; every rank
   * reconstructs its share of the boxes. Requires setGrids.
   */
  void
  verifySurface() const;

  /**
   * @brief Check that the interface triangles form a surface that is closed except on the domain boundary.
   * @details Every polyhedron has a neighbour across each edge of its interface unless that edge lies on the domain
   * boundary, so every interior edge of the triangulation must be used exactly twice. Vertices are welded at
   * s_weldSpacing times the finest spacing before edges are matched, which is far above the round-off two
   * cells produce for one shared vertex and far below any geometry. Gathers the triangles onto rank 0, where
   * an interior edge used once or more than twice is reported with its position and stops the run. Switched by
   * PolyhedralGeometryShop.sanity_check, on by default.
   * @param[in] a_facets This rank's interface triangles, nine reals each.
   */
  void
  sanityCheck(const Vector<Real>& a_facets) const;

  /**
   * @brief Whether a cut body holds too little solid to be a cut cell.
   * @details The counterpart, on the fluid side, of dropping a body whose fluid volume is below the threshold: a body
   * with neither solid volume nor boundary area is a regular cell. A cell whose interface lies in one of its
   * faces, wholly or in part, has no solid volume either but does have boundary area, and is not dust.
   * @param[in] a_body      Body of the cell.
   * @param[in] a_threshold Volume-fraction threshold; nothing is dust at zero.
   * @return True if the cell should be treated as regular.
   */
  static bool
  isDust(const PolyhedralEB::CutCellBody& a_body, const Real a_threshold) noexcept;

  /**
   * @brief Evaluate the implicit function at every node of a region.
   * @details One evaluation per node rather than per cell corner, since the cells meeting at a node share it and must
   * agree on its value. Values are read through snappedValue.
   * @param[in]  a_function   Implicit function
   * @param[out] a_nodeValues Node-centred values over a_region
   * @param[in]  a_region     Cell-centred region whose nodes are wanted
   * @param[in]  a_probLo     Lower-left corner of computational domain
   * @param[in]  a_dx         Grid resolution
   */
  static void
  fillNodeValues(const BaseIF&   a_function,
                 BaseFab<Real>&  a_nodeValues,
                 const Box&      a_region,
                 const RealVect& a_probLo,
                 const Real&     a_dx);

  /**
   * @brief Size the edge-crossing cache for a region and mark every edge as not yet bisected.
   * @details One fab per direction, over the edges of the region's cells, addressed by the node at the low end of the
   * edge as edgeCrossing addresses them.
   * @param[out] a_intercept Edge-centred cache, one fab per direction
   * @param[in]  a_region    Cell-centred region whose edges are wanted
   */
  static void
  defineIntercepts(BaseFab<Real> a_intercept[SpaceDim], const Box& a_region);

  /**
   * @brief Find where the implicit function changes sign along a cell edge.
   * @details Bisection between the two endpoints, looked up in a cache first. Every cut edge is reached from both of
   * its adjoining faces and from each neighbouring cell, always with the same endpoints, so the root is found
   * once and reused. That is what makes two cells agree on a shared aperture bit for bit.
   * @param[in]     a_function  Implicit function
   * @param[in,out] a_intercept Edge-centred cache, one component per direction
   * @param[in]     a_cell      Cell the edge is being read from
   * @param[in]     a_edge      Edge index within the cell
   * @param[in]     a_lo        Implicit function at the edge's low corner
   * @param[in]     a_hi        Implicit function at the edge's high corner
   * @param[in]     a_probLo    Lower-left corner of computational domain
   * @param[in]     a_dx        Grid resolution
   * @return Position along the edge, parameterised from its low corner.
   */
  static Real
  edgeCrossing(const BaseIF&   a_function,
               BaseFab<Real>   a_intercept[SpaceDim],
               const IntVect&  a_cell,
               const int       a_edge,
               const Real      a_lo,
               const Real      a_hi,
               const RealVect& a_probLo,
               const Real&     a_dx);

  /**
   * @brief Read one cell's corner values off the node values.
   * @details Corner c of the cell is the node offset by ((c >> d) & 1) in each direction d.
   * @param[out] a_surface    Surface whose corners are filled. Its crossings are left alone.
   * @param[in]  a_nodeValues Node-centred implicit function values
   * @param[in]  a_cell       The cell
   */
  static void
  fillCorners(PolyhedralEB::CutCellSurface& a_surface, const BaseFab<Real>& a_nodeValues, const IntVect& a_cell);

  /**
   * @brief Gather the corner values and edge crossings one cell's surface is built from.
   * @details The corners come from the shared node values and every edge whose ends disagree gets its crossing through
   * the shared cache, so every caller that reconstructs a cell -- the graph, the surface pass and the grid
   * builder's curvature test -- builds the same cell.
   * @param[in]     a_function   Implicit function
   * @param[in,out] a_intercept  Edge-centred crossing cache
   * @param[out]    a_surface    The values the cell's body is reconstructed from
   * @param[in]     a_nodeValues Node-centred implicit function values
   * @param[in]     a_cell       The cell
   * @param[in]     a_probLo     Lower-left corner of computational domain
   * @param[in]     a_dx         Grid resolution
   */
  static void
  buildSurface(const BaseIF&                 a_function,
               BaseFab<Real>                 a_intercept[SpaceDim],
               PolyhedralEB::CutCellSurface& a_surface,
               const BaseFab<Real>&          a_nodeValues,
               const IntVect&                a_cell,
               const RealVect&               a_probLo,
               const Real&                   a_dx);

  /**
   * @brief Bisect one cell edge for the crossing on it.
   * @details An endpoint at exactly zero is not taken to be the crossing: the function may be zero along a stretch of
   * the edge (a face of the geometry lying in the edge's plane), and the crossing is where it leaves that
   * stretch. The bisection walks through zeros, which the fluid predicate reads as solid, and a root that
   * lands within s_snapTolerance of an endpoint that is exactly zero is put exactly on it.
   * @param[in] a_function Implicit function.
   * @param[in] a_edgeIV   Node at the low end of the edge.
   * @param[in] a_dir      Direction the edge runs in.
   * @param[in] a_loValue  Implicit function at the low end.
   * @param[in] a_hiValue  Implicit function at the high end.
   * @param[in] a_probLo   Lower-left corner of the domain.
   * @param[in] a_dx       Grid spacing.
   * @return Position along the edge, in [0,1].
   */
  static Real
  edgeRoot(const BaseIF&   a_function,
           const IntVect&  a_edgeIV,
           const int       a_dir,
           const Real      a_loValue,
           const Real      a_hiValue,
           const RealVect& a_probLo,
           const Real      a_dx) noexcept;

  /**
   * @brief Fraction of a cell width within which an implicit-function value is read as exactly zero.
   */
  static constexpr Real s_snapTolerance = 1.0E-12;

  /**
   * @brief Distance along an edge, from an endpoint that is exactly zero, within which a root is put on it.
   * @details A surface tangent to a node reads as zero for a stretch either side of it under the value snap, of length
   * about sqrt(2 R tol / dx) for a radius of curvature R -- a few parts in 1e6 of an edge for any R that is
   * not itself below the cell width. A root that close to a corner the surface passes through is that
   * stretch's end, not a second crossing; a genuine second crossing within 1e-4 of a cell width of the first
   * would be a feature no cell resolves.
   */
  static constexpr Real s_rootSnap = 1.0E-4;

  /**
   * @brief Read the implicit function at a point, with a value within s_snapTolerance of a cell width of
   * zero read as exactly zero.
   * @details A face of the geometry lying in a node plane evaluates to zero there only up to round-off: the same face
   * comes back as plus zero, minus zero and a few parts in 1e19 at neighbouring nodes, and cells on either
   * side of it would then be classified by which of those they were handed. Reading everything within a fixed
   * fraction of the cell width of zero as zero puts every such node on one side by rule, in the corner values
   * and in the bisection alike.
   * @param[in] a_function Implicit function.
   * @param[in] a_point    Point to evaluate at.
   * @param[in] a_dx       Grid spacing that sets the scale.
   * @return The value, or exactly zero.
   */
  static Real
  snappedValue(const BaseIF& a_function, const RealVect& a_point, const Real a_dx) noexcept;

protected:
  /**
   * @brief Fraction of the finest grid spacing within which two vertices are the same vertex in sanityCheck.
   */
  static constexpr Real s_weldSpacing = 1.0E-8;

  /**
   * @brief The geometry whose grids the polyhedra are built over, or null before setGrids.
   */
  const ComputationalGeometry* m_compGeom;

  /**
   * @brief Graph of every level, indexed as the geometry indexes its levels; undefined where a level has no tiles.
   */
  Vector<RefCountedPtr<PolyhedralEBGraph>> m_graphs;

  /**
   * @brief Phase this generator serves on those grids.
   */
  phase::which_phase m_phase;

  /**
   * @brief Whether verifySurface writes the surface as STL. PolyhedralGeometryShop.write_stl, off by default.
   */
  bool m_writeSTL;

  /**
   * @brief Whether verifySurface runs sanityCheck. PolyhedralGeometryShop.sanity_check, on by default.
   */
  bool m_sanityCheck;

  /**
   * @brief Whether verifySurface reports its timings. PolyhedralGeometryShop.profile, off by default.
   */
  bool m_profile;

  /**
   * @brief Interface triangles of one level, from this rank's share of its irregular boxes, nine reals each.
   * @details Cells the next finer level carries are skipped. A cut cell with such a cell across one of its faces has
   * that face restricted to the finer cells' chords (restrictFace) and its interface rebuilt (closeInterface).
   * A coarse edge that carries no crossing while both finer edges under it do is a multi-valued coarse cell in
   * disguise and stops the run, as does a body that will not close.
   * @param[out] a_facets Triangle coordinates, appended to.
   * @param[in]  a_level  Level to collect from.
   */
  void
  collectFacets(Vector<Real>& a_facets, const int a_level) const;

  /**
   * @brief Write triangles gathered from every rank as one STL file.
   * @param[in] a_fileName File to write.
   * @param[in] a_name     Solid name inside the file.
   * @param[in] a_facets   This rank's triangles, nine reals each.
   */
  void
  writeSTL(const std::string& a_fileName, const std::string& a_name, const Vector<Real>& a_facets) const;

  /**
   * @brief Residual above which the stored moments and the ones PolyGeom recomputes disagree.
   */
  static constexpr Real s_divergenceTolerance = 1.0E-9;

  /**
   * @brief Whether a declined cell aborts the run rather than falling back.
   */
  bool m_strict;

  /**
   * @brief Volume fraction below which a cut cell carries no fluid worth keeping.
   */
  Real m_volumeThreshold;

  /**
   * @brief Write a cut cell's node from its body.
   * @details The arcs come from the covered set and the domain rather than from the moments, so the graph this produces
   * is the one GeometryShop would have produced. The apertures, the centroids and the volume fraction come
   * from the body.
   * @param[out] a_node             The node to fill
   * @param[in]  a_body             The cell's assembled polyhedron
   * @param[in]  a_regIrregCovered  Regular/Covered/Irregular cells
   * @param[in]  a_cell             The cell
   * @param[in]  a_domain           Domain
   */
  void
  fillNode(IrregNode&                       a_node,
           const PolyhedralEB::CutCellBody& a_body,
           const BaseFab<int>&              a_regIrregCovered,
           const IntVect&                   a_cell,
           const ProblemDomain&             a_domain) const;
};

#include <CD_NamespaceFooter.H>

#endif
