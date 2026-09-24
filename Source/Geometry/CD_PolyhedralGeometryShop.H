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
   * @details One PolyhedralEBGraph per level from the geometry's start level up, over that level's cut tiles, each
   * defined from this phase's implicit function with one ghost cell and then linked to the level above it so that
   * the coarse side of every level boundary knows which of its cells and faces the finer level describes. The levels
   * below the start level carry no graph. Collective. Requires setGrids.
   */
  void
  buildGraphs();

  /**
   * @brief The graph of a level, undefined below the start level. Stops the run before buildGraphs.
   * @param[in] a_level Level, 0 being the coarsest of the geometry's levels.
   * @return The graph.
   */
  const PolyhedralEBGraph&
  getGraph(const int a_level) const noexcept;

  /**
   * @brief Build the polyhedra over the grids, check that they form a surface, and write it if asked.
   * @details Every cut cell of every irregular box is reconstructed at its own level, cells the next finer level
   * carries are skipped, and a face shared with such a cell is restricted to the finer chords and the
   * interface closed. The interface triangles of every level are then checked by sanityCheck when
   * PolyhedralGeometryShop.sanity_check is set, and written to the geometry directory as an indexed surface mesh
   * when PolyhedralGeometryShop.write_surface is set. Collective; every rank reconstructs its share of the boxes
   * and writes its own share of the mesh. Requires setGrids.
   */
  void
  verifySurface() const;

  /**
   * @brief Check that the interface triangles form a surface that is closed except on the domain boundary.
   * @details A triangle edge can be shared only by cells that touch, so every edge of a cell's interface must be
   * used exactly twice among the interface triangles of the cell and of its 3^D neighbourhood: on the cell's own
   * level through the ghost surfaces the graph holds, and across a level boundary through the finer cells behind a
   * face the finer level describes, reconstructed for the purpose. An edge in the domain boundary is used once, and
   * an edge in a face the coarser level describes is checked from that side. Every rank checks its own cells,
   * vertices are welded at s_weldSpacing times the finest spacing, the counts are summed, and a violation stops the
   * run. A regular cell sharing a face with a covered one stops the run as well: the four nodes of that face
   * belong to both cells and cannot be fluid for one and solid for the other, so the pair can only come from a
   * classification that disagrees with itself, and there would be no cut cell to carry the interface between
   * them. That test is dimension-independent and runs on every graph, including the cells beyond the tiles, which
   * is where the edge of the tiled region meets what no tile carries. Switched by
   * PolyhedralGeometryShop.sanity_check, off by default. Requires buildGraphs.
   */
  void
  sanityCheck() const;

  /**
   * @brief Copy every graph onto a layout that splits each of its tiles into octants with the ranks shuffled, and then
   * onto the tiles themselves with the ranks shuffled; run sanityCheck over the copies each time, copy them back onto
   * the original layouts and require equality bit for bit.
   * @details A stand-in, until the index space is generated over the graphs, for what a consumer on another
   * layout does: take the graph's data through a Copier and read it there. Switched by
   * PolyhedralGeometryShop.test_copy. Collective. Requires buildGraphs.
   */
  void
  testGraphCopy() const;

  /**
   * @brief Whether a cut body holds nothing at all and is therefore a regular cell.
   * @details A node that lands exactly on the interface is solid by rule, so the cells around it classify as cut
   * from their corners while their bodies hold no solid and no interface whatsoever -- the reconstruction finds
   * an empty polyhedron. Those cells are regular, and dropping them discards nothing, since there is no surface
   * in them for a neighbour to meet. Emptiness is exact rather than thresholded: a body with any solid in it,
   * however little, keeps its cell, because its interface is what closes the surface against its neighbours, and
   * a rule that discarded it by size would leave a hole the size of what it discarded. A cell whose interface
   * lies in one of its faces has no solid volume either, but it does have boundary area, and it is not dust.
   * @param[in] a_body Body of the cell.
   * @return True if the cell should be treated as regular.
   */
  static bool
  isDust(const PolyhedralEB::CutCellBody& a_body) noexcept;

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
   * @brief Fraction of a cell width the implicit function is read at either side of a point, to tell a surface
   * through the point from a zero no surface passes through.
   * @details Large enough that a function of unit gradient reads well clear of s_snapTolerance beside the point, and
   * small enough that the reading belongs to the point rather than to its neighbourhood.
   */
  static constexpr Real s_probeSpacing = 1.0E-6;

  /**
   * @brief Read the implicit function at a point, with a value within s_snapTolerance of a cell width of
   * zero read as exactly zero, and a zero that is a tangency read from beside the point instead.
   * @details A face of the geometry lying in a node plane evaluates to zero there only up to round-off: the same face
   * comes back as plus zero, minus zero and a few parts in 1e19 at neighbouring nodes, and cells on either
   * side of it would then be classified by which of those they were handed. Reading everything within a fixed
   * fraction of the cell width of zero as zero puts every such node on one side by rule, in the corner values
   * and in the bisection alike. A zero that no surface passes through is then read through resolveTangency, so
   * that nothing downstream sees it.
   * @param[in] a_function Implicit function.
   * @param[in] a_point    Point to evaluate at.
   * @param[in] a_dx       Grid spacing that sets the scale.
   * @return The value, or exactly zero.
   */
  static Real
  snappedValue(const BaseIF& a_function, const RealVect& a_point, const Real a_dx) noexcept;

  /**
   * @brief Read the implicit function beside a point where it is zero, when no surface passes through the point.
   * @details A composite function can be zero over a whole region that is not its surface: a difference of two bodies
   * sharing a face evaluates to the maximum of zero and something negative over the whole of that face. Taken
   * at face value such a zero is material by the fluid rule, and it invents an interface where the geometry
   * has none and displaces the interfaces of whatever else the cell holds. It is told apart from a surface by
   * the function's sign a probe spacing off the point, along each axis in turn. A surface through the point
   * changes the sign across it along any axis the surface is not parallel to, and since a surface has a normal
   * there is always such an axis, so one sign change is enough to leave the zero alone. The zero is read
   * through only on the evidence of an axis whose two readings are both away from zero and of the same sign,
   * which is what a function that touches zero without a surface being there looks like; the largest such
   * reading stands in for the zero. One reading of zero is not evidence of anything -- it is the point's own
   * zero again, which is what an axis along a face of the geometry returns, and what both axes along the
   * faces meeting at a convex edge return.
   * @param[in] a_function Implicit function.
   * @param[in] a_point    Point the function is zero at.
   * @param[in] a_dx       Grid spacing that sets the probe spacing.
   * @return The reading beside the point, or exactly zero if a surface passes through it.
   */
  static Real
  resolveTangency(const BaseIF& a_function, const RealVect& a_point, const Real a_dx) noexcept;

protected:
  /**
   * @brief Fraction of the finest grid spacing within which two vertices are the same vertex in sanityCheck.
   */
  static constexpr Real s_weldSpacing = 1.0E-6;

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
   * @brief Whether verifySurface writes the surface. PolyhedralGeometryShop.write_surface, off by default.
   */
  bool m_writeSurface;

  /**
   * @brief Whether verifySurface runs sanityCheck. PolyhedralGeometryShop.sanity_check.
   * @details On by default in a build that keeps its assertions and off in one that does not, since the check
   *          tells the generator nothing it uses and only reports. Either default can be overridden from the
   *          inputs.
   */
  bool m_sanityCheck;

  /**
   * @brief Whether verifySurface reports its timings. PolyhedralGeometryShop.profile, off by default.
   */
  bool m_profile;

  /**
   * @brief Whether verifySurface runs testGraphCopy. PolyhedralGeometryShop.test_copy, off by default.
   */
  bool m_testCopy;

  /**
   * @brief Whether the member functions that run once per level or per phase announce themselves in pout.
   * PolyhedralGeometryShop.verbose, off by default. The functions that run per cell never announce themselves.
   */
  bool m_verbose;

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
   * @brief The closure check over a given set of graphs, one per level as m_graphs holds them.
   * @param[in] a_graphs Graphs, undefined entries skipped.
   */
  void
  sanityCheck(const Vector<RefCountedPtr<PolyhedralEBGraph>>& a_graphs) const;

#if CH_SPACEDIM == 3
  /**
   * @brief Build the body of a cut cell of a graph, restricted on every face the finer level describes.
   * @details The body is defined from the cell's surface, and every face onto a refined neighbour is replaced by the
   * chords of the cell's own children, reconstructed from the implicit function at the finer spacing, after the
   * check that no coarse edge on that face is crossed twice at the finer spacing; the interface is closed again
   * afterwards. A body that will not close, a face that cannot be restricted, or a doubly crossed edge stops the
   * run.
   * @param[out] a_body     The body.
   * @param[in]  a_graph    Graph of the cell's level.
   * @param[in]  a_surfaces Surfaces of the cell's box, from the graph.
   * @param[in]  a_refined  Refined mask of the cell's box, from the graph.
   * @param[in]  a_cell     The cell.
   */
  void
  defineBody(PolyhedralEB::CutCellBody&                  a_body,
             const PolyhedralEBGraph&                    a_graph,
             const IVSFAB<PolyhedralEB::CutCellSurface>& a_surfaces,
             const BaseFab<signed char>&                 a_refined,
             const IntVect&                              a_cell) const;
#endif

  /**
   * @brief Write the surface, the boxes it was built over and the boxes that split, as one HDF5 file per phase.
   * @details One group per level, holding the level's triangles as an indexed mesh -- the vertices once each and a
   *          triple of indices per triangle -- rather than as a soup of coordinates, which is what a viewer wants
   *          and what an XDMF description can point at. The boxes of the level and the boxes that split on it go
   *          in the same group, so that the surface can be read against the grids it was built on without a second
   *          file. Vertices are identified the way sanityCheck identifies them, by distance rather than on a
   *          lattice, and the distance used is written alongside them. An XDMF description of the file is written
   *          next to it. Collective: every rank sends its triangles to the master, which writes.
   * @param[in] a_fileName File to write, without a suffix.
   * @param[in] a_facets   This rank's triangles of each level, nine reals each.
   */
  void
  writeSurface(const std::string& a_fileName, const Vector<Vector<Real>>& a_facets) const;

  /**
   * @brief Identify the vertices of a triangle soup, giving each one once and each triangle three indices.
   * @details Two positions closer than the tolerance are one vertex. They are found by sorting on the first
   *          coordinate and scanning back over the window the tolerance allows, which is what makes this
   *          linear-ish rather than quadratic.
   * @param[in]  a_facets       Triangles, nine reals each.
   * @param[in]  a_tolerance    Distance within which two positions are the same vertex.
   * @param[out] a_vertices     The vertices, three reals each.
   * @param[out] a_connectivity Three vertex indices per triangle.
   */
  static void
  indexFacets(const Vector<Real>& a_facets,
              const Real          a_tolerance,
              std::vector<Real>&  a_vertices,
              std::vector<int>&   a_connectivity);

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
