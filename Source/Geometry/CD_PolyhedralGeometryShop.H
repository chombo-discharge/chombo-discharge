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
   * @brief What to do with a cell this generator cannot describe at the refinement the mesh will ask for.
   * @details A cell whose fluid falls into more than one piece when it is cut for a refinement cannot be
   * carried single valued at that ratio, and there is no reformulation of the interior that avoids it: the
   * only representation guaranteed to keep every descendant single valued is a planar interface, and a planar
   * interface cannot meet the face chords a cell shares with its neighbours. So the cell is a choice rather
   * than a bug, and the choice belongs to whoever set up the run.
   */
  enum class UnresolvedAction
  {
    /**
     * @brief Report the cells and carry on, leaving it to the user to raise the geometric resolution.
     */
    Warn,

    /**
     * @brief Turn the cell solid, so that neither it nor anything cut from it is ever multi valued.
     */
    Cover,

    /**
     * @brief Keep the pieces in one body, taking every moment from their sum and moving only the volume centroid.
     * @details For a cell whose interface is one sheet twisted into a saddle, so that cutting it for a refinement
     * leaves a child holding its fluid in two pieces. The pieces are joined a cell away, so the saddle is resolved
     * geometry and is merged rather than covered. A cell holding more than one sheet is still filled.
     *
     * Conserved: the volume fraction, the apertures, the boundary area vector and the boundary centroid are the
     * exact sums of the pieces', satisfy the divergence identity, and coarsen back to the parent exactly. Face
     * centroids are left as the exact first moment of the open face, even where two patches put it on covered face
     * between them, since the flux stencils read it as a quadrature point.
     *
     * Not conserved: the scalar boundary area, the length of the summed area vector, is smaller than the area of the
     * interface itself when the pieces face different ways; the interface's own area is kept beside it. And the
     * volume centroid, which is the one moment moved.
     *
     * The repair: the volume centroid is the exact first moment of the merged fluid, a weighted average of the
     * pieces' centroids, and it can fall in the solid between them. Where it lies in the fluid it is kept. Where it
     * does not, it is replaced by the nearest point of the fluid, moved slightly inside it. That point no longer
     * sums back to the parent's first moment, so coarsening this body's centroid misses the parent's by the distance
     * moved, weighted by this body's volume.
     */
    Collapse
  };

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
   * face the finer level describes, reconstructed for the purpose. An edge in the domain boundary is used once, an
   * edge in a face the coarser level describes is checked from that side, and an edge in a face onto a cell the graph
   * filled has no partner to meet, since that face is closed. Every rank checks its own cells, vertices are welded at
   * s_weldSpacing times the finest spacing, the counts are summed, and a violation stops the run. A regular cell
   * sharing a face with a covered one stops the run as well: the four nodes of that face belong to both cells and
   * cannot be fluid for one and solid for the other, so the pair can only come from a classification that disagrees
   * with itself, and there would be no cut cell to carry the interface between them. A filled cell is covered whatever
   * its nodes read, which is why the graph keeps a neighbour of one as a cut cell holding the whole of itself rather
   * than a regular one; the test is what says it did. That test is dimension-independent and runs on every graph,
   * including the cells beyond the tiles, which is where the edge of the tiled region meets what no tile carries.
   * Switched by PolyhedralGeometryShop.sanity_check, on by default in a build that keeps its assertions. Requires
   * buildGraphs.
   */
  void
  sanityCheck() const;

  /**
   * @brief Report the cells that cannot be cut for a refinement without falling into more than one piece.
   * @details Asked of every cut cell the level itself describes, at the ratio m_maxRefinement names, and of the
   * body the consumers build -- which on a refinement boundary is the one whose faces were taken from the finer
   * level, since stitching lengthens the interface loop and it is that loop the question is about. A cell whose
   * interface is planar is exempt without being cut, because its fluid is then convex and convexity survives
   * every ratio.
   *
   * What happens to a failing cell is m_unresolvedAction's to decide. Only Warn is wired: the cells are named
   * and the run continues, leaving it to whoever set it up to raise the geometric resolution. Runs regardless
   * of m_sanityCheck, since it reports on the geometry the run will actually use rather than on this
   * generator's own consistency. Requires buildGraphs.
   */
  void
  reportUnresolvedRefinement() const;

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
   * @brief State of a cell of the finer level that the finer level does not hold, in a copy of its states.
   */
  static constexpr int s_absent = 2;

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
   * tells the generator nothing it uses and only reports. Either default can be overridden from the inputs.
   */
  bool m_sanityCheck;

  /**
   * @brief Whether verifySurface reports its timings. PolyhedralGeometryShop.profile, off by default.
   */
  bool m_profile;

  /**
   * @brief How far below the geometry's own grids a cell must stay single valued.
   * PolyhedralGeometryShop.max_refinement, one -- that is, off -- by default.
   * @details A power of two. One asks for nothing and skips the check entirely; two asks that the cells one
   * refinement below the geometry grid are each a single piece of fluid, four that the cells two below are,
   * and so on. The property is not inherited -- a cell whose children are all single valued can still have a
   * grandchild that is not -- so this has to name the ratio the mesh will actually use rather than being
   * derived from one level.
   *
   * Off by default because the test costs a subdivision of every cut cell whose interface is not planar, which
   * is two to three times what building the geometry costs in the first place, and because nothing downstream
   * reads a subdivided cell yet. It becomes worth its price once it does.
   */
  int m_maxRefinement;

  /**
   * @brief What becomes of a cell that fails that test. PolyhedralGeometryShop.unresolved_cells, warn by
   * default.
   */
  UnresolvedAction m_unresolvedAction;

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
   * @brief Stitch every level of a set of graphs to the level above it, from the graphs alone.
   * @details In two dimensions no face is stitched and every level's overrides stay empty. Collective.
   * @param[in] a_graphs Graphs, one per level, linked to one another; undefined entries are skipped.
   */
  void
  stitchGraphs(const Vector<RefCountedPtr<PolyhedralEBGraph>>& a_graphs) const;

  /**
   * @brief The closure check over a given set of graphs, one per level as m_graphs holds them.
   * @param[in] a_graphs Graphs, undefined entries skipped.
   */
  void
  sanityCheck(const Vector<RefCountedPtr<PolyhedralEBGraph>>& a_graphs) const;

#if CH_SPACEDIM == 3
  /**
   * @brief Whether the graph filled a cell: covered, with corners that disagree, and not carried by the finer level.
   * @details The one kind of covered neighbour a cut cell's face has to be closed onto. A covered cell whose corners
   * agree closes the face by its corners already, and a cell the finer level carries is described up there.
   * @param[in] a_states      Cell states around the cell.
   * @param[in] a_refined     Refined mask around the cell.
   * @param[in] a_withSurface Cells around the cell that hold a surface.
   * @param[in] a_surfaces    Their surfaces.
   * @param[in] a_cell        The cell.
   * @return True if the cell was filled.
   */
  static bool
  isFilled(const BaseFab<signed char>&                 a_states,
           const BaseFab<signed char>&                 a_refined,
           const IntVectSet&                           a_withSurface,
           const IVSFAB<PolyhedralEB::CutCellSurface>& a_surfaces,
           const IntVect&                              a_cell);

  /**
   * @brief Which of a cell's edges the finer level describes.
   * @details An edge is shared by the cells that meet along it, and the crossing on it has to be one point for all of
   * them, so the question is asked of the edge: it is refined if any cell meeting it is.
   * @param[out] a_edgeRefined One flag per edge of the cell.
   * @param[in]  a_refined     Refined mask around the cell.
   * @param[in]  a_domain      Domain box of the cell's level.
   * @param[in]  a_cell        The cell.
   * @return True if any edge is refined.
   */
  static bool
  edgesRefined(bool                        a_edgeRefined[PolyhedralEB::CutCellSurface::s_numEdges],
               const BaseFab<signed char>& a_refined,
               const Box&                  a_domain,
               const IntVect&              a_cell);

  /**
   * @brief Stitch a cut cell's body to the finer level and to the covered cells around it.
   * @details Every refined edge takes its crossing from the half the finer level put it in, which rewrites a_surface.
   * The body is then defined from it -- as the whole cell where every corner is fluid -- every face onto a refined
   * neighbour is replaced by the chords of the children on that face, after the check that no coarse edge on it is
   * crossed twice at the finer spacing, and the interface is closed again, which also recomputes the moments. A face
   * onto a filled cell (a_filled) is closed instead, through CutCellBody::snapFace, and the footprint it gives up
   * becomes interface lying in the plane of that face. Every face so changed is recorded, in the order it was changed,
   * when overrides are given. A body that will not close, a face that cannot be restricted or closed, or a doubly
   * crossed edge stops the run.
   * @param[out]    a_body        The body.
   * @param[in,out] a_surface     The cell's surface; its refined edges' crossings are rewritten.
   * @param[in]     a_children    The cell's eight children in corner order, which must be given when any edge is
   * refined and describe at least the cell's refined faces and edges; null otherwise.
   * @param[in]     a_edgeRefined Which edges are refined, from edgesRefined.
   * @param[in]     a_filled      One where a cell around the cell was filled, from isFilled.
   * @param[in]     a_refined     Refined mask around the cell.
   * @param[in]     a_fineFilled  One where a cell of the finer level was filled, over the finer cells across the
   * cell's faces; a quadrant of a restricted face whose finer cell was filled is closed. Null where nothing is known of
   * the finer level's fills.
   * @param[in]     a_domain      Domain box of the cell's level.
   * @param[in]     a_cell        The cell.
   * @param[in,out] a_overrides   Overrides to record the changed faces in, with the cell begun; null records nothing.
   * @return Number of quadrants closed onto filled finer cells.
   */
  static int
  stitchBody(PolyhedralEB::CutCellBody&          a_body,
             PolyhedralEB::CutCellSurface&       a_surface,
             const PolyhedralEB::CutCellSurface* a_children,
             const bool*                         a_edgeRefined,
             const BaseFab<signed char>&         a_filled,
             const BaseFab<signed char>&         a_refined,
             const BaseFab<signed char>*         a_fineFilled,
             const Box&                          a_domain,
             const IntVect&                      a_cell,
             PolyhedralEB::CutCellFaceOverrides* a_overrides);

  /**
   * @brief Stitch a cut cell's body with its children reconstructed from the implicit function.
   * @details The reference the stitched graph is checked against under sanity_check: the same stitching as the stored
   * faces were made with, from children built at the finer spacing from the function rather than from the finer
   * level's stored surfaces. Nothing else calls the function to build a body.
   * @param[out] a_body    The body.
   * @param[in]  a_surface The cell's surface as the graph first built it.
   * @param[in]  a_graph   Graph of the cell's level.
   * @param[in]  a_filled  One where a cell around the cell was filled, from isFilled.
   * @param[in]  a_refined Refined mask around the cell.
   * @param[in]  a_cell    The cell.
   */
  void
  defineBodyFromFunction(PolyhedralEB::CutCellBody&          a_body,
                         const PolyhedralEB::CutCellSurface& a_surface,
                         const PolyhedralEBGraph&            a_graph,
                         const BaseFab<signed char>&         a_filled,
                         const BaseFab<signed char>&         a_refined,
                         const IntVect&                      a_cell) const;

  /**
   * @brief Build the body of a cut cell of a stitched graph from its stored surface and face overrides.
   * @details Defines the body from the surface, replaces every restricted face with its stored polygons, closes the
   * interface over them, and closes every closed face, in the order they were recorded -- which rebuilds the body
   * stitchBody made, bit for bit, from the cell's own data alone.
   * @param[out] a_body      The body.
   * @param[in]  a_surfaces  Surfaces of the cell's box, from the graph.
   * @param[in]  a_overrides Face overrides of the cell's box, from the graph.
   * @param[in]  a_cell      The cell.
   */
  static void
  defineBody(PolyhedralEB::CutCellBody&                  a_body,
             const IVSFAB<PolyhedralEB::CutCellSurface>& a_surfaces,
             const PolyhedralEB::CutCellFaceOverrides&   a_overrides,
             const IntVect&                              a_cell);

  /**
   * @brief Read the nodes and edges of a region of the finer level off its stored cells.
   * @details A node or edge belongs to every cell around it, so it is written from each cell of the region that holds a
   * surface, and two cells that disagree about one stop the run. A node that no such cell touches lies among cells
   * whose corners agree, and takes their sign: minus one where they are regular, one where they are covered. An edge
   * that no such cell touches carries no crossing. Cells the finer level does not hold write nothing.
   * @param[out] a_nodes     Node values over the nodes of a_region.
   * @param[out] a_crossings Crossings over the edges of a_region, addressed as defineIntercepts addresses them.
   * @param[in]  a_region    Cells of the finer level.
   * @param[in]  a_states    Cell states of the finer level over the region, s_absent where it holds none.
   * @param[in]  a_mask      One where a cell holds a surface.
   * @param[in]  a_surfaces  The surfaces.
   */
  static void
  readFineCells(BaseFab<Real>&                              a_nodes,
                BaseFab<Real>                               a_crossings[SpaceDim],
                const Box&                                  a_region,
                const BaseFab<signed char>&                 a_states,
                const BaseFab<signed char>&                 a_mask,
                const IVSFAB<PolyhedralEB::CutCellSurface>& a_surfaces);

  /**
   * @brief One cell's surface, read off node values and edge crossings that are already known.
   * @param[out] a_surface   The surface.
   * @param[in]  a_nodes     Node values.
   * @param[in]  a_crossings Edge crossings, addressed as defineIntercepts addresses them.
   * @param[in]  a_cell      The cell.
   */
  static void
  assembleSurface(PolyhedralEB::CutCellSurface& a_surface,
                  const BaseFab<Real>&          a_nodes,
                  const BaseFab<Real>           a_crossings[SpaceDim],
                  const IntVect&                a_cell);

  /**
   * @brief Copy a finer level's cells onto a coarser level's layout refined, with a ring around every box.
   * @details The cell states, a mask of the cells that hold a surface, the sets of those cells and their surfaces --
   * everything a cell of the finer level is rebuilt from, apart from faces it restricts to a level above it. Cells
   * the finer level does not hold read s_absent. Collective.
   * @param[in]  a_fine     The finer level.
   * @param[in]  a_grids    The coarser level's layout.
   * @param[in]  a_ring     Width of the ring around every box, in cells of the coarser level.
   * @param[out] a_states   Cell states.
   * @param[out] a_mask     One where a cell holds a surface.
   * @param[out] a_sets     The cells that hold a surface, per box.
   * @param[out] a_surfaces Their surfaces.
   */
  static void
  copyFinerCells(const PolyhedralEBGraph&                         a_fine,
                 const DisjointBoxLayout&                         a_grids,
                 const int                                        a_ring,
                 LevelData<BaseFab<signed char>>&                 a_states,
                 LevelData<BaseFab<signed char>>&                 a_mask,
                 LayoutData<IntVectSet>&                          a_sets,
                 LevelData<IVSFAB<PolyhedralEB::CutCellSurface>>& a_surfaces);

  /**
   * @brief Stitch one level: take its refined edges' crossings and its seam faces from the finer level's stored cells,
   * close its faces onto filled cells, write the crossings into its surfaces and the faces into its overrides.
   * @details The finer level's cell states and surfaces over the refined coarse boxes and a ring around them are
   * copied onto this level's layout refined, and every cut cell of every box and its ghost ring is stitched from them.
   * A quadrant of a seam face across which the finer level filled a cell is closed. A regular or covered cell beside
   * the finer level is made cut where the finer level describes the face between them otherwise than the cell's
   * corners do -- the node at the middle of the face on the other side, or a filled cell across a quadrant of a regular
   * one -- and stitched like any other; a node at the middle of such a face's edge that disagrees with the edge's ends
   * is an edge the finer level crosses twice, which is counted and reported rather than acted on. Every stitched
   * surface carries its role. Under sanity_check each body is also built from the stored data, and from the implicit
   * function where its children can say the same, and they must agree bit for bit. Collective.
   * @param[in,out] a_coarse The level.
   * @param[in]     a_fine   The finer level, or null where there is none.
   */
  void
  stitchLevel(PolyhedralEBGraph& a_coarse, const PolyhedralEBGraph* a_fine) const;
#endif

  /**
   * @brief Write the surface, the boxes it was built over and the boxes that split, as one HDF5 file per phase.
   * @details One group per level, holding the level's triangles as an indexed mesh -- the vertices once each and a
   * triple of indices per triangle -- rather than as a soup of coordinates, which is what a viewer wants and what
   * an XDMF description can point at. The boxes of the level and the boxes that split on it go in the same group,
   * so that the surface can be read against the grids it was built on without a second file. The distance within
   * which two positions are one vertex is written alongside them. An XDMF description of the file is written next
   * to it. Collective: every rank writes its own stretch of every dataset, so a vertex on a boundary between two
   * ranks is written by each, since vertices are identified within a rank and not across them.
   * @param[in] a_fileName File to write, without a suffix.
   * @param[in] a_facets   This rank's triangles of each level, nine reals each.
   */
  void
  writeSurface(const std::string& a_fileName, const Vector<Vector<Real>>& a_facets) const;

  /**
   * @brief Identify the vertices of a triangle soup, giving each one once and each triangle three indices.
   * @details Two positions closer than the tolerance are one vertex. Candidates are taken from the lattice cell of
   * the tolerance a position falls in and the twenty-six around it, which is what keeps this out of the quadratic
   * scan a surface meeting the mesh squarely would otherwise provoke, with every candidate still measured against
   * so that two positions either side of a lattice wall are not told apart.
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
