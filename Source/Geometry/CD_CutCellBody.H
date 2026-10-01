/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CutCellBody.H
 * @brief  Declaration of a closed polyhedron bounding the fluid in one cut cell.
 * @author Robert Marskar
 */

#ifndef CD_CUTCELLBODY_H
#define CD_CUTCELLBODY_H

// Std includes
#include <cstdint>

// Chombo includes
#include <RealVect.H>
#include <LoHiSide.H>

// Our includes
#include <CD_CutCellFaceOverrides.H>
#include <CD_CutCellSurface.H>
#include <CD_PolyhedralEBUtils.H>
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

/**
 * @brief The closed polyhedron bounding the fluid in one cut cell, and its moments.
 * @details The body is assembled from two kinds of polygon. Each cell face contributes the
 * fluid part of that face, found by walking the face's circuit of corners and edges and closing
 * the region off with a straight chord between the crossings. The interface contributes
 * triangles spanning the loop those crossings form, which is not planar in general and is what
 * closes the body.
 *
 * The moments are then exact integrals of that polyhedron rather than quadratures of the
 * implicit function, which is what makes the divergence identity
 * kappa = (1/SpaceDim)*[sum(alpha)/2 - a_B*(x_B . n)] hold to rounding, and what makes the
 * eight children of a subdivided body sum back to it exactly.
 *
 * The price is the chord: where Chombo integrates the curved arc the interface really follows,
 * this integrates the straight line between its endpoints. In two dimensions the interface in a
 * cut cell genuinely is that line and nothing is lost.
 */
class CutCellBody
{
public:
  /**
   * @brief How a cell relates to the interface.
   */
  enum class Kind
  {
    Regular,
    Covered,
    Cut
  };

  /**
   * @brief Constructor. Leaves the body empty; call define before reading any moment.
   */
  CutCellBody() noexcept;

  /**
   * @brief Destructor.
   */
  ~CutCellBody() noexcept = default;

  /**
   * @brief Build the body from the values the surface is reconstructed from.
   * @details Fails, rather than returning an unusable body, when the polygons do not close or
   * when the resulting volume leaves [0,1]. Both are checked rather than assumed, because a
   * folded or self-intersecting patch can pass every individual moment test.
   * @param[in] a_surface Corner values and edge crossings in the cell's own frame.
   * @return True if the body closed and its moments are usable.
   */
  bool
  define(const CutCellSurface& a_surface) noexcept;

  /**
   * @brief Build the body of a cut cell from its surface and its face overrides, as the polyhedral graph stores them.
   * @details The body is defined from the surface -- as the whole cell where every corner is fluid -- and the cell's
   * face overrides are applied in the order they were recorded: in three dimensions every face restricted to the finer
   * level, the interface closed over them, then every closed face; in two, every closed half, then every closed face.
   * That is the order the graph's stitching made them in, so the body comes out as the graph builds it, bit for bit.
   * A cell with no entry in the overrides is its surface's body.
   * @param[in] a_surface   The cell's surface.
   * @param[in] a_overrides Face overrides holding the cell's entry, if it has one.
   * @param[in] a_cell      The cell, as the overrides know it.
   * @return True if every step closed.
   */
  bool
  define(const CutCellSurface& a_surface, const CutCellFaceOverrides& a_overrides, const IntVect& a_cell) noexcept;

  /**
   * @brief How the cell relates to the interface.
   * @details The cell is cut as soon as two corners disagree under isFluid, regular if every corner
   * is fluid and covered otherwise. An exact zero is solid, so a cell whose face lies in the
   * interface is cut from the fluid side, full, with that face as its interface.
   *
   * A cell the surface enters as more than one sheet is still called cut here. Whether such a cell should be
   * filled instead -- see numSheets for why nothing it can hold describes the feature that put two sheets in
   * it -- is a decision for whoever knows whether a finer level carries the cell: one that does is left to it,
   * and one that does not is filled, which is what PolyhedralEBGraph::findUnresolvedCells decides.
   * @param[in] a_surface Corner values in the cell's own frame.
   * @return Regular, Covered or Cut.
   */
  static Kind
  classify(const CutCellSurface& a_surface) noexcept;

  /**
   * @brief How many separate sheets of interface the surface enters the cell as.
   * @details The crossings of a cut cell close into cycles, and each cycle is spanned by its own patch. One
   * cycle is the ordinary case: one patch, and the fluid either side of it. More than one means the cell holds
   * a feature thinner than itself -- a plate the nodes only see where it happens to swallow a corner, say --
   * and neither the corner values nor the crossings carry enough to say what that feature is. In two dimensions
   * the crossings pair into chords, one around each run of solid corners.
   *
   * A sheet whose solid side is nothing but corners at exactly zero is not counted. Such a corner lies on the
   * surface, and the fluid rule calls it solid, so the edges leaving it for fluid corners carry crossings -- but
   * every one of them lies on the corner itself, and the sheet they bound is a point, or a segment where two such
   * corners share an edge: the surface touches the cell there and does not enter it. A surface tangent to a node
   * plane along a line through nodes touches every cell along the line this way, often at two corners diagonal
   * across a face, which read as two sheets and would have a cell that holds nothing but fluid filled.
   * @param[in] a_surface Corner values and edge crossings in the cell's own frame.
   * @return Number of sheets that hold area, zero where the corners agree, or -1 where the crossings do not form
   * cycles.
   */
  static int
  numSheets(const CutCellSurface& a_surface) noexcept;

  /**
   * @brief Whether the fluid this body holds is one connected region.
   * @details Asked of the polygons rather than of the fluid: they bound a closed surface, so two polygons that
   * share an edge bound the same piece of fluid, and the fluid is in one piece exactly when walking that
   * sharing reaches every polygon from any one of them. A cell the graph builds should never hold two pieces --
   * the interface is anchored to the cell's edges, so it cannot close on itself inside the cell and leave a
   * cavity, and a cell the surface enters as more than one sheet is filled before any body is built. A body cut
   * from a coarser one is different: it was not built from an implicit function, and a single-sheet parent
   * whose interface is twisted can hand a child two pieces. That is what this is asked for, and it is also
   * what decides which polygons bound which piece, and so what opens onto a neighbouring cell.
   * @return True if every polygon is reachable from any other by shared edges, or if the body holds none.
   */
  bool
  isConnected() const noexcept;

  /**
   * @brief Whether every polygon of the interface lies in one plane.
   * @details A planar interface leaves the fluid the cell intersected with a half-space, which is convex. A
   * convex region meets every axis-aligned box in a convex, and therefore connected, set -- so a body that
   * answers true here has single-valued children at every refinement ratio, and their children in turn, without
   * any of them being built. It is the one property that settles the question for all ratios at once:
   * everything else is an interpolant of the same crossings and carries no such guarantee. A body whose
   * interface is a single triangle is planar by construction, which is most of the cut cells of a smooth
   * surface.
   * @return True if the interface polygons share a plane, or if the body has no interface at all.
   */
  bool
  interfaceIsPlanar() const noexcept;

  /**
   * @brief Cut this body into the bodies of the 2^SpaceDim cells one refinement finer that fill its cell.
   * @details A genuine subdivision of the polyhedron, not a reconstruction of it at the finer spacing: each
   * child is this body intersected with its own octant, so the children partition this body and their moments
   * sum back to its own exactly rather than to within the curvature over the cell. That is what makes a cell
   * cut from its parent join the parent's own description -- the two share the faces between them by
   * construction, which two independent reconstructions of the same region do not.
   *
   * The children come back in quadrant order: bit d of a child's index is the side it lies on in direction d,
   * zero low and one high, matching faceCorners and restrictFace. Each is expressed in its own cell's frame,
   * so a length here is twice one there. A child that the fluid misses entirely comes back empty, which reads
   * as covered; one the fluid fills comes back whole.
   *
   * A vertex of a child is not a crossing of the child's own edges, so the provenance those carry is dropped;
   * closeInterface already handles polygons in that state, since a restricted face leaves them the same way.
   * @param[out] a_children The 2^SpaceDim bodies, in quadrant order.
   * @return False if a child's fluid region is not connected, if a child would need more polygons or vertices
   * than it can hold, or if the children do not add back up to this body in every moment partitions checks.
   */
  bool
  subdivide(CutCellBody* a_children) const noexcept;

  /**
   * @brief Whether cutting this body for a refinement asks for a child whose fluid is in more than one piece.
   * @details The single-valued property is not inherited: a body whose eight children are all single valued can
   * still have a grandchild that is not, so asking one level down says nothing about two. This asks to the ratio
   * the mesh will actually use, and it is the test to put in front of any cell that a finer grid may subdivide,
   * whether or not the cell sits on a refinement boundary.
   *
   * A planar interface is answered from interfaceIsPlanar without building anything, since convexity settles
   * every ratio. Otherwise the body is cut in halves repeatedly, log2(a_refRat) times, and every body produced
   * on the way down is asked -- not only those at the final ratio. That is deliberate: the intermediate bodies
   * are what a recursive subdivision would actually build, so one of them being multi-valued is a failure even
   * when the cells the mesh ends up with are not. A direct a_refRat-way clip of this body is not the same
   * computation, because each round re-expresses its children in their own frames and closes their interfaces
   * there; when subdivide learns to clip by a ratio in one step, this should follow it.
   *
   * A body that cannot be cut at all -- one whose children would need more polygons or vertices than they can
   * hold -- is not reported here. That is a separate failure, and subdivide already returns it.
   * @param[in] a_refRat Refinement ratio to test, a power of two and at least two.
   * @return True if any body produced on the way to that ratio holds fluid in more than one piece.
   */
  bool
  hasMultiValuedChildren(const int a_refRat) const noexcept;

  /**
   * @brief Volume fraction.
   * @return Fraction of the cell occupied by fluid.
   */
  Real
  volumeFraction() const noexcept;

  /**
   * @brief Volume centroid, relative to the cell centre, in units of the cell size.
   * @return Centroid of the fluid, relative to the cell centre.
   */
  const RealVect&
  volumeCentroid() const noexcept;

  /**
   * @brief Area fraction of one cell face.
   * @param[in] a_dir  Coordinate direction of the face normal.
   * @param[in] a_side Which of the two faces in that direction.
   * @return Fraction of that face open to flux.
   */
  Real
  areaFraction(const int a_dir, const Side::LoHiSide a_side) const noexcept;

  /**
   * @brief Centroid of one cell face, relative to the face centre, in units of the cell size.
   * @param[in] a_dir  Coordinate direction of the face normal.
   * @param[in] a_side Which of the two faces in that direction.
   * @return Centroid of the open part of that face, relative to the face centre.
   */
  const RealVect&
  faceCentroid(const int a_dir, const Side::LoHiSide a_side) const noexcept;

  /**
   * @brief Magnitude of the interface's area vector.
   * @details This is the quantity that makes sum(alpha_hi - alpha_lo) equal a_B*n exactly, and
   * it is what PolyGeom::bndryArea derives from the apertures. It under-counts a folded patch,
   * for which see trueBoundaryArea.
   * @return Magnitude of the interface's area vector.
   */
  Real
  boundaryArea() const noexcept;

  /**
   * @brief Area of the interface itself, as the sum of its triangles' areas.
   * @details Equal to boundaryArea only for a planar patch. This is the weight under which a
   * boundary centroid partitions, so it is the one a coarsening rule has to use.
   * @return Area of the interface itself.
   */
  Real
  trueBoundaryArea() const noexcept;

  /**
   * @brief Unit normal to the interface, pointing into the fluid.
   * @return Unit normal to the interface, pointing into the fluid.
   */
  const RealVect&
  normal() const noexcept;

  /**
   * @brief Interface centroid, relative to the cell centre, in units of the cell size.
   * @return Centroid of the interface, relative to the cell centre.
   */
  const RealVect&
  boundaryCentroid() const noexcept;

  /**
   * @brief How far the assembled polygons are from enclosing a volume.
   * @details The sum of the outward area vectors over every polygon, which vanishes for a
   * closed surface. Non-zero means the face polygons and the interface patch disagree about
   * where the body's boundary is, and no moment taken from it can be trusted.
   * @return Length of the summed outward area vector, which vanishes for a closed body.
   */
  Real
  closureResidual() const noexcept;

#if CH_SPACEDIM == 3
  /**
   * @brief Append this body's interface, as triangles in physical coordinates.
   * @details The interface is the embedded boundary inside this cell: the polygons the body holds that belong to no
   * cell face. Each is written as a fan of triangles, nine reals per triangle, so that the result can be
   * gathered across ranks and written as one surface. Nothing is appended by a body that is not cut.
   * Coordinates are formed from the cell index rather than the cell centre, so that a vertex two cells share
   * -- a point on their common face, or a crossing on a common edge -- is the same double from both.
   * @param[in,out] a_facets Coordinates to append to.
   * @param[in]     a_cell   Index of this cell.
   * @param[in]     a_probLo Lower-left corner of the domain.
   * @param[in]     a_dx     Grid spacing of the level this cell belongs to.
   */
  void
  appendInterfaceFacets(Vector<Real>&   a_facets,
                        const IntVect&  a_cell,
                        const RealVect& a_probLo,
                        const Real      a_dx) const noexcept;

  /**
   * @brief Print every polygon this body holds, in the cell's own frame, for inspection.
   * @param[in,out] a_out Stream to print to.
   */
  void
  printPolygons(std::ostream& a_out) const noexcept;

  /**
   * @brief Replace one face's chord with the chords of the children that cover it.
   * @details The face is described at the resolution of the level below rather than this cell's, which is what makes it
   * agree with cells on the other side of a refinement boundary. The interface is dropped with it, since it
   * was built to meet the chord being replaced; call closeInterface once every face that needs restricting has
   * been, which is also what recomputes the moments.
   *
   * A quadrant marked closed contributes nothing: the finer cell across it was filled, so that part of the face opens
   * onto nothing, and what the quadrant would have held becomes interface when the interface is closed, as a face
   * closed through snapFace does.
   * @param[in] a_children       Surfaces of the children covering this face, in quadrant order.
   * @param[in] a_closedQuadrant One flag per quadrant, true where the quadrant is closed; null where none is.
   * @param[in] a_dir            Direction of the face.
   * @param[in] a_side           Side of the cell the face is on.
   * @return True if the children's chords merged into a face this body can hold.
   */
  bool
  restrictFace(const CutCellSurface* a_children,
               const bool*           a_closedQuadrant,
               const int             a_dir,
               const int             a_side) noexcept;

  /**
   * @brief Replace one face's polygons with a face entry's, and drop the interface.
   * @details What restrictFace does, with the merged polygons read from the entry instead of walked from the children:
   * the face's polygons and the interface go, and the entry's polygons are appended in the order they are stored.
   * Call closeInterface once every face that needs replacing has been.
   * @param[in] a_overrides The overrides.
   * @param[in] a_entry     Face entry to take the polygons from.
   * @return True if the polygons fitted in the body's storage.
   */
  bool
  replaceFace(const CutCellFaceOverrides& a_overrides, const int a_entry) noexcept;

  /**
   * @brief Rebuild the interface from the face polygons this body holds.
   * @details Every face-polygon edge that no other face polygon carries in reverse bounds the interface, and walking
   * those edges gives its loops. Used instead of the crossing loops once a face has been restricted, since a
   * restricted face has vertices that are not crossings of this cell's edges at all. The moments are recomputed
   * from the closed body.
   * @return True if every loop closed and the result fitted in the body's storage.
   */
  bool
  closeInterface() noexcept;

  /**
   * @brief Close the body over the edges no polygon carries in reverse, with a patch lying in a given face.
   * @details What closeInterface does, with the face the patch belongs to given rather than assumed to be the
   * interface: cutting the body with a plane leaves it open along that plane, and the patch closing it lies in
   * a face of the cell rather than in the interface. A patch in a cell face is oriented by that face, since
   * the edges bounding the opening do not say which way round it should be.
   * @param[in] a_face Face the closing patch lies in, or -1 for the interface.
   * @return True if every loop closed and the result fitted in the body's storage.
   */
  bool
  closeBoundary(const int a_face) noexcept;

  /**
   * @brief Split every polygon edge at the vertices of other polygons that lie on it.
   * @details Restricting a face rebuilds it at the finer level's resolution, so an edge the face shares
   * with a face that was left alone can be one segment on one side and several on the other, or carry a
   * vertex on one side that the other has not got. The two sides then
   * compare unequal and closeInterface reads both as boundary, which puts a stretch of a cell face into
   * the interface and leaves the interface itself unclosed. Inserting each side's vertices into the other's
   * edge makes the two agree segment for segment; the inserted vertices lie on the edge they are inserted
   * into, so no polygon changes shape.
   * @return False if a polygon would need more vertices than it can hold.
   */
  bool
  weldTJunctions() noexcept;

#endif

  /**
   * @brief Make the body the whole cell: every face a polygon, and no interface.
   * @details What a cell whose corners are all fluid holds, written as polygons rather than taken as a degenerate body,
   * so that a face of it can be restricted to the finer level or closed onto a filled cell and the interface closed
   * over what changed.
   */
  void
  defineWhole() noexcept;

  /**
   * @brief Write one face's polygons, as this body holds them, into a set of face overrides.
   * @details Opens a face entry for the cell last begun and adds every polygon lying in that face. A face closed onto
   * a covered cell is recorded with no polygon.
   * @param[in]     a_face      Face index, 2 * dir + side.
   * @param[in]     a_reason    CutCellFaceOverrides::s_finer or CutCellFaceOverrides::s_closed.
   * @param[in,out] a_overrides The overrides, with the cell already begun.
   */
  void
  recordFace(const int a_face, const int a_reason, CutCellFaceOverrides& a_overrides) const;

  /**
   * @brief Make one face agree with a neighbour that holds no solid, or none that holds no fluid, and let what the
   * face loses become interface.
   * @details A neighbour whose kind was decided otherwise than from its own corners leaves this cell describing a
   * face the neighbour does not agree with: this cell's chord on the shared face ends against nothing, and the two
   * sides read different apertures for it. That happens when the volume threshold or the dust rule turns a cell
   * into a regular or a covered one from its own body, and when the graph fills a cell holding more than one
   * sheet, which is covered whatever its nodes read. The face takes the neighbour's answer -- the whole face, or none
   * of it -- and the footprint it gives up becomes part of the interface, which lies in the plane of that face. The
   * identity sum(alpha_hi - alpha_lo) = a_B*n survives it exactly, since the aperture and the boundary's area vector
   * gain the same vector. The moments are recomputed. In two dimensions only closing is done: the face's stretch of the
   * polygon becomes interface lying in the face.
   * @param[in] a_dir               Direction of the face.
   * @param[in] a_side              Side of the cell the face is on.
   * @param[in] a_neighbourIsFluid  True if the neighbour holds no solid, false if it holds no fluid.
   * @return True if the body closed again and its volume stayed in range.
   */
  bool
  snapFace(const int a_dir, const int a_side, const bool a_neighbourIsFluid) noexcept;

#if CH_SPACEDIM == 2
  /**
   * @brief Close one half of a face, and let what that half held become interface.
   * @details Two dimensions only, where a face is a segment and a face on a refinement boundary is covered by two
   * finer faces, one each side of its midpoint. Where the finer cell across one of them was filled, that half of the
   * face opens onto nothing: the polygon's stretch of the face is split at the midpoint and the part in that half
   * becomes interface lying in the face, as snapFace does for a whole face. The moments are recomputed.
   * @param[in] a_dir  Direction of the face.
   * @param[in] a_side Side of the cell the face is on.
   * @param[in] a_half Which half: 0 below the face's midpoint along the other direction, 1 above it.
   * @return True if the polygon had room for the vertex at the midpoint and the body still closes.
   */
  bool
  closeHalfFace(const int a_dir, const int a_side, const int a_half) noexcept;

#endif

  /**
   * @brief How many polygons this body holds.
   * @return Number of polygons.
   */
  int
  numPolygons() const noexcept;

  /**
   * @brief Largest number of vertices any of this body's polygons holds.
   * @return Vertex count of the widest polygon.
   */
  int
  widestPolygon() const noexcept;

  /**
   * @brief Residual in the identity the moments must satisfy among themselves.
   * @details The divergence theorem applied to a constant field, which holds for any closed
   * body whether or not its interface patch is planar. This is also the relation PolyGeom uses
   * to derive the boundary area and the normal from the apertures, so a non-zero residual means
   * the stored moments and the ones Chombo will recompute disagree. Cheap enough to evaluate on
   * every cut cell, and it is sensitive to errors that each moment's own range check passes.
   * @return Length of sum(alpha_hi - alpha_lo) - a_B*n, which vanishes for any closed body.
   */
  Real
  divergenceResidual() const noexcept;

  /**
   * @brief Whether another body holds the same polygons and the same moments, bit for bit.
   * @details Compares every polygon's face, vertices and vertex edges in order, then every moment.
   * @param[in] a_other The other body.
   * @return True if nothing differs.
   */
  bool
  identical(const CutCellBody& a_other) const noexcept;

  /**
   * @brief A 64-bit hash of everything identical compares: every polygon's face, vertices and vertex edges in order,
   * and every moment, bit for bit.
   * @details Two identical bodies hash the same; a sum of hashes over many bodies is what tells two ranks, or two
   * layouts, that they hold the same bodies without sending the bodies.
   * @return The hash.
   */
  std::uint64_t
  fingerprint() const noexcept;

protected:
  /**
   * @brief Number of cell faces.
   */
  static constexpr int s_numFaces = 2 * SpaceDim;

  /**
   * @brief Largest number of polygons a single cut cell's body can need.
   */
  static constexpr int s_maxPolygons = 20;

  /**
   * @brief Largest number of vertices one of those polygons can carry.
   * @details Fourteen was enough for a face carrying a single chord, where the widest polygon
   * measured across four geometries held ten. A face that carries the chords of the cells on the
   * other side of a seam instead of one of its own needs thirteen, measured over 16,640 cut faces,
   * which leaves a single vertex of headroom -- not a margin. Twenty restores one.
   *
   * Raising this costs memory on every body, since a body holds s_maxPolygons of them and the
   * vertex arrays dominate: a polygon goes from 456 to 648 bytes and a body from about 9.4 kB to
   * 13.3 kB. Bodies are transient, built per cut and dropped, so that is paid per thread during a
   * refinement rather than per cell of the index space.
   *
   * Lowering it is what does damage. Sixteen and eight refused between 156 and 5,226 cells per
   * geometry and lost conservation to 1.8e-2, which is why these two numbers are what they are.
   */
  static constexpr int s_maxVertices = 20;

  /**
   * @brief Distance by which a crossing is held off an edge's endpoints.
   * @details A corner exactly on the interface makes the surface tangent to the faces meeting
   * there, so the contour touches a face without cutting any of its edges and marching squares
   * cannot represent it. Displacing the crossing makes the combinatorics generic.
   */
  static constexpr Real s_edgeTolerance = 1.0E-12;

  /**
   * @brief Area below which a polygon carries no moment and is dropped.
   * @details Holding crossings off the edge endpoints leaves slivers of this order behind.
   * They contribute nothing, but a detached one reads as a second connected component, so they
   * are culled rather than carried.
   */
  static constexpr Real s_nullArea = 1.0E-18;
  /**
   * @brief One polygon of the body, in the cell's own frame.
   */
  struct Polygon
  {
    /**
     * @brief Vertices in circuit order, oriented so the area vector points out of the fluid.
     */
    RealVect m_vertex[s_maxVertices];

    /**
     * @brief Cell edge each vertex came from, or -1 for a cell corner.
     */
    int m_vertexEdge[s_maxVertices];

    /**
     * @brief Number of vertices in use.
     */
    int m_numVertices;

    /**
     * @brief Cell face this polygon lies in, or -1 if it is part of the interface.
     */
    int m_face;

    /**
     * @brief Cell face the segment leaving each vertex lies in, or -1 for the interface.
     * @details Two dimensions only, where the fluid region is one polygon whose segments lie in
     * different cell faces. Which face a segment belongs to is read off the circuit that
     * produced it rather than from where its endpoints sit: a sliver no wider than the distance
     * crossings are held off the edge endpoints has every vertex within any position tolerance
     * of the same face plane, and a positional test then assigns the whole polygon to it.
     */
    int m_segmentFace[s_maxVertices];
  };

  /**
   * @brief Moments of a cell the interface only touches, never enters.
   * @details The boundary of such a cell is one of its own faces, so the area, the normal and
   * the centroid all follow from which face is covered rather than from any polygon.
   * @param[in] a_surface Corner values in the cell's own frame.
   * @param[in] a_kind    Regular or Covered, as returned by classify.
   */
  void
  defineDegenerate(const CutCellSurface& a_surface, const Kind a_kind) noexcept;

  /**
   * @brief Assemble the polygons of a cell the interface passes through.
   * @param[in] a_surface Corner values and edge crossings in the cell's own frame.
   * @return True if the polygons were assembled; false if the face contours could not be paired
   * or the crossings did not form clean loops.
   */
  bool
  defineCut(const CutCellSurface& a_surface) noexcept;

#if CH_SPACEDIM == 3
  /**
   * @brief Merge coplanar polygons into the boundary of their union.
   * @details An edge shared by two of them is interior, so cancelling every edge that appears reversed among the others
   * leaves the union's boundary, and walking what remains gives it in order. Collinear vertices are then
   * dropped, which is what lets a sub-face's edge along a cell edge compare equal to the coarse face's edge
   * beside it. The union need not be connected, and a union of no area is an answer rather than a failure.
   * @param[in]  a_in      Polygons to merge, all in one plane.
   * @param[in]  a_num     How many.
   * @param[out] a_out     The union's boundary, one polygon per connected piece.
   * @param[in]  a_maxOut  How many pieces the caller can hold.
   * @param[out] a_numOut  How many pieces were written.
   * @return True if the union was walked and fitted.
   */
  bool
  mergeCoplanar(const Polygon* a_in, const int a_num, Polygon* a_out, const int a_maxOut, int& a_numOut) const noexcept;

  /**
   * @brief The fluid part of one cell face, as outward-oriented polygons.
   * @details Walks the face's circuit of corners and edges, keeping fluid corners and inserting
   * the crossing wherever an edge changes sign. On a face carrying four crossings the polygons
   * are read off the chord pairing rather than off the corner signs, so that the face cannot
   * disagree with the loop assembly, which trusts that pairing.
   * @param[in]  a_dir     Coordinate direction of the face normal.
   * @param[in]  a_side    0 for the low face in that direction, 1 for the high face.
   * @param[in]  a_surface Corner values and edge crossings in the cell's own frame.
   * @param[out] a_out     Room for up to two polygons.
   * @return Number of polygons written, or -1 if the face could not be walked.
   */
  int
  faceWalk(const int a_dir, const int a_side, const CutCellSurface& a_surface, Polygon* a_out) const noexcept;

  /**
   * @brief Reverse a face polygon if its area vector points into the fluid.
   * @param[in,out] a_polygon Polygon lying in the given cell face.
   * @param[in]     a_dir     Coordinate direction of the face normal.
   * @param[in]     a_side    0 for the low face in that direction, 1 for the high face.
   */
  void
  orientOutward(Polygon& a_polygon, const int a_dir, const int a_side) const noexcept;
#endif

  /**
   * @brief Whether a set of children partitions this body, in every moment that has to add up.
   * @details What subdivide has to deliver, checked rather than asserted so that a caller learns of it. First,
   * that every child holds its fluid in one piece: a child that does not is reported by returning false before
   * anything is summed, because it is a fact about a parent whose interface is twisted rather than a fault to
   * stop on, and hasMultiValuedChildren asks exactly that. Then four things: the children's volumes add to this
   * one's, the apertures of the children on a face of this cell add to that face's, the two children meeting at an
   * interior face agree about it -- which is what lets those sums cancel, and what a cell cut from this one relies on
   * to join its sibling -- and the children's interface area vectors add to this one's. Each child is also required
   * to satisfy the divergence identity on its own, since a child that does not is not a body whatever it sums to.
   *
   * A quantity spanning fewer dimensions than the cell is counted once too often per dimension it does not
   * span, so volumes carry 2^-SpaceDim and areas 2^-(SpaceDim-1).
   * @param[in] a_children The 2^SpaceDim bodies, in quadrant order.
   * @return True if every sum returns this body's own moment to rounding.
   */
  bool
  partitions(const CutCellBody* a_children) const noexcept;

  /**
   * @brief Take the moments of the assembled polygons.
   * @details Face apertures are accumulated signed about the outward normal, because where the
   * interface lies in a cell face it bounds a hole in that face and the aperture is the net
   * area open to flux.
   */
  void
  accumulateMoments() noexcept;

  /**
   * @brief The polygons bounding the fluid.
   */
  Polygon m_polygon[s_maxPolygons];

  /**
   * @brief Number of polygons in use.
   */
  int m_numPolygons;

  /**
   * @brief Volume fraction.
   */
  Real m_volumeFraction;

  /**
   * @brief Area fraction of each cell face, indexed 2*dir + side.
   */
  Real m_areaFraction[s_numFaces];

  /**
   * @brief Magnitude of the interface's area vector.
   */
  Real m_boundaryArea;

  /**
   * @brief Area of the interface itself.
   */
  Real m_trueBoundaryArea;

  /**
   * @brief Volume centroid.
   */
  RealVect m_volumeCentroid;

  /**
   * @brief Centroid of each cell face, indexed 2*dir + side.
   */
  RealVect m_faceCentroid[s_numFaces];

  /**
   * @brief Unit normal to the interface, into the fluid.
   */
  RealVect m_normal;

  /**
   * @brief Interface centroid.
   */
  RealVect m_boundaryCentroid;

  /**
   * @brief Sum of the outward area vectors over every polygon.
   */
  RealVect m_closure;
};

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#endif
