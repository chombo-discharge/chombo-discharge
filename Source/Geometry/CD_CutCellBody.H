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

// Chombo includes
#include <RealVect.H>
#include <LoHiSide.H>

// Our includes
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
   * @brief How the cell relates to the interface, from the corner values alone.
   * @details Follows GeometryShop::insideOutsideFromNodes: the cell is cut as soon as two
   * corners disagree in sign, with exact zeros compared through copysign so that a corner at
   * negative zero counts as fluid. A face lying exactly in the interface therefore leaves the
   * cell regular with that face covered, rather than cut.
   * @param[in] a_surface Corner values in the cell's own frame.
   * @return Regular, Covered or Cut.
   */
  static Kind
  classify(const CutCellSurface& a_surface) noexcept;

  /**
   * @brief Whether one side of the interface is present only as exact zeros.
   * @details A facet lying in a node plane leaves every corner on one side at exactly zero. The
   * signs then disagree and the corner test reads the cell as cut, but a side represented only
   * by exact zeros encloses nothing.
   * @param[in] a_surface   Corner values in the cell's own frame.
   * @param[in] a_fluidSide Which side to test.
   * @return True if that side is present and every one of its corners is exactly zero.
   */
  static bool
  touchesOnly(const CutCellSurface& a_surface, const bool a_fluidSide) noexcept;

  /**
   * @brief Whether the interface lies exactly in one of the cell's faces.
   * @details A cell the interface only grazes is full, and classify calls it regular for that
   * reason, but the face the interface lies in is closed: the fluid ends there. Such a cell is
   * not a regular cell, whose every face is open by definition, and it has to carry a body
   * saying so. Leaving it regular has the cut cell across that face read the face from the
   * interface, find nothing of it open, and drop its side of a face the full cell keeps.
   *
   * Only meaningful for a cell classify has called regular: the test is that every corner of
   * some face is on the solid side, which in such a cell means those corners are on the
   * interface.
   * @param[in] a_surface Corner values in the cell's own frame.
   * @return True if some face has all of its corners on the interface.
   */
  static bool
  interfaceLiesInFace(const CutCellSurface& a_surface) noexcept;

  /**
   * @brief How this body relates to the interface.
   * @details Settled when the body was built and carried from there, rather than worked out
   * again from corner values at whatever resolution the cell happens to sit on. A body built
   * from the cells that partition it takes theirs: cut if any of them is cut, and otherwise
   * whichever side they all lie on. That is what stops a cell being called whole because its
   * own corners agree, when the surface it was built from passes through it.
   * @return Regular, Covered or Cut.
   */
  Kind
  kind() const noexcept;

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
   * @brief How a cell relates to the interface, given how the cells partitioning it do.
   * @details Cut if any of them is cut, and cut also if some are whole and others covered, since
   * the interface has to pass between them. Otherwise whichever side they all lie on.
   * @param[in] a_kinds Kind of each of the cells partitioning the cell, in any order.
   * @return Regular, Covered or Cut.
   */
  static Kind
  coarsenKind(const Kind a_kinds[1 << SpaceDim]) noexcept;

  /**
   * @brief Whether one plane can stand for this body's interface.
   * @details False where the interface turns so far within the cell that its patches face away
   * from each other, so that the area carrying flux is not the area that is there. See
   * s_maxAreaSpread.
   *
   * Asked of a body rather than enforced by it: a coarse level cannot resolve everything and is
   * not meant to, so it is only on the level where the cells are claimed to be single valued
   * that a body failing this is a body that should not have been built.
   * @return True if the interface can be represented.
   */
  bool
  interfaceIsOneSided() const noexcept;

  /**
   * @brief Subdivide the body into the 2^SpaceDim cells that partition it.
   * @details Clips the polygons by the cell's own mid-planes and caps each opening that leaves,
   * so the children are literal sub-bodies of this one rather than reconstructions of it. Every
   * coarsening relation is then an identity about a partition rather than a condition to
   * satisfy, which is what makes the children's moments sum back to this body's exactly.
   * @param[out] a_children The children, each expressed in its own frame, ordered so that child
   * c sits on the high side of direction d when bit d of c is set.
   * @return True if every child closed; false if one of them could not be assembled.
   */
  bool
  refine(CutCellBody a_children[1 << SpaceDim]) const noexcept;

  /**
   * @brief Build this body's moments from the bodies of the cells that partition it.
   * @details The inverse of refine, and the operation that lets a cell whose geometry was
   * reconstructed at a finer resolution be carried up to a coarser one without being
   * reconstructed there a second time. Only the moments are taken: the children partition this
   * cell, so its volume, apertures and interface area are theirs summed and rescaled, and no
   * polyhedron has to be reassembled to read them off. Where the interface is a plane the result
   * is the same as reconstructing here; where it curves it is better, because it is the finer
   * description that was integrated.
   * @param[in] a_children The children, each in its own frame, ordered as refine orders them.
   * @return True if the children leave a body with a boundary that closes.
   */
  bool
  coarsen(const CutCellBody a_children[1 << SpaceDim]) noexcept;

  /**
   * @brief How far the assembled polygons are from enclosing a volume.
   * @details The sum of the outward area vectors over every polygon, which vanishes for a
   * closed surface. Non-zero means the face polygons and the interface patch disagree about
   * where the body's boundary is, and no moment taken from it can be trusted.
   * @return Length of the summed outward area vector, which vanishes for a closed body.
   */
  Real
  closureResidual() const noexcept;

  /**
   * @brief Append this body's interface, as triangles in physical coordinates.
   * @details The interface is the embedded boundary inside this cell: the polygons the body holds
   *          that belong to no cell face. Each is written as a fan of triangles, nine reals per
   *          triangle, so that the result can be gathered across ranks and written as one surface.
   *          Nothing is appended by a body that is not cut.
   * @param[in,out] a_facets     Coordinates to append to.
   * @param[in]     a_cellCentre Centre of this cell, in physical coordinates.
   * @param[in]     a_dx         Grid spacing of the level this cell belongs to.
   */
  void
  appendInterfaceFacets(Vector<Real>& a_facets, const RealVect& a_cellCentre, const Real a_dx) const noexcept;

  /**
   * @brief Replace one face's chord with the chords of the children that cover it.
   * @details The face is described at the resolution of the level below rather than this cell's,
   *          which is what makes it agree with cells on the other side of a refinement boundary.
   *          The interface is dropped with it, since it was built to meet the chord being replaced;
   *          call closeInterface once every face that needs restricting has been.
   * @param[in] a_children Surfaces of the children covering this face, in quadrant order.
   * @param[in] a_dir      Direction of the face.
   * @param[in] a_side     Side of the cell the face is on.
   * @return True if the children's chords merged into a face this body can hold.
   */
  bool
  restrictFace(const CutCellSurface* a_children, const int a_dir, const int a_side) noexcept;

  /**
   * @brief Rebuild the interface from the face polygons this body holds.
   * @details Every face-polygon edge that no other face polygon carries in reverse bounds the
   *          interface, and walking those edges gives its loops. Used instead of the crossing loops
   *          once a face has been restricted, since a restricted face has vertices that are not
   *          crossings of this cell's edges at all.
   * @return True if every loop closed and the result fitted in the body's storage.
   */
  bool
  closeInterface() noexcept;

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
   */
  static constexpr int s_maxVertices = 14;

  /**
   * @brief Distance by which a crossing is held off an edge's endpoints.
   * @details A corner exactly on the interface makes the surface tangent to the faces meeting
   * there, so the contour touches a face without cutting any of its edges and marching squares
   * cannot represent it. Displacing the crossing makes the combinatorics generic.
   */
  static constexpr Real s_edgeTolerance = 1.0E-12;

  /**
   * @brief How far the interface's area may exceed the magnitude of its area vector.
   * @details The ratio between the two is one over the cosine of the spread of the patch normals
   * about their mean, so it says how far the interface turns within the cell. Two patches meeting
   * at a right angle give the square root of two, three meeting at a corner give the square root
   * of three, and no convex feature gives more. Past two the interface has turned beyond a right
   * angle and is folding back on itself, which is a cell holding sheets that face away from each
   * other.
   *
   * Such a cell has no single plane, and EBData stores one plane. The area vectors cancel, so the
   * area that carries flux goes to nothing while the interface is really there, and the
   * divergence identity does not notice: an interface lying wholly inside a cell touches none of
   * its faces, so the identity reads zero against zero. Refusing is the only thing that catches
   * it. The cure is resolution, which is what the cell is refused in favour of.
   */
  static constexpr Real s_maxAreaSpread = 2.0;

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
   * @details An edge shared by two of them is interior, so cancelling every edge that appears
   *          reversed among the others leaves the union's boundary, and walking what remains gives
   *          it in order. Collinear vertices are then dropped, which is what lets a sub-face's edge
   *          along a cell edge compare equal to the coarse face's edge beside it. The union need
   *          not be connected, and a union of no area is an answer rather than a failure.
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
   * @brief Cut the body with an axis-aligned half-space, capping the opening that leaves.
   * @details The cap is assembled from the segments the cut leaves lying in the plane, and they
   * are consumed as segments rather than as vertices: a cross-section can have several
   * components and two of them may meet at a shared vertex, and marking that vertex used drops
   * one loop and splices the other.
   * @param[in]  a_dir        Coordinate direction of the cutting plane's normal.
   * @param[in]  a_coordinate Where the plane sits along that direction.
   * @param[in]  a_keepLow    Whether to keep the side below the plane.
   * @param[out] a_out        The cut body.
   * @return True if the result fitted in the body's fixed storage.
   */
  bool
  clip(const int a_dir, const Real a_coordinate, const bool a_keepLow, CutCellBody& a_out) const noexcept;

  /**
   * @brief Read the classification off the polygons this body is made of.
   * @details A body with no polygons encloses nothing and the cell is covered. One whose polygons
   * all lie in cell faces never meets the interface, so the cell is full. Anything else is cut.
   * @return Regular, Covered or Cut.
   */
  Kind
  kindFromPolygons() const noexcept;

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
   * @brief How this body relates to the interface, as settled when it was built.
   */
  Kind m_kind;

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
