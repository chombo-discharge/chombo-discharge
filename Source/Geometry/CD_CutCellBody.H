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
#include <CD_NamespaceHeader.H>

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
   * @brief Which side of the interface a value lies on.
   * @details One predicate for corners and crossings alike, matching the sign test in
   * GeometryShop::insideOutsideFromNodes. The side is read off the sign bit, so a node at
   * negative zero is fluid and one at positive zero is solid.
   * @param[in] a_value Implicit function value.
   * @return True if the value is on the fluid side.
   */
  static bool
  isFluid(const Real a_value) noexcept;

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
   * @brief Distance within which a vertex counts as lying in a cell face.
   */
  static constexpr Real s_planeTolerance = 1.0E-12;

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

  /**
   * @brief The fluid part of one cell face, as outward-oriented polygons.
   * @details Walks the face's circuit of corners and edges, keeping fluid corners and inserting
   * the crossing wherever an edge changes sign. On a face carrying four crossings the polygons
   * are read off the chord pairing rather than off the corner signs, so that the face cannot
   * disagree with the loop assembly, which trusts that pairing.
   * @param[in]  a_dir     Coordinate direction of the face normal.
   * @param[in]  a_side    Which of the two faces in that direction.
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
   * @param[in]     a_side    Which of the two faces in that direction.
   */
  void
  orientOutward(Polygon& a_polygon, const int a_dir, const int a_side) const noexcept;

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
   * @brief Volume centroid.
   */
  RealVect m_volumeCentroid;

  /**
   * @brief Area fraction of each cell face, indexed 2*dir + side.
   */
  Real m_areaFraction[s_numFaces];

  /**
   * @brief Centroid of each cell face, indexed 2*dir + side.
   */
  RealVect m_faceCentroid[s_numFaces];

  /**
   * @brief Magnitude of the interface's area vector.
   */
  Real m_boundaryArea;

  /**
   * @brief Area of the interface itself.
   */
  Real m_trueBoundaryArea;

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

#include <CD_NamespaceFooter.H>

#endif
