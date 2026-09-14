/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PolyhedralEBUtils.H
 * @brief  Declaration of the cell topology and polygon geometry the polyhedral cut cells are built from.
 * @author Robert Marskar
 */

#ifndef CD_POLYHEDRALEBUTILS_H
#define CD_POLYHEDRALEBUTILS_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_CutCellSurface.H>
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

/**
 * @brief Which side of the interface a value lies on.
 * @details One predicate for corners and crossings alike, matching the sign test in
 * GeometryShop::insideOutsideFromNodes. The side is read off the sign bit, so a node at
 * negative zero is fluid and one at positive zero is solid.
 * @param[in] a_value Implicit function value.
 * @return True if the value is on the fluid side.
 */
inline bool
isFluid(const Real a_value) noexcept;

/**
 * @brief Unit normal of the plane through a cut cell's edge crossings.
 * @details Total least squares through the crossing points: the normal is the direction the
 * points spread least in, which is the null space of their covariance, read off its adjugate
 * rather than from an eigenvalue solve. Three crossings span the plane exactly and more are
 * fitted. Only the positions of the sign changes enter, never how fast the implicit function
 * moves through them, so this survives implicit functions whose derivatives do not exist at
 * their own zero set. The orientation is taken from the corner values, pointing into the fluid.
 * @param[in] a_surface Corner values and edge crossings of the cell.
 * @return Unit normal in the cell's own frame, or zero where the crossings do not span a plane.
 */
inline RealVect
crossingNormal(const CutCellSurface& a_surface) noexcept;

/**
 * @brief Unit normal from the gradient of the multilinear interpolant of the corner values.
 * @details Cheaper than crossingNormal, needing no crossings at all, and likewise insensitive to
 * how the implicit function is scaled: a monotone reparameterisation rescales the gradient but
 * does not turn it. It is evaluated at the cell centre rather than on the surface, so it is the
 * less faithful of the two where the interface runs near a corner.
 * @param[in] a_surface Corner values of the cell.
 * @return Unit normal in the cell's own frame, pointing into the fluid, or zero where the
 * corner values are flat.
 */
inline RealVect
gradientNormal(const CutCellSurface& a_surface) noexcept;

/**
 * @brief Number of samples in the stencil the curvature is measured from.
 */
constexpr int s_curvatureStencilSize = (SpaceDim == 2) ? 9 : 27;

/**
 * @brief Largest principal curvature of the level set through a stencil of implicit function values.
 * @details The samples sit on a 3^SpaceDim lattice of spacing a_spacing, centred on the point the
 * curvature is wanted at, and are indexed by sum_d (o_d + 1) * 3^d over offsets o_d in {-1, 0, 1}.
 * Gradient and Hessian are central differences of those samples, and the principal curvatures
 * follow from them in closed form, so no eigenvalue solve is needed. Dividing through by the
 * gradient magnitude makes the result independent of how the implicit function is scaled, which
 * matters because only some of them are signed distance functions.
 *
 * The spacing is what the measurement is relative to. At a spacing of one cell width the result
 * is the bend of the surface across a cell, which stays finite where two surfaces meet in a
 * crease; measured at a spacing approaching zero the same crease diverges instead. The former is
 * what a refinement criterion wants, since a crease is not something more resolution can fix.
 *
 * Only valid where the implicit function is twice differentiable at its own zero set. That is a
 * real restriction rather than a formality: PerlinSphereSdf returns the square root of the
 * difference of squared radii, whose second derivative is unbounded at the surface, and the
 * curvature of a sphere measured through it diverges under refinement instead of converging to
 * the inverse radius. Where the implicit function cannot be trusted that far, take the normals
 * from crossingNormal and compare them between neighbouring cells instead; root finding only
 * asks where the function changes sign, never how fast.
 * @param[in] a_stencil Implicit function sampled on the lattice.
 * @param[in] a_spacing Lattice spacing.
 * @return Largest principal curvature in magnitude, or zero where the gradient vanishes.
 */
inline Real
maxPrincipalCurvature(const Real a_stencil[s_curvatureStencilSize], const Real a_spacing) noexcept;

/**
 * @brief Cell topology and polygon geometry the cut cell classes are assembled from.
 * @details Corners, edges and faces are indexed in the cell's own [-1/2, 1/2]^SpaceDim frame.
 * Corner c has offset ((c >> d) & 1) in direction d. Edges are numbered direction-major, with
 * s_numEdges / SpaceDim edges running along each direction: edge e runs along direction
 * e / (s_numEdges / SpaceDim), and the remainder holds the offsets of the edge's low corner in
 * the transverse directions, one bit per direction in increasing order. An edge is
 * parameterised from its low corner. Faces are indexed as 2 * dir + side with side 0 the low
 * face and side 1 the high face. Two neighbouring cells therefore address a shared corner,
 * edge or face identically.
 */
namespace detail {

#if CH_SPACEDIM == 3
/**
 * @brief The two directions transverse to each coordinate direction, in increasing order.
 * @details This is the order the transverse offsets are packed into an edge index and the order
 * the face corners and edges are walked in.
 */
constexpr int s_transverse[3][2] = {{1, 2}, {0, 2}, {0, 1}};
#endif

/**
 * @brief Direction a cell edge runs along.
 * @param[in] a_edge Edge index in [0, CutCellSurface::s_numEdges).
 * @return Coordinate direction of the edge.
 */
inline int
edgeDirection(const int a_edge) noexcept;

/**
 * @brief Corner offsets of a cell edge's low end.
 * @param[in]  a_edge   Edge index in [0, CutCellSurface::s_numEdges).
 * @param[out] a_offset Corner offsets, zero in the direction the edge runs along.
 */
inline void
edgeOrigin(const int a_edge, int a_offset[SpaceDim]) noexcept;

/**
 * @brief The corners a cell edge joins, low end first.
 * @param[in]  a_edge Edge index in [0, CutCellSurface::s_numEdges).
 * @param[out] a_lo   Corner at the edge's low end.
 * @param[out] a_hi   Corner at the edge's high end.
 */
inline void
edgeCorners(const int a_edge, int& a_lo, int& a_hi) noexcept;

/**
 * @brief Position of a cell corner in the cell's own frame.
 * @param[in] a_corner Corner index in [0, CutCellSurface::s_numCorners).
 * @return Corner position, each component -1/2 or 1/2.
 */
inline RealVect
cornerPosition(const int a_corner) noexcept;

/**
 * @brief The corners of a cell face, in circuit order around the face.
 * @details The circuit runs through the transverse offsets (0,0), (1,0), (1,1), (0,1). Its
 * handedness relative to the face's outward normal depends on the direction, so a polygon
 * assembled along it has to be oriented afterwards.
 * @param[in]  a_dir    Coordinate direction of the face normal.
 * @param[in]  a_side   0 for the low face in that direction, 1 for the high face.
 * @param[out] a_corner Face corners in circuit order.
 */
inline void
faceCorners(const int a_dir, const int a_side, int a_corner[1 << (SpaceDim - 1)]) noexcept;

/**
 * @brief Index of the edge running along a_dir whose low corner has the given offsets.
 * @details Only the offsets in the transverse directions take part; the offset along a_dir is
 * ignored, since the edge's low corner has offset zero there by definition.
 * @param[in] a_dir    Coordinate direction the edge runs along.
 * @param[in] a_offset Corner offsets of the edge's low end, each 0 or 1.
 * @return Edge index in [0, CutCellSurface::s_numEdges).
 */
inline int
edgeIndex(const int a_dir, const int a_offset[SpaceDim]) noexcept;

/**
 * @brief Position of an edge crossing in the cell's own frame.
 * @details The crossing is held at least a_tolerance away from either endpoint of the edge,
 * except a crossing recorded exactly at an endpoint: that one sits on a corner the interface
 * passes through and is left where it is, since displacing it would open a sliver of the
 * displacement's own width.
 * @param[in] a_surface   Corner values and edge crossings of the cell.
 * @param[in] a_edge      Edge index in [0, CutCellSurface::s_numEdges). The edge must carry a
 * crossing.
 * @param[in] a_tolerance Distance by which the crossing is held off the edge's endpoints, in
 * [0, 1/2).
 * @return Crossing position.
 */
inline RealVect
crossingPosition(const CutCellSurface& a_surface, const int a_edge, const Real a_tolerance) noexcept;

#if CH_SPACEDIM == 3
/**
 * @brief The edges of a cell face, edge i joining face corners i and i+1.
 * @param[in]  a_dir  Coordinate direction of the face normal.
 * @param[in]  a_side 0 for the low face in that direction, 1 for the high face.
 * @param[out] a_edge Face edges in the circuit order of faceCorners.
 */
inline void
faceEdges(const int a_dir, const int a_side, int a_edge[4]) noexcept;

/**
 * @brief Area vector, area and centroid of a planar polygon given in circuit order.
 * @details The area vector is the sum of the signed triangle area vectors of a fan from the
 * first vertex, so a non-convex polygon comes out right. A face contour that is a polyline
 * rather than a single chord leaves such a polygon, and summing unsigned triangle areas over
 * the fan over-counts it. The centroid is the area-weighted mean of the triangle centroids
 * with the same signs. Fewer than three vertices yield all-zero outputs.
 * @param[in]  a_vertex   Polygon vertices in circuit order.
 * @param[in]  a_num      Number of vertices.
 * @param[out] a_area     Magnitude of the area vector, which is the polygon area when the
 * polygon is planar and simple.
 * @param[out] a_vector   Polygon area vector, normal to the polygon with the handedness of the
 * circuit.
 * @param[out] a_centroid Polygon centroid, zero when the area is zero.
 */
inline void
polygonMoments(const RealVect* a_vertex,
               const int       a_num,
               Real&           a_area,
               RealVect&       a_vector,
               RealVect&       a_centroid) noexcept;

/**
 * @brief The chords on a cell face, each an unordered pair of edge indices.
 * @details A face with two crossings carries one chord. A face whose corners alternate fluid
 * and solid carries four crossings and there are two ways to join them; the bilinear
 * interpolant of the four corner values decides it. The saddle shares its side with exactly
 * one of the two diagonals: that diagonal meets through the middle, and the other is the pair
 * the chords cut off. Both cells adjoining the face hold the same four values, so they cannot
 * disagree. With four crossings each chord joins two consecutive edges of the face circuit.
 * @param[in]  a_dir     Coordinate direction of the face normal.
 * @param[in]  a_side    0 for the low face in that direction, 1 for the high face.
 * @param[in]  a_surface Corner values and edge crossings of the cell.
 * @param[out] a_pair    Each chord's two edge indices. Entries beyond the returned count are
 * left untouched.
 * @return Number of chords, 0, 1 or 2, or -1 if the face carries an odd number of crossings.
 */
inline int
facePairs(const int a_dir, const int a_side, const CutCellSurface& a_surface, int a_pair[2][2]) noexcept;

/**
 * @brief Order the crossings into closed loops, or fail if they do not form clean cycles.
 * @details Every crossing lies on exactly two faces and every face pairs its crossings up, so
 * each node has degree two and the graph is a disjoint union of cycles. Several cycles is
 * legal: the interface simply enters the cell as several sheets. Every loop returned has at
 * least three crossings.
 * @param[in]  a_surface Corner values and edge crossings of the cell.
 * @param[out] a_loop    Edge indices of every loop, one loop after the other.
 * @param[out] a_start   Where each loop begins in a_loop, with the total length last.
 * @return Number of loops, or -1 if some face carries an odd number of crossings, some
 * crossing is paired with other than two, the cell has fewer than three crossings in all, or a
 * loop closes on fewer than three.
 */
inline int
crossingLoops(const CutCellSurface& a_surface,
              int                   a_loop[CutCellSurface::s_numEdges],
              int                   a_start[CutCellSurface::s_numEdges + 1]) noexcept;
#endif
/**
 * @brief Distance within which a clipped vertex is placed exactly on the cutting plane.
 * @details A vertex kept within tolerance of the plane must be put exactly on it, or the cap
 * it belongs to stops being recognised as lying in a cell face and its area is charged to the
 * interface instead.
 */
constexpr Real s_clipTolerance = 1.0E-12;

/**
 * @brief Distance within which two cap vertices are taken to be the same one.
 * @details A vertex where two of a body's polygons meet is computed once along each polygon's
 * own edge, and those two intersections are ill conditioned where the edges are nearly parallel:
 * the answers agree to about a part in a thousand billion rather than bit for bit. Cap loops are
 * stitched by matching one segment's end to the next one's start, so the match has to be made at
 * that distance and not exactly.
 *
 * Set far below any cap edge that carries area -- a cell is one unit across in this frame -- and
 * far above the disagreement two intersections of the same point produce, which is what leaves a
 * band to choose from at all. Quantising instead, which is what this replaced, has no such band:
 * two answers either side of a rounding boundary go to different buckets however close they are,
 * the loop then fails to close, and the body silently loses a whole face.
 */
constexpr Real s_weldTolerance = 1.0E-10;

/**
 * @brief Whether two positions are the same vertex.
 * @param[in] a_a First position.
 * @param[in] a_b Second position.
 * @return True if they lie within s_weldTolerance of one another.
 */
inline bool
sameVertex(const RealVect& a_a, const RealVect& a_b) noexcept;

/**
 * @brief Order two points so an edge is always interpolated from the same end.
 * @details Two polygons sharing an edge walk it in opposite directions, and a + t(b - a) is
 * only algebraically equal to b + t'(a - b): the two differ in the last bits, which is enough
 * to put them in different quantisation buckets and break the cap loop that has to join them.
 * @param[in] a_a First point.
 * @param[in] a_b Second point.
 * @return True if the first point comes first.
 */
inline bool
lexLess(const RealVect& a_a, const RealVect& a_b) noexcept;

/**
 * @brief Area of a simple polygon, in whichever dimension it lives.
 * @details In three dimensions this is the length of the area vector; in two it is the shoelace
 * sum, where there is no area vector to take a length of. Used where only the size matters, so
 * that a clipped polygon can be tested for being worth keeping in either dimension.
 * @param[in] a_vertex Polygon vertices in circuit order.
 * @param[in] a_num    Number of vertices.
 * @return Polygon area.
 */
inline Real
polygonArea(const RealVect* a_vertex, const int a_num) noexcept;

} // namespace detail
} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#include <CD_PolyhedralEBUtilsImplem.H>

#endif
