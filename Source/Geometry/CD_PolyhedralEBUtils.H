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
 * @details Corners, edges and faces are indexed in the cell's own [-1/2, 1/2]^SpaceDim frame:
 * corner c has offset ((c >> d) & 1) in direction d, and an edge is parameterised from its low
 * corner. Two neighbouring cells therefore address a shared corner or edge identically.
 */
namespace detail {

#if CH_SPACEDIM == 3
/**
 * @brief The two directions transverse to each coordinate direction.
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
 * @return Corner position.
 */
inline RealVect
cornerPosition(const int a_corner) noexcept;

/**
 * @brief The corners of a cell face, in circuit order around the face.
 * @param[in]  a_dir    Coordinate direction of the face normal.
 * @param[in]  a_side   Which of the two faces in that direction.
 * @param[out] a_corner Face corners in circuit order.
 */
inline void
faceCorners(const int a_dir, const int a_side, int a_corner[1 << (SpaceDim - 1)]) noexcept;

/**
 * @brief Index of the edge running along a_dir whose low corner has the given offsets.
 * @param[in] a_dir    Coordinate direction the edge runs along.
 * @param[in] a_offset Corner offsets of the edge's low end.
 * @return Edge index in [0, CutCellSurface::s_numEdges).
 */
inline int
edgeIndex(const int a_dir, const int a_offset[SpaceDim]) noexcept;

/**
 * @brief The edges of a cell face, edge i joining face corners i and i+1.
 * @param[in]  a_dir  Coordinate direction of the face normal.
 * @param[in]  a_side Which of the two faces in that direction.
 * @param[out] a_edge Face edges in circuit order.
 */
inline void
faceEdges(const int a_dir, const int a_side, int a_edge[2 * (SpaceDim - 1)]) noexcept;

/**
 * @brief Position of an edge crossing in the cell's own frame.
 * @param[in] a_surface   Corner values and edge crossings of the cell.
 * @param[in] a_edge      Edge index in [0, CutCellSurface::s_numEdges).
 * @param[in] a_tolerance Distance by which the crossing is held off the edge's endpoints.
 * @return Crossing position.
 */
inline RealVect
crossingPosition(const CutCellSurface& a_surface, const int a_edge, const Real a_tolerance) noexcept;

/**
 * @brief Area, area vector and centroid of a simple planar polygon in circuit order.
 * @details Areas are taken signed about the polygon's own normal. A face contour that is a
 * polyline rather than a single chord leaves a non-convex polygon, and summing unsigned
 * triangle areas over a fan over-counts those.
 * @param[in]  a_vertex   Polygon vertices in circuit order.
 * @param[in]  a_num      Number of vertices.
 * @param[out] a_area     Polygon area.
 * @param[out] a_vector   Polygon area vector.
 * @param[out] a_centroid Polygon centroid.
 */
inline void
polygonMoments(const RealVect* a_vertex,
               const int       a_num,
               Real&           a_area,
               RealVect&       a_vector,
               RealVect&       a_centroid) noexcept;

#if CH_SPACEDIM == 3
/**
 * @brief The chords on a cell face, each an ordered pair of edge indices.
 * @details A face with two crossings carries one chord. A face whose corners alternate fluid
 * and solid carries four crossings and there are two ways to join them; the bilinear
 * interpolant of the four corner values decides it. The saddle shares its side with exactly
 * one of the two diagonals: that diagonal meets through the middle, and the other is the pair
 * the chords cut off. Both cells adjoining the face hold the same four values, so they cannot
 * disagree.
 * @param[in]  a_dir     Coordinate direction of the face normal.
 * @param[in]  a_side    Which of the two faces in that direction.
 * @param[in]  a_surface Corner values and edge crossings of the cell.
 * @param[out] a_pair    Each chord's two edge indices.
 * @return Number of chords, or -1 if the face carries an impossible number of crossings.
 */
inline int
facePairs(const int a_dir, const int a_side, const CutCellSurface& a_surface, int a_pair[2][2]) noexcept;

/**
 * @brief Order the crossings into closed loops, or fail if they do not form clean cycles.
 * @details Every crossing lies on exactly two faces and every face pairs its crossings up, so
 * each node has degree two and the graph is a disjoint union of cycles. Several cycles is
 * legal: the interface simply enters the cell as several sheets.
 * @param[in]  a_surface Corner values and edge crossings of the cell.
 * @param[out] a_loop    Edge indices of every loop, one loop after the other.
 * @param[out] a_start   Where each loop begins in a_loop, with the total length last.
 * @return Number of loops, or -1 if the crossings do not form clean cycles.
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
 * @brief A vertex's identity, quantised, for matching cap segments up.
 * @details Compared component by component rather than hashed: a hash collision would splice
 * two cap loops into one.
 */
struct VertexKey
{
  /**
   * @brief Quantised position.
   */
  long m_q[SpaceDim] = {0};

  /**
   * @brief Whether two vertices quantise to the same position.
   * @param[in] a_other Key to compare against.
   * @return True if every component matches.
   */
  bool
  operator==(const VertexKey& a_other) const noexcept
  {
    bool same = true;

    for (int d = 0; d < SpaceDim; d++) {
      same = same && (m_q[d] == a_other.m_q[d]);
    }

    return same;
  }

  /**
   * @brief Whether two vertices quantise to different positions.
   * @param[in] a_other Key to compare against.
   * @return True if any component differs.
   */
  bool
  operator!=(const VertexKey& a_other) const noexcept
  {
    return !(*this == a_other);
  }
};

/**
 * @brief Quantise a vertex to its identity.
 * @param[in] a_x Vertex position.
 * @return The vertex's key.
 */
inline VertexKey
vertexKey(const RealVect& a_x) noexcept;

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
