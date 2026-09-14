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
} // namespace detail
} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#include <CD_PolyhedralEBUtilsImplem.H>

#endif
