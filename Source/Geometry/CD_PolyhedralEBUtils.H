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
#include <BaseFab.H>
#include <BaseIF.H>
#include <Box.H>
#include <IntVect.H>
#include <RealVect.H>

// Our includes
#include <CD_CutCellSurface.H>
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

/**
 * @brief Which side of the interface a value lies on.
 * @details One predicate for corners and crossings alike: negative is fluid, and anything else,
 * an exact zero of either sign included, is solid. GeometryShop reads an exact zero by its sign
 * bit, which a negated distance function turns over, and then rescues the regular cell that
 * leaves beside a covered one by looking at its neighbour. A cell here is classified from its own
 * corners alone, so the zero has to fall on one side by rule: on the solid side, so that a fluid
 * cell whose face lies in the interface is cut with that face as its interface, and the cell on
 * the other side is covered.
 * @param[in] a_value Implicit function value.
 * @return True on the fluid side.
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
 * and solid carries four crossings and there are two ways to join them; the surface's face decision
 * (CutCellSurface::m_fluidJoined, read from the implicit function at the centre of the face) says
 * which diagonal meets through the middle, and the other is the pair the chords cut off. Both cells
 * adjoining the face read the same point, so they cannot disagree. With four crossings each chord
 * joins two consecutive edges of the face circuit.
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
 * @brief Distance within which two vertices are taken to be the same one.
 * @details A vertex where two polygons meet is computed once along each polygon's own edge, and
 * the two answers agree to about a part in a thousand billion rather than bit for bit. Loops are
 * stitched by matching one segment's end to the next one's start, so the match is made at that
 * distance and not exactly. Far below any edge that carries area -- a cell is one unit across in
 * this frame -- and far above the disagreement two computations of the same point produce.
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

} // namespace detail

class CutCellBody;

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
bool
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
void
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
void
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
Real
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
void
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
void
buildSurface(const BaseIF&                 a_function,
             BaseFab<Real>                 a_intercept[SpaceDim],
             PolyhedralEB::CutCellSurface& a_surface,
             const BaseFab<Real>&          a_nodeValues,
             const IntVect&                a_cell,
             const RealVect&               a_probLo,
             const Real&                   a_dx);

/**
 * @brief Decide, for every face of a cell whose corners alternate in sign, whether its fluid corners or its solid
 * ones meet across it.
 * @details Four alternating corners admit both readings, and the bilinear interpolant of the corners picks the one
 * its saddle falls on, which for a feature thinner than the face is the wrong one: the large values at two corners
 * outweigh the small ones of a thin plate at the other two, and the chords pierce the plate. The implicit function
 * at the centre of the face is read instead, with the value snap of a grid twice as fine, so that it is the very
 * node a finer level holds there. Every cell sharing the face reads the same point. A centre that reads exactly
 * zero falls back to the bilinear saddle. In two dimensions the face is the cell itself.
 * @param[in]     a_function Implicit function.
 * @param[in,out] a_surface  Surface whose corners are filled; its face decisions are set.
 * @param[in]     a_cell     The cell.
 * @param[in]     a_probLo   Lower-left corner of the domain.
 * @param[in]     a_dx       Grid spacing.
 */
void
decideFaces(const BaseIF&                 a_function,
            PolyhedralEB::CutCellSurface& a_surface,
            const IntVect&                a_cell,
            const RealVect&               a_probLo,
            const Real                    a_dx) noexcept;

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
Real
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
constexpr Real s_snapTolerance = 1.0E-12;

/**
 * @brief Distance along an edge, from an endpoint that is exactly zero, within which a root is put on it.
 * @details A surface tangent to a node reads as zero for a stretch either side of it under the value snap, of length
 * about sqrt(2 R tol / dx) for a radius of curvature R -- a few parts in 1e6 of an edge for any R that is
 * not itself below the cell width. A root that close to a corner the surface passes through is that
 * stretch's end, not a second crossing; a genuine second crossing within 1e-4 of a cell width of the first
 * would be a feature no cell resolves.
 */
constexpr Real s_rootSnap = 1.0E-4;

/**
 * @brief Fraction of a cell width the implicit function is read at either side of a point, to tell a surface
 * through the point from a zero no surface passes through.
 * @details Large enough that a function of unit gradient reads well clear of s_snapTolerance beside the point, and
 * small enough that the reading belongs to the point rather than to its neighbourhood.
 */
constexpr Real s_probeSpacing = 1.0E-6;

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
Real
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
Real
resolveTangency(const BaseIF& a_function, const RealVect& a_point, const Real a_dx) noexcept;
} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#include <CD_PolyhedralEBUtilsImplem.H>

#endif
