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

// Chombo includes
#include <BaseFab.H>
#include <IntVectSet.H>

// Our includes
#include <CD_ScanShop.H>
#include <CD_CutCellBody.H>
#include <CD_CutCellSurface.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Geometry generation whose moments are exact integrals of an explicit surface.
 * @details Takes the same topology as GeometryShop -- the same classification, the same graph
 * and the same covered-cell fix-ups -- and replaces only the quadrature. Where GeometryShop
 * integrates the arc the interface really follows, this assembles the closed polyhedron whose
 * faces are bounded by straight chords between the edge crossings and integrates that.
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
   * @details Classifies the cells, reuses GeometryShop's fix-up for the regular cells bordering
   * a covered one, and then builds each remaining irregular cell's node from the polyhedron
   * rather than from a quadrature of the implicit function.
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
   * @brief Whether a cut body holds too little solid to be a cut cell.
   * @details The counterpart, on the fluid side, of dropping a body whose fluid volume is below the
   *          threshold: a body with neither solid volume nor boundary area is a regular cell. A cell whose
   *          interface lies in one of its faces, wholly or in part, has no solid volume either but does
   *          have boundary area, and is not dust.
   * @param[in] a_body      Body of the cell.
   * @param[in] a_threshold Volume-fraction threshold; nothing is dust at zero.
   * @return True if the cell should be treated as regular.
   */
  static bool
  isDust(const PolyhedralEB::CutCellBody& a_body, const Real a_threshold) noexcept;

  /**
   * @brief Bisect one cell edge for the crossing on it.
   * @details An endpoint at exactly zero is not taken to be the crossing: the function may be zero along
   *          a stretch of the edge (a face of the geometry lying in the edge's plane), and the crossing
   *          is where it leaves that stretch. The bisection walks through zeros, which the fluid
   *          predicate reads as solid, and a root that lands within s_snapTolerance of an endpoint that
   *          is exactly zero is put exactly on it.
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
   * @details A surface tangent to a node reads as zero for a stretch either side of it under the value
   *          snap, of length about sqrt(2 R tol / dx) for a radius of curvature R -- a few parts in 1e6 of
   *          an edge for any R that is not itself below the cell width. A root that close to a corner the
   *          surface passes through is that stretch's end, not a second crossing; a genuine second
   *          crossing within 1e-4 of a cell width of the first would be a feature no cell resolves.
   */
  static constexpr Real s_rootSnap = 1.0E-4;

  /**
   * @brief Read the implicit function at a point, with a value within s_snapTolerance of a cell width of
   *        zero read as exactly zero.
   * @details A face of the geometry lying in a node plane evaluates to zero there only up to round-off:
   *          the same face comes back as plus zero, minus zero and a few parts in 1e19 at neighbouring
   *          nodes, and cells on either side of it would then be classified by which of those they were
   *          handed. Reading everything within a fixed fraction of the cell width of zero as zero puts
   *          every such node on one side by rule, in the corner values and in the bisection alike.
   * @param[in] a_function Implicit function.
   * @param[in] a_point    Point to evaluate at.
   * @param[in] a_dx       Grid spacing that sets the scale.
   * @return The value, or exactly zero.
   */
  static Real
  snappedValue(const BaseIF& a_function, const RealVect& a_point, const Real a_dx) noexcept;

protected:
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
   * @brief Evaluate the implicit function at every node of a region.
   * @details One evaluation per node rather than per cell corner, since the cells meeting at a
   * node share it and must agree on its value.
   * @param[out] a_nodeValues Node-centred values over a_region
   * @param[in]  a_region     Cell-centred region whose nodes are wanted
   * @param[in]  a_probLo     Lower-left corner of computational domain
   * @param[in]  a_dx         Grid resolution
   */
  void
  fillNodeValues(BaseFab<Real>& a_nodeValues, const Box& a_region, const RealVect& a_probLo, const Real& a_dx) const;

  /**
   * @brief Find where the implicit function changes sign along a cell edge.
   * @details Bisection between the two endpoints, looked up in a cache first. Every cut edge is
   * reached from both of its adjoining faces and from each neighbouring cell, always with the
   * same endpoints, so the root is found once and reused. That is what makes two cells agree on
   * a shared aperture bit for bit.
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
  edgeCrossing(BaseFab<Real>   a_intercept[SpaceDim],
               const IntVect&  a_cell,
               const int       a_edge,
               const Real      a_lo,
               const Real      a_hi,
               const RealVect& a_probLo,
               const Real&     a_dx) const;

  /**
   * @brief Read one cell's corner values off the node values.
   * @details Corner c of the cell is the node offset by ((c >> d) & 1) in each direction d.
   * @param[out] a_surface    Surface whose corners are filled. Its crossings are left alone.
   * @param[in]  a_nodeValues Node-centred implicit function values
   * @param[in]  a_cell       The cell
   */
  void
  fillCorners(PolyhedralEB::CutCellSurface& a_surface, const BaseFab<Real>& a_nodeValues, const IntVect& a_cell) const;

  /**
   * @brief Gather the corner values and edge crossings one cell's surface is built from.
   * @param[in,out] a_intercept  Edge-centred crossing cache
   * @param[out]    a_surface    The values the cell's body is reconstructed from
   * @param[in]     a_nodeValues Node-centred implicit function values
   * @param[in]     a_cell       The cell
   * @param[in]     a_probLo     Lower-left corner of computational domain
   * @param[in]     a_dx         Grid resolution
   */
  void
  buildSurface(BaseFab<Real>                 a_intercept[SpaceDim],
               PolyhedralEB::CutCellSurface& a_surface,
               const BaseFab<Real>&          a_nodeValues,
               const IntVect&                a_cell,
               const RealVect&               a_probLo,
               const Real&                   a_dx) const;

  /**
   * @brief Write a cut cell's node from its body.
   * @details The arcs come from the covered set and the domain rather than from the moments, so
   * the graph this produces is the one GeometryShop would have produced. The apertures, the
   * centroids and the volume fraction come from the body.
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
