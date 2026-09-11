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
   * @param[in] a_refinement   How many times the generated cells are a refinement of the cells
   * the surface is reconstructed on. One reconstructs on the generated cells themselves.
   */
  PolyhedralGeometryShop(const BaseIF&        a_localGeom,
                         const int            a_verbosity,
                         const Real           a_dx,
                         const RealVect&      a_probLo,
                         const ProblemDomain& a_finestDomain,
                         const ProblemDomain& a_scanLevel,
                         const int            a_ebGhost,
                         const Real           a_thrshdVoF,
                         const bool           a_strict,
                         const int            a_refinement);

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
   * @brief Refinement between the cells the surface is reconstructed on and the ones generated.
   * @details At one, each cell's surface is reconstructed on the cell itself. At two, it is
   * reconstructed on the cell's parent and cut in half along each direction, so the generated
   * cells are literal sub-bodies of a coarser reconstruction and their moments sum back to it
   * exactly. Higher powers of two cut repeatedly.
   */
  int m_refinement;

  /**
   * @brief Classify every cell of a region from the bodies of the cells they subdivide.
   * @details Only meaningful when the generated cells refine the ones the surface is
   * reconstructed on. The surface is reconstructed once per coarse cell and cut down, and each
   * generated cell takes its classification from the piece of the coarse body it received. A
   * coarse cell that never met the interface hands the same answer to all of its children
   * without any body being built at all.
   *
   * Classifying this way rather than from the generated cell's own corner values keeps one
   * classification for the geometry, made at the resolution where the surface was
   * reconstructed and single-valuedness was mandated, rather than one per resolution.
   * @param[out] a_regIrregCovered Regular/Covered/Irregular cells
   * @param[out] a_irregularCells  The irregular cells within the valid region
   * @param[in]  a_validRegion     Grid region
   * @param[in]  a_ghostRegion     Grid region, including ghost cells
   * @param[in]  a_probLo          Lower-left corner of computational domain
   * @param[in]  a_dx              Grid resolution of the generated cells
   */
  void
  classifyFromParents(BaseFab<int>&   a_regIrregCovered,
                      IntVectSet&     a_irregularCells,
                      const Box&      a_validRegion,
                      const Box&      a_ghostRegion,
                      const RealVect& a_probLo,
                      const Real&     a_dx) const;

  /**
   * @brief Build a cell's body by cutting the body of the coarser cell containing it.
   * @details The coarse cell's surface is reconstructed at the coarser resolution and cut,
   * rather than the cell's own surface being reconstructed at its own resolution, so that the
   * cells of a refinement partition the geometry the coarser level sees.
   * @param[out] a_body    The cell's body, in its own frame
   * @param[in]  a_cell    The cell
   * @param[in]  a_probLo  Lower-left corner of computational domain
   * @param[in]  a_dx      Grid resolution of the generated cells
   * @return True if the coarse body closed and the cut succeeded.
   */
  bool
  buildRefinedBody(PolyhedralEB::CutCellBody& a_body,
                   const IntVect&             a_cell,
                   const RealVect&            a_probLo,
                   const Real&                a_dx) const;

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
   * @param[in]     a_probLo    Lower-left corner of computational domain
   * @param[in]     a_dx        Grid resolution
   * @return Position along the edge, parameterised from its low corner.
   */
  Real
  edgeCrossing(BaseFab<Real>   a_intercept[SpaceDim],
               const IntVect&  a_cell,
               const int       a_edge,
               const Real      a_lo,
               const RealVect& a_probLo,
               const Real&     a_dx) const;

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
