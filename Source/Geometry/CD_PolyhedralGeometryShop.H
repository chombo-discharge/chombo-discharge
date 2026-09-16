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
   * @brief The embedded boundary inside one cell, as triangles in physical coordinates.
   * @details The cell is reconstructed the same way generation reconstructs it, from the same
   *          implicit function and the same bracketing, so the triangles are the ones the index
   *          space was built from. Nothing is appended for a cell that is not cut.
   * @param[in,out] a_facets Coordinates to append to, nine reals per triangle.
   * @param[in]     a_cell   Cell to reconstruct.
   * @param[in]     a_probLo Lower-left corner of the domain.
   * @param[in]     a_dx     Grid spacing of the level the cell belongs to.
   */
  void
  interfaceFacets(Vector<Real>&   a_facets,
                  const IntVect&  a_cell,
                  const RealVect& a_probLo,
                  const Real      a_dx) const noexcept;

  /**
   * @brief Say how much of each level to carry.
   * @details A level is generated only over these boxes; elsewhere it is left out, and the levels
   * below it take what lies there through coarsening instead. The boxes are indexed coarsest
   * first, as the pre-pass in ComputationalGeometry returns them, and a_coarsestDomain says which
   * of this object's levels the first of them belongs to.
   * @param[in] a_regions        Boxes each level is carried over.
   * @param[in] a_coarsestDomain Domain of the first entry in a_regions.
   * @param[in] a_buffer         Cells to carry beyond what the boxes ask for, on every level.
   */
  void
  setCoverage(const Vector<Vector<Box>>& a_regions, const ProblemDomain& a_coarsestDomain, const int a_buffer) noexcept;

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
   * @brief How many numbers this generator keeps per cut cell so that the cell can be refined.
   * @details None unless each cell's surface is reconstructed on the cell itself, since that is
   * the only path that records one.
   * @return Values per cut cell, or zero if none are kept.
   */
  int
  numSurfaceComponents() const override;

  /**
   * @brief Hand back what was kept for each cut cell of a box.
   * @param[out] a_cells  Cut cells of the box, in the order their values appear.
   * @param[out] a_values Values per entry of a_cells.
   * @param[in]  a_region Cells wanted.
   * @param[in]  a_dx     Grid resolution of the level the cells belong to.
   */
  void
  getSurfaces(Vector<IntVect>& a_cells, Vector<Real>& a_values, const Box& a_region, const Real& a_dx) const override;

  /**
   * @brief Resolution of the level a region has to be cut from.
   * @param[in] a_ghostRegion Region wanted, including its ghost cells, on its own level.
   * @param[in] a_dx          Grid resolution of that region.
   * @return Resolution of the level to cut from, or a negative number if none describes it.
   */
  Real
  refinedFillParentDx(const Box& a_ghostRegion, const Real& a_dx) const override;

  /**
   * @brief Fill the graph and the nodes of a box on a level that was never generated.
   * @details The cells are cut from the surfaces of the coarser level handed in, rather than
   * reconstructed at their own resolution, so they are literal sub-bodies of cells the index
   * space already holds and their moments sum back to them.
   *
   * Each parent is cut once and its whole refinement kept, rather than a path being cut per
   * output cell: the children a path throws away are exactly the cells asked for next, and at
   * depth that is the difference between paying for the surface and paying for the volume.
   * @param[out] a_regIrregCovered Regular/Covered/Irregular cells
   * @param[out] a_nodes           Nodes
   * @param[in]  a_validRegion     Grid region
   * @param[in]  a_ghostRegion     Grid region, including ghost cells
   * @param[in]  a_domain          Domain
   * @param[in]  a_probLo          Lower-left corner of computational domain
   * @param[in]  a_dx              Grid resolution of the cells being filled
   * @param[in]  a_parents         Surfaces of the coarser level's cut cells
   * @param[in]  a_parentDx        Grid resolution of that coarser level
   * @return True if the box could be cut from the parents handed in.
   */
  bool
  fillRefinedGraph(BaseFab<int>&          a_regIrregCovered,
                   Vector<IrregNode>&     a_nodes,
                   const Box&             a_validRegion,
                   const Box&             a_ghostRegion,
                   const ProblemDomain&   a_domain,
                   const RealVect&        a_probLo,
                   const Real&            a_dx,
                   const BaseIVFAB<Real>& a_parents,
                   const Real&            a_parentDx) const override;

protected:
  /**
   * @brief Residual above which the stored moments and the ones PolyGeom recomputes disagree.
   */
  static constexpr Real s_divergenceTolerance = 1.0E-9;

  /**
   * @brief Marks a cell no parent has answered for yet.
   */
  static constexpr int s_unclassified = -2;

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
   * @brief Boxes each level is carried over, coarsest first.
   */
  Vector<Vector<Box>> m_coverageRegions;

  /**
   * @brief The cut cells of each box, in the order their surfaces were recorded.
   * @details Indexed by level, coarsest first, then by box. Kept alongside m_surfaces rather than
   * in one container of pairs, so that both are plain Vectors inside the LayoutData the layout
   * already provides.
   */
  Vector<RefCountedPtr<LayoutData<Vector<IntVect>>>> m_surfaceCells;

  /**
   * @brief What each cut cell was reconstructed from, kept for refining it later.
   * @details A cell finer than any the index space carries can only be had by cutting its parent,
   * and cutting a parent needs the parent's body, which is rebuilt from this. Reconstructing the
   * surface instead would mean root finding again -- of order a hundred implicit-function
   * evaluations a crossing -- which is the cost this exists to avoid. At 160 bytes a cut cell it
   * is a fraction of the body it rebuilds.
   */
  Vector<RefCountedPtr<LayoutData<Vector<PolyhedralEB::CutCellSurface>>>> m_surfaces;

  /**
   * @brief Keep a level only where it is carried.
   * @details Asked once for the whole refinement of a coarse box, so a coarse box is carried in
   * full or not at all: coarsening reads a coarse cell's whole refinement at once and skips any
   * box it cannot read in full.
   * @param[in] a_box   Refinement of a coarse box, on that level.
   * @param[in] a_level Level the box belongs to, finest first.
   * @return True if the box is kept.
   */
  virtual bool
  retainBox(const Box& a_box, const int a_level) const noexcept override;

  /**
   * @brief Size the surface store for a level once its layout is known.
   * @details Called by EBISLevel after it makes the grids and before it asks for any graph, which
   * is the only point at which this object learns the layout the surfaces will be recorded on.
   * @param[in] a_dbl Layout of the level.
   * @param[in] a_dx  Grid spacing of the level, which is what says which level it is.
   */
  virtual void
  postMakeBoxLayout(const DisjointBoxLayout& a_dbl, const RealVect& a_dx) override;

  /**
   * @brief Which level of this hierarchy a spacing belongs to.
   * @param[in] a_dx Grid spacing.
   * @return Index into m_dx, or -1 if the spacing is not one of this hierarchy's.
   */
  int
  levelFromDx(const Real a_dx) const noexcept;

  /**
   * @brief Index into m_domains of the level the first entry of m_coverageRegions belongs to.
   * @details Negative when no boxes have been supplied, which carries every level in full.
   */
  int m_coverageLevel;

  /**
   * @brief Cells carried beyond what the coverage boxes ask for, on every level.
   * @details Zero carries exactly what was asked. Carrying wider costs cut cells now and saves
   * cutting them again when the simulation's grids move into them, so the right value is the one
   * that keeps the index space ahead of the grids for a few regrids. It is quantised by the box
   * size: boxes are retained whole, so this only bites where it reaches across a box boundary.
   */
  int m_coverageBuffer;

  /**
   * @brief Bisect one cell edge for the crossing on it.
   * @param[in] a_edgeIV  Node at the low end of the edge.
   * @param[in] a_dir     Direction the edge runs in.
   * @param[in] a_loValue Implicit function at the low end.
   * @param[in] a_probLo  Lower-left corner of the domain.
   * @param[in] a_dx      Grid spacing.
   * @return Position along the edge, in [0,1].
   */
  Real
  edgeRoot(const IntVect& a_edgeIV, const int a_dir, const Real a_loValue, const RealVect& a_probLo, const Real& a_dx)
    const noexcept;

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
   * @brief Finest generated level that describes the whole of a region.
   * @details Which level a cell is cut from is settled by the geometry and not by which rank is
   * asking. Two levels that were both generated are not refinements of one another -- each was
   * built at its own resolution -- so two neighbouring boxes cut from different ancestors would
   * not agree on the face between them, and taking whichever level the asking rank happened to
   * hold would make that depend on the load balance.
   *
   * Costs a walk over the level's boxes. Worth revisiting if the caller ever knows the answer
   * already, which the level being extended generally does.
   * @param[in] a_ghostRegion Region to be filled, on its own level.
   * @param[in] a_level       Level the region belongs to, finest first.
   * @return Level to cut from, finest first, or -1 if no level describes the whole region.
   */
  int
  cuttableParent(const Box& a_ghostRegion, const int a_level) const noexcept;

  /**
   * @brief Whether a cell is one the generator gives a body and keeps the surface of.
   * @details Single-sourced because a cell finer than the level it is cut from is only reachable
   * where its parent's surface was kept, and a rule that says which cells those are in one place
   * and which cells were kept in another is a rule that will drift.
   *
   * A cut cell, and also a cell the interface only grazes: that cell is full, which is why
   * classify calls it regular, but the face the interface lies in is closed and a regular cell
   * has every face open.
   * @param[in] a_surface The values the cell would be reconstructed from.
   * @return True if the cell is treated as irregular.
   */
  static bool
  isRecorded(const PolyhedralEB::CutCellSurface& a_surface) noexcept;

  /**
   * @brief Turn cells that carry no fluid worth keeping into covered ones.
   * @details A dropped cell has just become covered, so the faces its irregular neighbours point
   * back with no longer lead anywhere, and the regular ones among them are now full cells with a
   * covered face. Where a neighbour's node is is read off a_nodes rather than tracked by the
   * caller, since the nodes carry the cell they belong to.
   * @param[in,out] a_nodes           Nodes, appended to for the regular cells this converts
   * @param[in,out] a_regIrregCovered Regular/Covered/Irregular cells
   * @param[in]     a_droppedCells    Cells that were dropped
   * @param[in]     a_validRegion     Grid region
   * @param[in]     a_domain          Domain
   * @param[in]     a_dx              Grid resolution
   */
  void
  dropCells(Vector<IrregNode>&   a_nodes,
            BaseFab<int>&        a_regIrregCovered,
            const IntVectSet&    a_droppedCells,
            const Box&           a_validRegion,
            const ProblemDomain& a_domain,
            const Real&          a_dx) const;

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
