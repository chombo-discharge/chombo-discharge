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

// Std includes
#include <string>
#include <utility>
#include <vector>

// Chombo includes
#include <BaseFab.H>
#include <BaseIF.H>
#include <EBData.H>
#include <EBGraph.H>
#include <FArrayBox.H>
#include <IntVectSet.H>
#include <IrregNode.H>
#include <LevelData.H>

#ifdef CH_USE_HDF5
#include <hdf5.h>
#endif

// Our includes
#include <CD_AMRPolyhedralEBGraph.H>
#include <CD_ScanShop.H>
#include <CD_CutCellBody.H>
#include <CD_CutCellSurface.H>
#include <CD_GraphCopyHeader.H>
#include <CD_PolyhedralEBGraph.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_NamespaceHeader.H>

class ComputationalGeometry;

/**
 * @brief Geometry generation whose moments are exact integrals of an explicit surface.
 * @details Takes the same topology as GeometryShop -- the same classification, the same graph and the same covered-cell
 * fix-ups -- and replaces only the quadrature. Where GeometryShop integrates the arc the interface really
 * follows, this assembles the closed polyhedron whose faces are bounded by straight chords between the edge
 * crossings and integrates that.
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
   * @brief What to do with a cell this generator cannot describe at the refinement the mesh will ask for.
   * @details A cell whose fluid falls into more than one piece when it is cut for a refinement cannot be
   * carried single valued at that ratio, and there is no reformulation of the interior that avoids it: the
   * only representation guaranteed to keep every descendant single valued is a planar interface, and a planar
   * interface cannot meet the face chords a cell shares with its neighbours. So the cell is a choice rather
   * than a bug, and the choice belongs to whoever set up the run.
   */
  enum class UnresolvedAction
  {
    /**
     * @brief Report the cells and carry on, leaving it to the user to raise the geometric resolution.
     */
    Warn,

    /**
     * @brief Turn the cell solid, so that neither it nor anything cut from it is ever multi valued.
     */
    Cover,

    /**
     * @brief Keep the pieces in one body, taking every moment from their sum and moving only the volume centroid.
     * @details For a cell whose interface is one sheet twisted into a saddle, so that cutting it for a refinement
     * leaves a child holding its fluid in two pieces. The pieces are joined a cell away, so the saddle is resolved
     * geometry and is merged rather than covered. A cell holding more than one sheet is still filled.
     *
     * Conserved: the volume fraction, the apertures, the boundary area vector and the boundary centroid are the
     * exact sums of the pieces', satisfy the divergence identity, and coarsen back to the parent exactly. Face
     * centroids are left as the exact first moment of the open face, even where two patches put it on covered face
     * between them, since the flux stencils read it as a quadrature point.
     *
     * Not conserved: the scalar boundary area, the length of the summed area vector, is smaller than the area of the
     * interface itself when the pieces face different ways; the interface's own area is kept beside it. And the
     * volume centroid, which is the one moment moved.
     *
     * The repair: the volume centroid is the exact first moment of the merged fluid, a weighted average of the
     * pieces' centroids, and it can fall in the solid between them. Where it lies in the fluid it is kept. Where it
     * does not, it is replaced by the nearest point of the fluid, moved slightly inside it. That point no longer
     * sums back to the parent's first moment, so coarsening this body's centroid misses the parent's by the distance
     * moved, weighted by this body's volume.
     */
    Collapse
  };

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

  using ScanShop::InsideOutside;

  /**
   * @brief Whether a box is all regular, all covered or holds the surface, as EBISLevel asks it of a region.
   * @details Answered from the boxes the geometry's grids classify (ComputationalGeometry::classify), which the
   * graph's cut tiles come from, so the index space treats as uniform exactly what the graph does. The scan
   * ScanShop answers from samples the implicit function at each level's own spacing and reads it as a distance, so it
   * can call a coarse box covered while the coarsened finer cells hold fluid. EBISLevel fills a ghost graph from that
   * answer before it copies the coarsened graph over it, and then reads a regular cell beside a cell with no face
   * toward it. EBISLevel also passes a region of the level above or below on a domain, so the region is classified
   * as lying on the domain, one level finer and one level coarser, and is uniform only if all three agree. A domain
   * coarser than every level the grids describe is answered as irregular, which is always safe. Before setGrids,
   * ScanShop answers.
   * @param[in] a_region Region to classify.
   * @param[in] a_domain Domain of the region.
   * @param[in] a_probLo Lower-left corner of the domain.
   * @param[in] a_dx     Grid spacing.
   * @param[in] a_dit    Grid index.
   * @return Regular, covered or irregular.
   */
  GeometryService::InOut
  InsideOutside(const Box&           a_region,
                const ProblemDomain& a_domain,
                const RealVect&      a_probLo,
                const Real&          a_dx,
                const DataIndex&     a_dit) const override;

  /**
   * @brief Fill the graph and the nodes for one box.
   * @details Classifies the cells, reuses GeometryShop's fix-up for the regular cells bordering a covered one, and then
   * builds each remaining irregular cell's node from the polyhedron rather than from a quadrature of the
   * implicit function.
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
   * @brief Say which grids the polyhedra are built over, and which phase this generator serves on them.
   * @param[in] a_compGeom Builder holding the boxes and their classifications per level.
   * @param[in] a_phase Phase whose classifications and implicit function this generator uses.
   */
  void
  setGrids(const ComputationalGeometry& a_compGeom, const phase::which_phase a_phase) noexcept;

  /**
   * @brief Build the graph of every tiled level: the cut cells, their surfaces, and what every face opens onto.
   * @details One PolyhedralEBGraph per level from the geometry's start level up, over that level's cut tiles, each
   * defined from this phase's implicit function with one ghost cell and then linked to the level above it so that
   * the coarse side of every level boundary knows which of its cells and faces the finer level describes. The levels
   * below the start level carry no graph. Collective. Requires setGrids.
   */
  void
  buildGraphs();

  /**
   * @brief The graph of a level, undefined below the start level. Stops the run before buildGraphs.
   * @param[in] a_level Level, 0 being the coarsest of the geometry's levels.
   * @return The graph.
   */
  const PolyhedralEBGraph&
  getGraph(const int a_level) const noexcept;

  /**
   * @brief Compute the index-space data on the geometry's own grids, from the graph and nothing else.
   * @details One LevelData<EBGraph> and one LevelData<EBData> per level, over every box the geometry built on it, as
   * EBISLevel keeps them: the graph with one ghost cell, the moments with none. The graph is copied onto the grids
   * first, so that every box holds what it needs. Finest level first: a level's own cells take what the graph says --
   * regular, covered, or the moments of the stitched body -- and its cells under the next finer level are coarsened
   * from that level by Chombo's coarsening (coarsenLevel), then joined to the level's own (assembleLevel). Under
   * PolyhedralGeometryShop.test_copy every coarsened cell is checked against the sums over its children. Collective.
   * Requires buildGraphs.
   */
  void
  buildMoments();

  /**
   * @brief The index-space graph buildMoments made for a level.
   * @param[in] a_level Level, 0 being the coarsest of the geometry's levels.
   * @return The graph, over the level's boxes.
   */
  const LevelData<EBGraph>&
  getEBGraph(const int a_level) const;

  /**
   * @brief The index-space moments buildMoments made for a level.
   * @param[in] a_level Level, 0 being the coarsest of the geometry's levels.
   * @return The moments, over the level's boxes.
   */
  const LevelData<EBData>&
  getEBData(const int a_level) const;

  /**
   * @brief Fill any layout on any level with its index space: graph and first moments, ghost cells included.
   * @details The graph is copied onto the layout, and then each cell in the geometry's boxes on the level takes what
   * buildMoments made there, and each cell outside them is cut from the coarser cell it lies in. Nothing evaluates the
   * implicit function. Collective. Requires buildMoments.
   * @param[out] a_graph      The graph, defined here over the layout with a_graphGhost ghost cells.
   * @param[out] a_data       The moments, defined here over the layout with a_dataGhost ghost cells.
   * @param[in]  a_grids      The layout.
   * @param[in]  a_domain     Its domain: one of the geometry's, or one refined from the finest by two per level.
   * @param[in]  a_graphGhost Ghost width of the graph, above a_dataGhost.
   * @param[in]  a_dataGhost  Ghost width of the moments.
   */
  void
  fillLayout(LevelData<EBGraph>&      a_graph,
             LevelData<EBData>&       a_data,
             const DisjointBoxLayout& a_grids,
             const ProblemDomain&     a_domain,
             const int                a_graphGhost,
             const int                a_dataGhost) const;

  /**
   * @brief Build the polyhedra over the grids, check that they form a surface, and write it if asked.
   * @details Every cut cell of every irregular box is reconstructed at its own level, cells the next finer level
   * carries are skipped, and a face shared with such a cell is restricted to the finer chords and the
   * interface closed. The interface triangles of every level are then checked by sanityCheck when
   * PolyhedralGeometryShop.sanity_check is set, and written to the geometry directory as an indexed surface mesh
   * when PolyhedralGeometryShop.write_surface is set. Collective; every rank reconstructs its share of the boxes
   * and writes its own share of the mesh. Requires setGrids.
   */
  void
  verifySurface() const;

  /**
   * @brief Check that the interface triangles form a surface that is closed except on the domain boundary.
   * @details A triangle edge can be shared only by cells that touch, so every edge of a cell's interface must be
   * used exactly twice among the interface triangles of the cell and of its 3^D neighbourhood: on the cell's own
   * level through the ghost surfaces the graph holds, and across a level boundary through the finer cells behind a
   * face the finer level describes, reconstructed for the purpose. An edge in the domain boundary is used once, an
   * edge in a face the coarser level describes is checked from that side, and an edge in a face onto a cell the graph
   * filled has no partner to meet, since that face is closed. Every rank checks its own cells, vertices are welded at
   * s_weldSpacing times the level's spacing, the counts are summed, and a violation stops the run. A regular cell
   * sharing a face with a covered one stops the run as well: the four nodes of that face belong to both cells and
   * cannot be fluid for one and solid for the other, so the pair can only come from a classification that disagrees
   * with itself, and there would be no cut cell to carry the interface between them. A filled cell is covered whatever
   * its nodes read, which is why the graph keeps a neighbour of one as a cut cell holding the whole of itself rather
   * than a regular one; the test is what says it did. That test is dimension-independent and runs on every graph,
   * including the cells beyond the tiles, which is where the edge of the tiled region meets what no tile carries.
   * Switched by PolyhedralGeometryShop.sanity_check, on by default in a build that keeps its assertions. Requires
   * buildGraphs.
   */
  void
  sanityCheck() const;

  /**
   * @brief Report the cells that cannot be cut for a refinement without falling into more than one piece.
   * @details Asked of every cut cell the level itself describes, at the ratio m_maxRefinement names, and of the
   * body the consumers build -- which on a refinement boundary is the one whose faces were taken from the finer
   * level, since stitching lengthens the interface loop and it is that loop the question is about. A cell whose
   * interface is planar is exempt without being cut, because its fluid is then convex and convexity survives
   * every ratio.
   *
   * What happens to a failing cell is m_unresolvedAction's to decide. Only Warn is wired: the cells are named
   * and the run continues, leaving it to whoever set it up to raise the geometric resolution. Runs regardless
   * of m_sanityCheck, since it reports on the geometry the run will actually use rather than on this
   * generator's own consistency. Requires buildGraphs.
   */
  void
  reportUnresolvedRefinement() const;

  /**
   * @brief Copy every graph onto a layout that splits each of its tiles into octants with the ranks shuffled, and then
   * onto the tiles themselves with the ranks shuffled; run sanityCheck over the copies each time, copy them back onto
   * the original layouts and require equality bit for bit.
   * @details A stand-in, until the index space is generated over the graphs, for what a consumer on another
   * layout does: take the graph's data through a Copier and read it there. Switched by
   * PolyhedralGeometryShop.test_copy. Collective. Requires buildGraphs.
   */
  void
  testGraphCopy() const;

  /**
   * @brief The geometry level a destination domain lives on.
   * @details A stored domain's own level; a domain finer than every stored one is a refinement by two per level of
   * the finest stored domain, and lives that many levels past it. Any other domain stops the run.
   * @param[in] a_domain The domain.
   * @return Its geometry level.
   */
  int
  destinationLevel(const ProblemDomain& a_domain) const;

  /**
   * @brief The cells of a destination box a copy of the graph fills: the box grown by a ghost width, inside the
   * domain, less what the next finer destination level covers.
   * @param[out] a_valid  The cells.
   * @param[in]  a_box    The destination box.
   * @param[in]  a_domain Domain of the box's level.
   * @param[in]  a_finer  Layout of the next finer destination level, or null where there is none.
   * @param[in]  a_finerBoxes Its boxes, as boxArray gives them, or null where there is no finer level.
   * @param[in]  a_ratio  Refinement ratio to that level.
   * @param[in]  a_ghost  Ghost width the box is grown by.
   */
  static void
  validRegion(IntVectSet&              a_valid,
              const Box&               a_box,
              const ProblemDomain&     a_domain,
              const DisjointBoxLayout* a_finer,
              const Vector<Box>*       a_finerBoxes,
              const int                a_ratio,
              const int                a_ghost);

  /**
   * @brief The pieces of the graph a destination region needs: what each tiled level owns of it, tile by tile, down to
   * the cells whose footprint on the destination level lies in the region.
   * @details A cell of a level at or finer than the destination's is kept if the destination cell it lies in is in the
   * region; a cell of a coarser level if any destination cell it holds is.
   * @param[out] a_pieces One list per geometry level of (index into the level's cut tiles, box on that level) pairs.
   * @param[in]  a_valid  The region, on the destination level, from validRegion.
   * @param[in]  a_level  Geometry level of the destination, from destinationLevel.
   */
  void
  destinationPieces(std::vector<std::vector<std::pair<int, Box>>>& a_pieces,
                    const IntVectSet&                              a_valid,
                    const int                                      a_level) const;

  /**
   * @brief Plan a copy of the graph onto a destination hierarchy: every piece this rank sends, and every piece it
   * receives.
   * @details The receiver walks its own destination boxes and asks each one which pieces it needs. The sender walks
   * its own graph boxes, finds the destination boxes that could need them through the destination layouts' own spatial
   * index, and asks each of those the same question, keeping the pieces its box holds. Both answer from replicated
   * data alone, so the two sides agree with no communication. Every header's counts are left at zero.
   * @param[out] a_sends    Per destination rank, the pieces this rank sends there.
   * @param[out] a_receives Per source rank, the pieces this rank receives from there.
   * @param[in]  a_grids    Destination layout of every level, coarsest first.
   * @param[in]  a_domains  Domain of every destination level.
   * @param[in]  a_ghost    Ghost width every destination box is filled to.
   */
  void
  planCopy(std::vector<std::vector<PolyhedralEB::GraphCopyHeader>>& a_sends,
           std::vector<std::vector<PolyhedralEB::GraphCopyHeader>>& a_receives,
           const std::vector<DisjointBoxLayout>&                    a_grids,
           const std::vector<ProblemDomain>&                        a_domains,
           const int                                                a_ghost) const;

  /**
   * @brief Plan a copy onto a destination hierarchy, send every rank the headers it is to receive, and require them to
   * be exactly the ones it planned to receive. Collective.
   * @param[in] a_grids   Destination layout of every level, coarsest first.
   * @param[in] a_domains Domain of every destination level.
   * @param[in] a_ghost   Ghost width every destination box is filled to.
   * @return Number of pieces planned, over all ranks.
   */
  long long
  checkCopyPlan(const std::vector<DisjointBoxLayout>& a_grids,
                const std::vector<ProblemDomain>&     a_domains,
                const int                             a_ghost) const;

  /**
   * @brief Check the copy plan onto the geometry's own box hierarchy with the ranks shuffled, onto every level of it
   * alone, and onto a level past the finest. Collective. Switched by PolyhedralGeometryShop.test_copy.
   */
  void
  testCopyPlan() const;

  /**
   * @brief Copy the graphs onto a destination hierarchy.
   * @details Plans the copy and sends the pieces in two passes, each destination box keeping exactly the pieces it
   * receives. The first carries every cell's state and the surface of every cell that holds one, its role
   * included; the second carries the face overrides of every cell whose surface is overridden, which the receiver
   * knows from the first. Nothing on either side evaluates the implicit function. A piece that arrives unplanned, a
   * cell whose role is invalid, and an overridden cell without overrides all stop the run. Collective.
   * @param[out] a_dest    The graph on the destination hierarchy.
   * @param[in]  a_grids   Destination layout of every level, coarsest first.
   * @param[in]  a_domains Domain of every destination level.
   * @param[in]  a_ghost   Ghost width every destination box is filled to.
   */
  void
  copyGraphs(AMRPolyhedralEBGraph&                 a_dest,
             const std::vector<DisjointBoxLayout>& a_grids,
             const std::vector<ProblemDomain>&     a_domains,
             const int                             a_ghost) const;

  /**
   * @brief Copy every graph onto its own layout as a destination level of its own, and require every cell the level
   * owns to arrive as the graph holds it, bit for bit: state, surface, role and face overrides. Collective. Switched by
   * PolyhedralGeometryShop.test_copy.
   */
  void
  testCopyRoundTrip() const;

  /**
   * @brief Build the body of one cut cell of a geometry level, from what a destination box holds of that level.
   * @param[out] a_body        The body.
   * @param[in]  a_graph       The graph on the destination hierarchy.
   * @param[in]  a_level       Destination level.
   * @param[in]  a_dit         Destination box.
   * @param[in]  a_sourceLevel Geometry level of the cell.
   * @param[in]  a_cell        The cell, on that level.
   * @return True if the box holds the cell and it is cut; false, with the body untouched, otherwise.
   */
  static bool
  sourceBody(PolyhedralEB::CutCellBody&  a_body,
             const AMRPolyhedralEBGraph& a_graph,
             const int                   a_level,
             const DataIndex&            a_dit,
             const int                   a_sourceLevel,
             const IntVect&              a_cell);

  /**
   * @brief The bodies a cell of a destination level is made of, from the geometry levels that own it.
   * @details A cell owned at its own resolution is its own body. One owned by a finer level is the bodies of every cut
   * cell under it, each in its own frame and none merged, so the surface they bound is the finer level's exactly. One
   * owned by a coarser level is the cut ancestor subdivided down to the destination level, and the child that is this
   * cell, in this cell's frame; a child whose fluid falls into more than one piece is kept, and counted. A cell no
   * level owns, or owned by a level that holds no cut cell there, has no body.
   * @param[out]    a_bodies        The bodies.
   * @param[in,out] a_numSplit      Incremented for every subdivision on the way down that left a child whose fluid is
   * in more than one piece, or whose children did not partition their parent.
   * @param[in]     a_graph         The graph on the destination hierarchy.
   * @param[in]     a_level         Destination level.
   * @param[in]     a_dit           Destination box.
   * @param[in]     a_cell          The cell, on the destination level.
   */
  static void
  destinationBodies(std::vector<PolyhedralEB::CutCellBody>& a_bodies,
                    long long&                              a_numSplit,
                    const AMRPolyhedralEBGraph&             a_graph,
                    const int                               a_level,
                    const DataIndex&                        a_dit,
                    const IntVect&                          a_cell);

  /**
   * @brief Require every cut cell's body to survive a copy bit for bit.
   * @details Every cut cell every geometry level owns is built from the graph as it stands and hashed, bit for bit,
   * and the hashes are summed over all ranks. Then the graphs are copied onto the geometry's own box hierarchy with the
   * ranks shuffled, onto every level of it alone, and onto a level past the finest, and every valid cell of every
   * destination -- ghost cells left out -- rebuilds the source bodies it stands for from what it received: its own,
   * those of the finer cells under it, and, from the one destination cell at its low corner, a coarser ancestor's. Each
   * owned cell is then counted exactly once, and the sums and counts must be the source's. Every coarser cut cell is
   * also subdivided down to the destination level, and the child destinationBodies gives the destination cell at its
   * low corner must be the one found there. The hierarchy is also copied with every box split into octants, which must
   * change nothing.
   *
   * On every face between two valid cells owned by different levels, the two sides must read the same aperture to
   * 1e-10: the claim that a subdivided stitched ancestor meets the finer cells beside it.
   *
   * All of it -- the copies and every body rebuilt from them -- runs with the shop's implicit function replaced by one
   * that stops the run if it is evaluated. Collective. Switched by PolyhedralGeometryShop.test_copy.
   */
  void
  testCopyBodies() const;

  /**
   * @brief Check the index-space data buildMoments made against the graph it came from.
   * @details On every level built, each cell is what the graph says it is, and each VoF of a cut cell carries the
   * moments of the cell's stitched body exactly: volume fraction and centroid, boundary centroid, and per face the
   * aperture and centroid. Collective. Switched by PolyhedralGeometryShop.test_copy.
   */
  void
  testMoments() const;

  /**
   * @brief Check fillLayout on every level of the geometry and one past the finest, onto the level's boxes in octants
   * and onto the boxes of the level below refined.
   * @details The cells the geometry's level holds must arrive bit for bit, every ghost cell must equal the cell its own
   * box holds to round-off, every face must read the same aperture from its two sides, and every cell of the level
   * below must be the sum of the cells it holds, to round-off. Collective. Switched by
   * PolyhedralGeometryShop.test_copy.
   */
  void
  testFill() const;

  /**
   * @brief Fill a destination level of a copied graph with its index space: graph and first moments, ghost cells
   * included.
   * @details A cell in the geometry's boxes on the level takes what buildMoments made there -- its own body, or the
   * VoFs coarsened from the finer level -- and a cell outside them is cut from the coarser cell it lies in, from the
   * copied graph. A regular cell beside a multi-valued one is written as a full cut cell. A cut cell a finer level owns
   * outside the geometry's boxes on the level stops the run. Collective.
   * @param[out] a_graph      The graph, defined here over the level's layout with a_graphGhost ghost cells.
   * @param[out] a_data       The moments, defined here over the level's layout with a_dataGhost ghost cells.
   * @param[in]  a_copy       The graph copied onto a destination hierarchy, with a ghost width above a_graphGhost.
   * @param[in]  a_level      Destination level to fill.
   * @param[in]  a_graphGhost Ghost width of the graph, above a_dataGhost.
   * @param[in]  a_dataGhost  Ghost width of the moments.
   */
  void
  fillLevel(LevelData<EBGraph>&         a_graph,
            LevelData<EBData>&          a_data,
            const AMRPolyhedralEBGraph& a_copy,
            const int                   a_level,
            const int                   a_graphGhost,
            const int                   a_dataGhost) const;

  /**
   * @brief Coarsen a level of the index space onto its coarsened layout.
   * @details Reads the finished fine level only. The children of a coarse cell are all inside the fine boxes, which
   * are blocked by two. Each coarse cell's VoFs are the pieces its children's VoFs fall into when joined through the
   * fine faces between sibling cells, and each carries the sums of its fine VoFs' volumes and first moments, and its
   * boundary centroid weighted by their boundary areas. A coarse face is split by the VoF across that its fine faces
   * lead into: a piece of the coarse cell across when that cell is coarsened too, the one VoF of a cell the level owns,
   * or the domain boundary; its aperture is the scaled sum of those fine apertures, and its centroid their
   * aperture-weighted mean. Every sum runs in an order fixed by the fine level, so every box writes a coarse face with
   * the same bits. Coarse nodes are written for the valid cells only. Under PolyhedralGeometryShop.test_copy every
   * coarse cell's volume, first moment and apertures are compared with the sums over its children, and a cell differing
   * by more than round-off is counted.
   * @param[out] a_coarGraph  Coarsened graph, defined here over the coarsened fine layout.
   * @param[out] a_coarData   Coarsened moments, defined here over the coarsened fine layout.
   * @param[out] a_numInexact Incremented for every coarse cell that is not the sum of its children.
   * @param[in]  a_grids      The geometry's grids.
   * @param[in]  a_domains    Their domains.
   * @param[in]  a_fineLevel  Level to coarsen, already built.
   */
  void
  coarsenLevel(LevelData<EBGraph>&                   a_coarGraph,
               LevelData<EBData>&                    a_coarData,
               long long&                            a_numInexact,
               const std::vector<DisjointBoxLayout>& a_grids,
               const std::vector<ProblemDomain>&     a_domains,
               const int                             a_fineLevel) const;

  /**
   * @brief Build a level of the index space from its own cells and the cells coarsened from the level above.
   * @details The level's own cells take what the copied graph says. Its cells under the next finer level take the
   * coarsened VoFs, and an own cell beside one takes the faces between them from the coarsened side, which knows which
   * VoF each leads into. A regular cell, own or coarsened, beside a cell of more than one VoF is written as a full cut
   * cell, since a regular cell connects to one VoF only; this changes no moment.
   * @param[out] a_numPromoted     Incremented for every regular cell written as a full cut cell.
   * @param[in]  a_copy            The graph copied onto the geometry's grids.
   * @param[in]  a_coarsenedGraph  coarsenLevel's graph for the level above, or null on the finest level.
   * @param[in]  a_coarsenedData   coarsenLevel's moments for the level above, or null on the finest level.
   * @param[in]  a_grids           The geometry's grids.
   * @param[in]  a_domains         Their domains.
   * @param[in]  a_level           Level to build.
   */
  void
  assembleLevel(long long&                            a_numPromoted,
                const AMRPolyhedralEBGraph&           a_copy,
                const LevelData<EBGraph>*             a_coarsenedGraph,
                const LevelData<EBData>*              a_coarsenedData,
                const std::vector<DisjointBoxLayout>& a_grids,
                const std::vector<ProblemDomain>&     a_domains,
                const int                             a_level);

  /**
   * @brief A layout over given boxes with the ranks shuffled: they walk the rank list with a stride coprime to its
   * length, so that no rank holds what a load balance would have given it. For the copy tests.
   * @param[in] a_boxes  The boxes.
   * @param[in] a_domain Their domain.
   * @return The closed layout.
   */
  static DisjointBoxLayout
  shuffledLayout(const Vector<Box>& a_boxes, const ProblemDomain& a_domain);

  /**
   * @brief Write a cut cell's node from its body.
   * @details The arcs come from the covered set and the domain rather than from the moments, so the graph this produces
   * is the one GeometryShop would have produced. The apertures, the centroids and the volume fraction come
   * from the body. Reads nothing of the shop, so the copied graph writes its nodes with it too.
   * @param[out] a_node             The node to fill
   * @param[in]  a_body             The cell's assembled polyhedron
   * @param[in]  a_regIrregCovered  Regular/Covered/Irregular cells
   * @param[in]  a_cell             The cell
   * @param[in]  a_domain           Domain
   */
  static void
  fillNode(IrregNode&                       a_node,
           const PolyhedralEB::CutCellBody& a_body,
           const BaseFab<int>&              a_regIrregCovered,
           const IntVect&                   a_cell,
           const ProblemDomain&             a_domain);

protected:
  /**
   * @brief Vertices of one element of the written surface: a triangle in three dimensions, a segment in two.
   */
  static constexpr int s_elementVertices = SpaceDim;

  /**
   * @brief Width of the ring the graph is copied onto the geometry's grids with for buildMoments.
   * @details assembleLevel builds each box with one ghost cell, fillNodes writes nodes one cell short of the copied
   * ring, and a node's face toward a cut cell below it reads that cell's body, one cell further out. Nothing to do with
   * the ghost cells a solver asks of the index space.
   */
  static constexpr int s_momentRing = 2;

  /**
   * @brief Number of components cellSignature writes: the VoF count, the volume, the first moment and the apertures.
   */
  static constexpr int s_signatureComps = 2 + 3 * SpaceDim;

  /**
   * @brief Kind of a node of the written graph: a cut cell, at its fluid's centroid.
   */
  static constexpr int s_cellNode = 0;

  /**
   * @brief Kind of a node of the written graph: an open face of a cut cell, at the centroid of its open part.
   */
  static constexpr int s_faceNode = 1;

  /**
   * @brief What a field of the written graph reads on a node of the kind it does not belong to.
   */
  static constexpr int s_noValue = -1;

  /**
   * @brief Fraction of a level's grid spacing within which two vertices are the same vertex, in sanityCheck and in
   * the written surface.
   * @details Cells of one level compute a shared vertex from the same nodes and edge and agree bit for bit, and cells
   * rebuilt across a level boundary agree to about 1e-12 of a cell, so the radius sits just above that. A larger one
   * welds the distinct vertices of a feature thinner than it into one, and the check then reports its edges as used
   * more than twice.
   */
  static constexpr Real s_weldSpacing = 1.0E-10;

  /**
   * @brief The geometry whose grids the polyhedra are built over, or null before setGrids.
   */
  const ComputationalGeometry* m_compGeom;

  /**
   * @brief State of a cell of the finer level that the finer level does not hold, in a copy of its states.
   */
  static constexpr int s_absent = 2;

  /**
   * @brief Trade byte buffers among all ranks: counts first, then the payloads, into one flat buffer.
   * @details What this rank sends itself is copied rather than sent. Byte totals that do not fit in an int, which is
   * what MPI counts are, stop the run in every build. Collective.
   * @param[in]  a_send   Per destination rank, the bytes this rank sends there.
   * @param[out] a_recv   Everything this rank received, source rank after source rank.
   * @param[out] a_offset Per source rank, where its bytes begin in a_recv.
   * @param[out] a_count  Per source rank, how many bytes it sent.
   */
  static void
  exchangeBytes(const std::vector<std::vector<char>>& a_send,
                std::vector<char>&                    a_recv,
                std::vector<int>&                     a_offset,
                std::vector<int>&                     a_count);

  /**
   * @brief Graph of every level, indexed as the geometry indexes its levels; undefined where a level has no tiles.
   */
  Vector<RefCountedPtr<PolyhedralEBGraph>> m_graphs;

  /**
   * @brief Index-space graph of every level buildMoments built, over the geometry's boxes; null where it built none.
   */
  Vector<RefCountedPtr<LevelData<EBGraph>>> m_ebGraphs;

  /**
   * @brief Index-space moments of every level buildMoments built, over the geometry's boxes; null where it built none.
   */
  Vector<RefCountedPtr<LevelData<EBData>>> m_ebData;

  /**
   * @brief Phase this generator serves on those grids.
   */
  phase::which_phase m_phase;

  /**
   * @brief Whether verifySurface writes the surface. PolyhedralGeometryShop.write_surface, off by default.
   */
  bool m_writeSurface;

  /**
   * @brief Whether verifySurface writes the graph as a skeleton. PolyhedralGeometryShop.write_graph, off by default.
   */
  bool m_writeGraph;

  /**
   * @brief Whether verifySurface runs sanityCheck. PolyhedralGeometryShop.sanity_check.
   * @details On by default in a build that keeps its assertions and off in one that does not, since the check
   * tells the generator nothing it uses and only reports. Either default can be overridden from the inputs.
   */
  bool m_sanityCheck;

  /**
   * @brief Whether verifySurface reports its timings. PolyhedralGeometryShop.profile, off by default.
   */
  bool m_profile;

  /**
   * @brief How far below the geometry's own grids a cell must stay single valued.
   * PolyhedralGeometryShop.max_refinement, one -- that is, off -- by default.
   * @details A power of two. One asks for nothing and skips the check entirely; two asks that the cells one
   * refinement below the geometry grid are each a single piece of fluid, four that the cells two below are,
   * and so on. The property is not inherited -- a cell whose children are all single valued can still have a
   * grandchild that is not -- so this has to name the ratio the mesh will actually use rather than being
   * derived from one level.
   *
   * Off by default because the test costs a subdivision of every cut cell whose interface is not planar, which
   * is two to three times what building the geometry costs in the first place, and because nothing downstream
   * reads a subdivided cell yet. It becomes worth its price once it does.
   */
  int m_maxRefinement;

  /**
   * @brief What becomes of a cell that fails that test. PolyhedralGeometryShop.unresolved_cells, warn by
   * default.
   */
  UnresolvedAction m_unresolvedAction;

  /**
   * @brief Whether a cell whose fluid falls into more than one piece is always filled (true), or only when its fluid is
   * the larger part of it (false). PolyhedralGeometryShop.fill_split_fluid, always by default, or majority.
   */
  bool m_fillSplitFluid;

  /**
   * @brief Whether verifySurface runs testGraphCopy. PolyhedralGeometryShop.test_copy, off by default.
   */
  bool m_testCopy;

  /**
   * @brief Whether the member functions that run once per level or per phase announce themselves in pout.
   * PolyhedralGeometryShop.verbose, off by default. The functions that run per cell never announce themselves.
   */
  bool m_verbose;

  /**
   * @brief Interface triangles of one level, or segments in two dimensions, from this rank's share of its irregular
   * boxes, three reals per vertex.
   * @details Cells the next finer level carries are skipped. A cut cell with such a cell across one of its faces has
   * that face restricted to the finer cells' chords (restrictFace) and its interface rebuilt (closeInterface).
   * A coarse edge that carries no crossing while both finer edges under it do is a multi-valued coarse cell in
   * disguise and stops the run, as does a body that will not close.
   * @param[out] a_facets Element coordinates, appended to.
   * @param[in]  a_level  Level to collect from.
   */
  void
  collectFacets(Vector<Real>& a_facets, const int a_level) const;

  /**
   * @brief Stitch every level of a set of graphs to the level above it, from the graphs alone.
   * @details In two dimensions no face is stitched and every level's overrides stay empty. Collective.
   * @param[in] a_graphs Graphs, one per level, linked to one another; undefined entries are skipped.
   */
  void
  stitchGraphs(const Vector<RefCountedPtr<PolyhedralEBGraph>>& a_graphs) const;

  /**
   * @brief The closure check over a given set of graphs, one per level as m_graphs holds them.
   * @param[in] a_graphs Graphs, undefined entries skipped.
   */
  void
  sanityCheck(const Vector<RefCountedPtr<PolyhedralEBGraph>>& a_graphs) const;

  /**
   * @brief Which of a cell's edges the finer level describes.
   * @details An edge is shared by the cells that meet along it, and the crossing on it has to be one point for all of
   * them, so the question is asked of the edge: it is refined if any cell meeting it is.
   * @param[out] a_edgeRefined One flag per edge of the cell.
   * @param[in]  a_refined     Refined mask around the cell.
   * @param[in]  a_domain      Domain box of the cell's level.
   * @param[in]  a_cell        The cell.
   * @return True if any edge is refined.
   */
  static bool
  edgesRefined(bool                        a_edgeRefined[PolyhedralEB::CutCellSurface::s_numEdges],
               const BaseFab<signed char>& a_refined,
               const Box&                  a_domain,
               const IntVect&              a_cell);

  /**
   * @brief Stitch a cut cell's body to the finer level and to the covered cells around it.
   * @details Every refined edge takes its crossing from the half the finer level put it in, which rewrites a_surface.
   * The body is then defined from it -- as the whole cell where every corner is fluid -- every face onto a refined
   * neighbour is replaced by the chords of the children on that face, after the check that no coarse edge on it is
   * crossed twice at the finer spacing, and the interface is closed again, which also recomputes the moments. A face
   * onto a filled cell (a_filled) is closed instead, through CutCellBody::snapFace, and the footprint it gives up
   * becomes interface lying in the plane of that face. Every face so changed is recorded, in the order it was changed,
   * when overrides are given. A body that will not close, a face that cannot be restricted or closed, or a doubly
   * crossed edge stops the run.
   * @param[out]    a_body        The body.
   * @param[in,out] a_surface     The cell's surface; its refined edges' crossings are rewritten.
   * @param[in,out] a_overrides   Overrides to record the changed faces in, with the cell begun; null records nothing.
   * @param[in]     a_children    The cell's eight children in corner order, which must be given when any edge is
   * refined and describe at least the cell's refined faces and edges; null otherwise.
   * @param[in]     a_edgeRefined Which edges are refined, from edgesRefined.
   * @param[in]     a_filled      One where a cell around the cell was filled, from isFilled.
   * @param[in]     a_refined     Refined mask around the cell.
   * @param[in]     a_fineFilled  One where a cell of the finer level was filled, over the finer cells across the
   * cell's faces; a quadrant of a restricted face whose finer cell was filled is closed. Null where nothing is known of
   * the finer level's fills.
   * @param[in]     a_domain      Domain box of the cell's level.
   * @param[in]     a_cell        The cell.
   * @return Number of quadrants closed onto filled finer cells.
   */
  static int
  stitchBody(PolyhedralEB::CutCellBody&          a_body,
             PolyhedralEB::CutCellSurface&       a_surface,
             PolyhedralEB::CutCellFaceOverrides* a_overrides,
             const PolyhedralEB::CutCellSurface* a_children,
             const bool*                         a_edgeRefined,
             const BaseFab<signed char>&         a_filled,
             const BaseFab<signed char>&         a_refined,
             const BaseFab<signed char>*         a_fineFilled,
             const Box&                          a_domain,
             const IntVect&                      a_cell);

  /**
   * @brief Stitch a cut cell's body with its children reconstructed from the implicit function.
   * @details The reference the stitched graph is checked against under sanity_check: the same stitching as the stored
   * faces were made with, from children built at the finer spacing from the function rather than from the finer
   * level's stored surfaces. Nothing else calls the function to build a body.
   * @param[out] a_body    The body.
   * @param[in]  a_surface The cell's surface as the graph first built it.
   * @param[in]  a_graph   Graph of the cell's level.
   * @param[in]  a_filled  One where a cell around the cell was filled, from isFilled.
   * @param[in]  a_refined Refined mask around the cell.
   * @param[in]  a_cell    The cell.
   */
  void
  defineBodyFromFunction(PolyhedralEB::CutCellBody&          a_body,
                         const PolyhedralEB::CutCellSurface& a_surface,
                         const PolyhedralEBGraph&            a_graph,
                         const BaseFab<signed char>&         a_filled,
                         const BaseFab<signed char>&         a_refined,
                         const IntVect&                      a_cell) const;

  /**
   * @brief Build the body of a cut cell of a stitched graph from its stored surface and face overrides.
   * @details Defines the body from the surface, replaces every restricted face with its stored polygons, closes the
   * interface over them, and closes every closed face, in the order they were recorded -- which rebuilds the body
   * stitchBody made, bit for bit, from the cell's own data alone.
   * @param[out] a_body      The body.
   * @param[in]  a_surfaces  Surfaces of the cell's box, from the graph.
   * @param[in]  a_overrides Face overrides of the cell's box, from the graph.
   * @param[in]  a_cell      The cell.
   */
  static void
  defineBody(PolyhedralEB::CutCellBody&                  a_body,
             const IVSFAB<PolyhedralEB::CutCellSurface>& a_surfaces,
             const PolyhedralEB::CutCellFaceOverrides&   a_overrides,
             const IntVect&                              a_cell);

  /**
   * @brief Read the nodes and edges of a region of the finer level off its stored cells.
   * @details A node or edge belongs to every cell around it, so it is written from each cell of the region that holds a
   * surface, and two cells that disagree about one stop the run. A node that no such cell touches lies among cells
   * whose corners agree, and takes their sign: minus one where they are regular, one where they are covered. An edge
   * that no such cell touches carries no crossing. Cells the finer level does not hold write nothing.
   * @param[out] a_nodes       Node values over the nodes of a_region.
   * @param[out] a_crossings   Crossings over the edges of a_region, addressed as defineIntercepts addresses them.
   * @param[out] a_fluidJoined Face decisions of each cell of a_region that holds a surface, minus one elsewhere.
   * @param[in]  a_region      Cells of the finer level.
   * @param[in]  a_states      Cell states of the finer level over the region, s_absent where it holds none.
   * @param[in]  a_mask        One where a cell holds a surface.
   * @param[in]  a_surfaces    The surfaces.
   */
  static void
  readFineCells(BaseFab<Real>&                              a_nodes,
                BaseFab<Real>                               a_crossings[SpaceDim],
                BaseFab<int>&                               a_fluidJoined,
                const Box&                                  a_region,
                const BaseFab<signed char>&                 a_states,
                const BaseFab<signed char>&                 a_mask,
                const IVSFAB<PolyhedralEB::CutCellSurface>& a_surfaces);

  /**
   * @brief One cell's surface, read off node values, edge crossings and face decisions that are already known.
   * @details A face's decision is taken from the cell itself where it holds a surface, else from the cell across the
   * face, else from the bilinear saddle of the face's corners.
   * @param[out] a_surface     The surface.
   * @param[in]  a_nodes       Node values.
   * @param[in]  a_crossings   Edge crossings, addressed as defineIntercepts addresses them.
   * @param[in]  a_fluidJoined Face decisions per cell.
   * @param[in]  a_cell        The cell.
   */
  static void
  assembleSurface(PolyhedralEB::CutCellSurface& a_surface,
                  const BaseFab<Real>&          a_nodes,
                  const BaseFab<Real>           a_crossings[SpaceDim],
                  const BaseFab<int>&           a_fluidJoined,
                  const IntVect&                a_cell);

  /**
   * @brief Copy a finer level's cells onto a coarser level's layout refined, with a ring around every box.
   * @details The cell states, a mask of the cells that hold a surface, the sets of those cells and their surfaces --
   * everything a cell of the finer level is rebuilt from, apart from faces it restricts to a level above it. Cells
   * the finer level does not hold read s_absent. Collective.
   * @param[out] a_states   Cell states.
   * @param[out] a_mask     One where a cell holds a surface.
   * @param[out] a_sets     The cells that hold a surface, per box.
   * @param[out] a_surfaces Their surfaces.
   * @param[in]  a_fine     The finer level.
   * @param[in]  a_grids    The coarser level's layout.
   * @param[in]  a_ring     Width of the ring around every box, in cells of the coarser level.
   */
  static void
  copyFinerCells(LevelData<BaseFab<signed char>>&                 a_states,
                 LevelData<BaseFab<signed char>>&                 a_mask,
                 LayoutData<IntVectSet>&                          a_sets,
                 LevelData<IVSFAB<PolyhedralEB::CutCellSurface>>& a_surfaces,
                 const PolyhedralEBGraph&                         a_fine,
                 const DisjointBoxLayout&                         a_grids,
                 const int                                        a_ring);

  /**
   * @brief Stitch one level: take its refined edges' crossings and its seam faces from the finer level's stored cells,
   * close its faces onto filled cells, write the crossings into its surfaces and the faces into its overrides.
   * @details The finer level's cell states and surfaces over the refined coarse boxes and a ring around them are
   * copied onto this level's layout refined, and every cut cell of every box and its ghost ring is stitched from them.
   * A quadrant of a seam face across which the finer level filled a cell is closed. A regular or covered cell beside
   * the finer level is made cut where the finer level describes the face between them otherwise than the cell's
   * corners do -- the node at the middle of the face on the other side, or a filled cell across a quadrant of a regular
   * one -- and stitched like any other. A cell pierced at the middle of such a face should have been refined by
   * ComputationalGeometry, so each one is listed in pout and warned about. A node at the middle of such a face's edge
   * that disagrees with the edge's ends is an edge the finer level crosses twice, which is counted and reported rather
   * than acted on. Every stitched surface carries its role. Under sanity_check each body is also built from the stored
   * data, and from the implicit function where its children can say the same, and they must agree bit for bit.
   * Collective. A face onto a cell the coarser level filled is closed like a face onto a cell this level filled.
   * @param[in,out] a_coarse  The level.
   * @param[in]     a_fine    The finer level, or null where there is none.
   * @param[in]     a_coarser The coarser level, or null where there is none.
   */
  void
  stitchLevel(PolyhedralEBGraph& a_coarse, const PolyhedralEBGraph* a_fine, const PolyhedralEBGraph* a_coarser) const;

  /**
   * @brief Write the surface, the boxes it was built over and the boxes that split, as one HDF5 file per phase.
   * @details One group per level, holding the level's triangles, or segments in two dimensions, as an indexed mesh --
   * the vertices once each and s_elementVertices indices per element -- rather than as a soup of coordinates, which is
   * what a viewer wants and what an XDMF description can point at. The boxes of the level and the boxes that split on
   * it go in the same group, so that the surface can be read against the grids it was built on without a second file.
   * The distance within which two positions are one vertex is written alongside them. An XDMF description of the file
   * is written next to it. Collective: every rank writes its own stretch of every dataset, so a vertex on a boundary
   * between two ranks is written by each, since vertices are identified within a rank and not across them.
   * @param[in] a_fileName File to write, without a suffix.
   * @param[in] a_facets   This rank's elements of each level, three reals per vertex.
   */
  void
  writeSurface(const std::string& a_fileName, const Vector<Vector<Real>>& a_facets) const;

  /**
   * @brief Write the cut cells' connectivity as a skeleton a viewer can draw: a point per cut cell and per open face of
   * one, and a line from each cut cell to each of its open faces.
   * @details One group per level. A cell node is written for every cut cell of the level's tiles, at its fluid's
   * centroid, with its volume fraction, the role of its surface, its number of sheets and its rank. A face node is
   * written for every open face of one, at the centroid of the face's open part, with its aperture, its face state
   * and whether the stitch rewrote it; a face two cut cells of the level share is written once, by the cell on its low
   * side, on whichever rank holds it. The two kinds share one point mesh, told apart by a field
   * holding s_cellNode or s_faceNode, and a field reads s_noValue on the kind it does not belong to. An edge is written
   * from a cell node to the centroid of each of its open faces, so that the edges of two neighbouring cut cells meet in
   * the face between them, and a face onto a regular cell, another level or the domain boundary ends in a half; it
   * carries the face's fields too. Two XDMF descriptions are written next to the file, a_fileName followed by
   * _nodes.xmf with the nodes as points and a_fileName followed by _edges.xmf with the edges as segments, since a
   * reader may take only the first grid of a file.
   * @param[in] a_fileName File to write, without a suffix.
   */
  void
  writeGraph(const std::string& a_fileName) const;

  /**
   * @brief Identify the vertices of a triangle soup, giving each one once and each triangle three indices.
   * @details Two positions closer than the tolerance are one vertex. Candidates are taken from the lattice cell of
   * the tolerance a position falls in and the twenty-six around it, which is what keeps this out of the quadratic
   * scan a surface meeting the mesh squarely would otherwise provoke, with every candidate still measured against
   * so that two positions either side of a lattice wall are not told apart.
   * @param[out] a_vertices     The vertices, three reals each.
   * @param[out] a_connectivity Three vertex indices per triangle.
   * @param[in]  a_facets       Triangles, nine reals each.
   * @param[in]  a_tolerance    Distance within which two positions are the same vertex.
   */
  static void
  indexFacets(std::vector<Real>&  a_vertices,
              std::vector<int>&   a_connectivity,
              const Vector<Real>& a_facets,
              const Real          a_tolerance);

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

private:
  /**
   * @brief An implicit function that stops the run if it is ever evaluated.
   * @details Put in place of the shop's own function while a test copies the graph and rebuilds bodies from the copy,
   * which must be done from the graph's data alone.
   */
  class AbortingIF : public BaseIF
  {
  public:
    /**
     * @brief Stop the run.
     * @param[in] a_point Where the function was asked for.
     * @return Nothing; it does not return.
     */
    Real
    value(const RealVect& a_point) const override;

    /**
     * @brief A copy.
     * @return A new aborting function.
     */
    BaseIF*
    newImplicitFunction() const override;
  };

#ifdef CH_USE_HDF5
  /**
   * @brief Write one attribute, a scalar or a short array.
   * @param[in] a_where Group or file to attach it to.
   * @param[in] a_name  Name of the attribute.
   * @param[in] a_type  HDF5 type of its elements.
   * @param[in] a_count Number of elements, one for a scalar.
   * @param[in] a_data  The elements.
   */
  static void
  writeH5Attribute(const hid_t        a_where,
                   const std::string& a_name,
                   const hid_t        a_type,
                   const hsize_t      a_count,
                   const void*        a_data);

  /**
   * @brief Write one dataset, sized by what every rank holds together, with this rank's rows in its own stretch.
   * @details A rank with nothing to write selects nothing and takes part in the call all the same. It describes that
   * with the null dataspace: a simple dataspace of zero extent is not one, the call fails, and a rank that fails on its
   * way into a collective write leaves every other rank waiting in it.
   * @param[in] a_where    Group or file to put it in.
   * @param[in] a_name     Name of the dataset.
   * @param[in] a_type     HDF5 type of its elements.
   * @param[in] a_rows     Rows this rank writes.
   * @param[in] a_offset   Row this rank's stretch starts at.
   * @param[in] a_total    Rows over all ranks; nothing is written when there are none.
   * @param[in] a_columns  Columns per row; one makes the dataset one-dimensional.
   * @param[in] a_data     This rank's rows.
   * @param[in] a_transfer Transfer property list.
   */
  static void
  writeH5Slab(const hid_t        a_where,
              const std::string& a_name,
              const hid_t        a_type,
              const hsize_t      a_rows,
              const hsize_t      a_offset,
              const hsize_t      a_total,
              const hsize_t      a_columns,
              const void*        a_data,
              const hid_t        a_transfer);
#endif

  /**
   * @brief Append the nodes of a cell's VoFs as an index space holds them, so that a graph can be built from them
   * again.
   * @details Each VoF's arcs come from its faces -- the VoF across, or the domain boundary -- with their apertures and
   * centroids. Where the moments hold no data for the cell, the moments are left at zero: a graph is built from arcs
   * alone.
   * @param[in,out] a_nodes The nodes, appended to.
   * @param[in]     a_graph Graph holding the cell.
   * @param[in]     a_data  Moments over some part of the graph's region.
   * @param[in]     a_cell  The cell, irregular in the graph.
   * @param[in]     a_withMoments Whether the moments hold data for the cell.
   */
  static void
  appendNodes(Vector<IrregNode>& a_nodes,
              const EBGraph&     a_graph,
              const EBData&      a_data,
              const IntVect&     a_cell,
              const bool         a_withMoments);

  /**
   * @brief Make each regular cell beside a multi-valued one a cut cell of volume one.
   * @details Chombo's coarsening reads the face of a regular cell as leading to the first VoF of the cell across, so a
   * regular cell's face into any other VoF of a multi-valued neighbour is lost. EBISLevel::fixRegularNextToMultiValued
   * makes such cells full cut cells and divides the face evenly among the VoFs across; here each face into a VoF across
   * is the face that VoF holds toward the cell, with its aperture and centroid. The nodes of a cell must be contiguous,
   * and a cell at the region's edge sees only the neighbours inside it.
   * @param[in,out] a_flags  Regular (1), cut (0) or covered (-1) for the region's cells. The cells made cut are set to
   * 0.
   * @param[in,out] a_nodes  Nodes of the region's cut cells. A node is appended for each cell made cut.
   * @param[in]     a_region Region the flags describe.
   * @param[in]     a_valid  Part of the region whose cells made cut are counted.
   * @param[in]     a_domain Domain the region lies in.
   * @return The number of cells in a_valid made cut.
   */
  static long long
  cutRegularNextToMultiValued(BaseFab<int>&        a_flags,
                              Vector<IrregNode>&   a_nodes,
                              const Box&           a_region,
                              const Box&           a_valid,
                              const ProblemDomain& a_domain);

  /**
   * @brief Write the index space of each cell of a region as numbers, summed over the cell's VoFs: the number of VoFs,
   * the volume fraction, the first moment of volume and the aperture of every face.
   * @details Component 0 is the number of VoFs, 1 the volume fraction, 2 to 1 + SpaceDim the first moment about the
   * cell centre in units of the cell, and 2 + SpaceDim + 2 * dir + side the aperture on that side. A regular cell reads
   * as a single VoF of volume one with every face open, so a full cut cell reads the same.
   * @param[in,out] a_signature At least s_signatureComps components over the region.
   * @param[in]     a_graph     The graph.
   * @param[in]     a_data      The moments, over the region.
   * @param[in]     a_region    The cells to write.
   */
  static void
  cellSignature(FArrayBox& a_signature, const EBGraph& a_graph, const EBData& a_data, const Box& a_region);
};

#include <CD_NamespaceFooter.H>

#endif
