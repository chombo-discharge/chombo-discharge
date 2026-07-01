/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Realm.H
  @brief  Declaration of the Realm class. 
  @author Robert Marskar
*/

#ifndef CD_REALM_H
#define CD_REALM_H

// Std includes
#include <map>
#include <set>

// Our includes
#include <CD_PhaseRealm.H>
#include <CD_EBAMRParticleMesh.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_MFLevelGrid.H>
#include <CD_LevelTiles.H>
#include <CD_ParticleGhostMask.H>
#include <CD_PetscGrid.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Alias for cutting down on the typic of booleans defined over AMR grids.
*/
using AMRMask = Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>;

/**
  @brief Class for holding multifluid Realms, which are the Realms that we actually use. 
*/
class Realm
{
public:
  /**
    @brief Identifier for perimal realm
  */
  const static std::string Primal;

  /**
    @brief Identifier for perimal realm
  */
  const static std::string primal;

  /**
    @brief Default constructor. Must subsequently call define
  */
  Realm();

  /**
    @brief Disallowed copy ctor
    @param[in] a_other Other realm
  */
  Realm(const Realm& a_other) = delete;

  /**
    @brief Disallowed move ctor
    @param[in] a_other Other realm
  */
  Realm(const Realm&& a_other) = delete;

  /**
    @brief Disallowed copy assignment operator
    @param[in] a_other Other realm
  */
  Realm&
  operator=(const Realm& a_other) = delete;

  /**
    @brief Disallowed move assignment operator
    @param[in] a_other Other realm
  */
  Realm&
  operator=(const Realm&& a_other) = delete;

  /**
    @brief Destructor
  */
  ~Realm();

  /**
    @brief Full define function
    @param[in] a_grids            Grids on each level, i.e. boxes and processor distributions
    @param[in] a_domains          Problem domains
    @param[in] a_refRat           Refinement ratios
    @param[in] a_dx               Grid resolutions
    @param[in] a_probLo           Lower-left corner of physical domain
    @param[in] a_finestLevel      Finest grid level
    @param[in] a_minBlockSize   Grid blocking factor
    @param[in] a_ebGhost          Number of ghost cells used for geometric information when making operators
    @param[in] a_numGhost         Number of ghost cells in data holders
    @param[in] a_lsfGhost         Number of ghost cells used when computing the level set on the mesh
    @param[in] a_redistRad        Redistribution radius
    @param[in] a_mgInterpOrder    Multigrid ghost cell interpolation order
    @param[in] a_mgInterpRadius   Multigrid ghost cell interpolation radius
    @param[in] a_mgInterpWeight   Multigrid ghost cell interpolation weight (for least squares)
    @param[in] a_centroidInterpType  Stencil type for interpolation from cell centers to cell centroids
    @param[in] a_ebInterpType     Stencil type for interpolation from cell centers to eb centroids
    @param[in] a_baseif           Implicit/signed distance function 
    @param[in] a_mfis             EBIndexSpace discrete information
  */
  void
  define(const Vector<DisjointBoxLayout>&                           a_grids,
         const Vector<ProblemDomain>&                               a_domains,
         const Vector<int>&                                         a_refRat,
         const Vector<Real>&                                        a_dx,
         const RealVect&                                            a_probLo,
         const int                                                  a_finestLevel,
         const int                                                  a_minBlockSize,
         const int                                                  a_ebGhost,
         const int                                                  a_numGhost,
         const int                                                  a_lsfGhost,
         const int                                                  a_redistRad,
         const int                                                  a_mgInterpOrder,
         const int                                                  a_mgInterpRadius,
         const int                                                  a_mgInterpWeight,
         const CellCentroidInterpolation::Type                      a_centroidInterpType,
         const EBCentroidInterpolation::Type                        a_ebInterpType,
         const std::map<phase::which_phase, RefCountedPtr<BaseIF>>& a_baseif,
         const RefCountedPtr<MultiFluidIndexSpace>&                 a_mfis);

  /**
    @brief Set grid method
    @param[in] a_grids       New grids
    @param[in] a_finestLevel New finest grid level
    @details This sets m_grids to be a_grids, but does not perform a regrid of remaining members.
  */
  void
  setGrids(const Vector<DisjointBoxLayout>& a_grids, const int a_finestLevel);

  /**
    @brief Perform pre-regrid operations
  */
  void
  preRegrid();

  /**
    @brief Regrid method for EBAMR base
    @param[in] a_lmin Coarsest grid level that changed. 
    @details This regrids EBLevelGrid and iterators
  */
  void
  regridBase(const int a_lmin);

  /**
    @brief Regrid method for EBAMR operators.
    @param[in] a_lmin Coarsest grid level that changed. 
    @details This regrids all operators. 
  */
  void
  regridOperators(const int a_lmin);

  /**
    @brief Register an AMR operator
    @param[in] a_operator Operator name
    @param[in] a_phase Phase where operator is registered
    @details Issues run-time error if a_operator can't be recognized
  */
  void
  registerOperator(const std::string& a_operator, const phase::which_phase a_phase);

  /**
    @brief Query if an AMR operator has been registered
    @param[in] a_operator operator name
    @param[in] a_phase Phase where operator is registered
    @return True if operator has been registered and false otherwise
  */
  bool
  queryOperator(const std::string& a_operator, const phase::which_phase a_phase) const;

  /**
    @brief Register a mask
    @param[in] a_mask Mask name
    @param[in] a_buffer Buffer
  */
  void
  registerMask(const std::string& a_mask, const int a_buffer);

  /**
    @brief Query if a mask has been registered
    @param[in] a_mask Mask name
    @return true if mask has been registered
    @param[in] a_buffer Buffer
  */
  bool
  queryMask(const std::string& a_mask, const int a_buffer) const;

  /**
    @brief Get refinement ratios
    @return m_refinementRatios
  */
  const Vector<int>&
  getRefinementRatios() const;

  /**
    @brief Get grid resolutions
    @return m_dx
  */
  const Vector<Real>&
  getDx() const;

  /**
    @brief Get AMR grids
    @return m_grids
  */
  const Vector<DisjointBoxLayout>&
  getGrids() const;

  /**
    @brief Get problem domains
    @return Domains
  */
  const Vector<ProblemDomain>&
  getDomains() const;

  /**
    @brief Get MFLevelGrids
    @return M f level grid
  */
  Vector<RefCountedPtr<MFLevelGrid>>&
  getMFLevelGrid();

  /**
    @brief Get the coarsened grids
    @details The returned vector contains the coarsened grids. Entry on level lvl=0 contains the coarsening of level lvl=1 and so on. The top
    level is undefined because there was no finer grid from which we can coarsen.
    @return M f level grid co fi
  */
  Vector<RefCountedPtr<MFLevelGrid>>&
  getMFLevelGridCoFi();

  /**
    @brief Get the refined grids
    @details The returned vector contains the refined grids. Entry on level lvl=1 contains the refinement of level lvl=0 and so on. The bottom
    level is undefined because there was no coarser grid from which we can refine.
    @return M f level grid fi co
  */
  Vector<RefCountedPtr<MFLevelGrid>>&
  getMFLevelGridFiCo();

  /**
    @brief Get EBIndexSpace for particular phase
    @param[in] a_phase Phase
    @return E b index space
  */
  const RefCountedPtr<EBIndexSpace>&
  getEBIndexSpace(const phase::which_phase a_phase) const;

  /**
    @brief Get EBISLayout for particular phase
    @param[in] a_phase Phase
    @return E b i s layout
  */
  const Vector<EBISLayout>&
  getEBISLayout(const phase::which_phase a_phase) const;

  /**
    @brief Get EBLevelGrid for particular phase
    @param[in] a_phase Phase
    @return E b level grid
  */
  const Vector<RefCountedPtr<EBLevelGrid>>&
  getEBLevelGrid(const phase::which_phase a_phase) const;

  /**
    @brief Get coarsened fine EBLevelGrid for particular phase
    @param[in] a_phase Phase
    @return E b level grid co fi
  */
  const Vector<RefCountedPtr<EBLevelGrid>>&
  getEBLevelGridCoFi(const phase::which_phase a_phase) const;

  /**
    @brief Get VoFIterators covering all cut-cells (singly- and multiply-cut) for a particular phase.
    @details For a pass covering only multiply-cut cells, use getMultiCutVofIterator instead.
    @param[in] a_phase Phase (gas or solid)
    @return Per-level LayoutData of VoFIterators covering all irregular cells
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
  getVofIterator(const phase::which_phase a_phase) const;

  /**
    @brief Get VoFIterators restricted to multiply-cut cells for a particular phase
    @param[in] a_phase Phase
    @return Multi-cut vof iterator
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
  getMultiCutVofIterator(const phase::which_phase a_phase) const;

  /**
    @brief Get face iterators with FaceStop::SurroundingWithBoundary over cut-cell faces in each valid grid box for a particular phase.
    @details Includes domain-boundary faces. For cell-to-face averaging that excludes domain boundary
    faces, use getFaceIteratorNoBoundary. For faces that extend into the tangential ghost layer, use
    getFaceIteratorWithTangentialGhosts.
    @param[in] a_phase Phase (gas or solid)
    @return Per-level LayoutData of FaceIterator arrays (one per coordinate direction)
  */
  Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
  getFaceIterator(const phase::which_phase a_phase) const;

  /**
    @brief Get face iterators with FaceStop::SurroundingNoBoundary over cut-cell faces in each valid grid box for a particular phase.
    @details Domain boundary faces are excluded. Intended for cell-to-face averaging where boundary
    faces are fixed up separately. For iterators that include domain boundary faces, use getFaceIterator.
    @param[in] a_phase Phase (gas or solid)
    @return Per-level LayoutData of FaceIterator arrays (one per coordinate direction)
  */
  Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
  getFaceIteratorNoBoundary(const phase::which_phase a_phase) const;

  /**
    @brief Get per-direction face iterators covering cut-cell faces that reach into the tangential ghost layer.
    @details For the iterator in direction dir, the underlying IVS covers irregular cells in a box grown by 1
    in the two directions tangential to dir; the dir extent stays at the valid box. Domain boundary faces are
    excluded (FaceStop::SurroundingNoBoundary). Callers require at least 1 ghost layer in cell data.
    @param[in] a_phase Phase
    @return Tangential-ghost face iterator
  */
  Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
  getFaceIteratorWithTangentialGhosts(const phase::which_phase a_phase) const;

  /**
    @brief Get face iterators restricted to multiply-cut faces in each valid grid box.
    @details Only faces abutting multiply-cut cells (cells with more than one VoF) are visited.
    Faces are enumerated with FaceStop::SurroundingNoBoundary.
    @param[in] a_phase Phase
    @return Multi-cut face iterator
  */
  Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>&
  getMultiCutFaceIterator(const phase::which_phase a_phase) const;

  /**
    @brief Get objects for computing "non-conservative divergences"
    @param[in] a_phase Phase
    @return Non conservative divergence
  */
  const Vector<RefCountedPtr<EBNonConservativeDivergence>>&
  getNonConservativeDivergence(const phase::which_phase a_phase) const;

  /**
    @brief Get objects for computing interpolation of cell-centered data to cell centroids.
    @param[in] a_phase Phase
    @return Cell centroid interpolation
  */
  const Vector<RefCountedPtr<CellCentroidInterpolation>>&
  getCellCentroidInterpolation(const phase::which_phase a_phase) const;

  /**
    @brief Get objects for computing interpolation of cell-centered data to EB centroids
    @param[in] a_phase Phase
    @return E b centroid interpolation
  */
  const Vector<RefCountedPtr<EBCentroidInterpolation>>&
  getEBCentroidInterpolation(const phase::which_phase a_phase) const;

  /**
    @brief Get data coarsening utility
    @param[in] a_phase Phase
    @return Coarse average
  */
  Vector<RefCountedPtr<EBCoarAve>>&
  getCoarseAverage(const phase::which_phase a_phase);

  /**
    @brief Get the multigrid interpolator
    @param[in] a_phase Phase
    @return Multigrid interpolator
  */
  Vector<RefCountedPtr<EBMultigridInterpolator>>&
  getMultigridInterpolator(const phase::which_phase a_phase);

  /**
    @brief Get the particle mesh operator
    @param[in] a_phase Which phase
    @return Particle mesh
  */
  EBAMRParticleMesh&
  getParticleMesh(const phase::which_phase a_phase);

  /**
    @brief Get the particle mesh operator
    @param[in] a_phase Which phase
    @return Surface deposition
  */
  EBAMRSurfaceDeposition&
  getSurfaceDeposition(const phase::which_phase a_phase);

  /**
    @brief Get gradient operator
    @param[in] a_phase Phase
    @return Gradient op
  */
  const Vector<RefCountedPtr<EBGradient>>&
  getGradientOp(const phase::which_phase a_phase) const;

  /**
    @brief Get ghost cell interpolator
    @param[in] a_phase Phase
    @return Ghost cell interpolator
  */
  Vector<RefCountedPtr<EBGhostCellInterpolator>>&
  getGhostCellInterpolator(const phase::which_phase a_phase);

  /**
    @brief Get piecewise linear ghost cell interpolation utility
    @param[in] a_phase Phase
    @return Fine interp
  */
  Vector<RefCountedPtr<EBCoarseToFineInterp>>&
  getFineInterp(const phase::which_phase a_phase);

  /**
    @brief Get flux register utility
    @param[in] a_phase Phase
    @return Flux register
  */
  Vector<RefCountedPtr<EBReflux>>&
  getFluxRegister(const phase::which_phase a_phase);

  /**
    @brief Get the redistribution operators
    @param[in] a_phase Phase (gas or solid)
    @return Redistribution op
  */
  Vector<RefCountedPtr<EBFluxRedistribution>>&
  getRedistributionOp(const phase::which_phase a_phase);

  /**
    @brief Get level-set function
    @param[in] a_phase Phase
    @return Levelset
  */
  const EBAMRFAB&
  getLevelset(const phase::which_phase a_phase) const;

  /**
    @brief Get the regular-cell mask (value 1 in regular cells, 0 elsewhere) for a phase.
    @param[in] a_phase Phase
    @return Regular-cell mask
  */
  const EBAMRCellData&
  getRegularCells(const phase::which_phase a_phase) const;

  /**
    @brief Get the covered-cell mask (value 1 in covered cells, 0 elsewhere) for a phase.
    @param[in] a_phase Phase
    @return Covered-cell mask
  */
  const EBAMRCellData&
  getCoveredCells(const phase::which_phase a_phase) const;

  /**
    @brief Get the non-covered-cell mask (1 in regular/irregular cells, 0 in covered cells) for a phase.
    @param[in] a_phase Phase
    @return Non-covered-cell mask
  */
  const EBAMRCellData&
  getNotCoveredCells(const phase::which_phase a_phase) const;

  /**
    @brief Get the irregular-cell mask (value 1 in irregular cells, 0 elsewhere) for a phase.
    @param[in] a_phase Phase
    @return Irregular-cell mask
  */
  const EBAMRCellData&
  getIrregularCells(const phase::which_phase a_phase) const;

  /**
    @brief Get AMR mask
    @return Mask
    @param[in] a_mask Mask
    @param[in] a_buffer Buffer
  */
  const AMRMask&
  getMask(const std::string& a_mask, const int a_buffer) const;

  /**
    @brief Get all valid cells
    @return Valid cells
  */
  const AMRMask&
  getValidCells() const;

  /**
    @brief Register a particle ghost-mask width so it is (re)built at every regrid.
    @details Downstream code must register a width before it can fetch the corresponding masks. Widths
    are deduplicated, so several consumers may safely request the same width. No masks are built at all
    unless at least one width has been registered.
    @param[in] a_width Ghost width (cells). Must be > 0; a run-time error is thrown otherwise.
  */
  void
  registerParticleGhostMask(const int a_width) noexcept;

  /**
    @brief Get the same-level particle ghost mask for a registered width.
    @param[in] a_width Registered ghost width.
    @return The per-level same-level mask for a_width. Aborts if the width was not registered.
  */
  const AMRParticleGhostMask&
  getParticleGhostMask(const int a_width) const noexcept;

  /**
    @brief Get the fine-to-coarse particle ghost mask for a registered width.
    @param[in] a_width Registered ghost width.
    @return The per-level fine-to-coarse mask for a_width. Aborts if the width was not registered.
  */
  const AMRParticleGhostMask&
  getParticleGhostMaskFineToCoar(const int a_width) const noexcept;

  /**
    @brief Get the coarse-to-fine particle ghost mask for a registered width.
    @param[in] a_width Registered ghost width.
    @return The per-level coarse-to-fine mask for a_width. Aborts if the width was not registered.
  */
  const AMRParticleGhostMask&
  getParticleGhostMaskCoarToFine(const int a_width) const noexcept;

  /**
    @brief Get a trivial (empty) particle ghost mask that scatters no particles.
    @details Pass this to ParticleContainer::fillGhost for a direction that should not scatter. It holds
    no levels, so a fill over it is a no-op.
    @return A static empty mask.
  */
  const AMRParticleGhostMask&
  getTrivialParticleGhostMask() const noexcept;

  /**
    @brief Get the tiled space
    @return Level tiles
  */
  const Vector<RefCountedPtr<LevelTiles>>&
  getLevelTiles() const noexcept;

  /**
    @brief Result of a point->block query. See getLevelAndBox.
  */
  struct LevelAndBox
  {
    /** @brief AMR level containing the point, or -1 if not found. */
    int level;

    /** @brief Global box/grid index within that level's LevelTiles. */
    unsigned int boxIndex;

    /** @brief MPI rank owning that box. */
    int rank;

    /** @brief True if a covering box was found. */
    bool valid;
  };

  /**
    @brief Find the (finest) grid box covering a physical position.
    @details O(1) point->block query implemented on the Realm's LevelTiles: the position is mapped to a
    min_block_size tile by integer division, then looked up in the per-level tile->box hash, finest level
    first. Works for variable-sized/anisotropic boxes because every box is registered under all the min
    tiles it covers. This is the same mapping ParticleContainer uses to bin particles, exposed for any
    user that needs to know which box owns a point.
    @param[in] a_pos Physical position.
    @return {level, boxIndex, rank, valid}. valid is false (and the other fields are unspecified) when no
    box covers the point.
  */
  LevelAndBox
  getLevelAndBox(const RealVect& a_pos) const noexcept;

#ifdef CH_USE_PETSC
  /**
    @brief Get the PETSc grid object
  */
  const RefCountedPtr<PetscGrid>&
  getPetscGrid() const noexcept;
#endif

protected:
  /**
    @brief Realm defined or not
  */
  bool m_isDefined;

  /**
    @brief Realm verbosity. Used for debugging.
  */
  int m_verbosity;

  /**
    @brief Finest AMR level
  */
  int m_finestLevel;

  /**
    @brief Blocking factor
  */
  int m_minBlockSize;

  /**
    @brief Number of ghost cells
  */
  int m_numGhost;

  /**
    @brief Physical coordinates of lower-left corner in simulation domain
  */
  RealVect m_probLo;

  /**
    @brief Grid resolutions
  */
  Vector<Real> m_dx;

  /**
    @brief Refinement ratios
  */
  Vector<int> m_refinementRatios;

  /**
    @brief AMR grids
  */
  Vector<DisjointBoxLayout> m_grids;

  /**
    @brief Domains
  */
  Vector<ProblemDomain> m_domains;

  /**
    @brief MFLevelGrids
  */
  Vector<RefCountedPtr<MFLevelGrid>> m_mflg;

  /**
    @brief Coarsened grids. Entry lvl contains the coarsened version of grids on level lvl+1
  */
  Vector<RefCountedPtr<MFLevelGrid>> m_mflgCoFi;

  /**
    @brief Refined grids. Entry lvl+1 contains the refined version of grids on level lvl
  */
  Vector<RefCountedPtr<MFLevelGrid>> m_mflgFiCo;

  /**
    @brief Level "tiles" -- i.e. the grids viewed as constant-size boxes.
  */
  Vector<RefCountedPtr<LevelTiles>> m_levelTiles;

  /**
    @brief Index spaces
  */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /**
    @brief Mask which is true for all valid cells
  */
  AMRMask m_validCells;

  /**
    @brief Registered particle ghost-mask widths (>0). Persist across regrids; the masks below are
    (re)built for exactly these widths in defineParticleGhostMasks().
  */
  std::set<int> m_particleGhostMaskWidths;

  /**
    @brief Same-level particle ghost masks, keyed by ghost width. Per width: indexed by source level, and
    per source box maps a source cell to the same-level boxes whose ghosted region covers it.
  */
  std::map<int, AMRParticleGhostMask> m_particleGhostMask;

  /**
    @brief Fine-to-coarse particle ghost masks, keyed by ghost width. Per width: indexed by source level
    (>=1); per source (fine) box maps a source cell to the coarser-level boxes whose ghosted region covers it.
  */
  std::map<int, AMRParticleGhostMask> m_particleGhostMaskFineToCoar;

  /**
    @brief Coarse-to-fine particle ghost masks, keyed by ghost width. Per width: indexed by source level
    (<finest); per source (coarse) box maps a source cell to the finer-level boxes whose ghosted region covers it.
  */
  std::map<int, AMRParticleGhostMask> m_particleGhostMaskCoarToFine;

#ifdef CH_USE_PETSC
  /**
    @brief PETSc grid
  */
  RefCountedPtr<PetscGrid> m_petscGrid;
#endif

  /**
    @brief Phase realms. 
    @details Mutable because some of the operators are mutable. 
  */
  mutable std::map<phase::which_phase, RefCountedPtr<PhaseRealm>> m_realms;

  /**
    @brief Implicit/signed distance functions on each phase
  */
  std::map<phase::which_phase, RefCountedPtr<BaseIF>> m_baseif;

  /**
    @brief Masks
  */
  std::map<std::pair<std::string, int>, AMRMask> m_masks;

  /**
    @brief Define MFLevelGrid
    @param[in] a_lmin Coarsest level that changed during regrid
  */
  void
  defineMFLevelGrid(const int a_lmin);

  /**
    @brief Define masks
    @param[in] a_lmin Coarsest level that changed during regrid
  */
  void
  defineMasks(const int a_lmin);

  /**
    @brief Define outer particle halo mask.
    @details This is a mask which is true on the OUTSIDE of a refinement boundary.
    @param[in] a_lmin Coarsest level that changed during regrid
  */
  void
  defineOuterHaloMask(const int a_lmin);

  /**
    @brief Define the outer particle halo masks. This is the level-version. 
    @param[out] a_coarMask   Mask
    @param[in]  a_domainCoar Coarse domain
    @param[in]  a_domainFine Fine domain
    @param[in]  a_gridsCoar  Coarse grids
    @param[in]  a_gridsFine  Fine grids
    @param[in]  a_buffer     Buffer (i.e. width of "halo")
    @param[in]  a_refRat     Refinement ratio between the two levels
  */
  void
  defineOuterHaloMask(LevelData<BaseFab<bool>>& a_coarMask,
                      const ProblemDomain&      a_domainCoar,
                      const ProblemDomain&      a_domainFine,
                      const DisjointBoxLayout&  a_gridsCoar,
                      const DisjointBoxLayout&  a_gridsFine,
                      const int                 a_buffer,
                      const int                 a_refRat);

  /**
    @brief Define inner particle halo mask.
    @details This is a mask which is true on the INSIDE of a refinement boundary. The buffer size is always
    the number of COARSE grid cells. I.e., if one uses a refinement factor of 2, the first two strips of fine-grid
    cells are set to true and the others to false. 
    @param[in] a_lmin Coarsest level that changed during regrid
  */
  void
  defineInnerHaloMask(const int a_lmin);

  /**
    @brief Define a mask which is true on the outside (i.e., coarse side) of the refinement boundary
    @param[in] a_lmin Description
  */
  void
  defineOuterCFMask(const int a_lmin);

  /**
    @brief Define the outer cf region mask. This is the level-version. 
    @param[out] a_coarMask   Mask
    @param[in]  a_domainCoar Coarse domain
    @param[in]  a_domainFine Fine domain
    @param[in]  a_gridsCoar  Coarse grids
    @param[in]  a_gridsFine  Fine grids
    @param[in]  a_buffer     Buffer
    @param[in]  a_refRat     Refinement ratio between the two levels
  */
  void
  defineOuterCFMask(LevelData<BaseFab<bool>>& a_coarMask,
                    const ProblemDomain&      a_domainCoar,
                    const ProblemDomain&      a_domainFine,
                    const DisjointBoxLayout&  a_gridsCoar,
                    const DisjointBoxLayout&  a_gridsFine,
                    const int                 a_buffer,
                    const int                 a_refRat);

  /**
    @brief Define a mask which is true on the inside (i.e., fine side) of the refinement boundary
    @param[in] a_lmin Description
  */
  void
  defineInnerCFMask(const int a_lmin);

  /**
    @brief Define the coarse-fine interface. 
    @details This is a mask which is true on the immediate outside of a refinement boundary.
    @param[in] a_lmin Description
  */
  void
  defineCFIVS(const int a_lmin);

  /**
    @brief Define the valid cell region mask. 
  */
  void
  defineValidCells();

  /**
    @brief Define the tiled space
  */
  void
  defineLevelTiles() noexcept;

  /**
    @brief Define the same-level, fine-to-coarse, and coarse-to-fine particle ghost masks on every level.
    @details Allocates and fills the three per-level masks (m_particleGhostMask,
    m_particleGhostMaskFineToCoar, m_particleGhostMaskCoarToFine) by delegating to
    defineParticleGhostMaskSameLevel and defineParticleGhostMaskCrossLevel. Each mask records, per source
    cell, the target boxes (grid index + receiving rank) that a particle in that cell must be scattered
    to in that direction. Called from regridBase after the valid-cell and level-tile structures exist.
    The ghost width is m_numGhost. Periodic boundaries are not handled.
  */
  void
  defineParticleGhostMasks() noexcept;

  /**
    @brief Define the SAME-level particle ghost mask on one level.
    @details For each box, records (per cell) the abutting same-level boxes whose grown region (the box
    grown by a_ghost cells) covers the cell, together with the owning/receiving rank. Uses a
    NeighborIterator to enumerate the abutting boxes, so the build is fully box-local -- no coarsen/refine
    and no Copier -- and each box fills its own CSR in two passes (count, then packed fill).
    NeighborIterator returns grow-by-one (abutting) neighbours only; this loses nothing because every box
    is min_block_size-aligned and AmrMesh::sanityCheck enforces a_ghost < min_block_size, so two
    non-abutting same-level boxes are always at least min_block_size (> a_ghost) apart and can never reach
    into each other's ghost region. Periodic boundaries are not handled.
    @param[out] a_mask   Per-box CSR ghost mask over a_dbl to fill.
    @param[in]  a_dbl    This level's grids (both the source cells and the neighbour targets).
    @param[in]  a_domain Problem domain at this level's resolution.
    @param[in]  a_ghost  Ghost width (cells) defining each neighbour's grown region.
  */
  void
  defineParticleGhostMaskSameLevel(LayoutData<ParticleGhostMask>& a_mask,
                                   const DisjointBoxLayout&       a_dbl,
                                   const ProblemDomain&           a_domain,
                                   const int                      a_ghost) noexcept;

  /**
    @brief Define one CROSS-LEVEL particle ghost mask (targets on the next-coarser or next-finer level).
    @details The mask is two separable pieces per source cell. (1) TARGETS: every box on the target level
    whose ghosted region (grown by a_ghost cells at the target resolution) contains the cell. These are
    found by walking a ghostDefine Copier motion plan -- only overlapping box pairs are visited (never a
    LayoutIterator over all boxes) -- and are multi-valued, since a cell may scatter to several target
    boxes (e.g. diagonal coarse-fine cells). (2) SOURCE restriction to the coarse-fine halo, so the nested
    interior is not scattered: for a coarser target (fine->coarse) the inner halo is built geometrically
    from the coarsened-fine layout (the ring outside each footprint box, minus the abutting footprint
    boxes, grown back inward), independent of the coarse grid's box decomposition; for a finer target
    (coarse->fine) only VALID source cells (those not covered by the finer level) are kept. The target
    resolution is reached by coarsen()/refine() of this level (which preserve the DataIndex), and halo
    cells are mapped back with refine()/coarsen(). The ghost width is measured in DESTINATION
    (target-level) cells and is a MINIMUM: whole source cells overlapping the halo are recorded, i.e.
    over-communication across the interface is accepted. Periodic boundaries are not handled.
    @param[out] a_targets       Per-box CSR ghost mask over a_thisLayout to fill.
    @param[in]  a_thisLayout    This (source) level's grids.
    @param[in]  a_otherLayout   Target level's grids (coarser or finer than this level).
    @param[in]  a_validCells    Valid-cell mask for the halo restriction. For a finer target this is this
                                level's mask (source cells are kept only where valid); for a coarser
                                target the inner halo is derived geometrically and this argument is unused.
    @param[in]  a_otherDomain   Problem domain at the target level's resolution.
    @param[in]  a_refRat        Refinement ratio between this level and the target level.
    @param[in]  a_ghost         Ghost width in target-level (destination) cells.
    @param[in]  a_coarserTarget True if the target level is coarser than this level (fine->coarse).
  */
  void
  defineParticleGhostMaskCrossLevel(LayoutData<ParticleGhostMask>&  a_targets,
                                    const DisjointBoxLayout&        a_thisLayout,
                                    const DisjointBoxLayout&        a_otherLayout,
                                    const LevelData<BaseFab<bool>>& a_validCells,
                                    const ProblemDomain&            a_otherDomain,
                                    const int                       a_refRat,
                                    const int                       a_ghost,
                                    const bool                      a_coarserTarget) noexcept;

  /**
    @brief Define the PETSc composite grid
  */
  void
  definePetscGrid() noexcept;

  /**
    @brief Get a particular realm
    @param[in] a_phase Phase
    @return Realm
  */
  PhaseRealm&
  getRealm(const phase::which_phase a_phase);
};

#include <CD_NamespaceFooter.H>

#endif
