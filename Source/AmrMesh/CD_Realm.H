/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Realm.H
  @brief  Declaration of the Realm class. 
  @author Robert Marskar
*/

#ifndef CD_Realm_H
#define CD_Realm_H

// Our includes
#include <CD_PhaseRealm.H>
#include <CD_EBAMRParticleMesh.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_MFLevelGrid.H>
#include <CD_LevelTiles.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Identifier for particle halo masks
*/
static const std::string s_particle_halo = "particle_halo";

/*!
  @brief Identifier for particle valid region. 
*/
static const std::string s_particle_valid_region = "pvr";

/*!
  @brief Alias for cutting down on the typic of booleans defined over AMR grids. 
*/
using AMRMask = Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>;

/*!
  @brief Class for holding multifluid Realms, which are the Realms that we actually use. 
*/
class Realm
{
public:
  /*!
    @brief Identifier for perimal realm
  */
  const static std::string Primal;

  /*!
    @brief Identifier for perimal realm
  */
  const static std::string primal;

  /*!
    @brief Default constructor. Must subsequently call define
  */
  Realm();

  /*!
    @brief Disallowed copy ctor
    @param[in] a_other Other realm
  */
  Realm(const Realm& a_other) = delete;

  /*!
    @brief Disallowed move ctor
    @param[in] a_other Other realm
  */
  Realm(const Realm&& a_other) = delete;

  /*!
    @brief Disallowed copy assignment operator
    @param[in] a_other Other realm
  */
  Realm&
  operator=(const Realm& a_other) = delete;

  /*!
    @brief Disallowed move assignment operator
    @param[in] a_other Other realm
  */
  Realm&
  operator=(const Realm&& a_other) = delete;

  /*!
    @brief Destructor
  */
  ~Realm();

  /*!
    @brief Full define function
    @param[in] a_grids            Grids on each level, i.e. boxes and processor distributions
    @param[in] a_domains          Problem domains
    @param[in] a_refRat           Refinement ratios
    @param[in] a_dx               Grid resolutions
    @param[in] a_probLo           Lower-left corner of physical domain
    @param[in] a_finestLevel      Finest grid level
    @param[in] a_blockingFactor   Grid blocking factor
    @param[in] a_ebGhost          Number of ghost cells used for geometric information when making operators
    @param[in] a_numGhost         Number of ghost cells in data holders
    @param[in] a_lsfGhost         Number of ghost cells used when computing the level set on the mesh
    @param[in] a_redistRad        Redistribution radius
    @param[in] a_mgInterpOrder    Multigrid ghost cell interpolation order
    @param[in] a_mgInterpRadius   Multigrid ghost cell interpolation radius
    @param[in] a_mgInterpWeight   Multigrid ghost cell interpolation weight (for least squares)
    @param[in] a_centroidStencil  Stencil type for interpolation from cell centers to cell centroids
    @param[in] a_ebStencil        Stencil type for interpolation from cell centers to eb centroids
    @param[in] a_baseif           Implicit/signed distance function 
    @param[in] a_mfis             EBIndexSpace discrete information
  */
  void
  define(const Vector<DisjointBoxLayout>&                          a_grids,
         const Vector<ProblemDomain>&                              a_domains,
         const Vector<int>&                                        a_refRat,
         const Vector<Real>&                                       a_dx,
         const RealVect                                            a_probLo,
         const int                                                 a_finestLevel,
         const int                                                 a_blockingFactor,
         const int                                                 a_ebGhost,
         const int                                                 a_numGhost,
         const int                                                 a_lsfGhost,
         const int                                                 a_redistRad,
         const int                                                 a_mgInterpOrder,
         const int                                                 a_mgInterpRadius,
         const int                                                 a_mgInterpWeight,
         const CellCentroidInterpolation::Type                     a_centroidStencil,
         const IrregStencil::StencilType                           a_ebStencil,
         const std::map<phase::which_phase, RefCountedPtr<BaseIF>> a_baseif,
         const RefCountedPtr<MultiFluidIndexSpace>&                a_mfis);

  /*!
    @brief Set grid method
    @param[in] a_grids       New grids
    @param[in] a_finestLevel New finest grid level
    @details This sets m_grids to be a_grids, but does not perform a regrid of remaining members.
  */
  void
  setGrids(const Vector<DisjointBoxLayout>& a_grids, const int a_finestLevel);

  /*!
    @brief Perform pre-regrid operations
  */
  void
  preRegrid();

  /*!
    @brief Regrid method for EBAMR base
    @param[in] a_lmin Coarsest grid level that changed. 
    @details This regrids EBLevelGrid and iterators
  */
  void
  regridBase(const int a_lmin);

  /*!
    @brief Regrid method for EBAMR operators.
    @param[in] a_lmin Coarsest grid level that changed. 
    @details This regrids all operators. 
  */
  void
  regridOperators(const int a_lmin);

  /*!
    @brief Register an AMR operator
    @param[in] a_operator Operator name
    @param[in] a_phase Phase where operator is registered
    @details Issues run-time error if a_operator can't be recognized
  */
  void
  registerOperator(const std::string a_operator, const phase::which_phase a_phase);

  /*!
    @brief Query if an AMR operator has been registered
    @param[in] a_operator operator name
    @param[in] a_phase Phase where operator is registered
    @return True if operator has been registered and false otherwise
  */
  bool
  queryOperator(const std::string a_operator, const phase::which_phase a_phase) const;

  /*!
    @brief Register a mask
    @param[in] a_mask Mask name
    @param[in] a_phase Phase where mask is registered
  */
  void
  registerMask(const std::string a_mask, const int a_buffer);

  /*!
    @brief Query if a mask has been registered
    @param[in] a_mask Mask name
    @param[in] a_phase Phase where mask is registered
    @return true if mask has been registered
  */
  bool
  queryMask(const std::string a_mask, const int a_buffer) const;

  /*!
    @brief Get refinement ratios
    @return m_refinementRatios
  */
  const Vector<int>&
  getRefinementRatios() const;

  /*!
    @brief Get grid resolutions
    @return m_dx
  */
  const Vector<Real>&
  getDx() const;

  /*!
    @brief Get AMR grids
    @return m_grids
  */
  const Vector<DisjointBoxLayout>&
  getGrids() const;

  /*!
    @brief Get problem domains
  */
  const Vector<ProblemDomain>&
  getDomains() const;

  /*!
    @brief Get MFLevelGrids
  */
  Vector<RefCountedPtr<MFLevelGrid>>&
  getMFLevelGrid();

  /*!
    @brief Get EBIndexSpace for particular phase
    @param[in] a_phase Phase
  */
  const RefCountedPtr<EBIndexSpace>&
  getEBIndexSpace(const phase::which_phase a_phase) const;

  /*!
    @brief Get EBISLayout for particular phase
    @param[in] a_phase Phase
  */
  const Vector<EBISLayout>&
  getEBISLayout(const phase::which_phase a_phase) const;

  /*!
    @brief Get EBLevelGrid for particular phase
    @param[in] a_phase Phase
  */
  const Vector<RefCountedPtr<EBLevelGrid>>&
  getEBLevelGrid(const phase::which_phase a_phase) const;

  /*!
    @brief Get coarsened fine EBLevelGrid for particular phase
    @param[in] a_phase Phase
  */
  const Vector<RefCountedPtr<EBLevelGrid>>&
  getEBLevelGridCoFi(const phase::which_phase a_phase) const;

  /*!
    @brief Get VoFIterators for particular phase
    @param[in] a_phase Phase
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
  getVofIterator(const phase::which_phase a_phase) const;

  /*!
    @brief Get stencils for interpolating from cell center to eb centroids
    @param[in] a_phase Phase
  */
  const IrregAmrStencil<EbCentroidInterpolationStencil>&
  getEbCentroidInterpolationStencilStencils(const phase::which_phase a_phase) const;

  /*!
    @brief Get objects for computing interpolation of cell-centered data to cell centroids.
    @param[in] a_phase Phase
  */
  const Vector<RefCountedPtr<CellCentroidInterpolation>>&
  getCellCentroidInterpolation(const phase::which_phase a_phase) const;

  /*!
    @brief Get objects for computing "non-conservative divergences"
    @param[in] a_phase Phase
  */
  const Vector<RefCountedPtr<EBNonConservativeDivergence>>&
  getNonConservativeDivergence(const phase::which_phase a_phase) const;

  /*!
    @brief Get data coarsening utility
    @param[in] a_phase Phase
  */
  Vector<RefCountedPtr<EBCoarAve>>&
  getCoarseAverage(const phase::which_phase a_phase);

  /*!
    @brief Get the multigrid interpolator
    @param[in] a_phase Phase
  */
  Vector<RefCountedPtr<EBMultigridInterpolator>>&
  getMultigridInterpolator(const phase::which_phase a_phase);

  /*!
    @brief Get the particle mesh operator
    @param[in] a_phase Which phase
  */
  EBAMRParticleMesh&
  getParticleMesh(const phase::which_phase a_phase);

  /*!
    @brief Get the particle mesh operator
    @param[in] a_phase Which phase
  */
  EBAMRSurfaceDeposition&
  getSurfaceDeposition(const phase::which_phase a_phase);

  /*!
    @brief Get gradient operator
    @param[in] a_phase Phase
  */
  const Vector<RefCountedPtr<EBGradient>>&
  getGradientOp(const phase::which_phase a_phase) const;

  /*!
    @brief Get ghost cell interpolator
    @param[in] a_phase Phase
  */
  Vector<RefCountedPtr<EBGhostCellInterpolator>>&
  getGhostCellInterpolator(const phase::which_phase a_phase);

  /*!
    @brief Get piecewise linear ghost cell interpolation utility
    @param[in] a_phase Phase
  */
  Vector<RefCountedPtr<EBCoarseToFineInterp>>&
  getFineInterp(const phase::which_phase a_phase);

  /*!
    @brief Get flux register utility
    @param[in] a_phase Phase
  */
  Vector<RefCountedPtr<EBReflux>>&
  getFluxRegister(const phase::which_phase a_phase);

  /*!
    @brief Get the redistribution operators
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBFluxRedistribution>>&
  getRedistributionOp(const phase::which_phase a_phase);

  /*!
    @brief Get level-set function
    @param[in] a_phase Phase
  */
  const EBAMRFAB&
  getLevelset(const phase::which_phase a_phase) const;

  /*!
    @brief Get AMR mask
    @param[in] a_phase Phase
  */
  const AMRMask&
  getMask(const std::string a_mask, const int a_buffer) const;

  /*!
    @brief Get all valid cells
  */
  const AMRMask&
  getValidCells() const;

  /*!
    @brief Get the tiled space
  */
  const Vector<RefCountedPtr<LevelTiles>>&
  getLevelTiles() const noexcept;

protected:
  /*!
    @brief Realm defined or not
  */
  bool m_isDefined;

  /*!
    @brief Realm verbosity. Used for debugging.
  */
  int m_verbosity;

  /*!
    @brief Finest AMR level
  */
  int m_finestLevel;

  /*!
    @brief Blocking factor
  */
  int m_blockingFactor;

  /*!
    @brief Physical coordinates of lower-left corner in simulation domain
  */
  RealVect m_probLo;

  /*!
    @brief Grid resolutions
  */
  Vector<Real> m_dx;

  /*!
    @brief Refinement ratios
  */
  Vector<int> m_refinementRatios;

  /*!
    @brief AMR grids
  */
  Vector<DisjointBoxLayout> m_grids;

  /*!
    @brief Domains
  */
  Vector<ProblemDomain> m_domains;

  /*!
    @brief MFLevelGrids
  */
  Vector<RefCountedPtr<MFLevelGrid>> m_mflg;

  /*!
    @brief Level "tiles" -- i.e. the grids viewed as constant-size boxes.
  */
  Vector<RefCountedPtr<LevelTiles>> m_levelTiles;

  /*!
    @brief Index spaces
  */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /*!
    @brief Mask which is true for all valid cells
  */
  AMRMask m_validCells;

  /*!
    @brief Phase realms. 
    @details Mutable because some of the operators are mutable. 
  */
  mutable std::map<phase::which_phase, RefCountedPtr<PhaseRealm>> m_realms;

  /*!
    @brief Implicit/signed distance functions on each phase
  */
  std::map<phase::which_phase, RefCountedPtr<BaseIF>> m_baseif;

  /*!
    @brief Masks
  */
  std::map<std::pair<std::string, int>, AMRMask> m_masks;

  /*!
    @brief Define MFLevelGrid
    @param[in] a_lmin Coarsest level that changed during regrid
  */
  void
  defineMFLevelGrid(const int a_lmin);

  /*!
    @brief Define masks
    @param[in] a_lmin Coarsest level that changed during regrid
  */
  void
  defineMasks(const int a_lmin);

  /*!
    @brief Define particle halo amsks
    @param[in] a_lmin Coarsest level that changed during regrid
  */
  void
  defineHaloMasks(const int a_lmin);

  /*!
    @brief Define particle halo masks. This is the level-version. 
    @param[out] a_coarMask   Mask
    @param[in]  a_domainCoar Coarse domain
    @param[in]  a_domainFine Fine domain
    @param[in]  a_gridsCoar  Coarse grids
    @param[in]  a_gridsFine  Fine grids
    @param[in]  a_buffer     Buffer (i.e. width of "halo")
    @param[in]  a_refRat     Refinement ratio between the two levels
  */
  void
  defineHaloMask(LevelData<BaseFab<bool>>& a_coarMask,
                 const ProblemDomain&      a_domainCoar,
                 const ProblemDomain&      a_domainFine,
                 const DisjointBoxLayout&  a_gridsCoar,
                 const DisjointBoxLayout&  a_gridsFine,
                 const int                 a_buffer,
                 const int                 a_refRat);

  /*!
    @brief Define the valid cell region mask. 
    @param[in] a_lmin Coarsest level that changed during regrid. 
  */
  void
  defineValidCells();

  /*!
    @brief Define the tiled space
  */
  void
  defineLevelTiles() noexcept;

  /*!
    @brief Get a particular realm
    @param[in] a_phase Phase
  */
  PhaseRealm&
  getRealm(const phase::which_phase a_phase);
};

#include <CD_NamespaceFooter.H>

#endif
