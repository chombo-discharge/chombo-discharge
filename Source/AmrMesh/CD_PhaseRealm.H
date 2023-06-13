/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_PhaseRealm.H
  @brief  Declaration of the PhaseRealm class. 
  @author Robert Marskar
*/

#ifndef CD_PhaseRealm_H
#define CD_PhaseRealm_H

// Chombo includes
#include <DisjointBoxLayout.H>
#include <ProblemDomain.H>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_IrregAmrStencil.H>
#include <CD_LoadBalancing.H>
#include <CD_MFLevelGrid.H>
#include <CD_EBCoarAve.H>
#include <CD_EBReflux.H>
#include <CD_DomainFluxIFFAB.H>
#include <CD_DomainFluxIFFABFactory.H>
#include <CD_ParticleContainer.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_EBGradient.H>
#include <CD_EBRedistribution.H>
#include <CD_EBAMRParticleMesh.H>
#include <CD_EBAMRSurfaceDeposition.H>
#include <CD_EBCoarseToFineInterp.H>
#include <CD_EBGhostCellInterpolator.H>
#include <CD_NonConservativeDivergenceStencil.H>
#include <CD_EbCentroidInterpolationStencil.H>
#include <CD_CentroidInterpolationStencil.H>
#include <CD_NamespaceHeader.H>

// These are operator that can be defined.
static const std::string s_eb_gradient     = "eb_gradient";     // For computing gradients
static const std::string s_eb_irreg_interp = "eb_irreg_interp"; // For data recentering
static const std::string s_eb_coar_ave     = "eb_coar_ave";     // For coarsening
static const std::string s_eb_fill_patch   = "eb_fill_patch";   // For regridding data
static const std::string s_eb_fine_interp  = "eb_fine_interp";  // For linearly filling ghost cells
static const std::string s_eb_flux_reg     = "eb_flux_reg";     // For flux registeration
static const std::string s_eb_redist       = "eb_redist";       // For redistribution
static const std::string s_noncons_div     = "eb_non_cons_div"; // For computing non-conservative divergences
static const std::string s_eb_multigrid    = "eb_multigrid";    // For multigrid interpolation
static const std::string s_levelset        = "levelset";        // For putting level-set on mesh
static const std::string s_particle_mesh   = "particle_mesh";   // For doing particle-mesh operations.

/*!
  @brief Class that holds important things for doing AMR over a specific phase and processor distribution. 
  @details This class is essentially an EBAMR core-grid class, containing the usually EBAMR stuff within
  a single EBIndexSpace. In thus contains a set of grids (DisjointBoxLayout) and the usual EB stuff
  (EBLevelGrid, EBISLayout). In addition, it contains a bunch of operators that can be registered by the
  user (through the user interface AmrMesh). These operators live in this class, and they can be fetched
  e.g. by AmrMesh. 
  @note This class is a book-keeping class used by AmrMesh -- it is not meant for direct interaction. 
  @note To use run-time profiling/debugging, use PhaseRealm.profile=true or PhaseRealm.verbose=true
*/
class PhaseRealm
{
public:
  /*!
    @brief Default constructor. Must subsequently call define
  */
  PhaseRealm();

  /*!
    @brief Disallowed copy ctor
    @param[in] a_other Other realm
  */
  PhaseRealm(const PhaseRealm& a_other) = delete;

  /*!
    @brief Disallowed move ctor
    @param[in] a_other Other realm
  */
  PhaseRealm(const PhaseRealm&& a_other) = delete;

  /*!
    @brief Disallowed copy assignment operator
    @param[in] a_other Other realm
  */
  PhaseRealm&
  operator=(const PhaseRealm& a_other) = delete;

  /*!
    @brief Disallowed move assignment operator
    @param[in] a_other Other realm
  */
  PhaseRealm&
  operator=(const PhaseRealm&& a_other) = delete;

  /*!
    @brief Destructor
  */
  virtual ~PhaseRealm();

  /*!
    @brief Full define function
    @param[in] a_grids            Grids on each level, i.e. boxes and processor distributions
    @param[in] a_domains          Problem domains
    @param[in] a_refRat           Refinement ratios
    @param[in] a_dx               Grid resolutions
    @param[in] a_probLo           Lower-left corner of physical domain
    @param[in] a_finestLevel      Finest grid level
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
    @param[in] a_ebis             EBIndexSpace discrete information
  */
  void
  define(const Vector<DisjointBoxLayout>&   a_grids,
         const Vector<ProblemDomain>&       a_domains,
         const Vector<int>&                 a_refRat,
         const Vector<Real>&                a_dx,
         const RealVect                     a_probLo,
         const int                          a_finestLevel,
         const int                          a_ebGhost,
         const int                          a_numGhost,
         const int                          a_lsfGhost,
         const int                          a_redistRad,
         const int                          a_mgInterpOrder,
         const int                          a_mgInterpRadius,
         const int                          a_mgInterpWeight,
         const IrregStencil::StencilType    a_centroidStencil,
         const IrregStencil::StencilType    a_ebStencil,
         const RefCountedPtr<BaseIF>&       a_baseif,
         const RefCountedPtr<EBIndexSpace>& a_ebis);

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
    @details Issues run-time error if a_operator can't be recognized
  */
  void
  registerOperator(const std::string a_operator);

  /*!
    @brief Query if an AMR operator has been registered
    @param[in] a_operator operator name
    @return True if operator has been registered and false otherwise
  */
  bool
  queryOperator(const std::string a_operator) const;

  /*!
    @brief Return the EB index space
  */
  const RefCountedPtr<EBIndexSpace>&
  getEBIndexSpace() const;

  /*!
    @brief Get refinment ratios
    @return Returns m_refinementRatios
  */
  const Vector<int>&
  getRefinementRatios() const;

  /*!
    @brief Returns resolutions
    @return Returns m_dx
  */
  const Vector<Real>&
  getDx() const;

  /*!
    @brief Get AMR grids
    @return Immutable m_grids
  */
  const Vector<DisjointBoxLayout>&
  getGrids() const;

  /*!
    @brief Get problem domains
    @return Returns m_domains
  */
  const Vector<ProblemDomain>&
  getDomains() const;

  /*!
    @brief Get EBIS layout
    @return Immutable m_ebislayout
  */
  const Vector<EBISLayout>&
  getEBISLayout() const;

  /*!
    @brief Get EBLevelGrids
    @return Immutable m_eblg
  */
  const Vector<RefCountedPtr<EBLevelGrid>>&
  getEBLevelGrid() const;

  /*!
    @brief Get coarsened fine grids
    @return Immutable m_eblg
  */
  const Vector<RefCountedPtr<EBLevelGrid>>&
  getEBLevelGridCoFi() const;

  /*!
    @brief Return vof iterator for iterating over irregular cells in each grid patch
    @return m_vofIter (which is, in fact, mutable)
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
  getVofIterator() const;

  /*!
    @brief Get stencils for interpolating from cell center to centroids
  */
  const IrregAmrStencil<CentroidInterpolationStencil>&
  getCentroidInterpolationStencils() const;

  /*!
    @brief Get stencils for interpolating from cell center to eb centroids
  */
  const IrregAmrStencil<EbCentroidInterpolationStencil>&
  getEbCentroidInterpolationStencilStencils() const;

  /*!
    @brief Get stencils for computing "non-conservative divergences"
  */
  const IrregAmrStencil<NonConservativeDivergenceStencil>&
  getNonConservativeDivergenceStencils() const;

  /*!
    @brief Get gradient object
  */
  const Vector<RefCountedPtr<EBGradient>>&
  getGradientOp() const;

  /*!
    @brief Get particle-mesh operator
  */
  EBAMRParticleMesh&
  getParticleMesh() const;

  /*!
    @brief Get the surface deposition operator
  */
  EBAMRSurfaceDeposition&
  getSurfaceDeposition() const;

  /*!
    @brief Get coarsening operators
  */
  Vector<RefCountedPtr<EBCoarAve>>&
  getCoarseAverage() const;

  /*!
    @brief Get deposition operator
  */
  Vector<RefCountedPtr<EBCoarseFineParticleMesh>>&
  getEBCoarseFineParticleMesh() const;

  /*!
    @brief Get multigrid interpolator
  */
  Vector<RefCountedPtr<EBMultigridInterpolator>>&
  getMultigridInterpolator() const;

  /*!
    @brief Get ghost cell interpolation utility
  */
  Vector<RefCountedPtr<EBGhostCellInterpolator>>&
  getGhostCellInterpolator() const;

  /*!
    @brief Get piecewise linear ghost cell interpolation utility
  */
  Vector<RefCountedPtr<EBCoarseToFineInterp>>&
  getFineInterp() const;

  /*!
    @brief Get flux register utility
  */
  Vector<RefCountedPtr<EBReflux>>&
  getFluxRegister() const;

  /*!
    @brief Get the redistribution operators
    @param[in] a_phase Phase (gas or solid)
  */
  Vector<RefCountedPtr<EBRedistribution>>&
  getRedistributionOp() const;

  /*!
    @brief Get level-set function
  */
  const EBAMRFAB&
  getLevelset() const;

protected:
  /*!
    @brief True if things on this phase can be defined. False otherwise. Only used internally. 
  */
  bool m_isDefined;

  /*!
    @brief Relic of an ancient past. Always true. 
  */
  bool m_hasEbCf;

  /*!
    @brief Profile or not
  */
  bool m_profile;

  /*!
    @brief Verbose or not
  */
  bool m_verbose;

  /*!
    @brief Finest grid level
  */
  int m_finestLevel;

  /*!
    @brief Number of ghost cells in data holders
  */
  int m_numGhostCells;

  /*!
    @brief Number of ghost cells in EBLevelGrid/EBISLayout/EBISBox
  */
  int m_numEbGhostsCells;

  /*!
    @brief Number of ghost cells in level-set function (on the mesh)
  */
  int m_numLsfGhostCells;

  /*!
    @brief Redistribution radius
  */
  int m_redistributionRadius;

  /*!
    @brief Multigrid interpolator order
  */
  int m_multigridInterpolationOrder;

  /*!
    @brief Multigrid interpolator radius
  */
  int m_multigridInterpolationRadius;

  /*!
    @brief Multigrid interpolator weight (for least squares)
  */
  int m_multigridInterpolationWeight;

  /*!
    @brief Lower-left corner of computational domain
  */
  RealVect m_probLo;

  /*!
    @brief Stencil recentering type
  */
  IrregStencil::StencilType m_centroidStencilType;

  /*!
    @brief Stencil recentering type
  */
  IrregStencil::StencilType m_ebCentroidStencilType;

  /*!
    @brief Grid resolutions
  */
  Vector<Real> m_dx;

  /*!
    @brief Refinement ratios between levels
  */
  Vector<int> m_refinementRatios;

  /*!
    @brief EB index space
  */
  RefCountedPtr<EBIndexSpace> m_ebis;

  /*!
    @brief Implicit/SD function
  */
  RefCountedPtr<BaseIF> m_baseif;

  /*!
    @brief Operator map for checking which ones are registered
  */
  std::map<std::string, bool> m_operatorMap;

  /*!
    @brief AMR grids
  */
  Vector<DisjointBoxLayout> m_grids;

  /*!
    @brief Problem domains
  */
  Vector<ProblemDomain> m_domains;

  /*!
    @brief EBIS layouts
  */
  Vector<EBISLayout> m_ebisl;

  /*!
    @brief EB level grids
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblg;

  /*!
    @brief Coarsened fine-level EB grids. 
    @details m_eblgCoFi[lvl] holds the grids on level lvl that are a coarsening of the grids on lvl+1
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblgCoFi;

  /*!
    @brief Coarsened fine-level EB grids. 
    @details m_eblgFiCo[lvl] holds the refinement of the grids on level lvl-1
  */
  Vector<RefCountedPtr<EBLevelGrid>> m_eblgFiCo;

  /*!
    @brief Level-set function
  */
  EBAMRFAB m_levelset;

  /*! 
    @brief Cut-cell iterator
  */
  mutable Vector<RefCountedPtr<LayoutData<VoFIterator>>> m_vofIter;

  /*!
    @brief Coarsening operator
  */
  mutable Vector<RefCountedPtr<EBCoarAve>> m_coarAve;

  /*!
    @brief Multigrid interpolation utility
  */
  mutable Vector<RefCountedPtr<EBMultigridInterpolator>> m_multigridInterpolator;

  /*!
    @brief Ghost cell interpolation utility
  */
  mutable Vector<RefCountedPtr<EBGhostCellInterpolator>> m_ghostCellInterpolator;

  /*!
    @brief Regridding utility (for filling new grid patches)
  */
  mutable Vector<RefCountedPtr<EBCoarseToFineInterp>> m_ebFineInterp;

  /*!
    @brief Flux register utility
  */
  mutable Vector<RefCountedPtr<EBReflux>> m_ebReflux;

  /*!
    @brief Redistribution utilities
  */
  mutable Vector<RefCountedPtr<EBRedistribution>> m_redistributionOp;

  /*!
    @brief Gradient operator
  */
  mutable Vector<RefCountedPtr<EBGradient>> m_gradientOp;

  /*!
    @brief For doing particle-mesh operations
  */
  mutable EBAMRParticleMesh m_particleMesh;

  /*!
    @brief For doing particle deposition onto surfaces
  */
  mutable EBAMRSurfaceDeposition m_surfaceDeposition;

  /*!
    @brief Data recentering utility (center to centroid)
  */
  RefCountedPtr<IrregAmrStencil<CentroidInterpolationStencil>> m_centroidInterpolationStencil;

  /*!
    @brief Data recentering utility (center to eb centroid)
  */
  RefCountedPtr<IrregAmrStencil<EbCentroidInterpolationStencil>> m_ebCentroidInterpolationStencil;

  /*!
    @brief Stencils for non-conservative divergences
  */
  RefCountedPtr<IrregAmrStencil<NonConservativeDivergenceStencil>> m_NonConservativeDivergenceStencil;

  /*!
    @brief Define EBLevelGrids
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineEBLevelGrid(const int a_lmin);

  /*!
    @brief Define vof iterators
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineVofIterator(const int a_lmin);

  /*!
    @brief Define coarsening utility
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineEBCoarAve(const int a_lmin);

  /*!
    @brief Define multigrid interpolation utility
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineEBMultigrid(const int a_lmin);

  /*!
    @brief Define regrid utility
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineFillPatch(const int a_lmin);

  /*!
    @brief Define ghost cell interpolation utility
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineEBCoarseToFineInterp(const int a_lmin);

  /*!
    @brief Define multigrid injection utility
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineMultigridInjection(const int a_lmin);

  /*!
    @brief Define flux register utility
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineFluxReg(const int a_lmin, const int a_regsize);

  /*!
    @brief Define operator redistribution utility
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineRedistOper(const int a_lmin, const int a_regsize);

  /*!
    @brief Define operator redistribution utility
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineFineToCoarRedistOper(const int a_lmin, const int a_regsize);

  /*!
    @brief Define operator redistribution utility
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineCoarToFineRedistOper(const int a_lmin, const int a_regsize);

  /*!
    @brief Define operator redistribution utility
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineCoarToCoarRedistOper(const int a_lmin, const int a_regsize);

  /*!
    @brief Define particle mesh operator. 
  */
  void
  defineParticleMesh();

  /*!
    @brief Define gradient stencils
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineGradSten(const int a_lmin);

  /*!
    @brief Define data recentering stencils
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineIrregSten();

  /*!
    @brief Define non-conservative divergence stencils
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineNonConsDivSten();

  /*!
    @brief Put the level-set on the mesh
    @param[in] a_lmin Coarsest grid level that changes
  */
  void
  defineLevelSet(const int a_lmin, const int a_numGhost);
};

#include <CD_NamespaceFooter.H>

#endif
