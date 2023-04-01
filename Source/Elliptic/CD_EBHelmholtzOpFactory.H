/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBHelmholtzOpFactory.H
  @brief  Declaration of a factory class for making Poisson operators for multigrid. 
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzOpFactory_H
#define CD_EBHelmholtzOpFactory_H

// Chombo includes
#include <RealVect.H>
#include <EBLevelGrid.H>
#include <EBFluxRegister.H>
#include <EBCellFAB.H>
#include <EBFluxFAB.H>
#include <BaseIVFAB.H>
#include <ConductivityBaseDomainBC.H>
#include <BaseEBBC.H>

// Our includes
#include <CD_Location.H>
#include <CD_EBCoarAve.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_EBHelmholtzOp.H>
#include <CD_EBHelmholtzEBBCFactory.H>
#include <CD_EBHelmholtzDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making variable-coefficient Helmholtz operators. 
  @note  This factory is designed for making EBHelmholtzOps which have time-independent coefficients. If you need new coefficients you better set up multigrid again.
  @note  This factory is designed without time-dependence in BCs. Time-dependent BCs are still doable by letting the boundary condition classes carry a reference to 
         an externally updated time. 
*/
class EBHelmholtzOpFactory : public AMRLevelOpFactory<LevelData<EBCellFAB>>
{
public:
  // Various alias for cutting down on typing.
  using Smoother = EBHelmholtzOp::Smoother;

  using AMRMask          = Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>;
  using AmrLevelGrids    = Vector<RefCountedPtr<EBLevelGrid>>;
  using AmrInterpolators = Vector<RefCountedPtr<EBMultigridInterpolator>>;
  using AmrFluxRegisters = Vector<RefCountedPtr<EBFluxRegister>>;
  using AmrCoarseners    = Vector<RefCountedPtr<EBCoarAve>>;

  using AmrCellData = Vector<RefCountedPtr<LevelData<EBCellFAB>>>;
  using AmrFluxData = Vector<RefCountedPtr<LevelData<EBFluxFAB>>>;
  using AmrIrreData = Vector<RefCountedPtr<LevelData<BaseIVFAB<Real>>>>;

  using AmrResolutions = Vector<Real>;
  using AmrRefRatios   = Vector<int>;

  using DomainBCFactory = RefCountedPtr<EBHelmholtzDomainBCFactory>; // Will be replaced with better BC class...
  using EBBCFactory     = RefCountedPtr<EBHelmholtzEBBCFactory>;

  /*!
    @brief Disallowed constructor. Use the full constructor
  */
  EBHelmholtzOpFactory() = delete;

  /*!
    @brief Disallowed constructor. Use the full constructor
  */
  EBHelmholtzOpFactory(const EBHelmholtzOpFactory& a_otherFactory) = delete;

  /*!
    @brief Disallowed move constructor. Use the full constructor
  */
  EBHelmholtzOpFactory(const EBHelmholtzOpFactory&& a_otherFactory) = delete;

  /*!
    @brief Full constructor
    @param[in] a_dataLocation     Assumed data centering
    @param[in] a_alpha            alpha-coefficient in Helmholtz operator.
    @param[in] a_beta             beta-coefficient in Helmholtz operator.
    @param[in] a_probLo           Lower-left corner of domain
    @param[in] a_amrValidCells    Valid grid cells (true for cells not covered by a finer level)
    @param[in] a_amrLevelGrid     AMR grids
    @param[in] a_amrInterpolators Interpolator objects between AMR levels. 
    @param[in] a_amrCoarseners    Conservative coarseners between AMR levels. 
    @param[in] a_amrFluxRegisters Flux registers between AMR levels. 
    @param[in] a_amrResolutions   Grid resolutions for AMR levels. 
    @param[in] a_amrAcoef         A-coefficient in Helmholtz operator. 
    @param[in] a_amrBcoef         B-coefficient in Helmholtz operator. 
    @param[in] a_amrBcoefIrreg    B-coefficient in Helmholtz operator. This one is defined on EB faces. 
    @param[in] a_domainBCFactory  Factory class for making domain BC objects. 
    @param[in] a_ebbcFactory      Factory class for making BC objects for EB boundary conditions. 
    @param[in] a_ghostPhi         Number of ghost cells in solution vector. 
    @param[in] a_ghostRhs         Number of ghost cells in right-hand side. 
    @param[in] a_relaxationMethod Relaxation method. 
    @param[in] a_bottomDomain     Coarsest domain on which we run multigrid. Must be a coarsening of the AMR problem domains. 
    @param[in] a_deeperLevelGrids Optional object in case you want to pre-define the deeper multigrid levels. 
    @note a_deeperLevelGrids exists because the default behavior in this factory is to use direct coarsening for deeper AMR levels. 
    However, this can prevent reaching "deep enough" into the multigrid hierarchy if you use small boxes (e.g. 16^3). So, a_deeperLevelGrids 
    provide an option for using aggregation. The first entry a_deeperLevelGrids[0] should be a a factor 2 coarsening of the coarsest AMR level. 
    @note This operator is for time-independent a- and b-coefficients. 
  */
  EBHelmholtzOpFactory(const Location::Cell    a_dataLocation,
                       const Real&             a_alpha,
                       const Real&             a_beta,
                       const RealVect&         a_probLo,
                       const AMRMask&          a_amrValidCells,
                       const AmrLevelGrids&    a_amrLevelGrids,
                       const AmrInterpolators& a_amrInterpolators,
                       const AmrFluxRegisters& a_amrFluxRegisters,
                       const AmrCoarseners&    a_amrCoarseners,
                       const AmrRefRatios&     a_amrRefRatios,
                       const AmrResolutions&   a_amrResolutions,
                       const AmrCellData&      a_amrAcoef,
                       const AmrFluxData&      a_amrBcoef,
                       const AmrIrreData&      a_amrBcoefIrreg,
                       const DomainBCFactory&  a_domainBcFactory,
                       const EBBCFactory&      a_ebbcFactory,
                       const IntVect&          a_ghostPhi,
                       const IntVect&          a_ghostRHS,
                       const Smoother&         a_relaxationMethod,
                       const ProblemDomain&    a_bottomDomain,
                       const int&              a_mgBlockingFactor,
                       const AmrLevelGrids&    a_deeperLevelGrids = AmrLevelGrids());

  /*!
    @brief Destructor. Does nothing.
  */
  virtual ~EBHelmholtzOpFactory();

  /*!
    @brief Disallowed assignment operator
  */
  void
  operator=(const EBHelmholtzOpFactory& a_opin) = delete;

  /*!
    @brief Disallowed move assignment
  */
  void
  operator=(const EBHelmholtzOpFactory&& a_opin) = delete;

  /*!
    @brief Go through all MG levels and coarsen the coefficients from the finer levels
  */
  void
  coarsenCoefficientsMG();

  /*!
    @brief Create multigrid operator
    @param[in] a_fineDomain      Domain
    @param[in] a_depth           Depth. This specifies that the operator will be created at depth coarsen(a_fineDomain, 2^a_depth);
    @param[in] a_homogeneousOnly If true, only homogeneous boundary conditions will be needed. 
  */
  EBHelmholtzOp*
  MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly = true) override final;

  /*!
    @brief Create AMR operator for specified domain
    @param[in] a_domain Domain
  */
  EBHelmholtzOp*
  AMRnewOp(const ProblemDomain& a_domain) override final;

  /*!
    @brief Get refinement ratio to next finest level.
    @note Returns -1 when there are no finer levels. 
  */
  int
  refToFiner(const ProblemDomain& a_indexspace) const override final;

protected:
  /*!
    @brief Component number that is solved for
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Number of components that we solve for
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief Data location
  */
  Location::Cell m_dataLocation;

  /*!
    @brief Smoother
  */
  Smoother m_smoother;

  /*!
    @brief Number of AMR levels
  */
  int m_numAmrLevels;

  /*!
    @brief Number of ghost cells that are used. Need because of Chombo prolongation objects. 
  */
  IntVect m_ghostPhi;

  /*!
    @brief Number of ghost cells that are used. Need because of Chombo prolongation objects. 
  */
  IntVect m_ghostRhs;

  /*!
    @brief Operator alpha
  */
  Real m_alpha;

  /*!
    @brief Operator beta
  */
  Real m_beta;

  /*!
    @brief Lower-left corner of computational domain
  */
  RealVect m_probLo;

  /*!
    @brief Mask of cells that are valid on the AMR levels
  */
  AMRMask m_amrValidCells;

  // Things that pertain to AMR levels. The first entry corresponds to the coarsest AMR level.
  /*!
    @brief AMR grids
  */
  AmrLevelGrids m_amrLevelGrids;

  /*!
    @brief Ghost cell interpolations
  */
  AmrInterpolators m_amrInterpolators;

  /*!
    @brief Flux registers
  */
  AmrFluxRegisters m_amrFluxRegisters;

  /*!
    @brief Data coarseners
  */
  AmrCoarseners m_amrCoarseners;

  /*!
    @brief Refinement ratios
  */
  AmrRefRatios m_amrRefRatios;

  /*!
    @brief Resolutions one each level
  */
  AmrResolutions m_amrResolutions;

  /*!
    @brief Helmholtz A-coefficient
  */
  AmrCellData m_amrAcoef;

  /*!
    @brief Helmholtz B-coefficient
  */
  AmrFluxData m_amrBcoef;

  /*!
    @brief Helmholtz B-coefficient (on EB faces. 
  */
  AmrIrreData m_amrBcoefIrreg;

  /*!
    @brief Domain BC factory
  */
  DomainBCFactory m_domainBcFactory;

  /*!
    @brief EB BC factory. 
  */
  EBBCFactory m_ebBcFactory;

  /*!
    @brief Bottom domain, i.e. the coarsest domain which will be used in multigrid
  */
  ProblemDomain m_bottomDomain;

  /*!
    @brief Blocking factor for when we create intermediate and deep multigrid levels
  */
  int m_mgBlockingFactor;

  /*!
    @brief This is for using pre-defined grids for the deeper multigrid levels, i.e. for the levels that are coarsenings of m_amrLevelGrids[0]
  */
  AmrLevelGrids m_deeperLevelGrids;

  /*!
    @brief For checking if an AMR level has multigrid levels
  */
  std::vector<bool> m_hasMgLevels;

  /*!
    @brief Deeper grids. Always weird to write this but e.g. m_mgLevelGrids[0] corresponds to the the multigrid levels below amr level 0. 
  */
  Vector<AmrLevelGrids> m_mgLevelGrids;

  /*!
    @brief A-coefficient on deeper grids. Always weird to write this but e.g. m_mgLevelGrids[0] corresponds to the the multigrid levels below amr level 0. 
  */
  Vector<AmrCellData> m_mgAcoef;

  /*!
    @brief B-coefficient on deeper grids. Always weird to write this but e.g. m_mgLevelGrids[0] corresponds to the the multigrid levels below amr level 0. 
  */
  Vector<AmrFluxData> m_mgBcoef;

  /*!
    @brief B-coefficient on deeper grids. Always weird to write this but e.g. m_mgLevelGrids[0] corresponds to the the multigrid levels below amr level 0. 
  */
  Vector<AmrIrreData> m_mgBcoefIrreg;

  /*!
    @brief Function which defines the multigrid levels for this operator factory
  */
  void
  defineMultigridLevels();

  /*!
    @brief Check if a domain is coarser than the other
    @param[in] a_domainOne The first domain
    @param[in] a_domainTwo The second domain
    @return Returns true of a_domainOne has fewer grid points than a_domainTwo
  */
  bool
  isCoarser(const ProblemDomain& a_domainOne, const ProblemDomain& a_domainTwo) const;

  /*!
    @brief Check if a domain is finer than the other
    @param[in] a_domainOne The first domain
    @param[in] a_domainTwo The second domain
    @return Returns true of a_domainOne has more grid points than a_domainTwo
  */
  bool
  isFiner(const ProblemDomain& a_domainOne, const ProblemDomain& a_domainTwo) const;

  /*!
    @brief Construct coarsening of a grid level. 
    @param[out] a_coarseGrid     The coarse grid layout. Must be a pointer to an undefined EBLevelGrid on input
    @param[in]  a_fineGrid       The coarse grid layout. Must be a pointer to an undefined EBLevelGrid on input
    @param[in]  a_refRat         Refinement ratio
    @param[in]  a_blockingFactor Blocking factor to use for grid aggregation
    @return This will return a multigrid level (i.e. one that is completely overlapping) the fine level. If we can, we coarsen directly. 
  */
  bool
  getCoarserLayout(EBLevelGrid&       a_coarseGrid,
                   const EBLevelGrid& a_fineGrid,
                   const int          a_refRat,
                   const int          a_blockingFactor) const;

  /*!
    @brief Coarsen coefficients (conservatively)
    @param[out] a_coarAcoef      Coarse A-coefficient
    @param[out] a_coarBcoef      Coarse B-coefficient
    @param[out] a_coarBcoefIrreg Coarse B-coefficient on EB faces
    @param[in]  a_fineAcoef      Fine A-coefficient
    @param[in]  a_fineBcoef      Fine B-coefficient
    @param[in]  a_fineBcoefIrreg Fine B-coefficient on EB faces
    @param[in]  a_eblgCoar       Coarse grids
    @param[in]  a_eblgFine       Fine grids
    @param[in]  a_refRat         Coarsening factor
  */
  void
  coarsenCoefficients(LevelData<EBCellFAB>&             a_coarAcoef,
                      LevelData<EBFluxFAB>&             a_coarBcoef,
                      LevelData<BaseIVFAB<Real>>&       a_coarBcoefIrreg,
                      const LevelData<EBCellFAB>&       a_fineAcoef,
                      const LevelData<EBFluxFAB>&       a_fineBcoef,
                      const LevelData<BaseIVFAB<Real>>& a_fineBcoefIrreg,
                      const EBLevelGrid&                a_eblgCoar,
                      const EBLevelGrid&                a_eblgFine,
                      const int                         a_refRat);

  /*!
    @brief Find level corresponding to amr level
    @param[in] a_domain Problem domain. 
    @return Depth in m_amrLevelGrids corresponding to a_domain. 
    @note Run-time error if no level was found. 
  */
  int
  findAmrLevel(const ProblemDomain& a_domain) const;
};

#include <CD_NamespaceFooter.H>

#endif
