/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOpFactory.H
  @brief  Declaration of a factory class for making multifluid Helmholtz operators for multigrid. 
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzOpFactory_H
#define CD_MFHelmholtzOpFactory_H

// Chombo includes
#include <MFCellFAB.H>
#include <MFFluxFAB.H>

// Our includes
#include <CD_MultiFluidIndexSpace.H>
#include <CD_Location.H>
#include <CD_EBAMRData.H>
#include <CD_MFHelmholtzOp.H>
#include <CD_MFCoarAve.H>
#include <CD_MFFluxReg.H>
#include <CD_MFMultigridInterpolator.H>
#include <CD_MFLevelGrid.H>
#include <CD_MFBaseIVFAB.H>
#include <CD_MFHelmholtzEBBCFactory.H>
#include <CD_MFHelmholtzJumpBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making MFHelmholtzOp. 
  @note  This factory is designed for coupling EBHelmholtzOps across EBIS levels
  @note  This factory is designed without time-dependence in BCs. Time-dependent BCs are still doable by letting the boundary condition classes carry a reference to 
         an externally updated time. If you need new coefficients you will have to set up multigrid again. 
*/
class MFHelmholtzOpFactory : public AMRLevelOpFactory<LevelData<MFCellFAB>>
{
public:
  // Various alias for cutting down on typing.
  using MFIS             = RefCountedPtr<MultiFluidIndexSpace>;
  using AMRMask          = Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>;
  using AmrLevelGrids    = Vector<MFLevelGrid>;
  using AmrInterpolators = Vector<MFMultigridInterpolator>;
  using AmrFluxRegisters = Vector<MFFluxReg>;
  using AmrCoarseners    = Vector<MFCoarAve>;
  using AmrResolutions   = Vector<Real>;
  using AmrRefRatios     = Vector<int>;

  using AmrCellData = Vector<RefCountedPtr<LevelData<MFCellFAB>>>;
  using AmrFluxData = Vector<RefCountedPtr<LevelData<MFFluxFAB>>>;
  using AmrIrreData = Vector<RefCountedPtr<LevelData<MFBaseIVFAB>>>;

  using DomainBCFactory = RefCountedPtr<MFHelmholtzDomainBCFactory>;
  using EBBCFactory     = RefCountedPtr<MFHelmholtzEBBCFactory>;
  using JumpBCFactory   = RefCountedPtr<MFHelmholtzJumpBCFactory>;

  using Smoother = MFHelmholtzOp::Smoother;

  /*!
    @brief Disallowed constructor
  */
  MFHelmholtzOpFactory() = delete;

  /*!
    @brief Disallowed constructor
  */
  MFHelmholtzOpFactory(const MFHelmholtzOpFactory& a_otherFactory) = delete;

  /*!
    @brief Full constructor
    @param[in] a_mfis             Multi-fluid index space
    @param[in] a_dataLocation     Assumed data centering
    @param[in] a_alpha            alpha-coefficient in Helmholtz operator.
    @param[in] a_beta             beta-coefficient in Helmholtz operator.
    @param[in] a_probLo           Lower-left corner of domain
    @param[in] a_amrValidCells    Valid grid cells (true for cells not covered by a finer level)
    @param[in] a_amrLevelGrid     AMR grids
    @param[in] a_amrLevelGridFiCo Refined AMR grids
    @param[in] a_amrInterpolators Interpolator objects between AMR levels. 
    @param[in] a_amrCoarseners    Conservative coarseners between AMR levels. 
    @param[in] a_amrFluxRegisters Flux registers between AMR levels. 
    @param[in] a_amrResolutions   Grid resolutions for AMR levels. 
    @param[in] a_amrAcoef         A-coefficient in Helmholtz operator. 
    @param[in] a_amrBcoef         B-coefficient in Helmholtz operator. 
    @param[in] a_amrBcoefIrreg    B-coefficient in Helmholtz operator. This one is defined on EB faces. 
    @param[in] a_domainBCFactory  Factory class for making domain BC objects. 
    @param[in] a_ebbcFactory      Factory class for making BC objects for EB boundary conditions. 
    @param[in] a_jumpBCFactory    Factory class for making jump BC objects. 
    @param[in] a_ghostPhi         Number of ghost cells in solution vector. 
    @param[in] a_ghostRhs         Number of ghost cells in right-hand side. 
    @param[in] a_jumpOrder        Stencil order in multiphase cells
    @param[in] a_jumpWeight       Equation weighting in multiphase cells
    @param[in] a_relaxationMethod Relaxation method. 
    @param[in] a_bottomDomain     Coarsest domain on which we run multigrid. Must be a coarsening of the AMR problem domains. 
    @param[in] a_deeperLevelGrids Optional object in case you want to pre-define the deeper multigrid levels. 
  */
  MFHelmholtzOpFactory(const MFIS&                               a_mfis,
                       const Location::Cell                      a_dataLocation,
                       const Real&                               a_alpha,
                       const Real&                               a_beta,
                       const RealVect&                           a_probLo,
                       const AMRMask&                            a_amrValidCells,
                       const AmrLevelGrids&                      a_amrLevelGrids,
                       const Vector<RefCountedPtr<MFLevelGrid>>& a_amrLevelGridsFiCo,
                       const AmrInterpolators&                   a_amrInterpolators,
                       const AmrFluxRegisters&                   a_amrFluxRegisters,
                       const AmrCoarseners&                      a_amrCoarseners,
                       const AmrRefRatios&                       a_amrRefRatios,
                       const AmrResolutions&                     a_amrResolutions,
                       const AmrCellData&                        a_amrAcoef,
                       const AmrFluxData&                        a_amrBcoef,
                       const AmrIrreData&                        a_amrBcoefIrreg,
                       const DomainBCFactory&                    a_domainBcFactory,
                       const EBBCFactory&                        a_ebBcFactory,
                       const JumpBCFactory&                      a_jumpBcFactory,
                       const IntVect&                            a_ghostPhi,
                       const IntVect&                            a_ghostRhs,
                       const Smoother&                           a_smoother,
                       const ProblemDomain&                      a_bottomDomain,
                       const int&                                a_jumpOrder,
                       const int&                                a_jumpWeight,
                       const int&                                a_blockingFactor,
                       const AmrLevelGrids&                      a_deeperLevelGrids = AmrLevelGrids());

  /*!
    @brief Destructor. Does nothing. 
  */
  ~MFHelmholtzOpFactory();

  /*!
    @brief Set jump condition. This is passed to the operators by reference. 
    @param[in] a_sigma Jump (must be defined on gas phase)
    @param[in] a_scale Scaling factor
  */
  void
  setJump(const EBAMRIVData& a_sigma, const Real& a_scale);

  /*!
    @param[in] a_sigma Jump
    @param[in] a_scale Scaling factor
  */
  void
  setJump(const Real& a_sigma, const Real& a_scale);

  /*!
    @brief Get the BC jump factor. 
  */
  const EBAMRIVData&
  getSigma() const;

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
  MFHelmholtzOp*
  MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly = true) override final;

  /*!
    @brief Create AMR operator for specified domain
    @param[in] a_domain Domain
  */
  MFHelmholtzOp*
  AMRnewOp(const ProblemDomain& a_domain) override final;

  /*!
    @brief Get refinement ratio to next finest level.
    @note Returns -1 when there are no finer levels. 
  */
  int
  refToFiner(const ProblemDomain& a_indexspace) const override final;

protected:
  /*!
    @brief Component that this operator solves for
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Number of components that this operator solves for.
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief Main phase
  */
  static constexpr int m_mainPhase = 0;

  /*!
    @brief Index space
  */
  MFIS m_mfis;

  /*!
    @brief Data centering
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

  /*!
    @brief AMR grids
  */
  AmrLevelGrids m_amrLevelGrids;

  /*!
    @brief Refined AMR grids.
    @details Level lvl holds grids on level lvl that are a refinement of grids on level lvl-1
  */
  Vector<RefCountedPtr<MFLevelGrid>> m_amrLevelGridsFiCo;

  /*!
    @brief Multigrid interpolators
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
    @brief Grid resolutions
  */
  AmrResolutions m_amrResolutions;

  /*!
    @brief Operator A-coefficient
  */
  AmrCellData m_amrAcoef;

  /*!
    @brief Operator B-coefficient
  */
  AmrFluxData m_amrBcoef;

  /*!
    @brief Operator B-coefficient
  */
  AmrIrreData m_amrBcoefIrreg;

  /*!
    @brief Domain BC factory
  */
  DomainBCFactory m_domainBcFactory;

  /*!
    @brief EB BC factory
  */
  EBBCFactory m_ebBcFactory;

  /*!
    @brief Jump BC factory
  */
  JumpBCFactory m_jumpBcFactory;

  /*!
    @brief Jump on multiphase cells. 
  */
  EBAMRIVData m_amrJump;

  /*!
    @brief Coarsest level where we run multigrid
  */
  ProblemDomain m_bottomDomain;

  /*!
    @brief Blocking factor for grid aggregration
  */
  int m_mgBlockingFactor;

  /*!
    @brief Stencil order in jump cells
  */
  int m_jumpOrder;

  /*!
    @brief Weighting for least squares reconstruction in jump cells
  */
  int m_jumpWeight;

  /*!
    @brief Deeper grids
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
    @brief Jump coefficient on deeper grids
  */
  Vector<EBAMRIVData> m_mgJump;

  /*!
    @brief Coarsneing operator on deeper grids (used for coarsening jump)
  */
  Vector<Vector<RefCountedPtr<EBCoarAve>>> m_mgAveOp;

  /*!
    @brief Function which defines the multigrid levels for this operator factory
  */
  void
  defineMultigridLevels();

  /*!
    @brief Define jump data
  */
  void
  defineJump();

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
    @return This will return a multigrid level (i.e. one that is completely overlapping) the fine level. If we can, we use aggregation with
    the blocking factor. Otherwise we try to coarsen directly. 
  */
  bool
  getCoarserLayout(MFLevelGrid&       a_coarseGrid,
                   const MFLevelGrid& a_fineGrid,
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
  coarsenCoefficients(LevelData<MFCellFAB>&         a_coarAcoef,
                      LevelData<MFFluxFAB>&         a_coarBcoef,
                      LevelData<MFBaseIVFAB>&       a_coarBcoefIrreg,
                      const LevelData<MFCellFAB>&   a_fineAcoef,
                      const LevelData<MFFluxFAB>&   a_fineBcoef,
                      const LevelData<MFBaseIVFAB>& a_fineBcoefIrreg,
                      const MFLevelGrid&            a_eblgCoar,
                      const MFLevelGrid&            a_eblgFine,
                      const int                     a_refRat);

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
