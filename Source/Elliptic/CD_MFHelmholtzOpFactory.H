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
#include <CD_EBAMRData.H>
#include <CD_MFHelmholtzOp.H>
#include <CD_MFCoarAve.H>
#include <CD_MFFluxReg.H>
#include <CD_MFMultigridInterpolator.H>
#include <CD_MFLevelGrid.H>
#include <CD_MFBaseIVFAB.H>
#include <CD_MFHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making MFHelmholtzOp. 
*/
class MFHelmholtzOpFactory : public AMRLevelOpFactory<LevelData<MFCellFAB> > {
public:

  using MFIS             = RefCountedPtr<MultiFluidIndexSpace>;
  using AmrLevelGrids    = Vector<MFLevelGrid>;
  using AmrInterpolators = Vector<MFMultigridInterpolator>;
  using AmrFluxRegisters = Vector<MFFluxReg>;
  using AmrCoarseners    = Vector<MFCoarAve>;
  using AmrResolutions   = Vector<Real>;
  using AmrRefRatios     = Vector<int>;

  using AmrCellData      = Vector<RefCountedPtr<LevelData<MFCellFAB> > >;
  using AmrFluxData      = Vector<RefCountedPtr<LevelData<MFFluxFAB> > >;
  using AmrIrreData      = Vector<RefCountedPtr<LevelData<MFBaseIVFAB> > >;

  using DomainBCFactory = RefCountedPtr<MFHelmholtzDomainBCFactory>;
  using EBBCFactory     = RefCountedPtr<MFHelmholtzEBBCFactory>;

  using RelaxType       = MFHelmholtzOp::RelaxationMethod;

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
  */
  MFHelmholtzOpFactory(const MFIS&             a_mfis,
		       const Real&             a_alpha,
		       const Real&             a_beta,
		       const RealVect&         a_probLo,
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
		       const EBBCFactory&      a_ebBcFactory,
		       const IntVect&          a_ghostPhi,
		       const IntVect&          a_ghostRhs,
		       const RelaxType&        a_relaxationMethod,
		       const ProblemDomain&    a_bottomDomain,
		       const int&              a_jumpOrder,
		       const int&              a_jumpWeight,
		       const int&              a_blockingFactor,
		       const AmrLevelGrids&    a_deeperLevelGrids = AmrLevelGrids());

  /*!
    @brief Destructor. Does nothing. 
  */
  ~MFHelmholtzOpFactory();

  /*!
    @brief Set jump condition
  */
  void setJump(const EBAMRIVData& a_sigma, const Real& a_scale);

    /*!
    @brief Set jump condition
  */
  void setJump(const Real& a_sigma, const Real& a_scale);

  /*!
    @brief Create multigrid operator
    @param[in] a_fineDomain      Domain
    @param[in] a_depth           Depth. This specifies that the operator will be created at depth coarsen(a_fineDomain, 2^a_depth);
    @param[in] a_homogeneousOnly If true, only homogeneous boundary conditions will be needed. 
  */
  MFHelmholtzOp* MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly = true) override final;

  /*!
    @brief Create AMR operator for specified domain
    @param[in] a_domain Domain
  */
  MFHelmholtzOp* AMRnewOp(const ProblemDomain& a_domain) override final;

  /*!
    @brief Get refinement ratio to next finest level.
    @note Returns -1 when there are no finer levels. 
  */
  int refToFiner(const ProblemDomain& a_indexspace) const override final;

protected:

  /*!
    @brief Component that this operator solves for
  */
  static constexpr int m_comp  = 0;

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

  // Relaxation method
  RelaxType m_relaxMethod;

  // Number of AMR levels
  int m_numAmrLevels;

  // Ghost cells
  IntVect m_ghostPhi;
  IntVect m_ghostRhs;

  // Alpha, beta. 
  Real m_alpha;
  Real m_beta;

  // Lower-left corner of domain
  RealVect m_probLo;

  // Things that pertain to AMR levels. The first entry corresponds to the coarsest AMR level. 
  AmrLevelGrids    m_amrLevelGrids;
  AmrInterpolators m_amrInterpolators;
  AmrFluxRegisters m_amrFluxRegisters;
  AmrCoarseners    m_amrCoarseners;
  AmrRefRatios     m_amrRefRatios;
  AmrResolutions   m_amrResolutions;

  AmrCellData      m_amrAcoef;
  AmrFluxData      m_amrBcoef;
  AmrIrreData      m_amrBcoefIrreg;

  // BC factories
  DomainBCFactory m_domainBcFactory;
  EBBCFactory     m_ebBcFactory;

  // Jump
  EBAMRIVData m_amrJump;

  // Bottom drop domain
  ProblemDomain m_bottomDomain;

  // Blocking factor for when we create intermediate and deep multigrid levels
  int m_mgBlockingFactor;

  // Order for jump BCs.
  int m_jumpOrder;
  int m_jumpWeight;

  // This is for using pre-defined grids for the deeper multigrid levels, i.e. for the levels that are coarsenings of m_amrLevelGrids[0]
  AmrLevelGrids    m_deeperLevelGrids;

  // For checking if an AMR levels has multigrid levels
  std::vector<bool>  m_hasMgLevels;

  // MG operator things. Always weird to write this but e.g. m_mgLevelGrids[0] corresponds to the the multigrid levels below amr level 0. 
  Vector<AmrLevelGrids> m_mgLevelGrids;
  Vector<AmrCellData>   m_mgAcoef;
  Vector<AmrFluxData>   m_mgBcoef;
  Vector<AmrIrreData>   m_mgBcoefIrreg;

  // Jump stuff on multigrid levels. Since the jump is only defined on one phase, so is the averaging operator. 
  Vector<EBAMRIVData>   m_mgJump;
  Vector<Vector<RefCountedPtr<EbCoarAve> > > m_mgAveOp;

  /*!
    @brief Function which defines the multigrid levels for this operator factory
  */
  void defineMultigridLevels();

  /*!
    @brief Define jump data
  */
  void defineJump();

  /*!
    @brief Check if a domain is coarser than the other
    @param[in] a_domainOne The first domain
    @param[in] a_domainTwo The second domain
    @return Returns true of a_domainOne has fewer grid points than a_domainTwo
  */
  bool isCoarser(const ProblemDomain& a_domainOne, const ProblemDomain& a_domainTwo) const;

  /*!
    @brief Check if a domain is finer than the other
    @param[in] a_domainOne The first domain
    @param[in] a_domainTwo The second domain
    @return Returns true of a_domainOne has more grid points than a_domainTwo
  */
  bool isFiner(const ProblemDomain& a_domainOne, const ProblemDomain& a_domainTwo) const;

  /*!
    @brief Construct coarsening of a grid level. 
    @param[out] a_coarseGrid     The coarse grid layout. Must be a pointer to an undefined EBLevelGrid on input
    @param[in]  a_fineGrid       The coarse grid layout. Must be a pointer to an undefined EBLevelGrid on input
    @param[in]  a_refRat         Refinement ratio
    @param[in]  a_blockingFactor Blocking factor to use for grid aggregation
    @return This will return a multigrid level (i.e. one that is completely overlapping) the fine level. If we can, we use aggregation with
    the blocking factor. Otherwise we try to coarsen directly. 
  */
  bool getCoarserLayout(MFLevelGrid& a_coarseGrid, const MFLevelGrid& a_fineGrid, const int a_refRat, const int a_blockingFactor) const;

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
  void coarsenCoefficients(LevelData<MFCellFAB>&         a_coarAcoef,
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
  int findAmrLevel(const ProblemDomain& a_domain) const;
};

#include <CD_NamespaceFooter.H>

#endif
