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
#include <EBCoarseAverage.H>
#include <EBCellFAB.H>
#include <EBFluxFAB.H>
#include <BaseIVFAB.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>

// Our includes
#include <CD_MultigridInterpolator.H>
#include <CD_EBHelmholtzOp.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making variable-coefficient Helmholtz operators. 
*/
class EBHelmholtzOpFactory : public AMRLevelOpFactory<LevelData<EBCellFAB> > {
public:

  /*!
    @brief Relaxation method for the operators
  */
  enum class RelaxationMethod {
    Jacobi,
    GSRB,
    GSRBFast,
  };

  using AmrLevelGrids    = Vector<RefCountedPtr<EBLevelGrid> >;
  using AmrInterpolators = Vector<RefCountedPtr<MultigridInterpolator> >;
  using AmrFluxRegisters = Vector<RefCountedPtr<EBFluxRegister> >;
  using AmrCoarseners    = Vector<RefCountedPtr<EBCoarseAverage> >;
  
  using AmrCellData      = Vector<RefCountedPtr<LevelData<EBCellFAB> > >;
  using AmrFluxData      = Vector<RefCountedPtr<LevelData<EBFluxFAB> > >;
  using AmrIrreData      = Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >;

  using AmrResolutions   = Vector<Real>;
  using AmrRefRatios     = Vector<int>;

  using DomainBcFactory  = RefCountedPtr<BaseDomainBCFactory>; // Will be replaced with better BC class...
  using EbBcFactory      = RefCountedPtr<BaseEBBCFactory>;     // Will be replaced with better BC class...

  /*!
    @brief Disallowed constructor. Use the full constructor
  */
  EBHelmholtzOpFactory() = delete;

  /*!
    @brief Disallowed constructor. Use the full constructor
  */
  EBHelmholtzOpFactory(const EBHelmholtzOpFactory& a_otherFactory) = delete;

  /*!
    @brief Full constructor
    @param[in] a_alpha            alpha-coefficient in Helmholtz operator.
    @param[in] a_beta             beta-coefficient in Helmholtz operator.
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
  */
  EBHelmholtzOpFactory(const Real&              a_alpha,
		       const Real&              a_beta,
		       const AmrLevelGrids&     a_amrLevelGrids,
		       const AmrInterpolators&  a_amrInterpolators,
		       const AmrFluxRegisters&  a_amrFluxRegisters,
		       const AmrCoarseners&     a_amrCoarseners,
		       const AmrRefRatios&      a_amrRefRatios,
		       const AmrResolutions&    a_amrResolutions,					   
		       const AmrCellData&       a_amrAcoef,
		       const AmrFluxData&       a_amrBcoef,
		       const AmrIrreData&       a_amrBcoefIrreg,
		       const DomainBcFactory&   a_domainBCFactory,
		       const EbBcFactory&       a_ebbcFactory,
		       const IntVect&           a_ghostPhi,
		       const IntVect&           a_ghostRHS,
		       const RelaxationMethod&  a_relaxationMethod,
		       const ProblemDomain&     a_bottomDomain,
		       const int&               a_mgBlockingFactor,
		       const AmrLevelGrids&     a_deeperLevelGrids = AmrLevelGrids());
  
  /*!
    @brief Destructor. Does nothing.
  */
  ~EBHelmholtzOpFactory();

  /*!
    @brief Disallowed assignment operator
  */
  void operator=(const EBHelmholtzOpFactory& a_opin) = delete;

  /*!
    @brief Create multigrid operator
    @param[in] a_fineDomain      Domain
    @param[in] a_depth           Depth. This specifies that the operator will be created at depth coarsen(a_fineDomain, 2^a_depth);
    @param[in] a_homogeneousOnly If true, only homogeneous boundary conditions will be needed. 
  */
  EBHelmholtzOp* MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly = true) override final;

  /*!
    @brief Create AMR operator for specified domain
    @param[in] a_domain Domain
  */
  EBHelmholtzOp* AMRnewOp(const ProblemDomain& a_domain) override final;

  /*!
    @brief Get refinement ratio to next finest level.
    @note Returns -1 when there are no finer levels. 
  */
  int refToFiner(const ProblemDomain& a_indexspace) const override final;

protected:

  // Relaxation method
  RelaxationMethod m_relaxMethod;

  // Number of AMR levels
  int m_numAmrLevels;

  // Ghost cells
  IntVect m_ghostPhi;
  IntVect m_ghostRhs;

  // Alpha, beta. 
  Real m_alpha;
  Real m_beta;

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
  DomainBcFactory m_domainBcFactory;
  EbBcFactory     m_ebBcFactory;

  // Bottom drop domain
  ProblemDomain m_bottomDomain;

  // Blocking factor for when we create intermediate and deep multigrid levels
  int m_mgBlockingFactor;

  // This is for using pre-defined grids for the deeper multigrid levels, i.e. for the levels that are coarsenings of m_amrLevelGrids[0]
  AmrLevelGrids    m_deeperLevelGrids;

  // For checking if an AMR levels has multigrid levels
  std::vector<bool>  m_hasMgLevels;

  // MG operator things. Always weird to write this but e.g. m_mgLevelGrids[0] corresponds to the the multigrid levels below amr level 0. 
  Vector<AmrLevelGrids> m_mgLevelGrids;
  Vector<AmrCellData>   m_mgAcoef;
  Vector<AmrFluxData>   m_mgBcoef;
  Vector<AmrIrreData>   m_mgBcoefIrreg;

  void defineMultigridLevels();
  bool isCoarser(const ProblemDomain& a_eblgOne, const ProblemDomain& a_eblgTwo) const;
  bool isFiner(const ProblemDomain& a_eblgOne, const ProblemDomain& a_eblgTwo) const;
  bool getCoarserLayout(EBLevelGrid& a_coarseGrid, const EBLevelGrid& a_fineGrid, const int a_refRat, const int a_blockingFactor) const;
};

#include <CD_NamespaceFooter.H>

#endif
