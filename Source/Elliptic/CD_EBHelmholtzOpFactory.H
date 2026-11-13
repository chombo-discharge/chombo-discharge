/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBHelmholtzOpFactory.H
 * @brief  Declaration of a factory class for making Poisson operators for multigrid.
 * @author Robert Marskar
 */

#ifndef CD_EBHELMHOLTZOPFACTORY_H
#define CD_EBHELMHOLTZOPFACTORY_H

// Chombo includes
#include <RealVect.H>
#include <EBLevelGrid.H>
#include <EBCellFAB.H>
#include <EBFluxFAB.H>
#include <BaseIVFAB.H>
#include <ConductivityBaseDomainBC.H>
#include <BaseEBBC.H>

// Our includes
#include <CD_Location.H>
#include <CD_EBCoarAve.H>
#include <CD_EBReflux.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_EBHelmholtzOp.H>
#include <CD_EBHelmholtzEBBCFactory.H>
#include <CD_EBHelmholtzDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Factory class for making variable-coefficient Helmholtz operators.
 * @note  This factory is designed for making EBHelmholtzOps which have time-independent coefficients. If you need new
 * coefficients you better set up multigrid again.
 * @note  This factory is designed without time-dependence in BCs. Time-dependent BCs are still doable by letting the
 * boundary condition classes carry a reference to an externally updated time.
 */
class EBHelmholtzOpFactory : public AMRLevelOpFactory<LevelData<EBCellFAB>>
{
public:
  // Various alias for cutting down on typing.

  /**
   * @brief Smoother type alias.
   */
  using Smoother = EBHelmholtzOp::Smoother;

  /**
   * @brief alias type alias.
   */
  using AmrLevelGrids = Vector<RefCountedPtr<EBLevelGrid>>;

  /**
   * @brief alias type alias.
   */
  using AmrMask = Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>;

  /**
   * @brief alias type alias.
   */
  using AmrInterpolators = Vector<RefCountedPtr<EBMultigridInterpolator>>;

  /**
   * @brief alias type alias.
   */
  using AmrFluxRegisters = Vector<RefCountedPtr<EBReflux>>;

  /**
   * @brief alias type alias.
   */
  using AmrCoarseners = Vector<RefCountedPtr<EBCoarAve>>;

  /**
   * @brief alias type alias.
   */
  using AmrCellData = Vector<RefCountedPtr<LevelData<EBCellFAB>>>;

  /**
   * @brief alias type alias.
   */
  using AmrFluxData = Vector<RefCountedPtr<LevelData<EBFluxFAB>>>;

  /**
   * @brief alias type alias.
   */
  using AmrIrreData = Vector<RefCountedPtr<LevelData<BaseIVFAB<Real>>>>;

  /**
   * @brief alias type alias.
   */
  using AmrResolutions = Vector<Real>;

  /**
   * @brief alias type alias.
   */
  using AmrRefRatios = Vector<int>;

  /**
   * @brief alias type alias.
   */
  using DomainBCFactory = RefCountedPtr<EBHelmholtzDomainBCFactory>; // Will be replaced with better BC class...

  /**
   * @brief alias type alias.
   */
  using EBBCFactory = RefCountedPtr<EBHelmholtzEBBCFactory>;

  /**
   * @brief Disallowed constructor. Use the full constructor
   */
  EBHelmholtzOpFactory() = delete;

  /**
   * @brief Disallowed constructor. Use the full constructor
   */
  EBHelmholtzOpFactory(const EBHelmholtzOpFactory& a_otherFactory) = delete;

  /**
   * @brief Disallowed move constructor. Use the full constructor
   */
  EBHelmholtzOpFactory(const EBHelmholtzOpFactory&& a_otherFactory) = delete;

  /**
   * @brief Full constructor
   * @param[in] a_dataLocation     Assumed data centering
   * @param[in] a_alpha            alpha-coefficient in Helmholtz operator.
   * @param[in] a_beta             beta-coefficient in Helmholtz operator.
   * @param[in] a_probLo           Lower-left corner of domain
   * @param[in] a_amrInterpolators Interpolator objects between AMR levels.
   * @param[in] a_amrCoarseners    Conservative coarseners between AMR levels.
   * @param[in] a_amrFluxRegisters Flux registers between AMR levels.
   * @param[in] a_amrResolutions   Grid resolutions for AMR levels.
   * @param[in] a_amrAcoef         A-coefficient in Helmholtz operator.
   * @param[in] a_amrBcoef         B-coefficient in Helmholtz operator.
   * @param[in] a_amrBcoefIrreg    B-coefficient in Helmholtz operator. This one is defined on EB faces.
   * @param[in] a_ebbcFactory      Factory class for making BC objects for EB boundary conditions.
   * @param[in] a_ghostPhi         Number of ghost cells in solution vector.
   * @param[in] a_smoother         Relaxation method.
   * @param[in] a_relaxFactor      Relaxation factor
   * @param[in] a_chebyOrder       Chebyshev polynomial degree (only used when a_smoother is Smoother::Chebyshev)
   * @param[in] a_chebyEigRatio    Chebyshev eigenvalue ratio lambda_max/lambda_min (only used when a_smoother is
   * Smoother::Chebyshev)
   * @param[in] a_rasInnerSweeps   Number of inner frozen-ghost sweeps per restricted-additive-Schwarz iteration (only
   * used when a_smoother is Smoother::RestrictedAdditiveSchwarz)
   * @param[in] a_mgLevelGrids     Multigrid levels under each AMR level, as built by the realm (see
   * PhaseRealm::getMultigridLevelGrids). The factory uses the entries that are not coarser than a_bottomDomain.
   * @param[in] a_mgLevelGridsCoFi Coarsened grids between consecutive multigrid levels, as built by the realm (see
   * PhaseRealm::getMultigridLevelGridsCoFi).
   * @param[in] a_mgCoarseners     Coarse averaging between consecutive multigrid levels, as built by the realm (see
   * PhaseRealm::getMultigridCoarseAverage).
   * @param[in] a_amrLevelGridsCoFi Coarsened grids between AMR levels. Entry lvl is AMR level lvl+1 coarsened to AMR
   * level lvl (see PhaseRealm::getEBLevelGridCoFi).
   * @param[in] a_bottomDomain     Coarsest domain on which we run multigrid. Must be a coarsening of the AMR problem
   * domains.
   * @param[in] a_refluxFree       If true, use the reflux-free AMR operator (skips coarse-fine flux refluxing)
   * @note This operator is for time-independent a- and b-coefficients.
   * @param[in] a_amrLevelGrids Amr level grids
   * @param[in] a_validCells Valid cells
   * @param[in] a_amrRefRatios Amr ref ratios
   * @param[in] a_domainBcFactory Domain bc factory
   * @param[in] a_ghostRHS Ghost r h s
   */
  EBHelmholtzOpFactory(Location::Cell               a_dataLocation,
                       const Real&                  a_alpha,
                       const Real&                  a_beta,
                       const RealVect&              a_probLo,
                       const AmrLevelGrids&         a_amrLevelGrids,
                       const AmrMask&               a_validCells,
                       const AmrInterpolators&      a_amrInterpolators,
                       const AmrFluxRegisters&      a_amrFluxRegisters,
                       const AmrCoarseners&         a_amrCoarseners,
                       const AmrRefRatios&          a_amrRefRatios,
                       const AmrResolutions&        a_amrResolutions,
                       const AmrCellData&           a_amrAcoef,
                       const AmrFluxData&           a_amrBcoef,
                       const AmrIrreData&           a_amrBcoefIrreg,
                       const DomainBCFactory&       a_domainBcFactory,
                       const EBBCFactory&           a_ebbcFactory,
                       const IntVect&               a_ghostPhi,
                       const IntVect&               a_ghostRHS,
                       const Smoother&              a_smoother,
                       const Real&                  a_relaxFactor,
                       const int&                   a_chebyOrder,
                       const Real&                  a_chebyEigRatio,
                       const int&                   a_rasInnerSweeps,
                       const Vector<AmrLevelGrids>& a_mgLevelGrids,
                       const Vector<AmrLevelGrids>& a_mgLevelGridsCoFi,
                       const Vector<AmrCoarseners>& a_mgCoarseners,
                       const AmrLevelGrids&         a_amrLevelGridsCoFi,
                       const ProblemDomain&         a_bottomDomain,
                       const bool                   a_refluxFree);

  /**
   * @brief Destructor. Does nothing.
   */
  ~EBHelmholtzOpFactory() override;

  /**
   * @brief Disallowed assignment operator
   */
  void
  operator=(const EBHelmholtzOpFactory& a_opin) = delete;

  /**
   * @brief Disallowed move assignment
   */
  void
  operator=(const EBHelmholtzOpFactory&& a_opin) = delete;

  /**
   * @brief Go through all MG levels and coarsen the coefficients from the finer levels
   */
  void
  coarsenCoefficientsMG();

  /**
   * @brief Create multigrid operator
   * @param[in] a_fineDomain      Domain
   * @param[in] a_depth           Depth. This specifies that the operator will be created at depth coarsen(a_fineDomain,
   * 2^a_depth);
   * @param[in] a_homogeneousOnly If true, only homogeneous boundary conditions will be needed.
   * @return Return value
   */
  EBHelmholtzOp*
  MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly = true) final;

  /**
   * @brief Create AMR operator for specified domain
   * @param[in] a_domain Domain
   * @return Return value
   */
  EBHelmholtzOp*
  AMRnewOp(const ProblemDomain& a_domain) final;

  /**
   * @brief Get refinement ratio to next finest level.
   * @note Returns -1 when there are no finer levels.
   * @return Return value
   * @param[in] a_domain Domain
   */
  int
  refToFiner(const ProblemDomain& a_domain) const final;

protected:
  /**
   * @brief Component number that is solved for
   */
  static constexpr int m_comp = 0;

  /**
   * @brief Number of components that we solve for
   */
  static constexpr int m_nComp = 1;

  /**
   * @brief Data location
   */
  Location::Cell m_dataLocation;

  /**
   * @brief Smoother
   */
  Smoother m_smoother;

  /**
   * @brief Chebyshev polynomial degree (only used when m_smoother is Smoother::Chebyshev)
   */
  int m_chebyOrder;

  /**
   * @brief Chebyshev eigenvalue ratio lambda_max/lambda_min (only used when m_smoother is Smoother::Chebyshev)
   */
  Real m_chebyEigRatio;

  /**
   * @brief Number of inner frozen-ghost sweeps per restricted-additive-Schwarz iteration (only used when m_smoother is
   * Smoother::RestrictedAdditiveSchwarz)
   */
  int m_rasInnerSweeps;

  /**
   * @brief Number of AMR levels
   */
  int m_numAmrLevels;

  /**
   * @brief Number of ghost cells that are used. Need because of Chombo prolongation objects.
   */
  IntVect m_ghostPhi;

  /**
   * @brief Number of ghost cells that are used. Need because of Chombo prolongation objects.
   */
  IntVect m_ghostRhs;

  /**
   * @brief Operator alpha
   */
  Real m_alpha;

  /**
   * @brief Operator beta
   */
  Real m_beta;

  /**
   * @brief Relaxation factor
   */
  Real m_relaxFactor;

  /**
   * @brief Use reflux-free AMR operator
   */
  bool m_refluxFree;

  /**
   * @brief Lower-left corner of computational domain
   */
  RealVect m_probLo;

  // Things that pertain to AMR levels. The first entry corresponds to the coarsest AMR level.

  /**
   * @brief AMR grids
   */
  AmrLevelGrids m_amrLevelGrids;

  /**
   * @brief Valid cells
   */
  AmrMask m_validCells;

  /**
   * @brief Ghost cell interpolations
   */
  AmrInterpolators m_amrInterpolators;

  /**
   * @brief Flux registers
   */
  AmrFluxRegisters m_amrFluxRegisters;

  /**
   * @brief Data coarseners
   */
  AmrCoarseners m_amrCoarseners;

  /**
   * @brief Refinement ratios
   */
  AmrRefRatios m_amrRefRatios;

  /**
   * @brief Resolutions one each level
   */
  AmrResolutions m_amrResolutions;

  /**
   * @brief Helmholtz A-coefficient
   */
  AmrCellData m_amrAcoef;

  /**
   * @brief Helmholtz B-coefficient
   */
  AmrFluxData m_amrBcoef;

  /**
   * @brief Helmholtz B-coefficient (on EB faces.
   */
  AmrIrreData m_amrBcoefIrreg;

  /**
   * @brief Domain BC factory
   */
  DomainBCFactory m_domainBcFactory;

  /**
   * @brief EB BC factory.
   */
  EBBCFactory m_ebBcFactory;

  /**
   * @brief Bottom domain, i.e. the coarsest domain which will be used in multigrid
   */
  ProblemDomain m_bottomDomain;

  /**
   * @brief Coarsened grids between AMR levels. Entry lvl is AMR level lvl+1 coarsened to AMR level lvl.
   */
  AmrLevelGrids m_amrLevelGridsCoFi;

  /**
   * @brief For checking if an AMR level has multigrid levels
   */
  std::vector<bool> m_hasMgLevels;

  /**
   * @brief Deeper grids. Always weird to write this but e.g. m_mgLevelGrids[0] corresponds to the the multigrid levels
   * below amr level 0.
   */
  Vector<AmrLevelGrids> m_mgLevelGrids;

  /**
   * @brief A-coefficient on deeper grids. Always weird to write this but e.g. m_mgLevelGrids[0] corresponds to the the
   * multigrid levels below amr level 0.
   */
  Vector<AmrCellData> m_mgAcoef;

  /**
   * @brief B-coefficient on deeper grids. Always weird to write this but e.g. m_mgLevelGrids[0] corresponds to the the
   * multigrid levels below amr level 0.
   */
  Vector<AmrFluxData> m_mgBcoef;

  /**
   * @brief B-coefficient on deeper grids. Always weird to write this but e.g. m_mgLevelGrids[0] corresponds to the the
   * multigrid levels below amr level 0.
   */
  Vector<AmrIrreData> m_mgBcoefIrreg;

  /**
   * @brief Coarse averaging between consecutive multigrid levels. m_mgCoarseners[lvl][k] averages from
   * m_mgLevelGrids[lvl][k] onto m_mgLevelGrids[lvl][k+1].
   */
  Vector<AmrCoarseners> m_mgCoarseners;

  /**
   * @brief Coarsened grids between consecutive multigrid levels. m_mgLevelGridsCoFi[lvl][k] is m_mgLevelGrids[lvl][k]
   * coarsened by 2.
   */
  Vector<AmrLevelGrids> m_mgLevelGridsCoFi;

  /**
   * @brief Define the multigrid levels for this operator factory, and the coefficients on them.
   * @details Takes the entries of the realm's multigrid levels that are not coarser than the bottom domain.
   * @param[in] a_mgLevelGrids     Multigrid levels under each AMR level, as built by the realm.
   * @param[in] a_mgLevelGridsCoFi Coarsened grids between consecutive multigrid levels, as built by the realm.
   * @param[in] a_mgCoarseners     Coarse averaging between consecutive multigrid levels, as built by the realm.
   */
  void
  defineMultigridLevels(const Vector<AmrLevelGrids>& a_mgLevelGrids,
                        const Vector<AmrLevelGrids>& a_mgLevelGridsCoFi,
                        const Vector<AmrCoarseners>& a_mgCoarseners);

  /**
   * @brief Check if a domain is coarser than the other
   * @param[in] a_domainOne The first domain
   * @param[in] B           The second domain
   * @return Returns true of a_domainOne has fewer grid points than B
   */
  static bool
  isCoarser(const ProblemDomain& a_domainOne, const ProblemDomain& B);

  /**
   * @brief Check if a domain is finer than the other
   * @param[in] a_domainOne The first domain
   * @param[in] B           The second domain
   * @return Returns true of a_domainOne has more grid points than B
   */
  static bool
  isFiner(const ProblemDomain& a_domainOne, const ProblemDomain& B);

  /**
   * @brief Find level corresponding to amr level
   * @param[in] a_domain Problem domain.
   * @return Depth in m_amrLevelGrids corresponding to a_domain.
   * @note Run-time error if no level was found.
   */
  int
  findAmrLevel(const ProblemDomain& a_domain) const;
};

#include <CD_NamespaceFooter.H>

#endif
