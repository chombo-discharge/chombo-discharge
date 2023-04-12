/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EddingtonSP1.H
  @brief  Declaration of a first order Eddington solver for radiative transfer
  @author Robert Marskar
*/

#ifndef CD_EddingtonSP1_H
#define CD_EddingtonSP1_H

// Std includes
#include <random>

// Chombo includes
#include <AMRMultiGrid.H>
#include <AMRTGA.H>
#include <EBBackwardEuler.H>
#include <EBSimpleSolver.H>
#include <BiCGStabSolver.H>
#include <GMRESSolver.H>

// Our includes
#include <CD_RtSolver.H>
#include <CD_EBHelmholtzOpFactory.H>
#include <CD_EddingtonSP1DomainBc.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Radiative tranfer equation solver in the SP1 (diffusion) approximation. 
*/
class EddingtonSP1 : public RtSolver
{
public:
  /*!
    @brief Default function for space-time dependence of domain boundary conditions. 
    @details This is a utility function 
    @param[in] a_position Position in space
    @param[in] a_time     Time
    @return Returns 1.0
  */
  static Real
  s_defaultDomainBcFunction(const RealVect a_position, const Real a_time);

  /*!
    @brief Weak constructor
  */
  EddingtonSP1();

  /*!
    @brief Disallowed copy construction
  */
  EddingtonSP1(const EddingtonSP1& a_other) = delete;

  /*!
    @brief Disallowed move construction
  */
  EddingtonSP1(const EddingtonSP1&& a_other) = delete;

  /*!
    @brief Destructor
  */
  virtual ~EddingtonSP1();

  /*!
    @brief Disallowed copy assignment
  */
  EddingtonSP1&
  operator=(const EddingtonSP1& a_other) = delete;

  /*!
    @brief Disallowed move assignment
  */
  EddingtonSP1&
  operator=(const EddingtonSP1&& a_other) = delete;

  /*!
    @brief Advance RTE onto state a_phi
    @param[in]    a_dt      Time step
    @param[inout] a_phi     RTE solution
    @param[in]    a_source  Source term
    @param[in]    a_zeroPhi Set phi to zero in initial guess for multigrid solve
  */
  virtual bool
  advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zeroPhi = false) override;

  /*!
    @brief Parse class options
  */
  virtual void
  parseOptions() override;

  /*!
    @brief Parse class options
  */
  virtual void
  parseRuntimeOptions() override;

  /*!
    @brief Set the boundary condition function on a domain side
    @param[in] a_dir      Coordinate direction.
    @param[in] a_side     Side (low/high)
    @param[in] a_function Boundary condition function.
    @details This sets a boundary condition for a particular domain side. The user must also specify how to use this BC in the input script. 
  */
  virtual void
  setDomainSideBcFunction(const int                               a_dir,
                          const Side::LoHiSide                    a_side,
                          const EddingtonSP1DomainBc::BcFunction& a_function);

  /*!
    @brief Allocate internal storage
  */
  virtual void
  allocate() override;

  /*!
    @brief Cache state
    @param[in] a_base           Coarsest level that did not change during regrid
    @param[in] a_oldFinestLevel Finest level before the regrid
  */
  virtual void
  preRegrid(const int a_base, const int a_oldFinestLevel) override;

  /*!
    @brief Deallocate internal storage
  */
  virtual void
  deallocate() override;

  /*!
    @brief Regrid function for this class
    @param[in] a_lmin           Coarsest level that changed during regrid
    @param[in] a_oldFinestLevel Finest level before the regrid
    @param[in] a_newFinestLevel Finest level after the regrid
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

  /*!
    @brief Register operators
  */
  virtual void
  registerOperators() override;

  /*!
    @brief Compute the boundary flux. For Eddington, the boundary flux is = c*phi/2
    @param[out] a_ebFlux Flux on the EBs
    @param[in]  a_phi    Isotropic part of the RTE solution
  */
  virtual void
  computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi) override;

  /*!
    @brief Compute the domain flux. For Eddington, the domain flux is = c*phi/2
    @param[out] a_domainFlux Flux on the domain sides
    @param[in]  a_phi        Isotropic part of the RTE solution
  */
  virtual void
  computeDomainFlux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_phi) override;

  /*!
    @brief Compute the flux. For Eddington, the flux is F = -c/(3*kappa)*grad(phi)
    @param[out] a_flux Cell-centered flux
    @param[in]  a_phi  Isotropic part of RTE solution
  */
  virtual void
  computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi) override;

  /*!
    @brief Get isotropic part. 
    @param[out] a_isotropic Isotropic part of RTE solution
    @param[in]  a_phi       RTE solution
    @details For Eddington, only the isotropic part is solved for so this routine just copies. 
  */
  virtual void
  computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi) override;

  /*!
    @brief Write plot file
  */
  virtual void
  writePlotFile() override;

  /*!
    @brief Write checkpoint data into HDF5 file
    @param[out] a_handle HDF5 file
    @param[in]  a_level  Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const override;
#endif

  /*!
    @brief Read checkpoint data from handle
    @param[out] a_handle HDF5 file
    @param[in]  a_level  Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level) override;
#endif

protected:
  /*!
    @brief alpha-coefficient for Helmholtz operator
  */
  static constexpr Real m_alpha = 1.0;

  /*!
    @brief Beta-coefficient for Helmholtz operator
  */
  static constexpr Real m_beta = -1.0;

  /*!
    @brief Enum class for supported bottom solvers in multigrid.
  */
  enum class BottomSolverType
  {
    Simple,
    BiCGStab,
    GMRES,
  };

  /*!
    @brief Enum for multigrid cycle types. 
  */
  enum class MultigridType
  {
    VCycle,
    WCycle,
  };

  /*!
    @brief Enum for boundary condition types on EBs
  */
  enum class EBBCType
  {
    Dirichlet,
    Neumann,
    Larsen
  };

  /*!
    @brief Relaxation type for gmg
  */
  EBHelmholtzOp::Smoother m_multigridRelaxMethod;

  /*!
    @brief GMG multigrid type
  */
  MultigridType m_multigridType;

  /*!
    @brief Needs setup
  */
  bool m_isSolverSetup;

  /*!
    @brief Use tga for transient solevs
  */
  bool m_useTGA;

  /*!
    @brief Use kappa scaling for source or not
  */
  bool m_kappaScale;

  /*!
    @brief Use slopes when regridding (or not)
  */
  bool m_regridSlopes;

  /*!
    @brief Verbosity for geometric multigrid
  */
  int m_multigridVerbosity;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_multigridPreSmooth;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_multigridPostSmooth;

  /*!
    @brief Number of smoothing before bottom solver
  */
  int m_multigridBottomSmooth;

  /*!
    @brief Maximum number of iterations
  */
  int m_multigridMaxIterations;

  /*!
    @brief Minimum number of iterations
  */
  int m_multigridMinIterations;

  /*!
    @brief Multigrid EBBC order (only relevant for Dirichlet)
  */
  int m_multigridBcOrder;

  /*!
    @brief Multigrid EBBC weight (only relevant for Dirichlet)
  */
  int m_multigridBcWeight;

  /*!
    @brief Bottom solver type
  */
  BottomSolverType m_bottomSolverType;

  /*!
    @brief Number of smoothing for bottom solver
  */
  int m_numSmoothingsForSimpleSolver;

  /*!
    @brief Set bottom drop depth
  */
  int m_minCellsBottom;

  /*!
    @brief 
  */
  Real m_multigridExitTolerance;

  /*!
    @brief 
  */
  Real m_multigridExitHang;

  /*!
    @brief Reflection coefficient
  */
  Real m_reflectCoefOne;

  /*!
    @brief Reflection coefficient
  */
  Real m_reflectCoefTwo;

  /*!
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB>>> m_multigridSolver;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB>> m_bicgstab;

  /*!
    @brief GMRES solver
  */
  GMRESSolver<LevelData<EBCellFAB>> m_gmres;

  /*!
    @brief TGA solver
  */
  RefCountedPtr<AMRTGA<LevelData<EBCellFAB>>> m_tgaSolver;

  /*!
    @brief Euler solver
  */
  RefCountedPtr<EBBackwardEuler> m_eulerSolver;

  /*!
    @brief Operator factory
  */
  RefCountedPtr<EBHelmholtzOpFactory> m_helmholtzOpFactory;

  /*!
    @brief multi-fluid simple solver
  */
  EBSimpleSolver m_simpleSolver;

  /*!
    @brief For regridding the source term. This is needed when doing a stationary solve. 
  */
  EBAMRCellData m_cacheSrc;

  /*!
    @brief Multigrid residue
  */
  EBAMRCellData m_resid;

  /*!
    @brief Always zero
  */
  EBAMRCellData m_zero;

  /*!
    @brief Scaled source term
  */
  EBAMRCellData m_scaledSource;

  /*!
    @brief a-coefficient
  */
  EBAMRCellData m_helmAco;

  /*!
    @brief b-coefficient
  */
  EBAMRFluxData m_helmBco;

  /*!
    @brief b-coefficient
  */
  EBAMRIVData m_helmBcoIrreg;

  /*!
    @brief Wrapper calss. 
  */
  EddingtonSP1DomainBc m_domainBc;

  /*!
    @brief Actual functions on domain edges (faces). 
  */
  std::map<EddingtonSP1DomainBc::DomainSide, EddingtonSP1DomainBc::BcFunction> m_domainBcFunctions;

  /*!
    @brief Associated boundary condition on the embedded boundary
  */
  std::pair<EBBCType, Real> m_ebbc;

  /*!
    @brief Set up geometric multigrid
  */
  virtual void
  setupSolver();

  /*!
    @brief Set multigrid coefficients
  */
  virtual void
  setHelmholtzCoefficients();

  /*!
    @brief Set EBHelmholtzOp A- and B-coefficients. 
    @details For the B-coefficient, this also fills one of the tangential ghost faces. 
    @param[out] a_helmAco      Helmholtz A-coefficient on on cell centers
    @param[out] a_helmBco      Helmholtz B-coefficient on on face centers
    @param[out] a_helmBcoIrreg Helmholtz B-coefficient on EB faces
    @param[in]  a_lvl          Grid level
    @param[in]  a_dit          Grid index
  */
  virtual void
  setHelmholtzCoefficientsBox(EBCellFAB&       a_helmAco,
                              EBFluxFAB&       a_helmBco,
                              BaseIVFAB<Real>& a_helmBcoIrreg,
                              const int        a_lvl,
                              const DataIndex& a_dit);

  /*!
    @brief Set the operator factory
  */
  virtual void
  setupHelmholtzFactory();

  /*!
    @brief Set the GMG solver
  */
  virtual void
  setupMultigrid();

  /*!
    @brief Make domain bc string
    @param[in] a_dir  Coordinate direction
    @param[in] a_side Side
    @details Used for pasing boundary conditions from input script. 
    @returns Returns string of type "bc.x.lo", "bc.y.hi" etc. 
  */
  virtual std::string
  makeBcString(const int a_dir, const Side::LoHiSide a_side) const;

  /*!
    @brief Returns BC type based on string.
    @param[in] a_str Boundary condition string. Must be "dirichlet <number>", "neumann <number>", "robin <number>", "dirichlet_custom", "neumann_custom", or "robin_custom".
  */
  virtual EddingtonSP1DomainBc::BcType
  parseBcString(const std::string a_str) const;

  /*!
    @brief Set up the TGA solver
  */
  virtual void
  setupTGA();

  /*!
    @brief Set up the Euler solver
  */
  virtual void
  setupEuler();

  /*!
    @brief Parse domain BC settings
  */
  virtual void
  parseDomainBC();

  /*!
    @brief Parse domain BC settings
  */
  virtual void
  parseEBBC();

  /*!
    @brief Parse whether or not this is a stationary solver
  */
  virtual void
  parseStationary();

  /*!
    @brief Parse plot variables
  */
  virtual void
  parsePlotVariables();

  /*!
    @brief Parse multigrid settings
  */
  virtual void
  parseMultigridSettings();

  /*!
    @brief Parse kappa-scaling or not
  */
  virtual void
  parseKappaScale();

  /*!
    @brief Parse reflection coefficients for Robin bcs
  */
  virtual void
  parseReflection();

  /*!
    @brief Parse whether or not to use slopes when regridding
  */
  virtual void
  parseRegridSlopes();

  /*!
    @brief Set default domain BC functions.
  */
  virtual void
  setDefaultDomainBcFunctions();
};

#include <CD_NamespaceFooter.H>

#endif
