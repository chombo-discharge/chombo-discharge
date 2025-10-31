/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_FieldSolverMultigrid.H
  @brief  Declaration of FieldSolverMultigrid
  @author Robert Marskar
*/

#ifndef CD_FieldSolverMultigrid_H
#define CD_FieldSolverMultigrid_H

// Chombo includes
#include <AMRMultiGrid.H>
#include <BiCGStabSolver.H>
#include <MFSimpleSolver.H>
#include <GMRESSolver.H>

// Our includes
#include <CD_FieldSolver.H>
#include <CD_MFHelmholtzOpFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Implementation of FieldSolver which uses multigrid. 
  @details This class uses MFHelmholtzOp for discretizing the Poisson equation to second order. Although MFHelmholtzOp solves the Helmholtz equation
  alpha*a(x)*phi(x) + beta*div[b(x)*grad(phi)] = rho, we just set alpha to zero and eat the performance penalty. This class runs with electrostatic
  domain and EB boundary conditions that fit in the EBHelmholtzOp/MFHelmholtzOp constructs. 
*/
class FieldSolverMultigrid : public FieldSolver
{
public:
  /*!
    @brief Weak constructor. 
    @details This sets default options. 
  */
  FieldSolverMultigrid();

  /*!
    @brief Disallowed copy constructor. 
    @param[in] a_other Other solver
  */
  FieldSolverMultigrid(const FieldSolverMultigrid& a_other) = delete;

  /*!
    @brief Disallowed move constructor.
    @param[in] a_other Other solver
  */
  FieldSolverMultigrid(const FieldSolverMultigrid&& a_other) = delete;

  /*!
    @brief Disallowed copy assignment operator.
    @param[in] a_other Other solver
  */
  FieldSolverMultigrid&
  operator=(const FieldSolverMultigrid& a_other) = delete;

  /*!
    @brief Disallowed move assignment operator.
    @param[in] a_other Other solver
  */
  FieldSolverMultigrid&
  operator=(const FieldSolverMultigrid&& a_other) = delete;

  /*!
    @brief Constructor (does nothing)
  */
  virtual ~FieldSolverMultigrid();

  /*!
    @brief Solves Poisson equation onto a_phi using a_rho and a_sigma as right-hand sides. 
    @param[inout] a_potential Potential
    @param[in] a_rho Space charge density
    @param[in] a_sigma Surface charge density. 
    @param[in] a_zeroPhi Set a_potential to zero first. 
    @return True if we found a solution and false otherwise. 
    @note a_sigma must be defined on the gas phase. 
  */
  virtual bool
  solve(MFAMRCellData&       a_potential,
        const MFAMRCellData& a_rho,
        const EBAMRIVData&   a_sigma,
        const bool           a_zeroPhi = false) override;

  /*!
    @brief Parse all class options from command-line or input script. 
  */
  virtual void
  parseOptions() override;

  /*!
    @brief Parse runtime options from command line or input script. 
  */
  virtual void
  parseRuntimeOptions() override;

  /*!
    @brief Compute the cell-centered electric field. 
    @param[out] a_electricField Cell-centered electric field
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential) const override;

  /*!
    @brief Compute the face-centered electric field. 
    @param[out] a_electricField Face-centered electric field
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(MFAMRFluxData& a_E, const MFAMRCellData& a_potential) const override;

  /*!
    @brief Compute the cell-centered electric field on a specific phase. 
    @param[out] a_electricField Cell-centered electric field
    @param[in]  a_phase         Phase
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(EBAMRCellData&           a_E,
                       const phase::which_phase a_phase,
                       const MFAMRCellData&     a_potential) const override;

  /*!
    @brief Compute the face-centered electric field on a specific phase. 
    @param[out] a_electricField Face-centered electric field
    @param[in]  a_phase         Phase
    @param[in]  a_potential     Cell-centered potential
  */
  virtual void
  computeElectricField(EBAMRFluxData&           a_E,
                       const phase::which_phase a_phase,
                       const MFAMRCellData&     a_potential) const override;

  /*!
    @brief Cache state before regridding
    @param[in] a_lbase          Coarsest level which changes during regrid. 
    @param[in] a_oldFinestLevel Finest AMR level before regrid. 
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel) override;

  /*!
    @brief Regrid method. 
    @details This calls the parent regrid function first. We only need the overwrite because we have to signal that multigrid must be set up again. 
    @param[in] a_lmin           Coarsest level allowed to change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level before the regrid. 
    @details This linearly interpolates (with limiters) m_potential to the new grids and recomputes the electric field (from the interpolated potential). 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

  /*!
    @brief Registers required operators for AMR for this class. 
  */
  virtual void
  registerOperators() override;

  /*!
    @brief Set up multigrid solver
  */
  virtual void
  setupSolver() override;

  /*!
    @brief Set new permittivities for the multigrid solver. 
    @param[in] a_permittivityCell Permittivity on cell center
    @param[in] a_permittivityFace Permittivity on faces
    @param[in] a_permittivityEB   Permittivity on EB
  */
  virtual void
  setSolverPermittivities(const MFAMRCellData& a_permittivityCell,
                          const MFAMRFluxData& a_permittivityFace,
                          const MFAMRIVData&   a_permittivityEB) override;

  /*!
    @brief Compute grid loads
    @param[in] a_dbl   Grids on input level
    @param[in] a_level Grid level
    @return Returns a vector of estimated loads. The indices in the vector correspond to the indices in DisjointBoxLayout::boxArray()
  */
  virtual Vector<long long>
  computeLoads(const DisjointBoxLayout& a_dbl, const int a_level) override;

protected:
  /*!
    @brief alpha-coefficient (for Helmholtz operator)
  */
  static constexpr Real m_alpha = 0.0;

  /*!
    @brief Beta-coefficient (for Helmholtz operator)
  */
  static constexpr Real m_beta = -1.0;

  /*!
    @brief Enum class for supported bottom solvers in multigrid.
  */
  enum class BottomSolverType
  {
    Simple,
    BiCGStab,
    GMRES
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
    @brief Jump BC type
  */
  enum class JumpBCType
  {
    Natural,
    SaturationCharge,
  };

  /*!
    @brief Relaxation type for gmg
  */
  MFHelmholtzOp::Smoother m_multigridRelaxMethod;

  /*!
    @brief GMG multigrid type
  */
  MultigridType m_multigridType;

  /*!
    @brief Bottom solver type
  */
  BottomSolverType m_bottomSolverType;

  /*!
    @brief JumpBC type
  */
  JumpBCType m_jumpBcType;

  /*!
    @brief Scale by kappa or not
  */
  bool m_kappaSource;

  /*!
    @brief Needs setup
  */
  bool m_isSolverSetup;

  /*!
    @brief Drop order in coarsened grid cells
  */
  bool m_multigridDropOrder;

  /*!
    @brief Use default multigrid settings
  */
  bool m_multigridUseDefaultSettings;

  /*!
    @brief Verbosity for geometric multigrid
  */
  int m_multigridVerbosity;

  /*!
    @brief Number of smoothing for multigrid
  */
  int m_multigridSmooth;

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
    @brief Number of smoothing for bottom solver
  */
  int m_numSmoothingsForSimpleSolver;

  /*!
    @brief Set bottom drop depth
  */
  int m_minCellsBottom;

  /*!
    @brief Domain drop order
  */
  int m_domainDropOrder;

  /*!
    @brief Set bc order 
  */
  int m_multigridBcOrder;

  /*!
    @brief Set bc weight (used for weighted least squares)
  */
  int m_multigridBcWeight;

  /*!
    @brief Set bc jump order 
  */
  int m_multigridJumpOrder;

  /*!
    @brief Set bc jump weight (used for weighted least squares)
  */
  int m_multigridJumpWeight;

  /*!
    @brief Multigrid smoothing steps for preconditioner
  */
  int m_multigridPreCondSmooth;

  /*!
    @brief Exit tolerance for multigrid. 
    @details Multigrid exits if L(phi) < tolerance*L(phi=0)
  */
  Real m_multigridExitTolerance;

  /*!
    @brief Exit hang for multigrid
    @details Multigrid exits if residue is not reduce by at least this factor. 
  */
  Real m_multigridExitHang;

  /*!
    @brief Multigrid relaxation factor
  */
  Real m_multigridRelaxFactor;

  /*!
    @brief Multigrid operator factory. 
  */
  RefCountedPtr<MFHelmholtzOpFactory> m_helmholtzOpFactory;

  /*!
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<MFCellFAB>>> m_multigridSolver;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<MFCellFAB>> m_bicgstab;

  /*!
    @brief GMRES solver
  */
  GMRESSolver<LevelData<MFCellFAB>> m_gmres;

  /*!
    @brief multi-fluid simple solver
  */
  MFSimpleSolver m_mfsolver;

  /*!
    @brief Parse multigrid settings
  */
  virtual void
  parseMultigridSettings();

  /*!
    @brief Parse m_kappaSource
  */
  virtual void
  parseKappaSource();

  /*!
    @brief Parse JumpBC type
  */
  virtual void
  parseJumpBC();

  /*!
    @brief Set the permittivities
    @details This sets m_permittivityCell, m_permittivityFace, and m_permittivityEB and fills
    necessary ghost faces/cells. 
  */
  virtual void
  setPermittivities() override;

  /*!
    @brief Set up the multigrid operator factory
  */
  virtual void
  setupHelmholtzFactory();

  /*!
    @brief Set up AMRMultiGrid
  */
  virtual void
  setupMultigrid();
};

#include <CD_NamespaceFooter.H>

#endif
