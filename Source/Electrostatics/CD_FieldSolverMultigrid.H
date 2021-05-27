/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   FieldSolverMultigrid.H
  @brief  Implementation of FieldSolver which uses geometric multigrid.
  @author Robert Marskar
*/

#ifndef CD_FieldSolverMultigrid_H
#define CD_FieldSolverMultigrid_H

// Chombo includes
#include <AMRMultiGrid.H>
#include <BiCGStabSolver.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <MFSimpleSolver.H>
#include <GMRESSolver.H>
#include <BaseBCFuncEval.H>

// Our includes
#include <CD_FieldSolver.H>
#include <CD_MfHelmholtzOpFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Implementation of FieldSolver which uses a two-phase geometric multigrid solver
*/
class FieldSolverMultigrid : public FieldSolver {
protected:

  /*!
    @brief Enum class for supported bottom solvers in multigrid.
  */
  enum class BottomSolver {
    Simple,
    BiCGStab,
    GMRES
  };

  /*!
    @brief Enum for supported relaxation methods in multigrid. 
  */
  enum class RelaxationMethod {
    Jacobi,
    GaussSeidel,
    GSRBFast
  };

  /*!
    @brief Enum for multigrid cycle types. 
  */
  enum class MultigridType{
    FAS,
    VCycle,
    FCycle,
  };
  
public:

  /*!
    @brief Added because of CD_FieldStepper, which uses the operator factory to generate computational loads.
  */
  RefCountedPtr<MfHelmholtzOpFactory>& getOperatorFactory(){
    return m_operatorFactory;
  }

  /*!
    @brief Constructor
  */
  FieldSolverMultigrid();

  /*!
    @brief Constructor
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
  virtual bool solve(MFAMRCellData&       a_potential,
		     const MFAMRCellData& a_rho,
		     const EBAMRIVData&   a_sigma,
		     const bool           a_zerophi = false) override;

  /*!
    @brief Parse all class options from command-line or input script. 
  */
  virtual void parseOptions() override;

  /*!
    @brief Parse runtime options from command line or input script. 
  */
  virtual void parseRuntimeOptions() override;

    /*!
    @brief Regrid method. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Previous finest AMR level
    @param[in] a_newFinestLevel New finest AMR level. 
    @details This linearly interpolates (with limiters) m_potential to the new grids and recomputes the electric field (from the interpolated potential). 
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

  /*!
    @brief Registers required operators for AMR for this class. 
  */
  virtual void registerOperators() override;

  /*!
    @brief Allocates internal storage. 
  */
  virtual void allocateInternals() override;

  /*!
    @brief Set up the geometric multigrid solver 
  */
  virtual void setupMultigrid();
  
  /*!
    @brief Set multigrid coefficients. 
    @details This sets the face-centered permittivities and the diagonal term (to zero). 
  */
  virtual void setMultigridCoefficients();

  /*!
    @brief Set up the multigrid operator factory
  */
  virtual void setupOperatorFactory();

  /*!
    @brief Set up AMRMultiGrid
  */
  virtual void setupMultigridSolver();

  /*!
    @brief Set whether multigrid will be set up again before the next solve. 
    @param[in] a_needsSetup If true, multigrid will be set up before the next solve. 
  */
  virtual void setNeedsMultigridSetup(const bool a_needsSetup);

    /*!
    @brief Get a-coefficient
  */
  virtual MFAMRCellData& getACoefficient();

  /*!
    @brief Get b-coefficient
  */
  virtual MFAMRFluxData& getBCoefficient();

  /*!
    @brief Get irregular b coefficient
  */
  virtual MFAMRIVData& getBCoefficientIrreg();
  
protected:

  /*!
    @brief Relaxation type for gmg
  */
  RelaxationMethod m_multigridRelaxMethod;

  /*!
    @brief GMG multigrid type
  */
  MultigridType m_multigridType;

  /*!
    @brief Scale by kappa or not
  */
  bool m_kappaSource;

  /*!
    @brief Needs setup
  */
  bool m_needsMultigridSetup;

  /*!
    @brief
  */
  bool m_hasDeeperMultigridLevels;

  /*!
    @brief Number of coarsening levels for multigrid
  */
  int m_numCoarseningsBeforeAggregation;

  /*!
    @brief Verbosity for geometric multigrid
  */
  int m_multigridVerbosity;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_multigridPreSmoothg;

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
    @brief Bottom solver type
  */
  //  int m_bottomSolver;

  BottomSolver m_bottomSolver;

  /*!
    @brief Number of smoothing for bottom solver
  */
  int m_numSmoothingsForSimpleSolver;

  /*!
    @brief Set bottom drop depth
  */
  int m_numCellsBottomDrop;

  /*!
    @brief Set bc order 
  */
  int m_multigridBcOrder;

  /*!
    @brief Length scale for the Poisson equation
  */
  Real m_lengthScale;

  /*!
    @brief Target residue. If they residue is greater than this, we have a problematic cell. 
  */
  Real m_convergedResidue;

  /*!
    @brief 
  */
  Real m_multigridTolerance;
  
  /*!
    @brief 
  */
  Real m_multigridHang;

  /*!
    @brief Multigrid operator factory. 
  */
  RefCountedPtr<MfHelmholtzOpFactory> m_operatorFactory;
  
  /*!
    @brief b-coefficient
  */
  MFAMRCellData m_aCoef;

  /*!
    @brief b-coefficient
  */
  MFAMRFluxData m_bCoefficient;

  /*!
    @brief b-coefficient
  */
  MFAMRIVData m_bCoefficientIrreg;

  /*!
    @brief Zero
  */
  MFAMRCellData m_zero;

  /*!
    @brief Scaled soruce
  */
  MFAMRCellData m_scaledSource;

  /*!
    @brief Sigma
  */
  EBAMRIVData m_scaledSigma;

  /*!
    @brief Geometric multigrid solver
  */
  AMRMultiGrid<LevelData<MFCellFAB> > m_multigridSolver;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<MFCellFAB> > m_bicgstab;

  /*!
    @brief GMRES solver
  */
  GMRESSolver<LevelData<MFCellFAB> > m_gmres;

  /*!
    @brief multi-fluid simple solver
  */
  MFSimpleSolver m_mfsolver;

  /*!
    @brief Coarsened MG grids
  */
  Vector<DisjointBoxLayout> m_deeperMultigridGrids;

  /*!
    @brief Coarsened problem domains
  */
  Vector<ProblemDomain> m_deeperMultigridDomains;

  /*!
    @brief Coarsened EBLevelGrids
  */
  Vector<Vector<RefCountedPtr<EBLevelGrid> > > m_mg_eblg;

  /*!
    @brief EBLevelgrids for the multigrid pre-coarsening
  */
  Vector<RefCountedPtr<MFLevelGrid> > m_mg_mflg;

  /*!
    @brief Set face-centered permittivities
  */
  virtual void setFacePermittivities(EBFluxFAB&                a_perm,
				     const Box&                a_box,
				     const RealVect&           a_origin,
				     const Real&               a_dx,
				     const Vector<Dielectric>& a_dielectrics);

  /*!
    @brief Set EB-centered permittivities
  */
  virtual void setEbPermittivities(BaseIVFAB<Real>&          a_perm,
				   const Box&                a_box,
				   const RealVect&           a_origin,
				   const Real&               a_dx,
				   const Vector<Dielectric>& a_dielectrics);

  /*!
    @brief Parse multigrid settings
  */
  virtual void parseMultigridSettings();

  /*!
    @brief Parse m_kappaSource
  */
  virtual void parseKappaSource();

    /*!
    @brief Set the permittivities
  */
  virtual void setPermittivities(const Vector<Dielectric>& a_dielectrics);
  
  /*!
    @brief Define the deeper multigrid levels (the ones that live below the AMR levels)
  */
  virtual void defineDeeperMultigridLevels();
  
};

#include <CD_NamespaceFooter.H>

#endif
