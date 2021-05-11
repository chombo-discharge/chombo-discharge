/*!
  @file FieldSolverMultigrid.H
  @brief Geometric multigrid multifluid poisson solver
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _FIELD_SOLVER_MULTIGRID_
#define _FIELD_SOLVER_MULTIGRID_

#include "CD_FieldSolver.H"
#include "mfconductivityopfactory.H"

#include <AMRMultiGrid.H>
#include <BiCGStabSolver.H>
#include <EBConductivityOpFactory.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <MFSimpleSolver.H>
#include <GMRESSolver.H>
#include <BaseBCFuncEval.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Implementation of FieldSolver which uses a two-phase geometric multigrid solver
*/
class FieldSolverMultigrid : public FieldSolver {
protected:

  /*!
    @brief Silly enum for switching between relaxation methods.
  */
  enum class relax{
    jacobi,
    gauss_seidel,
    gsrb_fast
  };

  /*!
    @brief Silly enum for switching between multigrid cycling
  */
  enum class amrmg{
    full,
    vcycle,
    fcycle,
  };
  
public:

  /*!
    @brief Added because of CD_FieldStepper, which uses the operator factory to generate computational loads.
  */
  RefCountedPtr<mfconductivityopfactory>& getOperatorFactory(){
    return m_opfact;
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
    @brief Solve Poisson onto state 
  */
  virtual bool solve(MFAMRCellData&       a_state,
		     const MFAMRCellData& a_rho,
		     const EBAMRIVData&   a_sigma,
		     const bool           a_zerophi = false) override;

  /*!
    @brief parse options
  */
  virtual void parseOptions() override;

  /*!
    @brief Parse runtime options
  */
  virtual void parseRuntimeOptions() override;

  /*!
    @brief Regrid function
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;

  /*!
    @brief Register operators for AMR
  */
  virtual void registerOperators() override;

  /*!
    @brief Set solver parameters. These are normally good for this class. 
  */
  virtual void setMultigridSettings(relax      a_relax_type  = relax::gsrb_fast,
					 amrmg      a_gmg_type    = amrmg::vcycle,
					 const int  a_verbosity   = 100,
					 const int  a_pre_smooth  = 16,
					 const int  a_post_smooth = 16,
					 const int  a_bot_smooth  = 16,
					 const int  a_max_iter    = 32,
					 const int  a_min_iter    = 5,
					 const Real a_eps         = 1.E-6,
					 const Real a_hang        = 0.2);

  /*!
    @brief Set the bottom solver.  
    @param[in] a_which_solver Solver type. 0 = smoother, 1 = BiCGStab
  */
  virtual void setBottomSolver(const int a_which_solver);

  /*!
    @brief Set the depth for bottom drop
  */
  virtual void setBottomDrop(const int a_bottom_drop);

  /*!
    @brief Set number of smoothing for the bottom solver (simplified solver only)
  */
  virtual void setSmoothingsForBottomSolver(const int a_smooth);

  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals() override;

  /*!
    @brief Post checkpoint stuff. 
  */
  virtual void postCheckpoint() override;

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

  /*!
    @brief Set up the geometric multigrid solver 
  */
  virtual void setupMultigrid();
  
  /*!
    @brief Set coefficients
  */
  virtual void setMultigridCoefficients();

  /*!
    @brief Set up the operator factory
  */
  virtual void setupOperatorFactory();

  /*!
    @brief Set up AMRMultiGrid
  */
  virtual void setupMultigridSolver();

  /*!
    @brief Check if solver needs setup
  */
  virtual void setNeedsMultigridSetup(const bool& a_needs_setup);

  /*!
    @brief Set the permittivities
  */
  virtual void setPermittivities(const Vector<dielectric>& a_dielectrics);
  
protected:

  /*!
    @brief Relaxation type for gmg
  */
  relax m_multigridRelaxMethod;

  /*!
    @brief GMG multigrid type
  */
  amrmg m_multigridType;

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
  int m_mg_coarsen;

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
  int m_bottomSolver;

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
    @brief Operator factory
  */
  RefCountedPtr<mfconductivityopfactory> m_opfact;
  
  /*!
    @brief b-coefficient
  */
  MFAMRCellData m_aCoefficient;

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
			     const Vector<dielectric>& a_dielectrics);

  /*!
    @brief Set EB-centered permittivities
  */
  virtual void setEbPermittivities(BaseIVFAB<Real>&          a_perm,
			   const Box&                a_box,
			   const RealVect&           a_origin,
			   const Real&               a_dx,
			   const Vector<dielectric>& a_dielectrics);

  /*!
    @brief Parse multigrid settings
  */
  virtual void parseMultigridSettings();

  /*!
    @brief Parse m_kappaSource
  */
  virtual void parseKappaSource(); 

  // Define mg levels
public:
  virtual void defineDeeperMultigridLevels();

};
#include "CD_NamespaceFooter.H"
#endif
