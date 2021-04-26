/*!
  @file poisson_multifluid_gmg.H
  @brief Geometric multigrid multifluid poisson solver
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _POISSON_MULTIFLUID_GMG_
#define _POISSON_MULTIFLUID_GMG_

#include "poisson_solver.H"
#include "mfconductivityopfactory.H"
#include "dirichlet_func.H"

#include <AMRMultiGrid.H>
#include <BiCGStabSolver.H>
#include <EBConductivityOpFactory.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <MFSimpleSolver.H>
#include <GMRESSolver.H>
#include <BaseBCFuncEval.H>

/*!
  @brief Implementation of poisson_solver which uses a two-phase geometric multigrid solver
*/
class poisson_multifluid_gmg : public poisson_solver {
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

  RefCountedPtr<mfconductivityopfactory>& get_opfact(){
    return m_opfact;
  }

  /*!
    @brief Constructor
  */
  poisson_multifluid_gmg();

  /*!
    @brief Constructor
  */
  virtual ~poisson_multifluid_gmg();

  /*!
    @brief Override this one for testing purposes. 
  */
  virtual bool solve(const bool a_zerophi = false);

  /*!
    @brief Solve Poisson onto state 
  */
  virtual bool solve(MFAMRCellData&       a_state,
		     const MFAMRCellData& a_rho,
		     const EBAMRIVData&   a_sigma,
		     const bool           a_zerophi = false);

  /*!
    @brief parse options
  */
  virtual void parse_options() override;

    /*!
    @brief Parse runtime options
  */
  virtual void parse_runtime_options() override;

  /*!
    @brief Automatic tuning of the Poisson solver
  */
  virtual void auto_tune();

  /*!
    @brief Set potential
  */
  virtual void set_potential(Real (*a_potential)(const Real a_time));

  /*!
    @brief Regrid function
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Register operators for AMR
  */
  virtual void register_operators();

  /*!
    @brief Set solver parameters. These are normally good for this class. 
  */
  virtual void set_gmg_solver_parameters(relax      a_relax_type  = relax::gsrb_fast,
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
  virtual void set_bottom_solver(const int a_which_solver);

  /*!
    @brief Set the depth for bottom drop
  */
  virtual void set_bottom_drop(const int a_bottom_drop);

  /*!
    @brief Set number of smoothing for the bottom solver (simplified solver only)
  */
  virtual void set_botsolver_smooth(const int a_smooth);

    /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals();

  /*!
    @brief Post checkpoint stuff. 
  */
  virtual void post_checkpoint() override;

  /*!
    @brief Get a-coefficient
  */
  virtual MFAMRCellData& get_aco();

  /*!
    @brief Get b-coefficient
  */
  virtual MFAMRFluxData& get_bco();

  /*!
    @brief Get irregular b coefficient
  */
  virtual MFAMRIVData& get_bco_irreg();

  /*!
    @brief Return require number of ghost cells. 
  */
  virtual int query_ghost() const;

  /*!
    @brief Set up the geometric multigrid solver 
  */
  virtual void setup_gmg();
  
  /*!
    @brief Set coefficients
  */
  virtual void set_coefficients();

  /*!
    @brief Set up the operator factory
  */
  virtual void setup_operator_factory();

  /*!
    @brief Set up AMRMultiGrid
  */
  virtual void setup_solver();

  /*!
    @brief Check if solver needs setup
  */
  virtual void set_needs_setup(const bool& a_needs_setup);

  /*!
    @brief Set the permittivities
  */
  virtual void set_permittivities(const Vector<dielectric>& a_dielectrics);
  
protected:

  /*!
    @brief Relaxation type for gmg
  */
  relax m_gmg_relax_type;

  /*!
    @brief GMG multigrid type
  */
  amrmg m_gmg_type;

  /*!
    @brief Returns one. Used for mapping boundary conditions
  */
  static Real s_constant_one(const RealVect a_pos);

  /*!
    @brief Scale by kappa or not
  */
  bool m_kappa_source;

  /*!
    @brief Needs setup
  */
  bool m_needs_setup;

  /*!
    @brief
  */
  bool m_has_mg_stuff;

  /*!
    @brief Number of coarsening levels for multigrid
  */
  int m_mg_coarsen;

  /*!
    @brief Verbosity for geometric multigrid
  */
  int m_gmg_verbosity;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_gmg_pre_smooth;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_gmg_post_smooth;

  /*!
    @brief Number of smoothing before bottom solver
  */
  int m_gmg_bot_smooth;

  /*!
    @brief Maximum number of iterations
  */
  int m_gmg_max_iter;

  /*!
    @brief Minimum number of iterations
  */
  int m_gmg_min_iter;

  /*!
    @brief Bottom solver type
  */
  int m_bottomsolver;

  /*!
    @brief Number of smoothing for bottom solver
  */
  int m_numsmooth;

  /*!
    @brief Set bottom drop depth
  */
  int m_bottom_drop;

  /*!
    @brief Set bc order 
  */
  int m_bc_order;

  /*!
    @brief Length scale for the Poisson equation
  */
  Real m_length_scale;

  /*!
    @brief Target residue. If they residue is greater than this, we have a problematic cell. 
  */
  Real m_converged_resid;

  /*!
    @brief 
  */
  Real m_gmg_eps;
  
  /*!
    @brief 
  */
  Real m_gmg_hang;

  /*!
    @brief 
  */
  Real m_gmg_norm_thresh;

  /*!
    @brief Operator factory
  */
  RefCountedPtr<mfconductivityopfactory> m_opfact;
  
  /*!
    @brief b-coefficient
  */
  MFAMRCellData m_aco;

  /*!
    @brief b-coefficient
  */
  MFAMRFluxData m_bco;

  /*!
    @brief b-coefficient
  */
  MFAMRIVData m_bco_irreg;

  /*!
    @brief Zero
  */
  MFAMRCellData m_zero;

  /*!
    @brief Scaled soruce
  */
  MFAMRCellData m_scaled_source;

  /*!
    @brief Sigma
  */
  EBAMRIVData m_scaled_sigma;

  /*!
    @brief Geometric multigrid solver
  */
  AMRMultiGrid<LevelData<MFCellFAB> > m_gmg_solver;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<MFCellFAB> > m_bicgstab;

  /*!
    @brief GMRES solver
  */
  GMRESSolver<LevelData<MFCellFAB> > m_gmres;

  /*!
    @brief Geometric multigrid solver
  */
  AMRMultiGrid<LevelData<EBCellFAB> > m_cond_solver;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB> > m_cond_bicgstab;

  /*!
    @brief multi-fluid simple solver
  */
  MFSimpleSolver m_mfsolver;

  /*!
    @brief Potential function
  */
  RefCountedPtr<dirichlet_func> m_bcfunc;

  /*!
    @brief Mapped potential functions
  */
  Vector<RefCountedPtr<BaseBCFuncEval> > m_wall_bcfunc;

  /*!
    @brief Coarsened MG grids
  */
  Vector<DisjointBoxLayout> m_mg_grids;

  /*!
    @brief Coarsened problem domains
  */
  Vector<ProblemDomain> m_mg_domains;

  /*!
    @brief Coarsened EBLevelGrids
  */
  Vector<Vector<RefCountedPtr<EBLevelGrid> > > m_mg_eblg;

  /*!
    @brief EBLevelgrids for the multigrid pre-coarsening
  */
  Vector<RefCountedPtr<MFLevelGrid> > m_mg_mflg;



  /*!
    @brief Need to override this in order to pass in the correct time to boundary conditions
  */
  virtual void set_time(const int a_step, const Real a_time, const Real a_dt);


  /*!
    @brief Set face-centered permittivities
  */
  virtual void set_face_perm(EBFluxFAB&                a_perm,
			     const Box&                a_box,
			     const RealVect&           a_origin,
			     const Real&               a_dx,
			     const Vector<dielectric>& a_dielectrics);

  /*!
    @brief Set EB-centered permittivities
  */
  virtual void set_eb_perm(BaseIVFAB<Real>&          a_perm,
			   const Box&                a_box,
			   const RealVect&           a_origin,
			   const Real&               a_dx,
			   const Vector<dielectric>& a_dielectrics);



  // Option parsing functions
  virtual void parse_autotune();     // Parse autotuning
  virtual void parse_domain_bc();    // Parse domain boundary conditions
  virtual void parse_plot_vars();    // Parse plot variables
  virtual void parse_gmg_settings(); // Parse geometric multigrid settings
  virtual void parse_kappa_source(); // Parse source scaling

  // Define mg levels
public:
  virtual void define_mg_levels();

};
#endif
