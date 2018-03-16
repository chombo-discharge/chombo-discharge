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
#include "nwomfconductivityopfactory.H"

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
  @details 
  
  Class options
  -------------

      poisson_multifluid.gmg_verbosity     = -1        # GMG verbosity
      poisson_multifluid.gmg_pre_smooth    = 16        # Number of relaxations in downsweep
      poisson_multifluid.gmg_post_smooth   = 16        # Number of relaxations in upsweep
      poisson_multifluid.gmg_bott_smooth   = 16        # NUmber of relaxations before dropping to bottom solver
      poisson_multifluid.gmg_min_iter      = 5         # Minimum number of iterations
      poisson_multifluid.gmg_max_iter      = 32        # Maximum number of iterations
      poisson_multifluid.gmg_tolerance     = 1.E-6     # Residue tolerance
      poisson_multifluid.gmg_hang          = 0.2       # Solver hang
      poisson_multifluid.gmg_bottom_drop   = 8         # Bottom drop
      poisson_multifluid.gmg_bottom_solver = simple    # Bottom solver type. Valid options are 'simple' and 'bicgstab'
      poisson_multifluid.gmg_bottom_relax  = 32        # Number of relaxations in bottom solve ('simple' solver only)
      
*/
class poisson_multifluid_gmg : public poisson_solver {
public:

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
    @brief Automatic tuning of the Poisson solver
  */
  virtual void auto_tune();

  /*!
    @brief Regrid function
  */
  virtual void regrid(const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Set solver parameters. These are normally good for this class. 
  */
  virtual void set_gmg_solver_parameters(relax::which_relax a_relax_type = relax::gsrb_fast,
					 amrmg::which_mg a_gmg_type      = amrmg::vcycle,
					 const int a_verbosity           = 100,
					 const int a_pre_smooth          = 16,
					 const int a_post_smooth         = 16,
					 const int a_bot_smooth          = 16,
					 const int a_max_iter            = 32,
					 const int a_min_iter            = 5,
					 const Real a_eps                = 1.E-6,
					 const Real a_hang               = 0.2);

  /*!
    @brief Set the bottom solver.  
    @param[in] a_which_solver Solver type. 0 = smoother, 1 = BiCGStab
  */
  virtual void set_bottom_solver(const int a_whichsolver);

  /*!
    @brief Set the depth for bottom drop
  */
  virtual void set_bottom_drop(const int a_bottom_drop);

  /*!
    @brief Set number of smoothing for the bottom solver (simplified solver only)
  */
  virtual void set_botsolver_smooth(const int a_smooth);

  /*!
    @brief Use nwo factory or not
  */
  virtual void set_nwo(const bool a_use_nwo);

  /*!
    @brief Return require number of ghost cells. 
  */
  virtual int query_ghost() const;
  
protected:

  /*!
    @brief Relaxation type for gmg
  */
  relax::which_relax m_gmg_relax_type;

  /*!
    @brief GMG multigrid type
  */
  amrmg::which_mg m_gmg_type;

  /*!
    @brief Needs setup
  */
  bool m_needs_setup;

  /*!
    @brief Use nwo or not
  */
  bool m_use_nwo;

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
    @brief Operator factory
  */
  RefCountedPtr<nwomfconductivityopfactory> m_nwo_opfact;
  
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
    @brief Geometric multigrid solver
  */
  AMRMultiGrid<LevelData<MFCellFAB> > m_gmg_solver;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<MFCellFAB> > m_bicgstab;

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
    @brief Set up the geometric multigrid solver 
  */
  virtual void setup_gmg();

  /*!
    @brief Set coefficients
  */
  virtual void set_coefficients();

  /*!
    @brief Set the permittivities
  */
  virtual void set_permittivities(const Vector<dielectric>& a_dielectrics);

  /*!
    @brief Set face-centered permittivities
  */
  virtual void set_face_perm(EBFluxFAB&                a_perm,
			     const RealVect&           a_origin,
			     const Real&               a_dx,
			     const Vector<dielectric>& a_dielectrics);

  /*!
    @brief Set EB-centered permittivities
  */
  virtual void set_eb_perm(BaseIVFAB<Real>&          a_perm,
			   const RealVect&           a_origin,
			   const Real&               a_dx,
			   const Vector<dielectric>& a_dielectrics);

  /*!
    @brief Set up the operator factory
  */
  virtual void setup_operator_factory();

  /*!
    @brief Set nwo operator factory
  */
  virtual void setup_nwo_operator_factory();

  /*!
    @brief Set up AMRMultiGrid
  */
  virtual void setup_solver();

private:

  /*!
    @brief Class for passing time and space dependent boundary conditions into the GMG boundary condition classes. 
  */
  class potential_func : public BaseBCFuncEval {
  public:
    
    potential_func(Real (*a_potential)(const Real a_time)){
      m_potential = a_potential;
    }

    virtual ~potential_func(){

    }

    virtual void set_time(const Real a_time){
      m_time = a_time;
    }

    virtual Real value(const RealVect& a_point, const int& a_comp) const {
      return m_potential(m_time);
    }

    virtual Real derivative(const RealVect& a_point, const int& a_comp, const int& a_dir) const {
      MayDay::Abort("poisson_multifluid_gmg::potential_func::derivative - this should not be called. How did you get here?");
    }
  protected:
    Real (*m_potential)(const Real a_time);

    Real m_time;
  };
};
#endif
