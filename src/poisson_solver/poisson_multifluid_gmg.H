/*!
  @file poisson_multifluid_gmg.H
  @brief Geometric multigrid multifluid poisson solver
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _POISSON_MULTIFLUID_GMG_
#define _POISSON_MULTIFLUID_GMG_

#include "poisson_solver.H"
#include "mf_helmholtz_opfactory.H"

#include <AMRMultiGrid.H>
#include <BiCGStabSolver.H>
#include <EBConductivityOpFactory.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>

/*!
  @brief Abstract Poisson solver class. Just an interface, so pretty lightweight stuff.
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
  virtual void solve();

  /*!
    @brief Solve Poisson onto state 
  */
  virtual void solve(MFAMRCellData& a_state, const MFAMRCellData& a_source);

    /*!
    @brief Set solver parameters. These are normally good for this class. 
  */
  virtual void set_gmg_solver_parameters(relax::which_relax a_relax_type = relax::gsrb_fast,
					 amrmg::which_mg a_gmg_type      = amrmg::vcycle,
					 const int a_verbosity           = 10,
					 const int a_pre_smooth          = 16,
					 const int a_post_smooth         = 16,
					 const int a_bot_smooth          = 16,
					 const int a_max_iter            = 32,
					 const Real a_eps                = 1.E-30,
					 const Real a_hang               = 0.2,
					 const Real a_norm_thresh        = 1.E-30);

  /*!
    @brief 2 ghost cells. That's what we need. 
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
  RefCountedPtr<mf_helmholtz_opfactory> m_opfact;

  /*!
    @brief ConductivityOp operator factory
  */
  RefCountedPtr<EBConductivityOpFactory> m_cond_op_fact;

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
    @brief Set up AMRMultiGrid
  */
  virtual void setup_solver();

  // Putting this in a function of its own
  virtual void base_tests();

  // Test function, will be removed from final interface
  virtual void do_ebcond_test();
};
#endif
