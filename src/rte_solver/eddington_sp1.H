/*!
  @file eddington_sp1.H
  @brief Declaration of a first order Eddington solver
  @author Robert Marskar
  @date Jan. 2018
*/

#ifndef _EDDINGTON_SP1_
#define _EDDINGTON_SP1_

#include "rte_solver.H"
#include "ebconductivityopfactory.H"

#include <AMRMultiGrid.H>
#include <BiCGStabSolver.H>
#include <EBSimpleSolver.H>
#include <EBConductivityOpFactory.H>

/*!
  @brief Radiative tranfer equation solver in the SP1 (diffusion) approximation
  @details This class solves the RTE in the form 
  \f[
  \frac{\partial E}{\partial t} + c\kappa E - \nabla\cdot\left[\frac{c}{\kappa}\nabla E\right] = \eta
  \f]
*/

class eddington_sp1 : public rte_solver {
public:

  /*!
    @brief Constructor
  */
  eddington_sp1();

  /*!
    @brief Destructor
  */
  virtual ~eddington_sp1();

  /*!
    @brief Set solver parameters. These are normally good for this class. 
  */
  virtual void set_gmg_solver_parameters(relax::which_relax a_relax_type = relax::gsrb_fast,
					 amrmg::which_mg a_gmg_type      = amrmg::vcycle,
					 const int a_verbosity           = 100,
					 const int a_pre_smooth          = 16,
					 const int a_post_smooth         = 16,
					 const int a_bot_smooth          = 16,
					 const int a_max_iter            = 16,
					 const int a_min_iter            = 5,
					 const Real a_eps                = 1.E-8,
					 const Real a_hang               = 0.2,
					 const Real a_norm_thresh        = 1.E-99);

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
    @brief Advance RTE onto state a_state
  */
  virtual void advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source);

  /*!
    @brief Compute the boundary flux given a state (a_state will be different for different RTE approximations)
  */
  virtual void compute_boundary_flux(EBAMRIVData& a_ebflux, const EBAMRCellData& a_state);

  /*!
    @brief Compute the flux
    @details For diffusive models, the flux will be something like grad(a_state). For higher-order models, the flux 
    will be contained in a_state (somehow);
  */
  virtual void compute_flux(EBAMRCellData& a_flux, const EBAMRCellData& a_state);

  /*!
    @brief Get isotropic part. 
  */
  virtual void get_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_state);

  /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals();

#ifdef CH_USE_HDF5
  /*!
    @brief Write plot file
  */
  virtual void write_plot_file(const int a_step = 0);
#endif

  /*!
    @brief Return required number of ghost cells
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
    @brief Geometric multigrid solver
  */
  AMRMultiGrid<LevelData<EBCellFAB> > m_gmg_solver;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB> > m_bicgstab;

  /*!
    @brief Operator factory
  */
  RefCountedPtr<ebconductivityopfactory> m_opfact;

  /*!
    @brief multi-fluid simple solver
  */
  EBSimpleSolver m_simple_solver;
  
  /*!
    @brief Residue for stationary solver
  */
  EBAMRCellData m_resid;

  /*!
    @brief a-coefficient
  */
  EBAMRCellData m_aco;

  /*!
    @brief b-coefficient
  */
  EBAMRFluxData m_bco;

  /*!
    @brief b-coefficient
  */
  EBAMRIVData m_bco_irreg;

  /*!
    @brief Set up geometric multigrid
  */
  virtual void setup_gmg();

  /*!
    @brief Set kappa
  */
  virtual void set_coefficients();

  /*!
    @brief Set the a-coefficient and b-coefficient
  */
  virtual void set_aco_and_bco();

  /*!
    @brief Set a-coefficient
  */
  virtual void set_aco(EBCellFAB& a_aco, const RealVect a_origin, const Real a_dx);

  /*!
    @brief Set a-coefficient
  */
  virtual void set_bco_face(EBFluxFAB& a_bco, const RealVect a_origin, const Real a_dx);

    /*!
    @brief Set a-coefficient
  */
  virtual void set_bco_eb(BaseIVFAB<Real>& a_bco, const RealVect a_origin, const Real a_dx);

  /*!
    @brief Set the operator factory
  */
  virtual void setup_operator_factory();

  /*!
    @brief Set the GMG solver
  */
  virtual void setup_multigrid();
};

#endif
