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
#include "larsen_coefs.H"
#include "robinconductivityebbcfactory.H"
#include "robinconductivitydomainbcfactory.H"

#include <AMRMultiGrid.H>
#include <AMRTGA.H>
#include <BiCGStabSolver.H>
#include <EBBackwardEuler.H>
#include <EBSimpleSolver.H>
#include <EBConductivityOpFactory.H>

/*!
  @brief Radiative tranfer equation solver in the SP1 (diffusion) approximation
  @details This class solves the RTE in the form 
  \f[
  \frac{\partial E}{\partial t} + c\kappa E - \nabla\cdot\left[\frac{c}{\kappa}\nabla E\right] = \eta
  \f]

  Class options
  -------------

      eddington_sp1.reflectivity      = 1.E-4     # Reflectivity
      eddington_sp1.gmg_verbosity     = -1        # GMG verbosity
      eddington_sp1.gmg_pre_smooth    = 16        # Number of relaxations in downsweep
      eddington_sp1.gmg_post_smooth   = 16        # Number of relaxations in upsweep
      eddington_sp1.gmg_bott_smooth   = 16        # NUmber of relaxations before dropping to bottom solver
      eddington_sp1.gmg_min_iter      = 5         # Minimum number of iterations
      eddington_sp1.gmg_max_iter      = 32        # Maximum number of iterations
      eddington_sp1.gmg_tolerance     = 1.E-6     # Residue tolerance
      eddington_sp1.gmg_hang          = 0.2       # Solver hang
      eddington_sp1.gmg_bottom_drop   = 8         # Bottom drop
      eddington_sp1.gmg_bottom_solver = simple    # Bottom solver type. Valid options are 'simple' and 'bicgstab'
      eddington_sp1.gmg_bottom_relax  = 32        # Number of relaxations in bottom solve ('simple' solver only)

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
    @brief Advance RTE onto state a_state
  */
  virtual bool advance(const Real a_dt, EBAMRCellData& a_state, const EBAMRCellData& a_source, const bool a_zerophi = false);
  
  /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals();
  
  /*!
    @brief Cache state
  */
  virtual void cache_state();

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocate_internals();

  /*!
    @brief Regrid function for this class
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
					 const int a_max_iter            = 16,
					 const int a_min_iter            = 5,
					 const Real a_eps                = 1.E-6,
					 const Real a_hang               = 0.2);

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
    @brief Use tga for temporal integration
  */
  virtual void set_tga(const bool a_use_tga);

  /*!
    @brief Set reflectivity
  */
  virtual void set_reflectivity(const Real a_reflectivity);

  /*!
    @brief Set reflection coefficients. These must be obtained by integrals over the reflectivity.
  */
  virtual void set_reflection_coefficients(const Real r1, const Real r2);

  /*!
    @brief Compute the boundary flux. For Eddington, the boundary flux is = c*phi/2
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
  virtual void compute_density(EBAMRCellData& a_isotropic, const EBAMRCellData& a_state);


#ifdef CH_USE_HDF5
  /*!
    @brief Write plot file
  */
  virtual void write_plot_file();
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
    @brief Use tga
  */
  bool m_use_tga;

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
    @brief Reflection coefficient
  */
  Real m_r1;

  /*!
    @brief Reflection coefficient
  */
  Real m_r2;

  /*!
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > m_gmg_solver;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB> > m_bicgstab;

  /*!
    @brief TGA solver
  */
  RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > > m_tgasolver;

  /*!
    @brief Euler solver
  */
  RefCountedPtr<EBBackwardEuler> m_eulersolver;

  /*!
    @brief Operator factory
  */
  RefCountedPtr<ebconductivityopfactory> m_opfact;

  /*!
    @brief Larsen coefficients for photon outflow bc
  */
  RefCountedPtr<larsen_coefs> m_robinco;

  /*!
    @brief Photon outflow boundary conditions
  */
  RefCountedPtr<robinconductivityebbcfactory> m_ebfact;

  /*!
    @brief Photon outflow boundary conditions
  */
  RefCountedPtr<robinconductivitydomainbcfactory> m_domfact;

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

  /*!
    @brief Set up the TGA solver
  */
  virtual void setup_tga();
  
  /*!
    @brief Set up the Euler solver
  */
  virtual void setup_euler();

};
#endif
