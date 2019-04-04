/*!
  @file   cdr_tga.H
  @brief  Extensions of cdr_solver which use ebconductivityop and TGA for diffusion parts
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _CDR_TGA_
#define _CDR_TGA_

#include "cdr_solver.H"

#include "ebconductivityopfactory.H"

#include <AMRMultiGrid.H>
#include <AMRTGA.H>
#include <BiCGStabSolver.H>
#include <EBBackwardEuler.H>
#include <EBSimpleSolver.H>

/*!
  @brief Extension class of cdr_solver that uses ebconductivityop for diffusion part. This class also has a multigrid
  solver for implicit diffusion

  @details 

  Class options
  -------------

      cdr_tga.use_tga           = true      # Use TGA or not. Valid options are "true" (use TGA) or "false" (use backward Euler)
      cdr_tga.gmg_verbosity     = -1        # GMG verbosity
      cdr_tga.gmg_pre_smooth    = 16        # Number of relaxations in downsweep
      cdr_tga.gmg_post_smooth   = 16        # Number of relaxations in upsweep
      cdr_tga.gmg_bott_smooth   = 16        # NUmber of relaxations before dropping to bottom solver
      cdr_tga.gmg_min_iter      = 5         # Minimum number of iterations
      cdr_tga.gmg_max_iter      = 32        # Maximum number of iterations
      cdr_tga.gmg_tolerance     = 1.E-6     # Residue tolerance
      cdr_tga.gmg_hang          = 0.2       # Solver hang
      cdr_tga.gmg_bottom_drop   = 8         # Bottom drop
      cdr_tga.gmg_bottom_solver = simple    # Bottom solver type. Valid options are 'simple' and 'bicgstab'
      cdr_tga.gmg_bottom_relax  = 32        # Number of relaxations in bottom solve ('simple' solver only)

*/
class cdr_tga : public cdr_solver {
public:

  /*!
    @brief Constructor
  */
  cdr_tga();

  /*!
    @brief Destructor
  */
  virtual ~cdr_tga();

  /*!
    @brief Advance state. Current is explicit Euler.
    @details Make the default RK2 for non-diffusive media and semi-implicit TGA for diffusive media. 
  */
  virtual void advance(EBAMRCellData& a_state, const Real& a_dt);

  /*!
    @brief Advective advance. This must be implemented by the solver. 
  */
  virtual void advance_advect(EBAMRCellData& a_state, const Real a_dt) = 0;

  /*!
    @brief Diffusion advance
  */
  virtual void advance_diffusion(EBAMRCellData& a_state, EBAMRCellData& a_error, const Real a_dt);

  /*!
    @brief Compute div(nv - D*grad(n))
  */
  virtual void compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt);

  /*!
    @brief Compute advective derivative
  */
  virtual void compute_divF(EBAMRCellData& a_divF, const EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_redist);

  /*!
    @brief Compute pure hyperbolic Euler advance using subcycling in time
  */
  virtual void eulerF_subcycle(EBAMRCellData& a_state, const Real a_dt, const bool a_redist) = 0;

  /*!
    @brief Compute diffusion term
  */
  virtual void compute_divD(EBAMRCellData& a_diffusive_term, const EBAMRCellData& a_state);

  /*!
    @brief Set solver parameters. These are normally good for this class. 
  */
  virtual void set_gmg_solver_parameters(relax::which_relax a_relax_type = relax::gsrb_fast,
					 amrmg::which_mg a_gmg_type      = amrmg::vcycle,
					 const int a_verbosity           = -1,
					 const int a_pre_smooth          = 8,
					 const int a_post_smooth         = 8,
					 const int a_bot_smooth          = 8,
					 const int a_max_iter            = 16,
					 const int a_min_iter            = 5,
					 const Real a_eps                = 1.E-8,
					 const Real a_hang               = 0.2);

  /*!
    @brief Set the bottom solver.  
    @param[in] a_whichsolver Solver type. 0 = smoother, 1 = BiCGStab
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
    @brief Use tga for temporal integration
  */
  virtual void set_tga(const bool a_use_tga);

  /*!
    @brief Query number of ghost cells
  */
  virtual int query_ghost() const = 0;

  /*!
    @brief Pure TGA advance without source terms
  */
  virtual void advance_tga(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt);

  /*!
    @brief Pure Euler advance with source terms
  */
  virtual void advance_tga(EBAMRCellData&       a_new_state,
			     const EBAMRCellData& a_old_state,
			     const EBAMRCellData& a_source,
			     const Real           a_dt);

  /*!
    @brief Pure Euler advance without source terms
  */
  virtual void advance_euler(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt);

  /*!
    @brief Pure Euler advance with source terms
  */
  virtual void advance_euler(EBAMRCellData&       a_new_state,
			     const EBAMRCellData& a_old_state,
			     const EBAMRCellData& a_source,
			     const Real           a_dt);

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
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > m_gmg_solver;

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
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB> > m_bicgstab;

  /*!
    @brief Simple solver
  */
  EBSimpleSolver m_simple_solver;

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
    @brief Use tga
  */
  bool m_use_tga;

  /*!
    @brief 
  */
  Real m_gmg_eps;
  
  /*!
    @brief 
  */
  Real m_gmg_hang;

  /*!
    @brief TGA advancement algorithm
  */
  virtual void advance_tga(EBAMRCellData& a_state, const Real a_dt);

  /*!
    @brief Advection-only extrapolation to faces
  */
  virtual void advect_to_faces(EBAMRFluxData& a_face_state, const EBAMRCellData& a_state, const Real a_extrap_dt) = 0;

  /*!
    @brief Set up geometric multigrid
  */
  virtual void setup_gmg();

  /*!
    @brief Setup the operator factory
  */
  virtual void setup_operator_factory();

  /*!
    @brief Setup multigrid
  */
  virtual void setup_multigrid();

  /*!
    @brief Setup the TGA solver
  */
  virtual void setup_tga();

  /*!
    @brief Setup the Euler solver
  */
  virtual void setup_euler();
};
#endif
