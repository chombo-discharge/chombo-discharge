/*!
  @file   splitstep_rk2_tga_trapz.H
  @brief  Declaration of a RK2/TGA/Trapezoidal split-step method for advancing plasma equations. 
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _SPLITSTEP_RK2_TGA_TRAPZ_
#define _SPLITSTEP_RK2_TGA_TRAPZ_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

/*!
  @brief Class for evolvinv plasma equations with a split-step RK2/TGA/Trapezoidal method
  @details This class uses RK2 advancements for the advection term, the implicit TGA method
  for the diffusive part, and the trapezoidal method for the source term part. 

  Class options
  -------------

  splitstep_rk2_tga_trapz.alpha              = 1.0 # Set alpha. 1.0 = Heuns method, 0.5 = Midpoint method
  splitstep_rk2.tga_trapz.max_newton_iter    = 10     # Maximum number of allowed Newton iterations
  splitstep_rk2.tga_trapz.solution_tolerance = 1E-6   # Solution tolerance for Newton iteration
  splitstep_rk2.tga_trapz.function_tolerance = 1E-6   # Function tolerance for Newton iteration
  splitstep_rk2_tga_trapz.turn_off_advection = false  # Turn off advection, debugging only
  splitstep_rk2_tga_trapz.turn_off_diffusion = false  # Turn off diffusion, debugging only
  splitstep_rk2_tga_trapz.turn_off_source    = false  # Turn off source term, debugging only
  splitstep_rk2_tga_trapz.turn_off_poisson   = false  # Turn off Poisson, debugging only
  splitstep_rk2_tga_trapz.turn_off_rte       = false  # Turn off RTE, debugging only
*/

class splitstep_rk2_tga_trapz : public time_stepper {
public:

  // Forward declare some classes

  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  /*!
    @brief Constructor
  */
  splitstep_rk2_tga_trapz();

  /*!
    @brief Destructor
  */
  virtual ~splitstep_rk2_tga_trapz();

  /*!
    @brief Advance method
  */
  virtual Real advance(const Real a_dt);

  /*!
    @brief Deallocate storage
  */
  virtual void deallocate_internals();

  /*!
    @brief Regrid internals
  */
  virtual void regrid_internals();

  /*!
    @brief Compute dt
  */
  virtual void compute_dt(Real& a_dt, time_code::which_code& a_timecode);

protected:

  /*!
    @brief Scratch storage
  */
  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;

  /*!
    @brief Scratch storage
  */
  RefCountedPtr<poisson_storage> m_poisson_scratch;

  /*!
    @brief Scratch storage
  */
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;

  /*!
    @brief Scratch storage
  */
  RefCountedPtr<sigma_storage> m_sigma_scratch;

  /*!
    @brief RK2 method. 
  */
  Real m_alpha;

  /*!
    @brief Newton iteration solution tolerance
  */
  Real m_tol_x;

  /*!
    @brief Newton iteration function tolerance
  */
  Real m_tol_f;

  /*!
    @brief Maximum number of allowed Newton iterations
  */
  Real m_max_iter;

  /*!
    @brief Finite difference step for computing the Jacobian
  */
  Real m_epsj;

  /*!
    @brief Semi-implicit or fully explicit coupling
  */
  bool m_simpi;

  /*!
    @brief Do advection advancement or not
  */
  bool m_do_advection;

  /*!
    @brief Do diffusion advancement or not
  */
  bool m_do_diffusion;

  /*!
    @brief Do source term advancement or not
  */
  bool m_do_source;

  /*!
    @brief Solve RTE equations
  */
  bool m_do_rte;

  /*!
    @brief Solve RTE equations
  */
  bool m_do_poisson;

  /*!
    @brief Get storage
  */
  virtual RefCountedPtr<cdr_storage>& get_cdr_storage(const cdr_iterator& a_solverit);

  /*!
    @brief Get storage
  */
  virtual RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

  /*!
    @brief Advection advance
  */
  virtual void advance_advection(const Real a_dt);

  /*!
    @brief Advect diffusion
  */
  virtual void advance_diffusion(const Real a_dt);

  /*!
    @brief Advance source
  */
  virtual void advance_sources(const Real a_dt);

  /*!
    @brief Semi-implicit Newton advance
  */
  virtual void advance_semi_implicit_newton(const Real a_dt);

  /*!
    @brief Allocate storage
  */
  virtual void allocate_cdr_storage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_poisson_storage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_rte_storage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_sigma_storage();

  /*!
    @brief Restrict dt
  */
  virtual Real restrict_dt();

  /*!
    @brief Backup solutions. Must do this in case the time step is rejected. 
  */
  virtual void cache_solutions();

  /*!
    @brief Compute step size for finite-difference evaluation of the Jacobian
  */
  virtual void compute_epsj();

  /*!
    @briefCompute the right-hand side for the trapezoidal discretization
  */
  virtual void compute_trapz_rhs(const Real a_dt); 

    /*!
    @brief Compute electric field at start of time step
  */
  virtual void compute_E_at_start_of_time_step();

  /*!
    @brief Compute sources for semi-implicit newton advance
  */
  virtual void compute_cdr_sources_for_simp_newt();

  /*!
    @brief Compute the semi-implicitly advanced electric field
  */
  virtual void compute_semi_implicit_potential();

  /*!
    @brief Solve the stationary RTE equations by using the semi-implicitly advanced cdr states
  */
  virtual void compute_semi_implicit_rte(const Real a_dt);

  /*!
    @brief Recompute the E-field by using the updated Newton values
  */
  virtual void recompute_newton_E();

  /*!
    @brief Recompute the RTE equations by using the updated Newton values
  */
  virtual void recompute_newton_rte();

  /*!
    @brief Perform point-wise Newton iteration by using the implicit trapezoidal rule
    @param[out] a_p      Newton solution
    @param[in]  a_rhs    Right-hand side of discretization (n - 0.5*dt*S)
    @param[in]  a_grad_E Gradient of electric field
    @param[in]  a_x      Previous iterate
    @param[in]  a_gradx  Gradient of the previous iterate
    @param[in]  a_E      Electric field
    @param[in]  a_gradE  Gradient of electric field
    @param[in]  a_rte    RTE densities
    @param[in]  a_pos    Position
    @param[in]  a_time   Time

  */
  virtual void newton_point_trapz(Vector<Real>&           a_p,
				  const Vector<Real>&     a_rhs,
				  const Vector<Real>&     a_x,
				  const Vector<RealVect>& a_gradx,
				  const RealVect&         a_E,
				  const RealVect&         a_grad_E,
				  const Vector<Real>&     a_rte,
				  const RealVect&         a_pos,
				  const Real&             a_time,
				  const Real&             a_dt);
  
};
#endif
