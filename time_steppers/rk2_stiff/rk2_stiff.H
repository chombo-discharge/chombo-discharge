/*!
  @file   rk2_stiff.H
  @brief  Declaration of a RK2/TGA/Trapezoidal split-step method for advancing plasma equations. 
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _RK2_STIFF_
#define _RK2_STIFF_

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

  rk2_stiff.alpha              = 1.0    # Set alpha. 1.0 = Heuns method, 0.5 = Midpoint method
  rk2_stiff.max_newton_iter    = 10     # Maximum number of allowed Newton iterations
  rk2_stiff.solution_tolerance = 1E-6   # Solution tolerance for Newton iteration
  rk2_stiff.function_tolerance = 1E-6   # Function tolerance for Newton iteration
  rk2_stiff.turn_off_advection = false  # Turn off advection, debugging only
  rk2_stiff.turn_off_diffusion = false  # Turn off diffusion, debugging only
  rk2_stiff.turn_off_source    = false  # Turn off source term, debugging only
  rk2_stiff.turn_off_poisson   = false  # Turn off Poisson, debugging only
  rk2_stiff.turn_off_rte       = false  # Turn off RTE, debugging only
*/

class rk2_stiff : public time_stepper {
public:

  // Forward declare some classes

  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  /*!
    @brief Constructor
  */
  rk2_stiff();

  /*!
    @brief Destructor
  */
  virtual ~rk2_stiff();

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
  Real m_dnj;

  /*!
    @brief Maximum value. Used for error estimation
  */
  Real m_nmax;

  /*!
    @brief Fudge factor for computing Jacobian using finite differences
   */
  Real m_EPS;

  /*!
    @brief Semi-implicit or fully explicit coupling
  */
  bool m_simpi;

  /*!
    @brief Do advection advancement or not
  */
  bool m_do_adv_diff;

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
    @brief Advance source
  */
  virtual bool advance_sources(const Real a_dt);

  /*!
    @brief Advection advance
  */
  virtual void advance_advection_diffusion(const Real a_dt);

  /*!
    @brief Semi-implicit Newton advance
  */
  virtual void explicit_euler_predict_newton(const Real a_dt);

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
    @brief Backup solutions. Must do this in case the time step is rejected. 
  */
  virtual void cache_solutions();

  /*!
    @brief Compute step size for finite-difference evaluation of the Jacobian
  */
  virtual void compute_dnj();

  /*!
    @brief Compute the right-hand side for the trapezoidal discretization
  */
  virtual void compute_trapz_rhs(const Real a_dt); 

    /*!
    @brief Compute electric field at start of time step
  */
  virtual void compute_E_into_scratch();

  /*!
    @brief Compute sources for semi-implicit newton advance
  */
  virtual void compute_cdr_sources_for_newton_pred();

  /*!
    @brief Recompute the E-field by using the updated Newton values
  */
  virtual void recompute_newton_E();

  /*!
    @brief Recompute the RTE equations by using the updated Newton values
  */
  virtual void recompute_newton_rte();

  /*!
   @brief Setup the newton iterates
  */
  virtual void setup_newton_iterates(Vector<EBAMRCellData*>& a_iterates);

  /*!
    @brief Set cdr source terms to zero
  */
  virtual void set_cdr_sources_to_zero();              

  /*!
    @brief Compute the cdr velocities. This assumes that E is updated and resides in m_poisson_scratch
  */
  virtual void compute_cdr_velo();

  /*!
    @brief Compute the extrapolated things on the EB. This assumes E is updated, and that velocities have been updated
  */
  virtual void compute_cdr_eb_states();

  /*!
    @brief Compute diffusion coefficients. This assumes that E has been updated
  */
  virtual void compute_cdr_diffco();

  /*!
    @brief Compute the boundary conditions on the EB. This assumes that we've extrapolated velocities, densities, and fluxes
    to the EB
  */
  virtual void compute_cdr_fluxes();

  /*!
    @brief Compute the sigma flux as the sum of ion fluxes. 
  */
  virtual void compute_sigma_flux();

  /*!
    @brief Do the k1-advance for the cdr equations. This computes k1 and places the temporary solution right in the solver
  */
  virtual void advance_cdr_k1(const Real a_dt);

  /*!
    @brief Do the k1 advance for the sigma equation. This computes k1 and places the temporary solution right in the solver
  */
  virtual void advance_sigma_k1(const Real a_dt);

  /*!
    @brief Compute the k2 advance for cdr equations. This computes k2 and places the advanced solution into the solver
  */
  virtual void advance_cdr_k2(const Real a_dt);

  /*!
    @brrief Compute the k2 advance for sigma. This computse k2 and places the advanced solution into the solver. 
  */
  virtual void advance_sigma_k2(const Real a_dt);

  /*!
    @brief Perform point-wise Newton iteration by using the implicit trapezoidal rule
    @param[out] a_p      Newton solution
    @param[in]  a_rhs    Right-hand side of discretization (n - 0.5*dt*S)
    @param[in]  a_grad_E Gradient of electric field
    @param[in]  a_x      Previous iterate
    @param[in]  a_gradx  Gradient of the previous iterate
    @param[in]  a_E      Electric field
    @param[in]  a_rte    RTE densities
    @param[in]  a_pos    Position
    @param[in]  a_time   Time
    @param[in]  a_dt     Time increment
    @return Equation error \f$\sum_i |F_i|\f$
  */
  virtual Real newton_point_trapz(Vector<Real>&           a_p,
				  const Vector<Real>&     a_rhs,
				  const Vector<Real>&     a_x,
				  const Vector<RealVect>& a_gradx,
				  const RealVect&         a_E,
				  const RealVect&         a_grad_E,
				  const Vector<Real>&     a_rte,
				  const RealVect&         a_pos,
				  const Real&             a_time,
				  const Real&             a_dt);

    /*!
    @brief Restrict dt
  */
  virtual Real restrict_dt();
  
};
#endif
