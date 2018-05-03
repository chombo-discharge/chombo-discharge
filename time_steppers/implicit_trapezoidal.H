/*!
  @file   implicit_trapezoidal.H
  @brief  Declaration of a split-step implicit method which uses the implicit trapezoidal method
  @author Robert Marskar
  @date   May. 2018
*/

#ifndef _IMPLICIT_TRAPEZOIDAL_
#define _IMPLICIT_TRAPEZOIDAL_


#include "time_stepper.H"

/*!
  @brief Implicit solver class for advancing the plasma equations with the implicit trapezoidal rule for the convection
  and source advancements, and the TGA scheme for the diffusion advance. 

  implicit_trapezoidal.max_pc_iter        = 10     # Maximum number of allowed predictor-corrector iterations for convection
  implicit_trapezoidal.max_newton_iter    = 10     # Maximum number of allowed Newton iterations for source term advance
  implicit_trapezoidal.solution_tolerance = 1E-6   # Solution increment tolerance. 
  implicit_trapezoidal.function_tolerance = 1E-6   # Equation tolerance. 
  implicit_trapezoidal.turn_off_advection = false  # Turn off advection, debugging only
  implicit_trapezoidal.turn_off_diffusion = false  # Turn off diffusion, debugging only
  implicit_trapezoidal.turn_off_source    = false  # Turn off source term, debugging only
  implicit_trapezoidal.turn_off_poisson   = false  # Turn off Poisson, debugging only
  implicit_trapezoidal.turn_off_rte       = false  # Turn off RTE, debugging only
*/
class implicit_trapezoidal : public time_stepper {
public:

  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;
  
  /*!
    @brief Constructor
  */
  implicit_trapezoidal();

  /*!
    @brief Destructor
  */
  virtual ~implicit_trapezoidal();

  /*!
    @brief Advance method
  */
  virtual Real advance(const Real a_dt);

  /*!
    @brief Restrict dt
  */
  virtual Real restrict_dt();

  /*!
    @brief Regrid internal storage
  */
  virtual void regrid_internals();

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocate_internals();

protected:

  /*!
    @brief Integration sequence
  */
  Vector<std::string> m_sequence;

  /*!
    @brief Scratch storage
  */
  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;

  /*!
    @brief Scratch storage
  */
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;

  /*!
    @brief Scratch storage
  */
  RefCountedPtr<poisson_storage> m_poisson_scratch;

  /*!
    @brief Scratch storage
  */
  RefCountedPtr<sigma_storage> m_sigma_scratch;

  /*!
    @brief Maximum number of allowed Newton iterations
  */
  int m_max_newton_iter;

  /*!
    @brief Maximum number of PC iterations
  */
  int m_max_pc_iter;

  /*!
    @brief Iteration solution tolerance
  */
  Real m_tol_x;

  /*!
    @brief Iteration function tolerance
  */
  Real m_tol_f;

  /*!
    @brief Finite difference step for computing the Jacobian
  */
  Real m_dnj;

  /*!
    @brief Maximum value. Used for error estimation
  */
  Real m_cdr_max;

  /*!
    @brief Fudge factor for computing Jacobian using finite differences
   */
  Real m_EPS;

  /*!
    @brief Reduction factor for CFL 
  */
  Real m_cfl_redu;

  /*!
    @brief Semi-implicit coupling or not
  */
  bool m_semi_implicit;

    /*!
    @brief Do advection advancement or not
  */
  bool m_do_convection;

  /*!
    @brief Do diffusion advancement or not
  */
  bool m_do_diffusion;

  /*!
    @brief Do reaction term advancement or not
  */
  bool m_do_reaction;

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
    @brief Cache states
  */
  virtual void store_states();

  /*!
    @brief Restore solver states
  */
  virtual void restore_states();

  /*!
    @brief Convection advance
  */
  virtual bool advance_convection(const Real a_dt);

  /*!
    @brief Diffusion advance
  */
  virtual bool advance_diffusion(const Real a_dt);

  /*!
    @brief Reaction advance
  */
  virtual bool advance_reaction(const Real a_dt);
};

#endif
