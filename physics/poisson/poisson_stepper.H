/*!
  @file   poisson_stepper.H
  @brief  time_stepper class for only solving the Poisson equation (with surface charge)
  @author Robert Marskar
  @data   March 2020
*/

#ifndef _POISSON_STEPPER_
#define _POISSON_STEPPER_

#include "time_stepper.H"
#include "poisson_solver.H"
#include "sigma_solver.H"

/*!
  @brief Class that just solves the poisson equation
*/
template <class T>
class poisson_stepper : public time_stepper{
public:

  static Real s_constant_one(const RealVect a_pos);
  static Real s_potential_one(const Real a_time);

  poisson_stepper(){}
  ~poisson_stepper(){}

  // Stuff below here is derived from time_stepper.
  // ----------------------------------------------
  void setup_solvers();
  void initial_data();

  // IO routines
  void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const{}
  void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl);
  void post_checkpoint_setup(){}
  int get_num_plot_vars() const;
  void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const;

  // Advance routines
  void compute_dt(Real& a_dt, time_code::which_code& a_timecode){a_dt = 0.0;}
  Real advance(const Real a_dt);
  void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt);
  void print_step_report(){}

  // Regrid routines
  bool need_to_regrid(){return false;}
  void cache_internals(){}
  void deallocate_internals(){}
  void cache_states();
  void deallocate_solver_internals(){}
  void regrid_solvers(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);
  void allocate_internals(){}
  void regrid_internals(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){}
  void post_regrid(){}

protected:

  RefCountedPtr<poisson_solver> m_poisson;
  RefCountedPtr<sigma_solver>   m_sigma;

  Real m_init_rho;
  Real m_init_sigma;
};

#include "poisson_stepperI.H"

#endif
