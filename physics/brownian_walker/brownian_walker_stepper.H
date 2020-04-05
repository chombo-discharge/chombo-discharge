/*!
  @file   brownian_walker_stepper.H
  @brief  time_stepper class for only solving the Brownian_Walker equation (with surface charge)
  @author Robert Marskar
  @data   March 2020
*/

#ifndef _BROWNIAN_WALKER_STEPPER_
#define _BROWNIAN_WALKER_STEPPER_

#include "time_stepper.H"
#include "ito_solver.H"

namespace physics {
  namespace brownian_walker {

    /*!
      @brief Class that just solves the brownian_walker equation
    */
    class brownian_walker_stepper : public time_stepper {
    public:

      brownian_walker_stepper();
      brownian_walker_stepper(RefCountedPtr<cdr_solver>& a_solver);
      ~brownian_walker_stepper();

      // Setup routines
      void initial_data();

      // IO routines
      void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const;
      void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl);
      void post_checkpoint_setup();
      int get_num_plot_vars() const;
      void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const;

      // Advance routines
      void compute_dt(Real& a_dt, time_code::which_code& a_timecode);

      void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt);
      void print_step_report();

      // Regrid routines
      bool need_to_regrid();
      void cache();
      void deallocate();

      // Pure routines for derived classes
      void setup_solvers() = 0;
      Real advance(const Real a_dt) = 0;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) = 0;

    protected:

      // CDR solver & species
      RefCountedPtr<ito_solver>  m_solver;
      RefCountedPtr<ito_species> m_species;
    };
  }
}
#endif
