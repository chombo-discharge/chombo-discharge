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

namespace physics {
  namespace poisson {

    /*!
      @brief Class that just solves the poisson equation
    */
    template <class T>
    class poisson_stepper : public time_stepper{
    public:

      static Real s_constant_one(const RealVect a_pos);
      static Real s_potential_one(const Real a_time);

      poisson_stepper(){m_verbosity = -1;}
      ~poisson_stepper(){}

      // Stuff below here is derived from time_stepper.
      // ----------------------------------------------
      void setup_solvers() override;
      void initial_data() override;

      // IO routines
      void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const override{}
      void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) override;
      void post_checkpoint_setup() override {}
      int get_num_plot_vars() const override;
      void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const override;

      // Advance routines
      void compute_dt(Real& a_dt, time_code::which_code& a_timecode) override {a_dt = 0.0;}
      Real advance(const Real a_dt) override;
      void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) override;
      void print_step_report() override {};

      // Regrid routines
      bool need_to_regrid() override {return false;}
      void register_operators() override;
      void cache() override;
      void deallocate() override;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;

    protected:

      RefCountedPtr<poisson_solver> m_poisson;
      RefCountedPtr<sigma_solver>   m_sigma;

      Real m_init_rho;
      Real m_init_sigma;
    };
  }
}

#include "poisson_stepperI.H"

#endif
