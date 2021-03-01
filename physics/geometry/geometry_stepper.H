/*!
  @file   geometry_stepper.H
  @brief  Declaration of a time_stepper class without any solver functionality, i.e. for geometries only. 
  @author Robert Marskar
  @date   March 2020
*/

#ifndef _GEOMETRY_STEPPER_
#define _GEOMETRY_STEPPER_

#include "time_stepper.H"

namespace physics {
  namespace geometry {
    
    /*!
      @brief Abstract class for performing time stepping
    */
    class geometry_stepper: public time_stepper {
    public:

      geometry_stepper();
      ~geometry_stepper();

      // Setup routines
      virtual void setup_solvers() override;
      virtual void allocate() override;
      virtual void initial_data() override;
      virtual void post_initialize() override;
      virtual void post_checkpoint_setup() override;

      // Registration routines
      virtual void register_realms() override;
      virtual void register_operators() override;

      // IO routines
      virtual void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const override;
      virtual void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) override;
      virtual void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const override;
      virtual int  get_num_plot_vars() const override;

      // Advance routines
      virtual void compute_dt(Real& a_dt, time_code& a_timecode) override;
      virtual Real advance(const Real a_dt) override;
      virtual void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) override;
      virtual void print_step_report() override;

      // Regrid routines
      virtual void deallocate() override;
      virtual void pre_regrid(const int a_lmin, const int a_old_finest_level) override;
      virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;
      virtual void post_regrid() override;

    };
  }
}

#endif
