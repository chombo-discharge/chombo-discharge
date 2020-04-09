/*!
  @file   advection_diffusion_stepper.H
  @brief  time_stepper class for only solving the Advection_Diffusion equation (with surface charge)
  @author Robert Marskar
  @data   March 2020
*/

#ifndef _ADVECTION_DIFFUSION_STEPPER_
#define _ADVECTION_DIFFUSION_STEPPER_

#include "time_stepper.H"
#include "cdr_solver.H"

namespace physics {
  namespace advection_diffusion {

    /*!
      @brief Class that just solves the advection_diffusion equation
    */
    class advection_diffusion_stepper : public time_stepper {
    public:

      advection_diffusion_stepper();
      advection_diffusion_stepper(RefCountedPtr<cdr_solver>& a_solver);
      ~advection_diffusion_stepper();

      // Setup routines
      void initial_data() override;

      // IO routines
      void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const override;
      void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) override;
      void post_checkpoint_setup() override;
      int get_num_plot_vars() const override;
      void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const override;

      // Advance routines
      void compute_dt(Real& a_dt, time_code::which_code& a_timecode) override;

      void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) override;
      void print_step_report() override;

      // Regrid routines
      bool need_to_regrid() override;
      void cache() override;
      void deallocate() override;

      // Pure routines for derived classes
      void setup_solvers() override;
      Real advance(const Real a_dt) override;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;
      

    protected:

      // CDR solver & species
      RefCountedPtr<cdr_solver>  m_solver;
      RefCountedPtr<cdr_species> m_species;

      // Diffusion coefficient and velocity
      Real m_cfl;
      Real m_diffco;
      Real m_omega;

      // Use stochastic diffusion
      bool m_fhd;

      // Integrator (yes, I'm too lazy to write an enum for holding this)
      int m_integrator;

      void set_velocity();
      void set_velocity(const int a_level);

    private:

      EBAMRCellData m_tmp;
      EBAMRCellData m_k1;
      EBAMRCellData m_k2;
    };
  }
}
#endif
