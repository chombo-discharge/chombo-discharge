/*!
  @file   rte_stepper.H
  @brief  time_stepper class for only solving the Rte equation (with surface charge)
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _RTE_STEPPER_
#define _RTE_STEPPER_

#include "time_stepper.H"
#include "rte_solver.H"

namespace physics {
  namespace rte {

    /*!
      @brief Class that just solves the rte equation
    */
    template <class T>
    class rte_stepper : public time_stepper {
    public:
      rte_stepper();
      ~rte_stepper();

      // Stuff below here is derived from time_stepper.
      // ----------------------------------------------
      void setup_solvers() override;
      void initial_data() override;

      // IO routines
      void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const override;
      void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) override;
      void post_checkpoint_setup() override;
      int get_num_plot_vars() const override;
      void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const override;

      // Advance routines
      void compute_dt(Real& a_dt, time_code::which_code& a_timecode) override;
      Real advance(const Real a_dt) override;
      void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) override;
      void print_step_report() override;

      // Regrid routines
      bool need_to_regrid() override;
      void register_operators() override;
      void cache() override;
      void deallocate() override;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;

    protected:

      // Solver and species
      RefCountedPtr<rte_solver>  m_solver;
      RefCountedPtr<rte_species> m_species;

      // Stuff for source term
      Real     m_blob_radius;    // Gaussian radius
      Real     m_blob_amplitude; // Gaussian amplitude
      RealVect m_blob_center;    // Center

      void set_gaussian_source();
    };
  }
}

#include "rte_stepperI.H"

#endif
