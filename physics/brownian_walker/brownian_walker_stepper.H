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
      brownian_walker_stepper(RefCountedPtr<ito_solver>& a_solver);
      ~brownian_walker_stepper();

      // Setup routines
      void initial_data() override;
      void allocate() override;

      // Load balancing
      void load_balance(Vector<Vector<int> >&            a_procs,
			Vector<Vector<Box> >&            a_boxes,
			const Vector<DisjointBoxLayout>& a_grids,
			const int                        a_lmin,
			const int                        a_finest_level);

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
      void register_operators() override;
      void pre_regrid(const int a_lbase, const int a_old_finest_level) override;
      void deallocate() override;

      // Pure routines for derived classes
      void setup_solvers() override;
      Real advance(const Real a_dt) override;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;

    protected:

      // Phase
      phase::which_phase m_phase;

      // Ito solver & species
      RefCountedPtr<ito_solver>  m_solver;
      RefCountedPtr<ito_species> m_species;

      // Diffusion coefficient and velocity
      Real m_diffco;
      Real m_omega;

      // Particles per cell
      int m_ppc;

      // Maximum hops and stuff like that
      Real m_max_cells_hop;

      void set_velocity();
      void set_velocity(const int a_level);
    };
  }
}
#endif
