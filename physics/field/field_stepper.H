/*!
  @file   field_stepper.H
  @brief  time_stepper class for only solving the Poisson equation (with surface charge)
  @author Robert Marskar
  @data   March 2020
*/

#ifndef _FIELD_STEPPER_
#define _FIELD_STEPPER_

#include "time_stepper.H"
#include "field_solver.H"
#include "sigma_solver.H"

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace poisson {

    /*!
      @brief Class that just solves the poisson equation
    */
    template <class T>
    class field_stepper : public time_stepper{
    public:

      static Real s_constant_one(const RealVect a_pos);
      static Real s_potential_one(const Real a_time);

      field_stepper();
      ~field_stepper();

      // Stuff below here is derived from time_stepper.
      // ----------------------------------------------
      void setup_solvers() override;
      void initial_data() override;
      void post_initialize() override;
      void allocate() override;

      // IO routines
      void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const override{}
      void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) override;
      void post_checkpoint_setup() override {}
      int get_num_plot_vars() const override;
      void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const override;

      // Advance routines
      void compute_dt(Real& a_dt, time_code& a_timecode) override {a_dt = 0.0;}
      Real advance(const Real a_dt) override;
      void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) override;
      void print_step_report() override {};

      // Registration routines
      void register_realms() override;
      void register_operators() override;

      // Regrid routines
      bool need_to_regrid() override {return false;}
      void pre_regrid(const int a_lbase, const int a_finest_level) override;
      void deallocate() override;
      bool load_balance_realm(const std::string a_realm) const override;
      void load_balance_boxes(Vector<Vector<int> >&            a_procs,
			      Vector<Vector<Box> >&            a_boxes,
			      const std::string                a_realm,
			      const Vector<DisjointBoxLayout>& a_grids,
			      const int                        a_lmin,
			      const int                        a_finest_level) override;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;
      void post_regrid() override;

      // Solve routine
      void solve_poisson();
      

    protected:

      bool m_load_balance;
      box_sorting m_boxsort;

      RefCountedPtr<field_solver> m_poisson;
      RefCountedPtr<sigma_solver>   m_sigma;

      std::string m_realm;

      Real m_init_rho;
      Real m_init_sigma;
    };
  }
}
#include "CD_NamespaceFooter.H"
#include "field_stepperI.H"

#endif
