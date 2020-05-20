/*!
  @file   ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
  @brief  Abstract class for integrating the Ito plasma equations
*/

#ifndef _ITO_PLASMA_
#define _ITO_PLASMA_

#include "time_stepper.H"
#include "ito_layout.H"
#include "rte_layout.H"
#include "mc_photo.H"
#include "poisson_solver.H"
#include "sigma_solver.H"

namespace physics {
  namespace ito_plasma {

    class ito_plasma_stepper : public time_stepper {
    public:

      // Ctor and dtor
      ito_plasma_stepper();
      ~ito_plasma_stepper();

      // Setup routines
      void setup_solvers();
      void initial_data();
      void post_checkpoint_setup();
  
      // IO routines
      void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) override;
      void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) override;
      void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const override;
      int  get_num_plot_vars() const override;
  
      // Advance routines
      void compute_dt(Real& a_dt, time_code::which_code& a_timecode) override;
      Real advance(const Real a_dt) override;
      void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) override;
      void print_step_report() override;

      // New regrid routines
      bool need_to_regrid();
      void register_operators() override; 
      void pre_regrid(const int a_lmin, const int a_old_finest_level) override;
      void deallocate() override;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;

    protected:

      /*!
	@brief Class name
      */
      std::string m_name;

      /*!
	@brief Ito solvers
      */
      RefCountedPtr<ito_layout> m_ito;

      /*!
	@brief Radiative transfer solvers
      */
      RefCountedPtr<rte_layout> m_rte;

      /*!
	@brief Poisson solver
      */
      RefCountedPtr<poisson_solver> m_poisson;

      /*!
	@brief Surface charge solver
      */
      RefCountedPtr<sigma_solver> m_sigma;
  
    }
  }
};
#endif
