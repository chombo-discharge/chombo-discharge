/*!
  @file   ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
  @brief  Abstract class for integrating the Ito plasma equations
*/

#ifndef _ITO_PLASMA_STEPPER_
#define _ITO_PLASMA_STEPPER_

#include "time_stepper.H"
#include "ito_plasma_physics.H"
#include "ito_layout.H"
#include "rte_layout.H"
#include "mc_photo.H"
#include "poisson_solver.H"
#include "sigma_solver.H"

namespace physics {
  namespace ito_plasma {

    class ito_plasma_stepper : public time_stepper {
    public:

      static Real s_constant_one(const RealVect a_pos);

      // Ctor and dtor
      ito_plasma_stepper();
      ito_plasma_stepper(RefCountedPtr<ito_plasma_physics>& a_physics);
      ~ito_plasma_stepper();

      // Setup routines
      virtual void parse_options() = 0;
      virtual void allocate_internals() = 0;
      void setup_solvers() override;
      void initial_data() override;
      void initial_sigma();
      void post_checkpoint_setup() override;
      void setup_ito();
      void setup_poisson();
      void setup_rte();
      void setup_sigma();
      void set_poisson_wall_func(const int a_dir, const Side::LoHiSide a_side, Real (*a_func)(const RealVect a_pos));
      void set_poisson_wall_func(Real (*a_func)(const RealVect a_pos));
  
      // IO routines
      void write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const override;
      void read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) override;
      void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const override;
      int  get_num_plot_vars() const override;
  
      // Advance routines
      virtual Real advance(const Real a_dt) = 0;
      void compute_dt(Real& a_dt, time_code::which_code& a_timecode) override;
      void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) override;
      void print_step_report() override;

      // New regrid routines
      void register_operators() override; 
      void pre_regrid(const int a_lmin, const int a_old_finest_level) override;
      void deallocate() override;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;

      // Internal routines
      void set_verbosity(const int a_verbosity);

      // Set routines
      void set_ito(RefCountedPtr<ito_layout<ito_solver> >& a_cdr);
      void set_poisson(RefCountedPtr<poisson_solver>& a_poisson);
      void set_rte(RefCountedPtr<rte_layout<mc_photo> >& a_rte);
      void set_potential(Real (*a_potential)(const Real a_time));

      // Get routines
      Real get_time() const;

      // Poisson and related routines
      void compute_E(MFAMRCellData& a_E,      const MFAMRCellData& a_potential);
      void compute_E(EBAMRCellData& a_E,      const phase::which_phase a_phase);
      void compute_E(EBAMRCellData& a_E,      const phase::which_phase a_phase, const MFAMRCellData& a_potential);
      void compute_E(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell);
      void compute_E(EBAMRIVData& a_E_eb,     const phase::which_phase a_phase, const EBAMRCellData& a_E_cell);
      
      void compute_rho();
      void compute_rho(EBAMRCellData& a_rho, const phase::which_phase a_phase);
      void compute_rho(MFAMRCellData&                                     a_rho,
		       const Vector<particle_container<ito_particle>* >&  a_densities,
		       const centering::which_center                      a_centering);
      
      bool solve_poisson();
      bool solve_poisson(MFAMRCellData&                                   a_potential,
			 MFAMRCellData&                                   a_rhs,
			 const Vector<particle_container<ito_particle>* > a_densities,
			 const EBAMRIVData&                               a_sigma,
			 const centering::which_center                    a_centering);

    protected:

      /*!
	@brief Class name
      */
      std::string m_name;

      /*!
	@brief Phase
      */
      phase::which_phase m_phase;
      

      /*!
	@brief Physics
      */
      RefCountedPtr<ito_plasma_physics> m_physics;

      /*!
	@brief Ito solvers
      */
      RefCountedPtr<ito_layout<ito_solver> > m_ito;

      /*!
	@brief Radiative transfer solvers
      */
      RefCountedPtr<rte_layout<mc_photo> > m_rte;

      /*!
	@brief Poisson solver
      */
      RefCountedPtr<poisson_solver> m_poisson;

      /*!
	@brief Surface charge solver
      */
      RefCountedPtr<sigma_solver> m_sigma;

      /*!
	@brief Potential curve
      */
      Real (*m_potential)(const Real a_time);

            /*!
	@brief Function pointer for Poisson BC on domain walls
      */
      Real (*m_wall_func_x_lo)(const RealVect a_pos);

      /*!
	@brief Function pointer for Poisson BC on domain walls
      */
      Real (*m_wall_func_x_hi)(const RealVect a_pos);

      /*!
	@brief Function pointer for Poisson BC on domain walls
      */
      Real (*m_wall_func_y_lo)(const RealVect a_pos);

      /*!
	@brief Function pointer for Poisson BC on domain walls
      */
      Real (*m_wall_func_y_hi)(const RealVect a_pos);

      /*!
	@brief Function pointer for Poisson BC on domain walls
      */
      Real (*m_wall_func_z_lo)(const RealVect a_pos);

      /*!
	@brief Function pointer for Poisson BC on domain walls
      */
      Real (*m_wall_func_z_hi)(const RealVect a_pos);
    };
  }
};
#endif
