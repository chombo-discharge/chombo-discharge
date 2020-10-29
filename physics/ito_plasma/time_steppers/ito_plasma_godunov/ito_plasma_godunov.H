/*!
  @file   ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
  @brief  Abstract class for integrating the Ito plasma equations
*/

#ifndef _ITO_PLASMA_GODUNOV_
#define _ITO_PLASMA_GODUNOV_

#include "ito_plasma_stepper.H"
#include "godunov_particle.H"

namespace physics {
  namespace ito_plasma {

    class ito_plasma_godunov : public ito_plasma_stepper {
    public:

      // Ctor and dtor
      ito_plasma_godunov() = delete;
      ito_plasma_godunov(RefCountedPtr<ito_plasma_physics>& a_physics);
      ~ito_plasma_godunov();

      // IO routines
      int get_num_plot_vars() const override;
      void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const override;
      void write_conductivity(EBAMRCellData& a_output, int& a_icomp) const;

      // Advance routines
      Real advance(const Real a_dt) override;
      void compute_dt(Real& a_dt, time_code& a_timecode) override;

      // Setup routines
      void allocate() override;
      void parse_options() override;
      void allocate_internals() override;

      // Regrid routines
      void pre_regrid(const int a_lmin, const int a_old_finest_level) override;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;


    protected:

      // Switching between algorithms. 
      enum class which_algorithm {
	euler_maruyama,
	trapezoidal,
      };

      // How to restrict dt. 
      enum class which_dt {
	advection,
	diffusion,
	advection_diffusion
      };

      // How to deposit particles. 
      enum class which_particles {
	all,
	all_mobile,
	all_diffusive,
	charged_mobile,
	charged_diffusive,
	all_mobile_or_diffusive,
	charged_and_mobile_or_diffusive
      };


      // These are things that are used to compute the conductivity. 
      EBAMRCellData m_scratch1;
      EBAMRCellData m_scratch2;
      EBAMRCellData m_conduct_cell;
      EBAMRCellData m_cache; // Used to store the cell-centered conductivity
      EBAMRFluxData m_conduct_face;
      EBAMRIVData   m_conduct_eb;

      bool m_profile;
      
      Real m_relax_factor; // Factor for relaxation time
      Real m_avg_cfl;      // Average CFL number
      Real m_prevDt;

      which_algorithm m_algorithm;
      which_dt        m_whichDt;

      Vector<particle_container<godunov_particle>* > m_conductivity_particles; // Used to store the particle that gave the sigma^k
      Vector<particle_container<godunov_particle>* > m_rho_dagger_particles;   // Used to store the particles that gave rho^dagger

      int m_merge_interval; // How often to merge superparticles

      // Routine for setting up extra particle storage
      void setup_runtime_storage();

      // Semi-implicit regrid routine. 
      void regrid_si(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

      // Global routines
      void set_old_positions();
      void rewind_particles();
      void rewind_diffusive_particles();
      void intersect_particles(const Real a_dt);
      void remap_godunov_particles(Vector<particle_container<godunov_particle>* >& a_particles, const which_particles a_which);
      void deposit_godunov_particles(const Vector<particle_container<godunov_particle>* >& a_particles, const which_particles a_which);


      // Routines for semi-implicit Poisson solution
      void compute_conductivity();
      void compute_face_conductivity();
      void setup_semi_implicit_poisson(const Real a_dt);

      // Routines for standard Poisson
      void setup_standard_poisson();

      // Deposition functions for scratch particles


      // Routines for semi-implicit regridding
      void copy_conductivity_particles();
      void copy_rho_dagger_particles();
      void compute_regrid_conductivity();
      void compute_regrid_rho();

      // Specific routines for Euler-Maruyama
      void advance_particles_euler_maruyama(const Real a_dt);
      void advect_particles_euler_maruyama(const Real a_dt);
      void diffuse_particles_euler_maruyama(const Real a_dt);
      void swap_euler_maruyama();

      // Specific routines for trapezoidal integration
      void advance_particles_trapezoidal(const Real a_dt);
    };
  }
}
#endif
