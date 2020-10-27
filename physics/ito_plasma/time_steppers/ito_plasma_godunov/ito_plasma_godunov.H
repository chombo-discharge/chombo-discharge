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
	euler,
	euler_maruyama,
	midpoint,
      };

      // How to restrict dt. 
      enum class which_dt {
	advection,
	diffusion,
	advection_diffusion
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
      void intersect_particles(const Real a_dt);

      // Routines for semi-implicit Poisson solution
      void compute_conductivity();
      void compute_face_conductivity();
      void setup_semi_implicit_poisson(const Real a_dt);

      // Routines for standard Poisson
      void setup_standard_poisson();

      // Routines for semi-implicit regridding
      void copy_conductivity_particles();
      void copy_rho_dagger_particles();
      void compute_regrid_conductivity();
      void compute_regrid_rho();

      // Specific routines for Euler method
      void advance_particles_euler(const Real a_dt);
      void advect_particles_euler(const Real a_dt);
      void diffuse_particles_euler(const Real a_dt);

      // Specific routines for Euler-Maruyama
      void advance_particles_euler_maruyama(const Real a_dt);
      void advect_particles_euler_maruyama(const Real a_dt);
      void diffuse_particles_euler_maruyama(const Real a_dt);
      void swap_euler_maruyama();

      // Specific routines for midpoint integration
      void advance_particles_midpoint(const Real a_dt);
      void diffuse_particles_midpoint1(const Real a_dt);
      void diffuse_particles_midpoint2(const Real a_dt);
      void advect_particles_midpoint1(const Real a_dt);
      void advect_particles_midpoint2(const Real a_dt);
      void swap_midpoint1();
      void swap_midpoint2();
    };
  }
}
#endif
