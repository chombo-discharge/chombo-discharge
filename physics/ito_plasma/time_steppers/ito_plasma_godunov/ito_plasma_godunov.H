/*!
  @file   ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
  @brief  Abstract class for integrating the Ito plasma equations
*/

#ifndef _ITO_PLASMA_GODUNOV_
#define _ITO_PLASMA_GODUNOV_

#include "ito_plasma_stepper.H"

namespace physics {
  namespace ito_plasma {

    class ito_plasma_godunov : public ito_plasma_stepper {
    public:

      // Ctor and dtor
      ito_plasma_godunov();
      ito_plasma_godunov(RefCountedPtr<ito_plasma_physics>& a_physics);
      ~ito_plasma_godunov();

      // Advance routines
      Real advance(const Real a_dt) override;
      void compute_dt(Real& a_dt, time_code& a_timecode) override;

      // Setup routines
      void parse_options() override;
      void allocate_internals() override;

      // Regrid routines
      void pre_regrid(const int a_lmin, const int a_old_finest_level) override;
      void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) override;

    protected:

      enum class which_algorithm {
	euler, 
	predictor_corrector,
	semi_implicit,
      };

      // These are things that are used to compute the conductivity. 
      EBAMRCellData m_scratch1;
      EBAMRCellData m_scratch2;
      EBAMRCellData m_conduct_cell;
      EBAMRFluxData m_conduct_face;
      EBAMRIVData   m_conduct_eb;

      bool m_use_old_dt;
      
      Real m_relax_factor; // Factor for relaxation time

      which_algorithm m_algorithm;

      int m_merge_interval; // How often to merge superparticles

      // Global routines
      void intersect_particles(const Real a_dt);
      void set_old_positions();

      void rewind_particles();
      void diffuse_particles_euler(const Real a_dt);
      void advect_particles_euler(const Real a_dt);
      void advect_particles_rk2(const Real a_dt);      

      // Specific routines for semi-implicit Euler
      void compute_conductivity();
      void swap_particle_positions();
      void setup_semi_implicit_poisson(const Real a_dt);
      void advect_particles_si(const Real a_dt);

      // Particle advancement routines. 
      void advance_particles_euler(const Real a_dt);
      void advance_particles_pc(const Real a_dt);
      void advance_particles_si(const Real a_dt);

    };
  }
}
#endif
