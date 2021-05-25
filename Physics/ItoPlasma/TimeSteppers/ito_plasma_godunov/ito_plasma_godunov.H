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

#include "CD_NamespaceHeader.H"
namespace Physics {
  namespace ItoPlasma {

    class ito_plasma_godunov : public ito_plasma_stepper {
    public:

      // Ctor and dtor
      ito_plasma_godunov() = delete;
      ito_plasma_godunov(RefCountedPtr<ItoPlasmaPhysics>& a_physics);
      ~ito_plasma_godunov();

      // IO routines
      int getNumberOfPlotVariables() const override;
      void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;
      void write_conductivity(EBAMRCellData& a_output, int& a_icomp) const;

      // Advance routines
      Real advance(const Real a_dt) override;
      void computeDt(Real& a_dt, TimeCode& a_timeCode) override;

      // Setup routines
      void allocate() override;
      void parseOptions() override;
      void parseRuntimeOptions() override;
      void allocateInternals() override;

      // Regrid routines
      void preRegrid(const int a_lmin, const int a_oldFinestLevel) override;
      void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;


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
	AdvectionDiffusion
      };

      // These are things that are used to compute the conductivity. 
      EBAMRCellData m_scratch1;
      EBAMRCellData m_scratch2;
      EBAMRCellData m_conduct_cell;
      EBAMRCellData m_cache; // Used to store the cell-centered conductivity
      EBAMRFluxData m_conduct_face;
      EBAMRIVData   m_conduct_eb;

      bool m_profile;
      bool m_filter_cond;
      
      Real m_relax_factor; // Factor for relaxation time
      Real m_avg_cfl;      // Average CFL number
      Real m_prevDt;
      Real m_eb_tolerance; // Tolerance for EB particle removal. 

      which_algorithm m_algorithm;
      which_dt        m_whichDt;

      Vector<ParticleContainer<godunov_particle>* > m_conductivity_particles; // Used to store the particle that gave the sigma^k
      Vector<ParticleContainer<godunov_particle>* > m_rho_dagger_particles;   // Used to store the particles that gave rho^dagger

      int m_merge_interval; // How often to merge superparticles

      // Routine for setting up extra particle storage
      void setup_runtime_storage();

      // Global routines
      void set_old_positions();
      void remap_godunov_particles(Vector<ParticleContainer<godunov_particle>* >& a_particles, const which_particles a_which);
      void deposit_godunov_particles(const Vector<ParticleContainer<godunov_particle>* >& a_particles, const which_particles a_which);
      void clear_godunov_particles(const Vector<ParticleContainer<godunov_particle>* >& a_particles, const which_particles a_which);

      // Routines for semi-implicit Poisson solution
      void compute_all_conductivities(const Vector<ParticleContainer<godunov_particle>* >& a_particles);
      void compute_cell_conductivity(EBAMRCellData& a_conductivity, const Vector<ParticleContainer<godunov_particle>* >& a_particles);
      void compute_face_conductivity();
      void setup_semi_implicit_poisson(const Real a_dt);

      // Routines for standard Poisson
      void setup_standard_poisson();

      // Routines for semi-implicit regridding
      void copy_conductivity_particles(Vector<ParticleContainer<godunov_particle>* >& a_conductivity_particles);
      void copy_rho_dagger_particles(Vector<ParticleContainer<godunov_particle>* >& a_rho_dagger_particles);
      void compute_regrid_conductivity();
      void compute_regrid_rho();

      // Specific routines for Euler-Maruyama
      void advanceParticles_euler_maruyama(const Real a_dt);
      void diffuse_particles_euler_maruyama(Vector<ParticleContainer<godunov_particle>* >& a_rho_dagger, const Real a_dt);
      void step_euler_maruyama(const Real a_dt);
      
      // Specific routines for trapezoidal integration
      void advanceParticles_trapezoidal(const Real a_dt);
      void pre_trapezoidal_predictor(Vector<ParticleContainer<godunov_particle>* >& a_rho_dagger, const Real a_dt);
      void pre_trapezoidal_corrector(Vector<ParticleContainer<godunov_particle>* >& a_rho_dagger, const Real a_dt);
      void trapezoidal_predictor(const Real a_dt);
      void trapezoidal_corrector(const Real a_dt);
    };
  }
}
#include "CD_NamespaceFooter.H"
#endif
