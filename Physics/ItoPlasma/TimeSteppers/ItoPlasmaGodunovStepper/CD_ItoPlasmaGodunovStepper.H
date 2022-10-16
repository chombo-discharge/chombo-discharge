/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaGodunovStepper.H
  @brief  Declaration of a class which uses a semi-implicit Godunov method for Ito plasma equations
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaGodunovStepper_H
#define CD_ItoPlasmaGodunovStepper_H

// Our includes
#include <CD_ItoPlasmaStepper.H>
#include <CD_PointParticle.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief Implementation of ItoPlasmaStepper that uses a semi-implicit split-step formalism for advancing the Ito-Poisson equations
    */
    class ItoPlasmaGodunovStepper : public ItoPlasmaStepper
    {
    public:
      /*!
	@brief Default constructor - disallowed
      */
      ItoPlasmaGodunovStepper() = delete;

      /*!
	@brief Full constructor. 
	@param[in] a_phyics Physics implementation. 
      */
      ItoPlasmaGodunovStepper(RefCountedPtr<ItoPlasmaPhysics>& a_physics);

      /*!
	@brief Destructor. Does nothing
      */
      virtual ~ItoPlasmaGodunovStepper();

      // Advance routines
      Real
      advance(const Real a_dt) override;
      Real
      computeDt() override;

      // Setup routines
      void
      allocate() override;
      void
      parseOptions() override;
      void
      parseRuntimeOptions() override;
      void
      allocateInternals() override;

      // Regrid routines
      void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;
      void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

    protected:
      // Switching between algorithms.
      enum class which_algorithm
      {
        euler_maruyama,
        trapezoidal,
      };

      // How to restrict dt.
      enum class which_dt
      {
        advection,
        diffusion,
        AdvectionDiffusion
      };

      // These are things that are used to compute the conductivity.
      EBAMRCellData m_scratch1;
      EBAMRCellData m_scratch2;

      EBAMRCellData m_cache; // Used to store the cell-centered conductivity

      bool m_profile;
      bool m_filter_cond;

      Real m_relax_factor; // Factor for relaxation time
      Real m_averageCFL;   // Average CFL number
      Real m_prevDt;
      Real m_eb_tolerance; // Tolerance for EB particle removal.

      which_algorithm m_algorithm;
      which_dt        m_whichDt;

      Vector<ParticleContainer<PointParticle>*>
        m_conductivity_particles; // Used to store the particle that gave the sigma^k
      Vector<ParticleContainer<PointParticle>*>
        m_rho_dagger_particles; // Used to store the particles that gave rho^dagger

      // Routine for setting up extra particle storage
      void
      setupRuntimeStorage();

      // Global routines
      void
      setOldPositions();
      void
      remapGodunovParticles(Vector<ParticleContainer<PointParticle>*>& a_particles, const SpeciesSubset a_which);
      void
      deposit_PointParticles(const Vector<ParticleContainer<PointParticle>*>& a_particles, const SpeciesSubset a_which);
      void
      clearGodunovParticles(const Vector<ParticleContainer<PointParticle>*>& a_particles, const SpeciesSubset a_which);

      // Routines for semi-implicit Poisson solution
      void
      computeAllConductivities(const Vector<ParticleContainer<PointParticle>*>& a_particles);
      void
      compute_cell_conductivity(EBAMRCellData&                                   a_conductivity,
                                const Vector<ParticleContainer<PointParticle>*>& a_particles);
      void
      compute_face_conductivity();
      void
      setupSemiImplicitPoisson(const Real a_dt);

      // Routines for standard Poisson
      void
      setupStandardPoisson();

      // Routines for semi-implicit regridding
      void
      copyConductivityParticles(Vector<ParticleContainer<PointParticle>*>& a_conductivity_particles);
      void
      copyRhoDaggerParticles(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger_particles);
      void
      computeRegridConductivity();
      void
      computeRegridRho();

      // Specific routines for Euler-Maruyama
      void
      advanceParticlesEulerMaruyama(const Real a_dt);
      void
      diffuseParticlesEulerMaruyama(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger, const Real a_dt);
      void
      stepEulerMaruyama(const Real a_dt);

      // Specific routines for trapezoidal integration
      void
      advanceParticlesTrapezoidal(const Real a_dt);
      void
      preTrapezoidalPredictor(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger, const Real a_dt);
      void
      preTrapezoidalCorrector(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger, const Real a_dt);
      void
      trapezoidalPredictor(const Real a_dt);
      void
      trapezoidalCorrector(const Real a_dt);
    };
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
