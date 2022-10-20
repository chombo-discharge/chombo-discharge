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
      @brief Implementation of ItoPlasmaStepper that uses a semi-implicit split-step formalism for advancing the Ito-Poisson-KMC system. 
    */
    class ItoPlasmaGodunovStepper : public ItoPlasmaStepper
    {
    public:
      /*!
	@brief Disallowed default constructor. Use the full constructor.
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

      /*!
	@brief Advance the Ito-Poisson-KMC system over a_dt
	@param[in] a_dt Time step
      */
      virtual Real
      advance(const Real a_dt) override;

      /*!
	@brief Allocate storage required for advancing the equations
      */
      virtual void
      allocate() override;

      /*!
	@brief Parse options
      */
      virtual void
      parseOptions() override;

      /*!
	@brief Parse run-time options
      */
      virtual void
      parseRuntimeOptions() override;

      /*!
	@brief Allocate internal storage
      */
      virtual void
      allocateInternals() override;

      /*!
	@brief Perform pre-regrid operations.
	@details This stores the particles and mesh data on the old grids. 
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@note This will store the "conductivity" and "rho dagger" particles. 
      */
      virtual void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;

      /*!
	@brief Regrid methods -- puts all data on the new mesh. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

    protected:
      /*!
	@brief Simple enum for distinguishing between algorithms
	@note Trapezoidal is a work-in-progress.
      */
      enum class WhichAlgorithm
      {
        EulerMaruyama,
        Trapezoidal
      };

      /*!
	@brief Average CFL number
      */
      Real m_averageCFL;

      /*!
	@brief Previous time step.
      */
      Real m_prevDt;

      /*!
	@brief Last 20 CFL numbers. Used for computing the average CFL.
      */
      std::deque<Real> m_cflNumbers;

      /*!
	@brief Which advancement algorithm to use.
      */
      WhichAlgorithm m_algorithm;

      /*!
	@brief Storage for simplified particles that gave us sigma^k
      */
      Vector<ParticleContainer<PointParticle>*> m_conductivityParticles;

      /*!
	@brief Storage for particles that gave rho^dagger
      */
      Vector<ParticleContainer<PointParticle>*> m_rhoDaggerParticles;

      /*!
	@brief Parse advancement algorithm
      */
      virtual void
      parseAlgorithm() noexcept;

      /*!
	@brief Set up extra particle storage for all Ito particles. 
      */
      virtual void
      setRuntimeParticleStorage() noexcept;

      /*!
	@brief Discard extra patricle storge for all Ito particles
      */
      virtual void
      resetRuntimeParticleStorage() noexcept;

      /*!
	@brief Set the starting positions for the ItoSolver particles
      */
      virtual void
      setOldPositions() noexcept;

      /*!
	@brief Remap the input point particles
	@param[inout] a_particles List of particle containers to remap. Indices must correspond to indices in the ItoSolvers
	@param[in]    a_subset    Particle subset to remap.
      */
      virtual void
      remapPointParticles(Vector<ParticleContainer<PointParticle>*>& a_particles, const SpeciesSubset a_subset) noexcept;

      /*!
	@brief Deposit the input point particles on the mesh
	@param[inout] a_particles List of particle containers to remap. Indices must correspond to indices in the ItoSolvers
	@param[in]    a_subset    Particle subset to remap. 
	@note This deposits into the ItoSolver m_phi data. 
      */     
      virtual void
      depositPointParticles(const Vector<ParticleContainer<PointParticle>*>& a_particles, const SpeciesSubset a_subset) noexcept;

      /*!
	@brief Clear the input particle data holders
	@param[inout] a_particles List of particle containers to remap. Indices must correspond to indices in the ItoSolvers
	@param[in]    a_subset    Particle subset to remap. 
      */
      virtual void
      clearPointParticles(const Vector<ParticleContainer<PointParticle>*>& a_particles, const SpeciesSubset a_subset) noexcept;

      /*!
	@brief Compute all conductivities (cell, face, and EB) from the input point particles
	@param[in] a_particles Particles to deposit. Particle weights must hold the mobility*weight
	@note This deposits the input particle weights onto m_conductivityCell, m_conductivityFace etc. 
      */
      virtual void
      computeConductivities(const Vector<ParticleContainer<PointParticle>*>& a_particles) noexcept;

      /*!
	@brief Compute the cell-centered conductivity
	@param[out] a_conductivityCell Cell-centered conductivity
	@param[in]  a_particles        Particles to deposit. Particle weights must hold the mobility*weight
      */
      virtual void
      computeCellConductivity(EBAMRCellData&                                   a_conductivityCell,
			      const Vector<ParticleContainer<PointParticle>*>& a_particles) noexcept;

      /*!
	@brief Compute the cell-centered conductivity
	@note Interpolates m_conductivityCell to m_conductivityFace
      */
      virtual void
      computeFaceConductivity() noexcept;

      /*!
	@brief Set up the semi-implicit Poisson solver
      */
      virtual void
      setupSemiImplicitPoisson(const Real a_dt) noexcept;

      /*!
	@brief Copy particles from the ItoSolver into PointParticles whose weight are ItoParticle::m_weight * ItoParticle::m_mobility
	@param[out] a_conductivityParticles Point-particles holding the mobility * weight
      */
      virtual void
      copyConductivityParticles(Vector<ParticleContainer<PointParticle>*>& a_conductivityParticles) noexcept;

      /*!
	@brief Copy particles from the ItoSolver (if they are charged) into PointParticles whose weight are ItoParticle::m_weight.
	@param[out] a_conductivityParticles Point-particles holding the mobility * weight     
      */
      virtual void
      copyRhoDaggerParticles(Vector<ParticleContainer<PointParticle>*>& a_rhoDaggerParticles) noexcept;
      
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
