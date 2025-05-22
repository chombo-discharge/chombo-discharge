/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCGodunovStepper.H
  @brief  Declaration of a class which uses a semi-implicit Godunov method for Ito plasma equations
  @author Robert Marskar
*/

#ifndef CD_ItoKMCGodunovStepper_H
#define CD_ItoKMCGodunovStepper_H

// Our includes
#include <CD_ItoKMCStepper.H>
#include <CD_Timer.H>
#include <CD_PointParticle.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Implementation of ItoKMCStepper that uses a semi-implicit split-step formalism for advancing the Ito-Poisson-KMC system. 
    */
    template <typename I = ItoSolver, typename C = CdrCTU, typename R = McPhoto, typename F = FieldSolverMultigrid>
    class ItoKMCGodunovStepper : public ItoKMCStepper<I, C, R, F>
    {
    public:
      /*!
	@brief Disallowed default constructor. Use the full constructor.
      */
      ItoKMCGodunovStepper() = delete;

      /*!
	@brief Full constructor. 
	@param[in] a_phyics Physics implementation. 
      */
      ItoKMCGodunovStepper(RefCountedPtr<ItoKMCPhysics>& a_physics);

      /*!
	@brief Destructor. Does nothing
      */
      virtual ~ItoKMCGodunovStepper();

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
      allocate() noexcept override;

      /*!
	@brief Parse options
      */
      virtual void
      parseOptions() noexcept override;

      /*!
	@brief Parse run-time options
      */
      virtual void
      parseRuntimeOptions() noexcept override;

      /*!
	@brief Perform pre-regrid operations.
	@details This stores the particles and mesh data on the old grids. 
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@note This will store the "conductivity" and "rho dagger" particles. 
      */
      virtual void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) noexcept override;

      /*!
	@brief Regrid methods -- puts all data on the new mesh. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept override;

#ifdef CH_USE_HDF5
      /*!  
	@brief Write checkpoint data to file
	@param[inout] a_handle HDF5 fil
	@param[in]    a_lvl    Grid level
      */
      virtual void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const noexcept override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read header data from checkpoint file.
	@param[inout] a_header HDF5 header.
      */
      virtual void
      readCheckpointHeader(HDF5HeaderData& a_header) noexcept override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Write header data to checkpoint file.
	@param[inout] a_header HDF5 header.
      */
      virtual void
      writeCheckpointHeader(HDF5HeaderData& a_header) const noexcept override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) noexcept override;
#endif

      /*!
	@brief Perform post-plot operations
      */
      virtual void
      postPlot() noexcept override;

    protected:
      /*!
	@brief Simple enum for distinguishing between algorithms
      */
      enum class WhichAlgorithm
      {
        EulerMaruyama,
      };

      /*!
	@brief If true, then the particles are checkpointed so we can regrid on checkpoint-restart.
      */
      bool m_writeCheckpointParticles;

      /*!
	@brief If true, then the HDF5 checkpoint file contained particles that we can read.
      */
      bool m_readCheckpointParticles;

      /*!
	@brief If true, then the class supports regrid-on-restart
	@details This is populated in the constructor or when reading the checkpoint file. 
      */
      bool m_canRegridOnRestart;

      /*!
	@brief For achieving a slightly smoother gradient in the conductivity near the EB.
      */
      bool m_extendConductivityEB;

      /*!
	@brief Smooth the conductivity or not.
      */
      bool m_smoothConductivity;

      /*!
	@brief If true, the algorithm will limit particle diffusion against the drift direction
      */
      bool m_limitParallelDiffusion;

      /*!
	@brief If true, particles will be emitted before the reactive step
      */
      bool m_emitSecondaryParticlesBeforeReactions;

      /*!
	@brief Which advancement algorithm to use.
      */
      WhichAlgorithm m_algorithm;

      /*!
	@brief Timer used for run-time logging of routines
      */
      Timer m_timer;

      /*!
	@brief Number of filterings when computing the semi-implicit space charge density
      */
      int m_filterNum;

      /*!
	@brief Number of filterings when computing the semi-implicit space charge density
      */
      int m_filterMaxStride;

      /*!
	@brief Factor for strided filter. Alpha=0.5 yields a bilinear filter
      */
      Real m_filterAlpha;

      /*!
	@brief Storage for simplified particles that gave us sigma^k
	@note Needed for regrids
      */
      Vector<RefCountedPtr<ParticleContainer<PointParticle>>> m_conductivityParticles;

      /*!
	@brief Storage for particles that fell inside the EB but should still contribute to the conductivity
      */
      Vector<RefCountedPtr<ParticleContainer<PointParticle>>> m_irregularParticles;

      /*!
	@brief Storage for particles that gave rho^dagger
	@note Needed for regrids
      */
      Vector<RefCountedPtr<ParticleContainer<PointParticle>>> m_rhoDaggerParticles;

      /*!
	@brief Storage for CDR densities used during the semi-implicit solve. 
	@note This is needed because the semi-implicit requires the space charge from
	before the reaction step. 
      */
      EBAMRCellData m_semiImplicitRhoCDR;

      /*!
	@brief Storage for conductivity term due to mobile CDR species. 
	@note This is needed because the semi-implicit requires the conductivity when solving
	the Poisson equation.
      */
      EBAMRCellData m_semiImplicitConductivityCDR;

      /*!
	@brief Scratch storage for CDR contribution to space charge density.
	@note Only used during regrids. 
      */
      EBAMRCellData m_scratchSemiImplicitRhoCDR;

      /*!
	@brief Scratch storage for CDR contribution to conductivity
	@note Only used during regrids. 
      */
      EBAMRCellData m_scratchSemiImplicitConductivityCDR;

      /*!
	@brief Storage for the semi-implicit approximation to the surface charge
      */
      EBAMRIVData m_semiImplicitSigma;

      /*!
	@brief Scratch for the semi-implicit approximation to the surface charge
	@note Only used during regrids
      */
      EBAMRIVData m_scratchSemiImplicitSigma;

      /*!
	@brief Storage for the finite-volume approximation of div(D*grad(phi)) for the CDR equations.
      */
      Vector<EBAMRCellData> m_cdrDivD;

      /*!
	@brief Allocate "internal" storage.
      */
      virtual void
      allocateInternals() noexcept override;

      /*!
	@brief Parse advancement algorithm
      */
      virtual void
      parseAlgorithm() noexcept;

      /*!
	@brief Parse filter settings
      */
      virtual void
      parseFiltering() noexcept;

      /*!
	@brief Parse checkpoint-restart functionality.
      */
      virtual void
      parseCheckpointParticles() noexcept;

      /*!
	@brief Parse when secondary particles are emitted
      */
      virtual void
      parseSecondaryEmissionSpecification() noexcept;

      /*!
	@brief Set the starting positions for the ItoSolver particles
      */
      virtual void
      setOldPositions() noexcept;

      /*!
	@brief Set an MPI barrier if using debug mode. 
	@details This defaults to ParallelOps::barrier() if m_debug=true
      */
      virtual void
      barrier() const noexcept;

      /*!
	@brief Remap the input point particles
	@param[inout] a_particles List of particle containers to remap. Indices must correspond to indices in the ItoSolvers
	@param[in]    a_subset    Particle subset to remap.
      */
      virtual void
      remapPointParticles(Vector<RefCountedPtr<ParticleContainer<PointParticle>>>& a_particles,
                          const SpeciesSubset                                      a_subset) noexcept;

      /*!
	@brief Deposit the input point particles on the mesh
	@param[inout] a_particles List of particle containers to remap. Indices must correspond to indices in the ItoSolvers
	@param[in]    a_subset    Particle subset to remap. 
	@note This deposits into the ItoSolver m_phi data. 
      */
      virtual void
      depositPointParticles(const Vector<RefCountedPtr<ParticleContainer<PointParticle>>>& a_particles,
                            const SpeciesSubset                                            a_subset) noexcept;

      /*!
	@brief Clear the input particle data holders
	@param[inout] a_particles List of particle containers to remap. Indices must correspond to indices in the ItoSolvers
	@param[in]    a_subset    Particle subset to remap. 
      */
      virtual void
      clearPointParticles(const Vector<RefCountedPtr<ParticleContainer<PointParticle>>>& a_particles,
                          const SpeciesSubset                                            a_subset) noexcept;

      /*!
	@brief Compute all conductivities (cell, face, and EB) from the input point particles
	@param[in] a_particles Particles to deposit. Particle weights must hold the mobility*weight
	@note This deposits the input particle weights onto m_conductivityCell, m_conductivityFace.
	@note Also adds in the conductive contribution from the CDR solvers. 
      */
      virtual void
      computeConductivities(const Vector<RefCountedPtr<ParticleContainer<PointParticle>>>& a_particles) noexcept;

      /*!
	@brief Compute the cell-centered conductivity
	@param[out] a_conductivityCell Cell-centered conductivity
	@param[in]  a_particles        Particles to deposit. Particle weights must hold the mobility*weight
	@note Also adds in the conductive contribution from the CDR solvers. 
      */
      virtual void
      computeCellConductivity(EBAMRCellData&                                                 a_conductivityCell,
                              const Vector<RefCountedPtr<ParticleContainer<PointParticle>>>& a_particles) noexcept;

      /*!
	@brief Compute the cell-centered conductivity
	@note Interpolates m_conductivityCell to m_conductivityFace
      */
      virtual void
      computeFaceConductivity() noexcept;

      /*!
	@brief Set up the space charge density for the regrid operation
      */
      virtual void
      computeSemiImplicitRho() noexcept;

      /*!
	@brief Set up the semi-implicit Poisson solver
      */
      virtual void
      setupSemiImplicitPoisson(const Real a_dt) noexcept;

      /*!
	@brief Solve the electrostatic problem
	@details This calls the other version with all lhs/rhs taken from the field solver. Space charge density is computed
	in this routine. 
	@return Returns true if the solver converged. 
      */
      virtual bool
      solvePoisson() noexcept override;

      /*!
	@brief Remove covered particles
	@param[inout] a_particles      Particles to remove
	@param[in]    a_representation Representation of the EB
	@param[in]    a_tolerance      EB tolerance
      */
      virtual void
      removeCoveredPointParticles(Vector<RefCountedPtr<ParticleContainer<PointParticle>>>& a_particles,
                                  const EBRepresentation                                   a_representation,
                                  const Real                                               a_tolerance) const noexcept;

      /*!
	@brief Copy particles from the ItoSolver into PointParticles whose weight are ItoParticle::m_weight * ItoParticle::m_mobility
	@param[out] a_conductivityParticles Point-particles holding the mobility * weight
      */
      virtual void
      copyConductivityParticles(
        Vector<RefCountedPtr<ParticleContainer<PointParticle>>>& a_conductivityParticles) noexcept;

      /*!
	@brief Advance the particles using the Euler-Maruyama scheme
	@param[in] a_dt Time step
      */
      virtual void
      advanceEulerMaruyama(const Real a_dt) noexcept;

      /*!
	@brief Perform the diffusive Ito advance in the Euler-Maruyama step.
	@details This computes the diffusion hop onto run-time enabled storage on the ItoParticles. We also store the
	"rho dagger" particles which correspond to the position of the charged particles after the diffusionhop. 
      */
      virtual void
      diffuseParticlesEulerMaruyama(Vector<RefCountedPtr<ParticleContainer<PointParticle>>>& a_rhoDaggerParticles,
                                    const Real                                               a_dt) noexcept;

      /*!
	@brief Compute the diffusion term for the CDR equations as well as the resulting CDR-contributions to 
	the space charge density.
	@details This does an explicit diffusion advance over a_dt for the diffusive CDR species. This also
	stores the resultant space charge on m_semiImplicitRhoCDR.
      */
      virtual void
      computeDiffusionTermCDR(EBAMRCellData& m_semiImplicitRhoCDR, const Real a_dt) noexcept;

      /*!
	@brief Step the particles according to the regular Euler-Maruyama scheme.
	@param[in] a_dt Time step. 
	@details This sets X^(k+1) = X^k + a_dt * V + sqrt(2*D*a_dt) * W
	@note The hop MUST have been stored on the particle run-time storage 
      */
      virtual void
      stepEulerMaruyamaParticles(const Real a_dt) noexcept;

      /*!
	@brief Step the CDR equations according to the regular Euler-Maruyama scheme.
	@param[in] a_dt Time step. 
	@details This sets phi^(k+1) = phi^k - a_dt * div(v * phi) + a_dt * div(D*grad(phi))
      */
      virtual void
      stepEulerMaruyamaCDR(const Real a_dt) noexcept;

      /*!
	@brief Utility function for plotting the ItoSolver particles. These are written in a particles folder
      */
      virtual void
      plotParticles() const noexcept;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCGodunovStepperImplem.H>

#endif
