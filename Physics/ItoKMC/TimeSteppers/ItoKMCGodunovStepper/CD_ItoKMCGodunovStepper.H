/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ItoKMCGodunovStepper.H
 * @brief  Declaration of a class which uses a semi-implicit Godunov method for Ito plasma equations
 * @author Robert Marskar
 */

#ifndef CD_ITOKMCGODUNOVSTEPPER_H
#define CD_ITOKMCGODUNOVSTEPPER_H

// Our includes
#include <CD_ItoKMCStepper.H>
#include <CD_Timer.H>
#include <CD_ParticleContainer.H>
#include <CD_EBIntersection.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
namespace ItoKMC {

/**
 * @brief Implementation of ItoKMCStepper that uses a semi-implicit split-step formalism for advancing the
 * Ito-Poisson-KMC system.
 */
template <typename I = ItoSolver, typename C = CdrCTU, typename R = McPhoto, typename F = FieldSolverGMG>
class ItoKMCGodunovStepper : public ItoKMCStepper<I, C, R, F>
{
public:
  /**
   * @brief Disallowed default constructor. Use the full constructor.
   */
  ItoKMCGodunovStepper() = delete;

  /**
   * @brief Full constructor.
   * @param[in] a_physics      Physics implementation.
   * @param[in] a_parseOptions If true, parse options immediately. Set to false when constructing through a derived
   * class that will re-parse with its own name.
   */
  ItoKMCGodunovStepper(RefCountedPtr<ItoKMCPhysics>& a_physics, bool a_parseOptions = true);

  /**
   * @brief Destructor. Does nothing
   */
  virtual ~ItoKMCGodunovStepper();

  /**
   * @brief Compute a time step used for the advance method.
   * @return Time step satisfying all solver CFL and physics restrictions.
   */
  virtual Real
  computeDt() override;

  /**
   * @brief Advance the Ito-Poisson-KMC system over a_dt.
   * @param[in] a_dt Time step.
   * @return Time step that was used.
   */
  virtual Real
  advance(const Real a_dt) override;

  /**
   * @brief Allocate storage required for advancing the equations
   */
  virtual void
  allocate() noexcept override;

  /**
   * @brief Parse options
   */
  virtual void
  parseOptions() noexcept override;

  /**
   * @brief Parse run-time options
   */
  virtual void
  parseRuntimeOptions() noexcept override;

  /**
   * @brief Perform pre-regrid operations.
   * @details This stores the particles and mesh data on the old grids.
   * @param[in] a_lmin           The coarsest level that changes
   * @param[in] a_oldFinestLevel The finest level before the regrid.
   * @note This will store the "conductivity" and "rho dagger" particles.
   */
  virtual void
  preRegrid(const int a_lmin, const int a_oldFinestLevel) noexcept override;

  /**
   * @brief Register operators used for the simulation
   */
  virtual void
  registerOperators() noexcept override;

  /**
   * @brief Regrid methods -- puts all data on the new mesh.
   * @param[in] a_lmin           The coarsest level that changed.
   * @param[in] a_oldFinestLevel The finest level before the regrid.
   * @param[in] a_newFinestLevel The finest level after the regrid.
   */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept override;

#ifdef CH_USE_HDF5
  /**
   * @brief Write checkpoint data to file.
   * @param[in,out] a_handle HDF5 file.
   * @param[in]    a_lvl    Grid level.
   */
  virtual void
  writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const noexcept override;
#endif

#ifdef CH_USE_HDF5
  /**
   * @brief Read header data from checkpoint file.
   * @param[in,out] a_header HDF5 header.
   */
  virtual void
  readCheckpointHeader(HDF5HeaderData& a_header) noexcept override;
#endif

#ifdef CH_USE_HDF5
  /**
   * @brief Write header data to checkpoint file.
   * @param[in,out] a_header HDF5 header.
   */
  virtual void
  writeCheckpointHeader(HDF5HeaderData& a_header) const noexcept override;
#endif

#ifdef CH_USE_HDF5
  /**
   * @brief Read checkpoint data from file
   * @param[in,out] a_handle HDF5 file
   * @param[in]    a_lvl    Grid level
   */
  virtual void
  readCheckpointData(HDF5Handle& a_handle, const int a_lvl) noexcept override;
#endif

  /**
   * @brief Perform post-plot operations
   */
  virtual void
  postPlot() noexcept override;

protected:
  /**
   * @brief Simple enum for distinguishing between algorithms
   */
  enum class WhichAlgorithm
  {
    EulerMaruyama,
  };

  /**
   * @brief What to do with a diffusion hop that would carry a particle, or its charge, inside the embedded boundary.
   * @details A deposit inside the solid falls in a covered cell, which carries no VoF and contributes nothing to rho,
   * so the charge leaves rho^dagger entirely. The geometry makes that one-way -- a particle in the fluid can hop into
   * the solid and there are none inside to hop back -- so the near-wall region drains at a rate set by
   * sqrt(2*D*dt)/dx rather than by any physics.
   * @details The settings differ in what they correct. Inside, Cancel and ReflectRho leave the particle's own
   * displacement alone, so a particle that hops into the solid stays there and is absorbed by the intersection test
   * later in the advance; only where its charge is deposited changes. ReflectIto corrects the displacement itself, so
   * the particle never reaches the surface by diffusion at all and the deposit follows it.
   * @note ReflectIto makes the wall reflecting for the diffusive part of the step while it stays absorbing for the
   * drift part, which is a split along the operator splitting rather than along anything physical. It is retained
   * because it is closer to the observed behaviour than the alternatives, not because the boundary condition is
   * settled.
   */
  enum class DiffusiveDeposit
  {
    Inside,     ///< Deposit from inside the solid, into a covered cell that contributes nothing to rho.
    Cancel,     ///< Deposit at the pre-hop position. Conserves charge, but biases density towards where it started.
    ReflectRho, ///< Reflect the deposit about the boundary. The particle's own displacement is untouched.
    ReflectIto, ///< Reflect the particle's displacement about the boundary. The deposit follows the corrected hop.
  };

  /**
   * @brief If true, then the particles are checkpointed so we can regrid on checkpoint-restart.
   */
  bool m_writeCheckpointParticles;

  /**
   * @brief If true, the stochastic diffusion hop enters rho^dagger along with the rest of the explicit displacement.
   * @details The hop has zero mean, so to leading order dropping it shifts the source of the semi-implicit Poisson
   * equation by dt*D*grad^2(n): the field no longer sees the charge spread over the step, which matters most where
   * diffusion competes with drift. What it buys is a quieter rho^dagger, since the deposited source no longer
   * carries the hop's shot noise. The lagged grad(D) drift is unaffected either way -- it is deterministic mean
   * transport and always enters rho^dagger.
   */
  bool m_rhoDaggerHop;

  /**
   * @brief If true, then the HDF5 checkpoint file contained particles that we can read.
   */
  bool m_readCheckpointParticles;

  /**
   * @brief If true, then the class supports regrid-on-restart
   * @details This is populated in the constructor or when reading the checkpoint file.
   */
  bool m_canRegridOnRestart;

  /**
   * @brief For achieving a slightly smoother gradient in the conductivity near the EB.
   */
  bool m_extendConductivityEB;

  /**
   * @brief If true, particles will be emitted before the reactive step
   */
  bool m_emitSecondaryParticlesBeforeReactions;

  /**
   * @brief What to do with a rho^dagger deposit that would land inside the embedded boundary. See DiffusiveDeposit.
   */
  DiffusiveDeposit m_diffusiveDeposit;

  /**
   * @brief Which advancement algorithm to use.
   */
  WhichAlgorithm m_algorithm;

  /**
   * @brief Timer used for run-time logging of routines
   */
  Timer m_timer;

  /**
   * @brief Number of filterings when computing the semi-implicit space charge density
   */
  int m_rhoFilterNum;

  /**
   * @brief Stride when filtering the space charge density.
   */
  int m_rhoFilterMaxStride;

  /**
   * @brief Factor for strided filter. Alpha=0.5 yields a bilinear filter
   */
  Real m_rhoFilterAlpha;

  /**
   * @brief Number of filterings of the conductivity.
   */
  int m_condFilterNum;

  /**
   * @brief Stride when filtering the conductivity.
   */
  int m_condFilterMaxStride;

  /**
   * @brief Factor for strided filter. Alpha=0.5 yields a bilinear filter
   */
  Real m_condFilterAlpha;

  /**
   * @brief Limit for maximum field abort
   */
  Real m_maxFieldAbort;

  /**
   * @brief Time-centering of the electric field that the chemistry is evaluated at.
   * @details The reactive substep is evaluated at (1-theta)*E^k + theta*E^(k+1), where theta is this
   * parameter and E^k/E^(k+1) are the fields from the start and end of the transport step. theta = 0.5
   * removes the leading-order bias from the semi-implicit screening of the field, theta = 1 evaluates at
   * the end-of-step field, and theta = 0 at the field the step started from. Restricted to [0,1].
   */
  Real m_reactiveFieldCentering;

  /**
   * @brief Storage for simplified particles that gave us sigma^k
   * @note Needed for regrids
   */
  Vector<RefCountedPtr<ParticleContainer<NoPayload>>> m_conductivityParticles;

  /**
   * @brief Storage for particles that fell inside the EB but should still contribute to the conductivity
   */
  Vector<RefCountedPtr<ParticleContainer<NoPayload>>> m_irregularParticles;

  /**
   * @brief Storage for particles that gave rho^dagger
   * @note Needed for regrids
   */
  Vector<RefCountedPtr<ParticleContainer<NoPayload>>> m_rhoDaggerParticles;

  /**
   * @brief Storage for CDR densities used during the semi-implicit solve.
   * @note This is needed because the semi-implicit requires the space charge from
   * before the reaction step.
   */
  EBAMRCellData m_semiImplicitRhoCDR;

  /**
   * @brief Storage for conductivity term due to mobile CDR species.
   * @note This is needed because the semi-implicit requires the conductivity when solving
   * the Poisson equation.
   */
  EBAMRCellData m_semiImplicitConductivityCDR;

  /**
   * @brief Scratch storage for CDR contribution to space charge density.
   * @note Only used during regrids.
   */
  EBAMRCellData m_scratchSemiImplicitRhoCDR;

  /**
   * @brief Scratch storage for CDR contribution to conductivity
   * @note Only used during regrids.
   */
  EBAMRCellData m_scratchSemiImplicitConductivityCDR;

  /**
   * @brief Storage for the finite-volume approximation of div(D*grad(phi)) for the CDR equations.
   */
  Vector<EBAMRCellData> m_cdrDivD;

  /**
   * @brief Storage for the electric field that the chemistry is evaluated at.
   * @details This holds E^k from the top of advance() until the transport step is done, and is then turned
   * into the interpolant (1-theta)*E^k + theta*E^(k+1) that is handed to advanceReactionNetwork.
   */
  EBAMRCellData m_reactiveElectricField;

  /**
   * @brief Allocate "internal" storage.
   */
  virtual void
  allocateInternals() noexcept override;

  /**
   * @brief Parse advancement algorithm
   */
  virtual void
  parseAlgorithm() noexcept;

  /**
   * @brief Parse filter settings
   */
  virtual void
  parseFiltering() noexcept;

  /**
   * @brief Parse checkpoint-restart functionality.
   */
  virtual void
  parseCheckpointParticles() noexcept;

  /**
   * @brief Parse when secondary particles are emitted
   */
  virtual void
  parseSecondaryEmissionSpecification() noexcept;

  /**
   * @brief Parse what happens to a diffusion hop that would cross the embedded boundary.
   */
  virtual void
  parseDiffusiveDeposit() noexcept;

  /**
   * @brief Parse whether the stochastic diffusion hop enters rho^dagger.
   */
  virtual void
  parseRhoDaggerHop() noexcept;

  /**
   * @brief Reflect a rho^dagger deposit displacement that crosses the embedded boundary back into the fluid.
   * @details Locates the crossing point with the solver's own intersection routine, takes the surface normal from the
   * gradient of the implicit function there, and mirrors the remainder of the displacement about it. The crossing
   * search always brackets from a_pos, which is known to lie in the fluid, so the bracket is valid on every attempt
   * and no search is ever started from a point sitting on the surface.
   * @details Reflection can fail to place the particle in the fluid -- at a concave corner, or in a feature thinner
   * than the hop, the mirrored point is inside the solid again. A second reflection is attempted, and if that also
   * fails the hop is cancelled, which is always safe because the pre-hop position is known to be in the fluid.
   * @note Whether a hop passing clean THROUGH thin geometry is caught is decided by the intersection routine and not
   * here. Bisection marches the path in intervals of a_bisectStep and sees such a feature only when it is thicker
   * than one interval; the caller's endpoint sign test never sees one at all. Sphere-tracing the particle would cover
   * the case for a genuine signed distance function, since the march is driven by the distance to the surface rather
   * than by a fixed step.
   * @param[in] a_pos             Pre-hop position, in the fluid.
   * @param[in] a_disp            Deposit displacement that crosses the boundary. This is not the particle's own
   * displacement -- the trajectory is never corrected, only where its charge is deposited.
   * @param[in] a_fOld            Implicit function at a_pos, so the fluid's sign is known without re-evaluating.
   * @param[in] a_dx              Grid resolution, used to size the gradient stencil and the ray-cast tolerance.
   * @param[in] a_intersectionAlg Intersection algorithm, taken from the solver that owns the particle.
   * @param[in] a_bisectStep      Bisection path-march length, taken from the solver that owns the particle.
   * @return The reflected displacement, or zero if reflection could not place the particle in the fluid.
   */
  virtual RealVect
  reflectDiffusionHop(const RealVect&      a_pos,
                      const RealVect&      a_disp,
                      const Real           a_fOld,
                      const Real           a_dx,
                      const EBIntersection a_intersectionAlg,
                      const Real           a_bisectStep) const noexcept;

  /**
   * @brief Parse the time-centering of the electric field used for the reactive substep.
   * @note Aborts if the requested centering falls outside [0,1].
   */
  virtual void
  parseReactiveFieldCentering() noexcept;

  /**
   * @brief Set the starting positions for the ItoSolver particles
   */
  virtual void
  setOldPositions() noexcept;

  /**
   * @brief Set an MPI barrier if using debug mode.
   * @details This defaults to ParallelOps::barrier() if m_debug=true
   */
  virtual void
  barrier() const noexcept;

  /**
   * @brief Remap the input point particles
   * @param[in,out] a_particles List of particle containers to remap. Indices must correspond to indices in the
   * ItoSolvers
   * @param[in]    a_subset    Particle subset to remap.
   */
  virtual void
  remapPointParticles(Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles,
                      const SpeciesSubset                                  a_subset) noexcept;

  /**
   * @brief Deposit the input point particles on the mesh
   * @param[in,out] a_particles List of particle containers to remap. Indices must correspond to indices in the
   * ItoSolvers
   * @param[in]    a_subset    Particle subset to remap.
   * @note This deposits into the ItoSolver m_phi data.
   */
  virtual void
  depositPointParticles(const Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles,
                        const SpeciesSubset                                        a_subset) noexcept;

  /**
   * @brief Clear the input particle data holders
   * @param[in,out] a_particles List of particle containers to remap. Indices must correspond to indices in the
   * ItoSolvers
   * @param[in]    a_subset    Particle subset to remap.
   */
  virtual void
  clearPointParticles(const Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles,
                      const SpeciesSubset                                        a_subset) noexcept;

  /**
   * @brief Compute the CDR contribution to the semi-implicit conductivity, i.e. sum(|Z| * phi * mu).
   * @details The result is stored on m_semiImplicitConductivityCDR, which is carried across regrids. It is
   * not scaled by Units::Qe -- that scaling is applied once to the total conductivity in
   * computeCellConductivity -- nor by dt/eps0, which setupSemiImplicitPoisson applies.
   * @note Reads m_cdrMobilities, so it may only be called when those hold the mobilities for the current
   * grids. That is not the case inside regrid(), where allocateInternals() has just reallocated them.
   */
  virtual void
  computeCdrConductivity() noexcept;

  /**
   * @brief Compute all conductivities (cell, face, and EB) from the input point particles
   * @param[in] a_particles Particles to deposit. Particle weights must hold the mobility*weight
   * @param[in] a_useStoredCdrConductivity If true, take the CDR contribution from
   * m_semiImplicitConductivityCDR as it stands; if false, recompute it first through computeCdrConductivity.
   * @note This deposits the input particle weights onto m_conductivityCell, m_conductivityFace.
   * @note Also adds in the conductive contribution from the CDR solvers.
   */
  virtual void
  computeConductivities(const Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles,
                        const bool                                                 a_useStoredCdrConductivity) noexcept;

  /**
   * @brief Compute the cell-centered conductivity
   * @param[out] a_conductivityCell Cell-centered conductivity
   * @param[in]  a_particles        Particles to deposit. Particle weights must hold the mobility*weight
   * @param[in]  a_useStoredCdrConductivity If true, take the CDR contribution from
   * m_semiImplicitConductivityCDR as it stands; if false, recompute it first through computeCdrConductivity.
   * @note Also adds in the conductive contribution from the CDR solvers.
   */
  virtual void
  computeCellConductivity(EBAMRCellData&                                             a_conductivityCell,
                          const Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles,
                          const bool a_useStoredCdrConductivity) noexcept;

  /**
   * @brief Compute the cell-centered conductivity
   * @note Interpolates m_conductivityCell to m_conductivityFace
   */
  virtual void
  computeFaceConductivity() noexcept;

  /**
   * @brief Set up the space charge density for the regrid operation
   */
  virtual void
  computeSemiImplicitRho() noexcept;

  /**
   * @brief Set up the semi-implicit Poisson solver
   *         @param[in] a_dt Description
   */
  virtual void
  setupSemiImplicitPoisson(const Real a_dt) noexcept;

  /**
   * @brief Solve the electrostatic problem
   * @details This calls the other version with all lhs/rhs taken from the field solver. Space charge density is
   * computed in this routine.
   * @return Returns true if the solver converged.
   */
  virtual bool
  solvePoisson() noexcept override;

  /**
   * @brief Remove covered particles
   * @param[in,out] a_particles      Particles to remove
   * @param[in]    a_representation Representation of the EB
   * @param[in]    a_tolerance      EB tolerance
   */
  virtual void
  removeCoveredPointParticles(Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles,
                              const EBRepresentation                               a_representation,
                              const Real                                           a_tolerance) const noexcept;

  /**
   * @brief Copy particles from the ItoSolver into PointParticles whose weight are ItoParticle::m_weight *
   * ItoParticle::m_mobility
   * @param[out] a_conductivityParticles Point-particles holding the mobility * weight
   */
  virtual void
  copyConductivityParticles(Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_conductivityParticles) noexcept;

  /**
   * @brief Advance the particles using the Euler-Maruyama scheme
   * @param[in] a_dt Time step
   */
  virtual void
  advanceEulerMaruyama(const Real a_dt) noexcept;

  /**
   * @brief Perform the explicit half of the Ito advance in the Euler-Maruyama step.
   * @details This computes each particle's explicit displacement -- the stochastic diffusion hop plus the lagged
   * a_dt*grad(D) drift correction -- and stores it on the ItoParticle scratch vector columns. It also fills the
   * "rho dagger" particles, which sit at the charged particles' positions displaced by that same amount. The
   * grad(D) term belongs here rather than in the particle velocities because it is independent of E^(k+1): the
   * rho dagger deposit accounts for it exactly, leaving the implicit half proportional to E^(k+1) as the
   * semi-implicit Poisson operator requires.
   * @param[in,out] a_rhoDaggerParticles Point particles at the displaced positions, used for the space charge.
   * @param[in]     a_dt                 Time step.
   */
  virtual void
  diffuseParticlesEulerMaruyama(Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_rhoDaggerParticles,
                                const Real                                           a_dt) noexcept;

  /**
   * @brief Compute the diffusion term for the CDR equations as well as the resulting CDR-contributions to
   * the space charge density.
   * @details This does an explicit diffusion advance over a_dt for the diffusive CDR species. This also
   *         @param[in,out] m_semiImplicitRhoCDR Description
   *         @param[in] a_dt Description
   * stores the resultant space charge on m_semiImplicitRhoCDR.
   */
  virtual void
  computeDiffusionTermCDR(EBAMRCellData& m_semiImplicitRhoCDR, const Real a_dt) noexcept;

  /**
   * @brief Step the particles according to the regular Euler-Maruyama scheme.
   * @param[in] a_dt Time step.
   * @details This sets X^(k+1) = X^k + a_dt*V + dX, where V = mu^k * E^(k+1) is the implicit half of the update
   * and dX is the explicit displacement computed by diffuseParticlesEulerMaruyama, i.e. the diffusion hop plus
   * the a_dt*grad(D) drift correction.
   * @note The explicit displacement MUST have been stored on the particle run-time storage
   */
  virtual void
  stepEulerMaruyamaParticles(const Real a_dt) noexcept;

  /**
   * @brief Step the CDR equations according to the regular Euler-Maruyama scheme.
   * @param[in] a_dt Time step.
   * @details This sets phi^(k+1) = phi^k - a_dt * div(v * phi) + a_dt * div(D*grad(phi))
   */
  virtual void
  stepEulerMaruyamaCDR(const Real a_dt) noexcept;

  /**
   * @brief Utility function for plotting the ItoSolver particles. These are written in a particles folder
   */
  virtual void
  plotParticles() const noexcept;
};
} // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCGodunovStepperImplem.H>

#endif
