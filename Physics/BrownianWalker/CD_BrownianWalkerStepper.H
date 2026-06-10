/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_BrownianWalkerStepper.H
   @brief  Declaration of the Physics::BrownianWalker::BrownianWalkerStepper class
   @author Robert Marskar
*/

#ifndef CD_BROWNIANWALKERSTEPPER_H
#define CD_BROWNIANWALKERSTEPPER_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_ItoSolver.H>
#include <CD_ItoLayout.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace BrownianWalker {

    /**
       @brief TimeStepper implementation for advancing Brownian walker particles.
       @details This class drives a single ItoSolver using an Euler-Maruyama scheme:
       @f[ \mathbf{X}^{n+1} = \mathbf{X}^n + \mathbf{V}\,\Delta t + \sqrt{2D\,\Delta t}\,\mathbf{N}(0,1) @f]
       where @f$\mathbf{V}@f$ is the advection velocity, @f$D@f$ is the diffusion coefficient,
       and @f$\mathbf{N}(0,1)@f$ is a standard Gaussian random vector.

       The velocity field is a solid-body rotation with angular velocity @f$\omega@f$, and both the
       mobility and diffusion coefficient are spatially uniform and set directly on each particle.

       This module demonstrates two particle-based load balancing strategies:
       - Mesh-based: estimates the particle count on new grids from mesh-deposited data (faster).
       - Particle-based: directly remaps and counts particles on the new grids (more accurate).
    */
    class BrownianWalkerStepper : public TimeStepper
    {
    public:
      /**
         @brief Default constructor. Parses all options from the BrownianWalker ParmParse block.
      */
      BrownianWalkerStepper();

      /**
         @brief Full constructor. Delegates to the default constructor and stores the provided solver.
         @param[in] a_solver Ito solver to use for particle advancement.
      */
      BrownianWalkerStepper(RefCountedPtr<ItoSolver>& a_solver);

      /**
         @brief Destructor.
      */
      virtual ~BrownianWalkerStepper();

      /**
         @brief Fill the solver with initial particles, make superparticles, and set flow fields.
      */
      void
      initialData() override;

      /**
         @brief Allocate internal storage in the solver.
      */
      void
      allocate() override;

      /**
         @brief Post-initialization hook. No operations are performed.
      */
      void
      postInitialize() override;

      /**
         @brief Return whether the given realm should be load balanced by this stepper.
         @param[in] a_realm Realm identifier to query.
         @return True if a_realm matches m_realm and load balancing is enabled.
      */
      bool
      loadBalanceThisRealm(const std::string a_realm) const override;

      /**
         @brief Return per-patch computational loads for checkpointing and plotting.
         @details Uses the number of particles in each grid patch as the load metric.
         @param[in] a_realm Realm identifier.
         @param[in] a_level AMR level.
         @return Vector of particle counts, one entry per grid patch on a_level.
      */
      Vector<long int>
      getCheckpointLoads(const std::string a_realm, const int a_level) const override;

      /**
         @brief Load balance grid patches for the specified realm.
         @details Dispatches to either loadBalanceBoxesMesh or loadBalanceBoxesParticles
         depending on m_whichLoadBalance.
         @param[out] a_procs       MPI rank assigned to each grid patch.
         @param[out] a_boxes       Grid patches on each level.
         @param[in]  a_realm       Realm identifier.
         @param[in]  a_grids       Current grids before load balancing.
         @param[in]  a_lmin        Coarsest level that changed.
         @param[in]  a_finestLevel Finest level after regridding.
      */
      void
      loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                       Vector<Vector<Box>>&             a_boxes,
                       const std::string                a_realm,
                       const Vector<DisjointBoxLayout>& a_grids,
                       const int                        a_lmin,
                       const int                        a_finestLevel) override;

#ifdef CH_USE_HDF5
      /**
         @brief Write solver particle data to an HDF5 checkpoint file.
         @param[in,out] a_handle HDF5 file handle.
         @param[in]    a_lvl    AMR level to write.
      */
      void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /**
         @brief Read solver particle data from an HDF5 checkpoint file.
         @param[in,out] a_handle HDF5 file handle.
         @param[in]    a_lvl    AMR level to read.
      */
      void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /**
         @brief Post-checkpoint setup. Remaps particles, makes superparticles, and restores
         advection and diffusion fields after a restart.
      */
      void
      postCheckpointSetup() override;

      /**
         @brief Return the number of output variables contributed by this stepper.
         @return Number of plot variables (delegated to the solver).
      */
      int
      getNumberOfPlotVariables() const override;

      /**
         @brief Return the names of the output variables contributed by this stepper.
         @return Variable name list (delegated to the solver).
      */
      Vector<std::string>
      getPlotVariableNames() const override;

      /**
         @brief Write plot data into the output holder.
         @param[in,out] a_output      Output data holder.
         @param[in,out] a_icomp       Starting component index in a_output; incremented on return.
         @param[in]    a_outputRealm Realm that owns a_output.
         @param[in]    a_level       AMR level to write.
      */
      void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_outputRealm,
                    const int             a_level) const override;

      /**
         @brief Compute the next time step as m_cfl * ItoSolver::computeDt().
         @return Time step size in seconds.
      */
      Real
      computeDt() override;

      /**
         @brief Update internal time/step counters and forward them to the solver.
         @param[in] a_step Current time step number.
         @param[in] a_time Current simulation time in seconds.
         @param[in] a_dt   Time step that was used.
      */
      void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /**
         @brief Print a step report showing the local and global particle counts.
      */
      void
      printStepReport() override;

      /**
         @brief Register the operators required by the solver.
      */
      void
      registerOperators() override;

      /**
         @brief Register the realm(s) used in the simulation.
      */
      void
      registerRealms() override;

      /**
         @brief Parse runtime options from ParmParse. Re-reads verbosity, ppc, CFL, and load balancing flag.
      */
      void
      parseRuntimeOptions() override;

      /**
         @brief Always returns false; regridding is triggered externally by Driver.
         @return False.
      */
      bool
      needToRegrid() override;

      /**
         @brief Pre-regrid hook. Deposits particle counts onto the mesh and puts the solver in pre-regrid mode.
         @param[in] a_lbase          Coarsest level that will change.
         @param[in] a_oldFinestLevel Finest level before regridding.
      */
      void
      preRegrid(const int a_lbase, const int a_oldFinestLevel) override;

      /**
         @brief Instantiate the species and configure the solver with realm, phase, and AMR data.
      */
      void
      setupSolvers() override;

      /**
         @brief Advance the particle system by one Euler-Maruyama step.
         @details Performs: advection hop, diffusion hop, remap, EB absorption, superparticle
         merging, field update, and mesh deposition.
         @param[in] a_dt Time step size in seconds.
         @return The time step actually used (equal to a_dt).
      */
      Real
      advance(const Real a_dt) override;

      /**
         @brief Regrid the solver onto the new mesh and remake superparticles.
         @param[in] a_lmin           Coarsest level that changed.
         @param[in] a_oldFinestLevel Finest level before regridding.
         @param[in] a_newFinestLevel Finest level after regridding.
      */
      void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /**
         @brief Post-regrid hook. Clears temporary mesh storage and restores flow fields.
      */
      void
      postRegrid() override;

    protected:
      /**
         @brief Supported load balancing strategies.
      */
      enum LoadBalancingMethod
      {
        Mesh,    ///< Estimate particle counts from mesh-deposited data (fast).
        Particle ///< Count particles directly by remapping to the new grids (accurate).
      };

      /**
         @brief Selected load balancing strategy.
      */
      LoadBalancingMethod m_whichLoadBalance;

      /**
         @brief Name of the realm on which the particles live.
      */
      std::string m_realm;

      /**
         @brief Computational phase on which the particles live.
      */
      phase::which_phase m_phase;

      /**
         @brief If true, the application performs particle-based load balancing at each regrid.
      */
      bool m_loadBalance;

      /**
         @brief The Ito solver used to advance the particle system.
      */
      RefCountedPtr<ItoSolver> m_solver;

      /**
         @brief Species definition passed to the Ito solver.
      */
      RefCountedPtr<ItoSpecies> m_species;

      /**
         @brief Temporary storage for the number of particles per cell, used during load balancing.
      */
      EBAMRCellData m_regridPPC;

      /**
         @brief Spatially uniform diffusion coefficient for all particles.
      */
      Real m_diffCo;

      /**
         @brief Spatially uniform particle mobility.
      */
      Real m_mobility;

      /**
         @brief Angular velocity of the solid-body rotation advection field.
      */
      Real m_omega;

      /**
         @brief Target number of computational particles per cell after superparticle merging.
      */
      int m_ppc;

      /**
         @brief CFL-like scaling factor applied to the solver's intrinsic time step estimate.
      */
      Real m_cfl;

      /**
         @brief Set the velocity and diffusion fields on the mesh if the solver is mobile/diffusive.
      */
      void
      setAdvectionDiffusion();

      /**
         @brief Compute and assign the solid-body rotation velocity field to the solver.
      */
      void
      setVelocity();

      /**
         @brief Set the mesh diffusion function on the solver.
         @details Sets a sentinel large value since actual per-particle diffusion is assigned directly.
      */
      void
      setDiffusion();

      /**
         @brief Sort particles by cell, merge to m_ppc per cell, then re-sort by patch.
      */
      void
      makeSuperParticles();

      /**
         @brief Load balance using particle counts estimated from mesh-deposited data.
         @details Uses the particle-per-cell data deposited during preRegrid to estimate the
         post-regrid particle distribution without remapping particles. Faster than
         loadBalanceBoxesParticles but may be less accurate after superparticle merging.
         @param[out] a_procs       MPI rank assigned to each grid patch.
         @param[out] a_boxes       Grid patches on each level.
         @param[in]  a_realm       Realm identifier.
         @param[in]  a_grids       Current grids.
         @param[in]  a_lmin        Coarsest level that changed.
         @param[in]  a_finestLevel Finest level after regridding.
      */
      void
      loadBalanceBoxesMesh(Vector<Vector<int>>&             a_procs,
                           Vector<Vector<Box>>&             a_boxes,
                           const std::string&               a_realm,
                           const Vector<DisjointBoxLayout>& a_grids,
                           const int                        a_lmin,
                           const int                        a_finestLevel);

      /**
         @brief Load balance by directly remapping particles to the new grids and counting them.
         @details More accurate than loadBalanceBoxesMesh, but requires remapping the particle
         data twice (once here and once in regrid), making it slower at large particle counts.
         @param[out] a_procs       MPI rank assigned to each grid patch.
         @param[out] a_boxes       Grid patches on each level.
         @param[in]  a_realm       Realm identifier.
         @param[in]  a_grids       Current grids.
         @param[in]  a_lmin        Coarsest level that changed.
         @param[in]  a_finestLevel Finest level after regridding.
      */
      void
      loadBalanceBoxesParticles(Vector<Vector<int>>&             a_procs,
                                Vector<Vector<Box>>&             a_boxes,
                                const std::string&               a_realm,
                                const Vector<DisjointBoxLayout>& a_grids,
                                const int                        a_lmin,
                                const int                        a_finestLevel);
    };
  } // namespace BrownianWalker
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
