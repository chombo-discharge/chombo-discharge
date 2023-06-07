/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BrownianWalkerStepper.H
  @brief  TimeStepper class for only solving the Brownian_Walker equation (with surface charge)
  @author Robert Marskar
*/

#ifndef CD_BrownianWalkerStepper_H
#define CD_BrownianWalkerStepper_H

// Our includes
#include "CD_TimeStepper.H"
#include <CD_ItoSolver.H>
#include <CD_ItoLayout.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace BrownianWalker {

    /*!
      @brief Implementation of TimeStepper which advances an ItoSolver. 
      @details This class is fairly simple -- we want to advance particles using an Ito kernel dX = V*dt + sqrt(2*D*dt)*N0 where N0 is a Gaussian random number. The particle
      velocity in the ItoSolver is V = mu*v(Xp) where mu is a particle mobility and v is a "velocity field". The ItoSolver has various way of representing the velocity, but
      here we consider a constant mobility. The diffusion coefficient is also constant (the value is fetched from the input script). 

      Note that while the ItoSolver CAN interpolate the particle diffusion coefficient from a mesh quantity -- this physics module sets it directly. 

      This physics module shows how to load balance an application introspectively based on the number of particles assigned to each rank. We support two methods,
      either by estimating the number of particles on the new mesh, or counting them directly. 
    */
    class BrownianWalkerStepper : public TimeStepper
    {
    public:
      /*!
	@brief Default constructor -- parses basic options
      */
      BrownianWalkerStepper();

      /*!
	@brief Full constructor. Calls the basic constructor and sets solver.
      */
      BrownianWalkerStepper(RefCountedPtr<ItoSolver>& a_solver);

      /*!
	@brief Destructor (does nothing)
      */
      virtual ~BrownianWalkerStepper();

      /*!
	@brief Fill solver with initial data.
      */
      void
      initialData() override;

      /*!
	@brief Allocate storage for solvers and TimeStepper.
      */
      void
      allocate() override;

      /*!
	@brief Perform post-initialization routines. This does nothing. 
      */
      void
      postInitialize() override;

      /*!
	@brief Check if the input realm should be load balanced.
	@param[in] a_realm Realm to be load balanced.
	@details This returns true if a_realm==m_realm and m_loadBalance==true.
      */
      bool
      loadBalanceThisRealm(const std::string a_realm) const override;

      /*!
	@brief Get computational loads to be checkpointed. This is used by Driver both for plotting computational loads AND for checkpointing the loads for restart
	purposes. 
	@details This uses the number of computational particles as a load metric. 
	@param[in] a_realm Realm
	@param[in] a_level Grid level
      */
      Vector<long int>
      getCheckpointLoads(const std::string a_realm, const int a_level) const override;

      /*!
	@brief Load balance grid boxes for a specific realm. 
	@details This will count the number of particles in each box and then load balance them. 
	@param[out] a_procs       MPI ranks owning the various grid boxes. 
	@param[out] a_boxes       Grid boxes on every level (obtain them with a_grids[lvl].boxArray())
	@param[in]  a_realm       Realm identifier
	@param[in]  a_grids       Original grids
	@param[in]  a_lmin        Coarsest grid level that changed
	@param[in]  a_finestLevel New finest grid level
	@note This should only be called if a_realm==m_realm and m_loadBalance==true
      */
      void
      loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                       Vector<Vector<Box>>&             a_boxes,
                       const std::string                a_realm,
                       const Vector<DisjointBoxLayout>& a_grids,
                       const int                        a_lmin,
                       const int                        a_finestLevel) override;

#ifdef CH_USE_HDF5
      /*!  
	@brief Write checkpoint data to file
	@param[inout] a_handle HDF5 fil
	@param[in]    a_lvl    Grid level
      */
      void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /*!
	@brief Perform post-initialization routines.
	@details This remaps the particles and makes new superparticles. It also sets the velocity fields. 
      */
      void
      postCheckpointSetup() override;

      /*!
	@brief Get the number of plot variables for this time stepper. 
	@details Only the solver plots to file so this just returns the number of plot variables specified for the solver
      */
      int
      getNumberOfPlotVariables() const override;

      /*!
	@brief Get plot variable names
      */
      Vector<std::string>
      getPlotVariableNames() const override;

      /*!
	@brief Write plot data to output holder. 
	@param[inout] a_output            Output data holder.
	@param[inout] a_icomp             Starting component in a_output to begin at. 
	@param[in]    a_level             Grid level
      */
      void
      writePlotData(LevelData<EBCellFAB>& a_output, int& a_icomp, const int a_level) const override;

      /*!
	@brief Compute a time step to be used by Driver. 
	@details This returns m_cfl*ItoSolver::computeDt (i.e. m_cfl is a CFL-like number)
      */
      Real
      computeDt() override;

      /*!
	@brief Synchronize solver times
	@param[in] a_step Time step
	@param[in] a_time Time (in seconds)
	@param[in] a_dt   Time step that was used. 
      */
      void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;
      /*!
	@brief Print a step report -- this prints the local and global number of particles in the simulation. 
	@details This is called by Driver after time step. 
      */
      void
      printStepReport() override;

      /*!
	@brief Register operators -- we only need the ones the solver needs to have the solver do this.
      */
      void
      registerOperators() override;

      /*!
	@brief Register realms to be used in the simulation. 
	@details The primal realm is always register -- but you can also register a separate particle realm. 
      */
      void
      registerRealms() override;

      /*!
	@brief Parse runtime options
      */
      void
      parseRuntimeOptions() override;

      /*!
	@brief Always returns false.
      */
      bool
      needToRegrid() override;

      /*!
	@brief Perform pre-regrid operations.
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      void
      preRegrid(const int a_lbase, const int a_oldFinestLevel) override;

      /*!
	@brief Set up the Ito solver. This instantiates the solver and parses solver options. 
      */
      void
      setupSolvers() override;

      /*!
	@brief Advancement method. 
	@return a_dt
      */
      Real
      advance(const Real a_dt) override;

      /*!
	@brief Regrid method. This just regrids the solver. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Post-regrid operations. 
      */
      void
      postRegrid() override;

    protected:
      /*!
	@brief For classifying load balancing methods.
      */
      enum LoadBalancingMethod
      {
        Mesh,
        Particle
      };

      /*!
	@brief Load balancing method
      */
      LoadBalancingMethod m_whichLoadBalance;

      /*!
	@brief Realm where the particles live
      */
      std::string m_realm;

      /*!
	@brief Phase where the particles live
      */
      phase::which_phase m_phase;

      /*!
	@brief If true, application will be load balanced
      */
      bool m_loadBalance;

      /*!
	@brief Reference to underlying solver
      */
      RefCountedPtr<ItoSolver> m_solver;

      /*!
	@brief Reference to underlying species
      */
      RefCountedPtr<ItoSpecies> m_species;

      /*!
	@brief Number of particles per cell. Used only during regrids
      */
      EBAMRCellData m_regridPPC;

      /*!
	@brief Fluid diffusion coefficient
      */
      Real m_diffCo;

      /*!
	@brief Particle mobility coefficients
      */
      Real m_mobility;

      /*!
	@brief Angular velocity for advection field
      */
      Real m_omega;

      /*!
	@brief Number of particles per cell
      */
      int m_ppc;

      /*! 
	@brief Effective CFL number for particle hops
      */
      Real m_cfl;

      /*!
	@brief Set advection and diffusion fields
      */
      void
      setAdvectionDiffusion();

      /*!
	@brief Set the Brownian walker advection field
      */
      void
      setVelocity();

      /*!
	@brief Set the ItoSolver "diffusion field"
      */
      void
      setDiffusion();

      /*!
	@brief Make super particles
      */
      void
      makeSuperParticles();

      /*!
	@brief Method that shows how to load balance the application using the number of particles per cell stored on the mesh. 
	@details This routine will load balance based on an estimated number of particles per cell. This number is computed ON THE MESH, i.e. it does
	not directly count the particles. Rather, we happen to know that ItoSolver deposited his particle numbers to the mesh before in the preRegrid, and we
	use that information to estimate the number of particles per patch in the new grids. 

	This is different from the other method which assign particles to the new grids and then counts them. Because this routine does not remap particles
	(which can be expensive), this routine is typically faster than loadBalanceBoxesParticles. 
	loadBalanceBoxesParticles 
	@param[out] a_procs       MPI ranks owning the various grid boxes. 
	@param[out] a_boxes       Grid boxes on every level (obtain them with a_grids[lvl].boxArray())
	@param[in]  a_realm       Realm identifier
	@param[in]  a_grids       Original grids
	@param[in]  a_lmin        Coarsest grid level that changed
	@param[in]  a_finestLevel New finest grid level
      */
      void
      loadBalanceBoxesMesh(Vector<Vector<int>>&             a_procs,
                           Vector<Vector<Box>>&             a_boxes,
                           const std::string                a_realm,
                           const Vector<DisjointBoxLayout>& a_grids,
                           const int                        a_lmin,
                           const int                        a_finestLevel);

      /*!
	@brief Method that shows how to load balance the application using the number of particles per cell stored on the mesh. 
	@details This routine will load balance by directly counting the particles on the new mesh. The result of that is that we regrid a particle data holder twice. 

	Because this routine ends up remapping particles twice (once here, and once in BrownianWalkerStepper::regrid), this routine is usually slower than
	loadBalanceBoxesMesh when used at large particle counts. In addition, loadBalanceBoxesMesh does it's best at estimating the final number of computatational particles
	after particle merging/splitting, while this routine only counts the original computational particles before the regrid. 
	@param[out] a_procs       MPI ranks owning the various grid boxes. 
	@param[out] a_boxes       Grid boxes on every level (obtain them with a_grids[lvl].boxArray())
	@param[in]  a_realm       Realm identifier
	@param[in]  a_grids       Original grids
	@param[in]  a_lmin        Coarsest grid level that changed
	@param[in]  a_finestLevel New finest grid level
      */
      void
      loadBalanceBoxesParticles(Vector<Vector<int>>&             a_procs,
                                Vector<Vector<Box>>&             a_boxes,
                                const std::string                a_realm,
                                const Vector<DisjointBoxLayout>& a_grids,
                                const int                        a_lmin,
                                const int                        a_finestLevel);
    };
  } // namespace BrownianWalker
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
