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
    */
    class BrownianWalkerStepper : public TimeStepper {
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
      ~BrownianWalkerStepper();

      /*!
	@brief Fill solver with initial data.
      */
      void initialData() override;

      /*!
	@brief Allocate storage for solvers and TimeStepper.
      */
      void allocate() override;

      /*!
	@brief Perform post-initialization routines. This does nothing. 
      */
      void postInitialize() override;

      /*!
	@brief Check if the input realm should be load balanced.
	@param[in] a_realm Realm to be load balanced.
	@details This returns true if a_realm==m_realm and m_loadBalance==true.
      */
      bool loadBalanceThisRealm(const std::string a_realm) const override;

      /*!
	@brief Load balance grid boxes for a specific realm. 
	@param[out] a_procs       MPI ranks owning the various grid boxes. 
	@param[out] a_boxes       Grid boxes on every level (obtain them with a_grids[lvl].boxArray())
	@param[in]  a_realm       Realm identifier
	@param[in]  a_grids       Original grids
	@param[in]  a_lmin        Coarsest grid level that changed
	@param[in]  a_finestLevel New finest grid level
	@note This should only be called if a_realm==m_realm and m_loadBalance==true
      */      
      void loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
			    Vector<Vector<Box> >&            a_boxes,
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
      void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif
      
#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */      
      void readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /*!
	@brief Perform post-initialization routines.
	@details This remaps the particles and makes new superparticles. It also sets the velocity fields. 
      */      
      void postCheckpointSetup() override;

      /*!
	@brief Get the number of plot variables for this time stepper. 
	@details Only the solver plots to file so this just returns the number of plot variables specified for the solver
      */
      int getNumberOfPlotVariables() const override;

      /*!
	@brief Write plot data to file. Only the solver writes plot data so this just calls the solver methods. 
	@param[inout] a_output            Output data holder.
	@param[in]    a_plotVariableNames Plot variable names.
	@param[inout] a_icomp             Starting component in a_output to begin at. 
      */      
      void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;

      /*!
	@brief Compute a time step to be used by Driver. This returns m_maxCellsHop*ItoSolver::computeDt (i.e. m_maxCellsHop is a CFL-like number)
	@param[out] a_dt       Time step to be used.
	@param[out] a_timeCode Identifier for how the time step was computed. 
      */
      void computeDt(Real& a_dt, TimeCode& a_timeCode) override;

      /*!
	@brief Synchronize solver times
	@param[in] a_step Time step
	@param[in] a_time Time (in seconds)
	@param[in] a_dt   Time step that was used. 
      */
      void synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;
      /*!
	@brief Print a step report -- this prints the local and global number of particles in the simulation. 
	@details This is called by Driver after time step. 
      */      
      void printStepReport() override;
      
      /*!
	@brief Register operators -- we only need the ones the solver needs to have the solver do this.
      */
      void registerOperators() override;

      /*!
	@brief Register realms to be used in the simulation. 
	@details The primal realm is always register -- but you can also register a separate particle realm. 
      */      
      void registerRealms() override;

      /*!
	@brief Parse runtime options
      */
      void parseRuntimeOptions() override;

      /*!
	@brief Always returns false.
      */
      bool needToRegrid() override;

      /*!
	@brief Perform pre-regrid operations.
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      void preRegrid(const int a_lbase, const int a_oldFinestLevel) override;

      /*!
	@brief Set up the Ito solver. This instantiates the solver and parses solver options. 
      */
      void setupSolvers() override;

      /*!
	@brief Advancement method. 
	@return a_dt
      */
      Real advance(const Real a_dt) override;

      /*!
	@brief Regrid method. This just regrids the solver. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */      
      void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Post-regrid operations. This does nothing (because we set velocities and diffusion coefficients in the regrid routine.
      */
      void postRegrid() override;

    protected:

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
      RefCountedPtr<ItoSolver>  m_solver;

      /*!
	@brief Reference to underlying species
      */      
      RefCountedPtr<ItoSpecies> m_species;

      /*!
	@brief Fluid diffusion coefficient
      */
      Real m_diffCo;

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
      Real m_maxCellsHop;

      /*!
	@brief Set the Brownian walker advection field
      */
      void setVelocity();
      
      /*!
	@brief Set the Brownian walker advection field (level version)
	@param[in] a_level Grid level
      */      
      void setVelocity(const int a_level);
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
