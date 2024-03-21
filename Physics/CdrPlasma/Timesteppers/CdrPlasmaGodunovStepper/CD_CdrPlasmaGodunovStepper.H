/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaGodunovStepper.H
  @brief  Declaration of the Godunov time stepper method
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaGodunovStepper_H
#define CD_CdrPlasmaGodunovStepper_H

// Std includes
#include <memory>

// Our includes
#include <CD_AmrMesh.H>
#include <CD_CdrPlasmaStepper.H>
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>
#include <CD_Timer.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Class for evolving plasma equations using a split step Godunov method. 
      @details This class can use semi-implicit formulations for the transport. 
    */
    class CdrPlasmaGodunovStepper : public CdrPlasmaStepper
    {
    public:
      // Forward declared utility classes for handling the extra memory when advancing the solvers.
      class CdrStorage;
      class FieldStorage;
      class RtStorage;
      class SigmaStorage;

      /*!
	@brief Disallowed constructor -- use strong construction
      */
      CdrPlasmaGodunovStepper() = delete;

      /*!
	@brief Full constructor
	@param[in] a_physics Reference to the plasma physics implementation.
      */
      CdrPlasmaGodunovStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics);

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaGodunovStepper();

      /*!
	@brief Parse startup options
      */
      void
      parseOptions() override;

      /*!
	@brief Parse run-time adjustable parameters
      */
      void
      parseRuntimeOptions() override;

      /*!
	@brief Implementation of the advance method
	@details This will switch between the various available implementations. 
	@param[in] a_dt Time step to try
	@return Returns a_dt (no adjustable time step for this class)
      */
      Real
      advance(const Real a_dt) override;

      /*!
	@brief Pre-regrid method. Used for storing important variables for the solvers
	@details This is overriden because the semi-implicit scheme requires us to store some extra variables on the old grids. 
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      void
      preRegrid(const int a_lbase, const int a_oldFinestLevel) override;

      /*!
	@brief Regrid method. 
	@details This overrides the parent class regrid method when using a semi-implicit scheme.
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Perform post-regrid operations. 
	@details This releases transient memory needed during the regrid. 
      */
      void
      postRegrid() override;

      /*!
	@brief Compute the time step -- this will be different for the different supported algorithms. 
      */
      Real
      computeDt() override;

      /*!
	@brief Run post-checkpoint setup operations. 
	@details The override is only relevant for the semi-implicit scheme because the field needs to be computed from a different equation. Fortunately,
	we've stored the necessary quantities (conductivity and space charge) to HDF5 so we can just set up a standard semi-implicit Poisson solve here. 
      */
      void
      postCheckpointSetup() override;

    protected:
      /*!
	@brief Enum for representing the coupling between the electric field and charge transport
      */
      enum class FieldCoupling
      {
        Explicit,
        SemiImplicit
      };

      /*!
	@brief Enum for switching between various advection solvers.
      */
      enum class AdvectionSolver
      {
        Euler,
        RK2,
        MUSCL
      };

      /*!
	@brief Enum for distinguishing various ways we handle diffusion. 
      */
      enum class DiffusionAlgorithm
      {
        Explicit,
        Implicit,
        Automatic,
      };

      /*!
	@brief For figuring out which transport algorithm we use. 
      */
      FieldCoupling m_fieldCoupling;

      /*!
	@brief For figuring out how we handle diffusion.
      */
      DiffusionAlgorithm m_diffusionAlgorithm;

      /*!
	@brief Advection solver
      */
      AdvectionSolver m_advectionSolver;

      /*!
	@brief Timer for run-time profiling
      */
      std::unique_ptr<Timer> m_timer;

      /*!
	@brief Scratch storage for the CDR solvers
      */
      Vector<RefCountedPtr<CdrStorage>> m_cdrScratch;

      /*!
	@brief Scratch storage for the RTE solvers
      */
      Vector<RefCountedPtr<RtStorage>> m_rteScratch;

      /*!
	@brief Scratch storage for the field solver
      */
      RefCountedPtr<FieldStorage> m_fieldScratch;

      /*!
	@brief Scratch storage for the surface charge solver
      */
      RefCountedPtr<SigmaStorage> m_sigmaScratch;

      /*!
	@brief Use for holding the "space charge" when doing semi-implicit solves. 
	@details This is different from the actual space charge because it also includes a diffusive contribution
      */
      EBAMRCellData m_semiImplicitRho;

      /*!
	@brief Used for storing conductivity*dt/eps0 on the cell center.
      */
      EBAMRCellData m_conductivityFactorCell;

      /*!
	@brief Used for storing conductivity*dt/eps0 on the face centers.
      */
      EBAMRFluxData m_conductivityFactorFace;

      /*!
	@brief Used for storing conductivity*dt/eps0 on the EB centers.
      */
      EBAMRIVData m_conductivityFactorEB;

      /*!
	@brief Scratch storage for storing m_semiImplicitRho on the old grids during regrid operations.
      */
      EBAMRCellData m_scratchSemiImplicitRho;

      /*!
	@brief Scratch storage for storing m_conductivityFactorCell on the old grids during regrid operations.
      */
      EBAMRCellData m_scratchConductivity;

      /*!
	@brief A special flag for letting the CDR solvers know that we want to compute advective derivatives at half time steps
      */
      bool m_extrapAdvect;

      /*!
	@brief Enable for debugging this class. 
      */
      bool m_debug;

      /*!
	@brief Enable performance profiling
      */
      bool m_profile;

      /*!
	@brief True if we floor the CDR densities. Due to reactions/EBs, they are not generally non-negative
      */
      bool m_floor;

      /*!
	@brief If true, add a stochastic diffusion flux. 
      */
      bool m_fhd;

      /*!
	@brief Regrid slopes or not
      */
      bool m_regridSlopes;

      /*!
	@brief Use compensation step in filter or not
      */
      bool m_filterCompensate;

      /*!
	@brief Diffusion threshold factor
      */
      Real m_implicitDiffusionThreshold;

      /*!
	@brief Integration order for diffusion equations
      */
      int m_diffusionOrder;

      /*!
	@brief Number of filterings for space charge
      */
      int m_numFilterRho;

      /*!
	@brief If true, we are using implicit diffusion. Otherwise it's explicit diffusion
      */
      std::vector<bool> m_useImplicitDiffusion;

      /*!
	@brief Function for getting the transient storage assocaited with a particular CDR solver
	@param[in] a_solverIt Solver iterator
	@return This returns the entry in m_cdrScratch corresponding to the solver. 
      */
      RefCountedPtr<CdrStorage>&
      getCdrStorage(const CdrIterator<CdrSolver>& a_solverIt);

      /*!
	@brief Function for getting the transient storage assocaited with a particular RTE solver
	@param[in] a_solverIt Solver iterator
	@return This returns the entry in m_rteScratch corresponding to the solver. 
      */
      RefCountedPtr<RtStorage>&
      getRtStorage(const RtIterator<RtSolver>& a_solverIt);

      /*!
	@brief For regridding internal storage. Does nothing. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void
      regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Allocate internal storage.
	@details This will set up extra memory needed in the semi-implicit coupling. 
      */
      void
      allocateInternals() override;

      /*!
	@brief Allocate internal storage.
	@details Set up some scratch storage needed in the time advancement algorithms.
      */
      void
      allocateScratch();

      /*!
	@brief Deallocate transient memory
      */
      void
      deallocateInternals() override;

      /*!
	@brief Allocate internal storage.
	@details Set up some scratch storage needed in the algorithms
      */
      void
      deallocateScratch();

      /*!
	@brief Compute electric field into scratch storage.
	@details This computes the electric field into the scratch storages in m_fieldScratch
      */
      void
      computeElectricFieldIntoScratch();

      /*!
	@brief Compute gradients of the CDR densities into scratch storage
      */
      void
      computeCdrGradients();

      /*!
	@brief Extrapolate cell-centered states with the source term
	@details This will compute cell-centered states (in the scratch storage) as n = n(t) + 0.5*a_dt * S(t)
	@param[in] a_dt Time step size. 
      */
      void
      extrapolateWithSourceTerm(const Real a_dt);

      /*!
	@brief Extrapolate the cell-centered states to the EB. 
      */
      void
      extrapolateCdrToEB();

      /*!
	@brief Compute the CDR boundary condition fluxes on the EBs.
	@details Prior to calling this routine, one must extrapolated the CDR states to the EB (calling extrapolateCdrToEB)
      */
      void
      computeCdrFluxesEB();

      /*!
	@brief Extrapolate the CDR cell-centered densities to domain edges/face
      */
      void
      extrapolateCdrToDomain();

      /*!
	@brief Compute the CDR boundary condition fluxes on the domain faces. 
	@details Prior to calling this routine, one must extrapolated the CDR states to the domain faces (calling extrapolateCdrToDomain)
      */
      void
      computeCdrDomainFluxes();

      /*!
	@brief Compute the surface flux for the surface charge solver. 
      */
      void
      computeSigmaFlux();

      /*!
	@brief Advance the transport problem.
	@param[in] a_dt Time step.
	@note This will switch between transport algorithms and call the other methods. 
      */
      void
      advanceTransport(const Real a_dt);

      /*!
	@brief Advance the transport problem using the explicit or explicit-implicit Euler rule
	@param[in] a_dt Time step.
      */
      void
      advanceTransportExplicitField(const Real a_dt);

      /*!
	@brief Advance the transport problem using a semi-implicit formulation for the electric field. 
	@details Under the hood, this is still an Euler solve.
	@param[in] a_dt Time step. 
      */
      void
      advanceTransportSemiImplicit(const Real a_dt);

      /*!
	@brief Advance the reactive problem over a time step dt.
	@param[in] a_dt Time step
	@note Only does the CDR solvers.
      */
      void
      advanceCdrReactions(const Real a_dt);

      /*!
	@brief Advance the radiative transfer problem. 
	@param[in] a_dt Time step
	@note This should be called after advanceReactions because that routine sets up the source terms
      */
      void
      advanceRadiativeTransfer(const Real a_dt);

      /*!
	@brief Print how much mass was injected into the system. 
	@details Occasionally, mass injection can be necessary due to negative densities arising the reactive problem, or due to redistribution near the
	embedded boundaries (which does not preserve non-negativeness). 
	@param[in] a_message Base message
	@param[in] a_solver  Which CDR solver to compute the mass for. 
      */
      void
      floorMass(EBAMRCellData& a_data, const std::string a_message, const RefCountedPtr<CdrSolver>& a_solver) const;

      /*!
	@brief Perform post-step operations
      */
      void
      postStep();

      /*!
	@brief Compute CDR drift velocities
	@details This is just like the parent method, except that we use the electric field stored in the scratch storage. 
	@param[in] a_time Time step
      */
      void
      computeCdrDriftVelocities(const Real a_time);

      /*!
	@brief Compute CDR diffusion coefficients
	@details This is just like the parent method, except that we use the electric field stored in the scratch storage. 
	@param[in] a_time Time step
      */
      void
      computeCdrDiffusionCoefficients(const Real a_time);

      /*!
	@brief Compute source terms for the CDR and RTE equations
	@param[in] a_dt Time step
      */
      void
      computeSourceTerms(const Real a_dt);

      /*!
	@brief Parse centering for the advective integrator.
      */
      void
      parseAdvection();

      /*!
	@brief Parse how we handle diffusion.
      */
      void
      parseDiffusion();

      /*!
	@brief Parse the transport algorithm
      */
      void
      parseField();

      /*!
	@brief Parse whether or not we should floor the mass
      */
      void
      parseFloor();

      /*!
	@brief Parse debug mode
      */
      void
      parseDebug();

      /*!
	@brief Parse profiling mode
      */
      void
      parseProfile();

      /*!
	@brief Parse inclusion of random diffusion flux.
      */
      void
      parseFHD();

      /*!
	@brief Parse whether or not to use slopes when regridding
      */
      void
      parseRegridSlopes();

      /*!
	@brief Parse filtering settings
      */
      void
      parseFiltering();

      /*!
	@brief Solve the semi-implicit Poisson equation
	@note This is just like the parent method for solving it, except that we use m_semiImplicitRho for the space charge density. Note that
	the user should first have set up the semi-implicit Poisson equation by computing the conductivities and calling CdrPlasmasetupSemiImplicitPoisson.
      */
      bool
      solveSemiImplicitPoisson();

#ifdef CH_USE_HDF5
      /*!  
	@brief Write checkpoint data to file.
	@details Overriden because semi-implicit scheme requires additioanl stuff in HDF5 files. 
	@param[inout] a_handle HDF5 fil
	@param[in]    a_lvl    Grid level
      */
      virtual void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file. Used for restarts. 
	@details Overriden because semi-implicit scheme requires additioanl stuff in HDF5 files. 
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
