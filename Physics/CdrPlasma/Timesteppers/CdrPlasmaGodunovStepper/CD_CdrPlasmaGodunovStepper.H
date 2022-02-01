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

// Our includes
#include <CD_AmrMesh.H>
#include <CD_CdrPlasmaStepper.H>
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {
  
    /*!
      @brief Class for evolving plasma equations using a split step Godunov method. 
      @details This class can use semi-implicit formulations for the transport. 
    */
    class CdrPlasmaGodunovStepper : public CdrPlasmaStepper {
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
      void parseOptions() override;

      /*!
	@brief Parse run-time adjustable parameters
      */
      void parseRuntimeOptions() override;      

      /*!
	@brief Implementation of the advance method
	@details This will switch between the various available implementations. 
	@param[in] a_dt Time step to try
	@return Returns a_dt (no adjustable time step for this class)
      */
      Real advance(const Real a_dt) override;

      /*!
	@brief Pre-regrid method. Used for storing important variables for the solvers
	@details This is overriden because the semi-implicit scheme requires us to store some extra variables on the old grids. 
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      void preRegrid(const int a_lbase, const int a_oldFinestLevel) override;      

      /*!
	@brief Regrid method. 
	@details This overrides the parent class regrid method when using a semi-implicit scheme.
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Perform post-regrid operations. 
	@details This releases transient memory needed during the regrid. 
      */
      void postRegrid() override;

      /*!
	@brief For regridding internal storage. Does nothing. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Allocate internal storage.
	@details This will set up the memory for handling all the various advances. 
      */
      void allocateInternals() override;

      /*!
	@brief Deallocate transient memory
      */
      void deallocateInternals() override;

      /*!
	@brief Compute the time step -- this will be different for the different supported algorithms. 
	@param[out] a_dt       Time step to be used.
	@param[out] a_timeCode Identifier for how the time step was restricted.
      */
      void computeDt(Real& a_dt, TimeCode& a_timeCode) override;

      /*!
	@brief Run post-checkpoint setup operations. 
	@details The override is only relevant for the semi-implicit scheme because the field needs to be computed from a different equation. Fortunately,
	we've stored the necessary quantities (conductivity and space charge) to HDF5 so we can just set up a standard semi-implicit Poisson solve here. 
      */
      void postCheckpointSetup() override;

    protected:

      /*!
	@brief Enum for distinguishing various transport algorithms
      */
      enum class TransportAlgorithm {
	SemiImplicit,
	Euler,
	RK2,
      };

      /*!
	@brief Enum for distinguishing various ways we handle diffusion. 
      */      
      enum class WhichDiffusionAlgorithm{
	Explicit,
	Implicit,
	Automatic,
      };

      /*!
	@brief Scratch storage for the CDR solvers
      */
      Vector<RefCountedPtr<CdrStorage> > m_cdrScratch;

      /*!
	@brief Scratch storage for the RTE solvers
      */      
      Vector<RefCountedPtr<RtStorage> > m_rteScratch;

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


      // CODE-WALK BEGINS HERE      

      bool m_extrapAdvect;
      bool m_debug;
      bool m_floor;
      bool m_fhd;
      bool m_implicitDiffusion;

      WhichDiffusionAlgorithm m_WhichDiffusionAlgorithm;
      TransportAlgorithm m_transportAlgorithm;

      // Get functions for storage
      RefCountedPtr<CdrStorage>& getCdrStorage(const CdrIterator<CdrSolver>& a_solverit);
      RefCountedPtr<RtStorage>& getRtStorage(const RtIterator<RtSolver>& a_solverit);

      Real restrictDt();

      // Aux functions
      void computeElectricFieldIntoScratch();
      void computeCdrGradients();
      void computeCdrEbStates();
      void computeCdrEbFluxes();
      void computeCdrDomainStates();
      void computeCdrDomainFluxes();
      void computeSigmaFlux();
      void computeReactionNetwork(const Real a_dt);

      void advanceTransport(const Real a_dt);
      void advanceTransportEuler(const Real a_dt);
      void advanceTransportRK2(const Real a_dt);
      void advanceTransportSemiImplicit(const Real a_dt);
      void advanceRadiativeTransfer(const Real a_dt);
      void postStep();
      void extrapolateSourceTerm(const Real a_dt);
    
      void computeCdrDriftVelocities(const Real a_time);
      void computeCdrDiffusionCoefficients(const Real a_time);

      // Various option parsing
      void parseDiffusion();
      void parseAdvection();
      void parseTransport();
      void parseFloor();
      void parseDebug();
      void parseFHD();

      // Routine for semi-implicit Poisson solve.
      bool solveSemiImplicitPoisson();

      // Need to override these because the semi-implicit scheme dump the conductivity to disk. 
#ifdef CH_USE_HDF5
      virtual void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif
#ifdef CH_USE_HDF5      
      virtual void readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif      
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
