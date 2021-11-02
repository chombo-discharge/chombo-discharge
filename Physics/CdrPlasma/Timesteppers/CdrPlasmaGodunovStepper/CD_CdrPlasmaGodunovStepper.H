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
    */
    class CdrPlasmaGodunovStepper : public CdrPlasmaStepper {
    public:
  
      // Forward declare nested classes. I don't want to clutter this file.
      class CdrStorage;
      class FieldStorage;
      class RtStorage;
      class SigmaStorage;

      CdrPlasmaGodunovStepper();
      CdrPlasmaGodunovStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics);
      ~CdrPlasmaGodunovStepper();

      // Pure functions that must be overwritten
      Real advance(const Real a_dt) override;
      void parseOptions() override;
      void parseRuntimeOptions() override;
      void init();
      void regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;
      void allocateInternals() override;
      void deallocateInternals() override;
      void computeDt(Real& a_dt, TimeCode& a_timeCode) override;
      bool needToRegrid() override;

    protected:

      enum class WhichTransportAlgorithm {
	SemiImplicit,
	Euler,
	RK2,
      };

      enum class WhichDiffusionAlgorithm{
	Explicit,
	Implicit,
	Automatic,
      };

      // Scratch storage
      Vector<RefCountedPtr<CdrStorage> > m_cdr_scratch;
      Vector<RefCountedPtr<RtStorage> >  m_rte_scratch;
      RefCountedPtr<FieldStorage>        m_fieldSolver_scratch;
      RefCountedPtr<SigmaStorage>        m_sigma_scratch;

      // Extra scratch storage when using semi-implicti scheme.
      EBAMRCellData m_conductivityCell;    // Used for holding the cell-centered conductivity.
      EBAMRFluxData m_conductivityFace;    // Used for holding the face-centered conductivity.
      EBAMRIVData   m_conductivityEB;      // Used for holding the EB-centered conductivity.

      EBAMRCellData m_scratchConductivity; // Scratch storage for holding the conductivity during regrids. 

      bool m_extrap_advect;
      bool m_debug;
      bool m_floor;
      bool m_fhd;
      bool m_implicit_diffusion;

      WhichDiffusionAlgorithm m_WhichDiffusionAlgorithm;
      WhichTransportAlgorithm m_WhichTransportAlgorithm;

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
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
