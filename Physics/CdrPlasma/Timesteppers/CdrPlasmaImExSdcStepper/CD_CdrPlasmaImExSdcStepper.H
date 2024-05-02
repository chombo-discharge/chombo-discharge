/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaImExSdcStepper.H
  @brief  Declaration of a semi-implicit spectral deferred correction class for Physics/CdrPlasma
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaImExSdcStepper_H
#define CD_CdrPlasmaImExSdcStepper_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_CdrPlasmaStepper.H>
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Class for evolving plasma equations using implicit-explicit spectral deferred corrections.
    */
    class CdrPlasmaImExSdcStepper : public CdrPlasmaStepper
    {
    public:
      // Forward declare utility classes for handlnig the extra memory when advancing the solvers.
      class CdrStorage;
      class FieldStorage;
      class RtStorage;
      class SigmaStorage;

    public:
      /*!
	@brief Disallowed constructor -- use the full constructor
      */
      CdrPlasmaImExSdcStepper() = delete;

      /*!
	@brief Full constructor.
	@param[in] a_physics Reference to plasma physics object
      */
      CdrPlasmaImExSdcStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics);

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaImExSdcStepper();

      /*!
	@brief Advance method. Implements the ImEx SDC algorithm
	@param[in] a_dt Tried time step
	@return Time step used (can be different from a_dt if using adaptive stepping).
      */
      Real
      advance(const Real a_dt) override;

      /*!
	@brief Allocation method for setting up internal storage.
	@details This will set up the various storages for the SDC terms.
      */
      void
      allocateInternals() override;

      /*!
	@brief De-allocation method for deleting run-time storage. 
      */
      void
      deallocateInternals() override;

      /*!
	@brief For regridding internal storage. Does nothing. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void
      regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Parse startup options
      */
      void
      parseOptions() override;

      /*!
	@brief Parse run-time options.
      */
      void
      parseRuntimeOptions() override;

    protected:
      /*!
	@brief Scratch storage for the CDR solvers.
      */
      Vector<RefCountedPtr<CdrStorage>> m_cdrScratch;

      /*!
	@brief Scratch storage for the radiative transfer solvers.
      */
      Vector<RefCountedPtr<RtStorage>> m_rteScratch;

      /*!
	@brief Scratch storage for the field solver.
      */
      RefCountedPtr<FieldStorage> m_fieldScratch;

      /*!
	@brief Scratch storage for the surface charge solver.
      */
      RefCountedPtr<SigmaStorage> m_sigmaScratch;

      /*!
	@brief Maximum number of retries when using adaptive time steps.
      */
      int m_maxRetries;

      /*!
	@brief Number of subintervals in the SDC discretization.
      */
      int m_p;

      /*!
	@brief Number of correction for the SDC discretization.
      */
      int m_k;

      /*!
	@brief Error norm to use when using adaptive time steps
      */
      int m_errorNorm;

      /*!
	@brief Minimum number of SDC corrections.
      */
      int m_minCorr;

      /*!
	@brief The species error index -- allows us to restrict error evaluation and adaptive times stepping to a specified species.
	@details If m_errorIdx < 0 we do all species.
      */
      int m_errorIdx;

      /*!
	@brief Maximum growth factor for the time step when using adaptive stepping
      */
      Real m_maxDtGrowth;
      Real m_minCFL;
      Real m_maxCFL;
      Real m_errThresh;
      Real m_safety;
      Real m_decreaseSafety;
      Real m_newDt; // dt for next time step
      Real m_extrapDt;

      bool m_extrapAdvect;
      bool m_adaptiveDt;
      bool m_printReport;
      bool m_haveError;
      bool m_useTGA; // Second order diffusion

      // Basically only debugging stuff
      bool m_doAdvectionSource;
      bool m_doDiffusion;
      bool m_doPoisson;
      bool m_doRTE;
      bool m_computeD;
      bool m_computeV;
      bool m_computeS;
      bool m_consistentE;
      bool m_consistentRTE;
      bool m_haveDtErr;
      bool m_profileSteps;

      // Error things
      Vector<Real> m_cdrError;
      Real         m_sigmaError;
      Real         m_maxError;
      Real         m_preError;

      // Stuff for Lobatto intervals, Vandermonde matrix, integration
      // weights, Lobatto nodes etc.
      Vector<Vector<Real>> m_qmj;
      Vector<Vector<Real>> m_vandermonde;
      Vector<Real>         m_nodes;
      Vector<Real>         m_tm;
      Vector<Real>         m_dtm;

      std::string m_whichNodes;

      // Get functions for storage
      RefCountedPtr<CdrStorage>&
      getCdrStorage(const CdrIterator<CdrSolver>& a_solverit);
      RefCountedPtr<RtStorage>&
      getRtStorage(const RtIterator<RtSolver>& a_solverit);

      // Allocation functions
      void
      allocateCdrStorage();
      void
      allocateFieldStorage();
      void
      allocateRtStorage();
      void
      allocateSigmaStorage();

      // Get functions for substep data.
      Vector<EBAMRCellData*>
      getCdrSolversPhiK(const int a_m);
      EBAMRIVData&
      getSigmaSolverK(const int a_m);

      // Routines for intervals and quadrature
      void
      setupQuadratureNodes(const int a_p);
      void
      setupUniformNodes(const int a_p);
      void
      setupLobattoNodes(const int a_p);
      void
      setupChebyshevNodes(const int a_p);
      void
      setupQmj(const int a_p);
      void
      setupSubintervals(const Real a_time, const Real a_dt);
      void
      quad(EBAMRCellData& a_quad, const Vector<EBAMRCellData>& a_integrand, const int a_m);
      void
      quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m);

      // Routines for copying data at tm[0] and tm[p]
      void
      copyCdrToPhiM0();
      void
      copySigmaToM0();

      void
      copyPhiPToCdr();
      void
      copySigmaPToSigma();

      void
      integrateRtTransient(const Real a_dt);
      void
      integrateRtStationary();

      // These are the standard IMEX_SDC integration routines
      void
      integrate(const Real a_dt, const Real a_time, const bool a_lagged_terms);
      void
      integrateAdvectionReaction(const Real a_dt, const int a_m, const bool a_lagged_terms);
      void
      integrateAdvection(const Real a_dt, const int a_m, const bool a_lagged_terms);
      void
      integrateDiffusion(const Real a_dt, const int a_m, const bool a_lagged_terms);

      /*!
	@brief Get errors for CDR equations
      */
      Vector<EBAMRCellData*>
      getCdrErrors();

      // Make integrands quadrature-ready
      void
      reconcileIntegrands();

      // IMEX_SDC predictor routines
      void
      computeFD0();

      // IMEX_SDC corrector routines
      void
      initializeErrors();
      void
      finalizeErrors();

      // Step size control routines
      void
      computeNewDt(bool& a_accept_step, const Real a_dt, const int a_num_corrections);
      void
      adaptiveReport(const Real a_first_dt,
                     const Real a_dt,
                     const Real a_new_dt,
                     const int  a_corr,
                     const int  a_rej,
                     const Real a_max_err);

      // Poisson, RTE, and diffusion coefficient updates
      void
      updateField(); // This uses the solver states
      void
      updateField(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma);
      void
      updateDiffusionCoefficients();

      // Store & restore functions
      void
      storeSolvers();
      void
      restoreSolvers();

      // Aux functions
      void
      computeReactionNetwork(const int a_m, const Real a_time, const Real a_dt);
      void
      computeElectricFieldIntoScratch();
      void
      computeCdrEbStates();
      void
      computeCdrEbStates(const Vector<EBAMRCellData*>& a_phis);
      void
      computeCdrDomainStates();
      void
      computeCdrDomainStates(const Vector<EBAMRCellData*>& a_phis);
      void
      computeCdrGradients();
      void
      computeCdrGradients(const Vector<EBAMRCellData*>& a_phis);
      void
      computeCdrFluxes(const Real a_time);
      void
      computeCdrFluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
      void
      computeCdrDomainFluxes(const Real a_time);
      void
      computeCdrDomainFluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
      void
      computeCdrVelo(const Real a_time);
      void
      computeCdrVelo(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
      Real
      computeDt() override;
      void
      computeSigmaFlux();

      //
      Real
      getMaxError();
      Real
      getMaxNodeDistance();

      // Profiling
      void
      writeStepProfile(const Real a_dt,
                       const Real a_error,
                       const int  a_substeps,
                       const int  a_corrections,
                       const int  a_rejections);

      // Parse optiosn that are not really used
      void
      parseNodes();
      void
      parseDiffusionCoupling();
      void
      parseAdaptiveOptions();
      void
      parseDebugOptions();
      void
      parseAdvectionOptions();
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
