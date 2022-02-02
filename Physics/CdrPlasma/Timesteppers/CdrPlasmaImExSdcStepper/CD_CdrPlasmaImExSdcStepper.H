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
      @brief Class for evolving plasma equations with a multirate method using semi-implicit spectral deferred corrections
    */
    class CdrPlasmaImExSdcStepper : public CdrPlasmaStepper {
    public:
  
      // Forward declare nested classes. I don't want to clutter this file.
      class CdrStorage;
      class FieldStorage;
      class RtStorage;
      class SigmaStorage;

  
    public:

      CdrPlasmaImExSdcStepper();
      CdrPlasmaImExSdcStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics);
      ~CdrPlasmaImExSdcStepper();
  
      Real advance(const Real a_dt);

      void init();
      void allocateInternals() override;
      void deallocateInternals() override;
      void regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;
      void parseOptions() override;
      void parseRuntimeOptions() override;
      bool needToRegrid() override;

      Vector<EBAMRCellData*> getCdrErrors();
  
    protected:

      Vector<RefCountedPtr<CdrStorage> > m_cdr_scratch;
      Vector<RefCountedPtr<RtStorage> > m_rte_scratch;
      RefCountedPtr<FieldStorage>      m_fieldSolver_scratch;
      RefCountedPtr<SigmaStorage>        m_sigma_scratch;

      // Extra storage outside of solver storage
      EBAMRCellData scratch1;
      EBAMRCellData scratchD;

      int m_max_retries;
      int m_p;
      int m_k;
      int m_error_norm;
      int m_min_corr;
      int m_error_idx;

      Real m_accum_cfl;
      Real m_regrid_cfl;
      Real m_max_growth;
      Real m_minCFL;
      Real m_maxCFL;
      Real m_cycleCFL;
      Real m_err_thresh;
      Real m_safety;
      Real m_decrease_safe;
      Real m_new_dt; // dt for next time step
      Real m_extrap_dt;

      bool m_extrap_advect;
      bool m_adaptive_dt;
      bool m_print_report;
      bool m_have_err; 
      bool m_useTGA; // Second order diffusion
  

      // Basically only debugging stuff
      bool m_do_advec_src;  
      bool m_do_diffusion;  
      bool m_do_poisson;
      bool m_do_rte;
      bool m_compute_D;
      bool m_compute_v;
      bool m_compute_S;
      bool m_consistent_E;
      bool m_consistent_rte;
      bool m_have_dt_err;
      bool m_profile_steps;

      // Error things
      Vector<Real> m_cdr_error;
      Real m_sigma_error;
      Real m_max_error;
      Real m_pre_error;

      // Stuff for Lobatto intervals, Vandermonde matrix, integration
      // weights, Lobatto nodes etc. 
      Vector<Vector<Real> > m_qmj;
      Vector<Vector<Real> > m_vandermonde;
      Vector<Real> m_nodes;
      Vector<Real> m_tm;
      Vector<Real> m_dtm;

      std::string m_which_nodes;

      // Get functions for storage
      RefCountedPtr<CdrStorage>& getCdrStorage(const CdrIterator<CdrSolver>& a_solverit);
      RefCountedPtr<RtStorage>& getRtStorage(const RtIterator<RtSolver>& a_solverit);

      // Allocation functions
      void allocateCdrStorage();
      void allocateFieldStorage();
      void allocateRtStorage();
      void allocateSigmaStorage();

      // Get functions for substep data. 
      Vector<EBAMRCellData*> getCdrSolversPhiK(const int a_m);
      EBAMRIVData& getSigmaSolverK(const int a_m);

      // Routines for intervals and quadrature
      void setupQuadratureNodes(const int a_p);
      void setupUniformNodes(const int a_p);
      void setupLobattoNodes(const int a_p);
      void setupChebyshevNodes(const int a_p);
      void setupQmj(const int a_p);
      void setupSubintervals(const Real a_time, const Real a_dt);
      void quad(EBAMRCellData& a_quad, const Vector<EBAMRCellData>& a_integrand, const int a_m);
      void quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m);

      // Routines for copying data at tm[0] and tm[p]
      void copyCdrToPhiM0();
      void copySigmaToM0();

      void copyPhiPToCdr();
      void copySigmaPToSigma();

      void integrateRtTransient(const Real a_dt);
      void integrateRtStationary();

      // These are the standard IMEX_SDC integration routines
      void integrate(const Real a_dt, const Real a_time, const bool a_lagged_terms);
      void integrateAdvectionReaction(const Real a_dt, const int a_m, const bool a_lagged_terms);
      void integrateAdvection(const Real a_dt, const int a_m, const bool a_lagged_terms);
      void integrateDiffusion(const Real a_dt, const int a_m, const bool a_lagged_terms);

      // Make integrands quadrature-ready
      void reconcileIntegrands();

      // IMEX_SDC predictor routines
      void computeFD0();

      // IMEX_SDC corrector routines
      void initializeErrors();
      void finalizeErrors();

      // Step size control routines
      void computeNewDt(bool& a_accept_step, const Real a_dt, const int a_num_corrections);
      void adaptiveReport(const Real a_first_dt, const Real a_dt, const Real a_new_dt, const int a_corr, const int a_rej, const Real a_max_err);

      // Poisson, RTE, and diffusion coefficient updates
      void updateField(); // This uses the solver states
      void updateField(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma);
      void updateDiffusionCoefficients();

      // Store & restore functions
      void storeSolvers();
      void restoreSolvers();

      // Aux functions
      void computeReactionNetwork(const int a_m, const Real a_time, const Real a_dt);
      void computeElectricFieldIntoScratch();
      void computeCdrEbStates();
      void computeCdrEbStates(const Vector<EBAMRCellData*>& a_phis);
      void computeCdrDomainStates();
      void computeCdrDomainStates(const Vector<EBAMRCellData*>& a_phis);
      void computeCdrGradients();
      void computeCdrGradients(const Vector<EBAMRCellData*>& a_phis);
      void computeCdrFluxes(const Real a_time);
      void computeCdrFluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
      void computeCdrDomainFluxes(const Real a_time);
      void computeCdrDomainFluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
      void computeCdrVelo(const Real a_time);
      void computeCdrVelo(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
      void computeDt(Real& a_dt, TimeCode& a_timeCode);
      void computeSigmaFlux();

      //
      Real getMaxError();
      Real getMaxNodeDistance();

      // Profiling
      void writeStepProfile(const Real a_dt,
			    const Real a_error,
			    const int  a_substeps,
			    const int  a_corrections,
			    const int  a_rejections);

      // Parse optiosn that are not really used
      void parseNodes();
      void parseDiffusionCoupling();
      void parseAdaptiveOptions();
      void parseDebugOptions();
      void parseAdvectionOptions();
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
