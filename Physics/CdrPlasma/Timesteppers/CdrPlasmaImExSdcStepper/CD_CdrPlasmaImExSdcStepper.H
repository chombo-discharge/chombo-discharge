/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_CdrPlasmaImExSdcStepper.H
   @brief  Declaration of Physics::CdrPlasma::CdrPlasmaImExSdcStepper
   @author Robert Marskar
*/

#ifndef CD_CDRPLASMAIMEXSDCSTEPPER_H
#define CD_CDRPLASMAIMEXSDCSTEPPER_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_CdrPlasmaStepper.H>
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /**
       @brief CdrPlasmaStepper subclass that advances plasma equations using implicit-explicit spectral deferred corrections (IMEX-SDC).
       @details This class implements a high-order IMEX time integration in which advection and
       reactions are treated explicitly while diffusion is treated implicitly. The method uses
       p+1 quadrature nodes on each time step (Gauss-Lobatto, Chebyshev, or uniform) to construct
       a high-order quadrature rule. Starting from an initial first-order predictor, k correction
       sweeps improve the solution towards the p-th order collocation solution.

       Adaptive time stepping is supported via a PI step-size controller that compares the solution
       at adjacent correction levels. The error norm (L-infinity or L2), the species used for error
       control, and the safety factors are all runtime-configurable.

       Scratch storage for the SDC sub-step data is managed by inner classes CdrStorage,
       FieldStorage, RtStorage, and SigmaStorage.
    */
    class CdrPlasmaImExSdcStepper : public CdrPlasmaStepper
    {
    public:
      /**
         @brief Forward declaration of per-CDR-solver scratch storage.
      */
      class CdrStorage;

      /**
         @brief Forward declaration of field solver scratch storage.
      */
      class FieldStorage;

      /**
         @brief Forward declaration of per-RTE-solver scratch storage.
      */
      class RtStorage;

      /**
         @brief Forward declaration of surface charge solver scratch storage.
      */
      class SigmaStorage;

      /**
         @brief Disallowed default constructor.
      */
      CdrPlasmaImExSdcStepper() = delete;

      /**
         @brief Full constructor.
         @param[in] a_physics Plasma kinetics interface.
      */
      CdrPlasmaImExSdcStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics);

      /**
         @brief Destructor.
      */
      virtual ~CdrPlasmaImExSdcStepper();

      /**
         @brief Advance equations by one time step using the IMEX-SDC algorithm.
         @param[in] a_dt Requested time step.
         @return Actual time step used. May differ from a_dt when adaptive stepping is enabled.
      */
      Real
      advance(const Real a_dt) override;

      /**
         @brief Allocate internal SDC scratch storage.
         @details Allocates phi^k, F_AR^k, F_D^k, and F^k arrays for all species at all SDC nodes,
         plus EB/domain state and flux arrays needed during the integration.
      */
      void
      allocateInternals() override;

      /**
         @brief Deallocate all internal SDC scratch storage.
      */
      void
      deallocateInternals() override;

      /**
         @brief Regrid hook for internal storage. No operations are performed.
         @param[in] a_lmin           Coarsest level that changed.
         @param[in] a_oldFinestLevel Finest level before the regrid.
         @param[in] a_newFinestLevel Finest level after the regrid.
      */
      void
      regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /**
         @brief Parse all startup options from ParmParse.
      */
      void
      parseOptions() override;

      /**
         @brief Parse runtime-adjustable options from ParmParse.
      */
      void
      parseRuntimeOptions() override;

    protected:
      /**
         @brief Per-CDR-solver SDC scratch storage.
      */
      Vector<RefCountedPtr<CdrStorage>> m_cdrScratch;

      /**
         @brief Per-RTE-solver SDC scratch storage.
      */
      Vector<RefCountedPtr<RtStorage>> m_rteScratch;

      /**
         @brief Field solver SDC scratch storage.
      */
      RefCountedPtr<FieldStorage> m_fieldScratch;

      /**
         @brief Surface charge solver SDC scratch storage.
      */
      RefCountedPtr<SigmaStorage> m_sigmaScratch;

      /**
         @brief Maximum number of step rejections before accepting regardless.
      */
      int m_maxRetries;

      /**
         @brief Number of SDC quadrature subintervals. There are p+1 nodes covering [t^n, t^{n+1}].
      */
      int m_p;

      /**
         @brief Number of SDC correction sweeps to perform per time step.
      */
      int m_k;

      /**
         @brief Error norm selector: 0 for L-infinity, 2 for L2.
      */
      int m_errorNorm;

      /**
         @brief Minimum number of corrections that must be completed before accepting a step.
      */
      int m_minCorr;

      /**
         @brief CDR species index used to drive the adaptive step-size controller.
         @details If negative, all species are used and the maximum error is taken.
      */
      int m_errorIdx;

      /**
         @brief Maximum growth factor applied to the time step by the PI controller.
      */
      Real m_maxDtGrowth;

      /**
         @brief Minimum CFL number allowed during adaptive stepping.
      */
      Real m_minCFL;

      /**
         @brief Maximum CFL number allowed during adaptive stepping.
      */
      Real m_maxCFL;

      /**
         @brief Error threshold below which a step is accepted.
      */
      Real m_errThresh;

      /**
         @brief Safety factor applied to the step-size controller's growth estimate.
      */
      Real m_safety;

      /**
         @brief Additional safety factor applied when decreasing the time step.
      */
      Real m_decreaseSafety;

      /**
         @brief Tentative time step for the next advance call, set by computeNewDt().
      */
      Real m_newDt;

      /**
         @brief Time step used for half-time extrapolation of advection states.
      */
      Real m_extrapDt;

      /**
         @brief If true, extrapolate CDR states to the half-time before computing advective fluxes.
      */
      bool m_extrapAdvect;

      /**
         @brief If true, use the PI controller to adaptively select the time step.
      */
      bool m_adaptiveDt;

      /**
         @brief If true, print a report of the adaptive stepping outcome after each step.
      */
      bool m_printReport;

      /**
         @brief True once an error estimate has been computed in the current step.
      */
      bool m_haveError;

      /**
         @brief If true, use the TGA (Crank-Nicholson/RK2) second-order implicit diffusion scheme.
      */
      bool m_useTGA;

      /**
         @brief Debug flag: if false, skip the advection and reaction integration.
      */
      bool m_doAdvectionSource;

      /**
         @brief Debug flag: if false, skip the implicit diffusion integration.
      */
      bool m_doDiffusion;

      /**
         @brief Debug flag: if false, skip the Poisson solve.
      */
      bool m_doPoisson;

      /**
         @brief Debug flag: if false, skip the RTE advance.
      */
      bool m_doRTE;

      /**
         @brief Debug flag: if false, skip diffusion coefficient updates.
      */
      bool m_computeD;

      /**
         @brief Debug flag: if false, skip drift velocity updates.
      */
      bool m_computeV;

      /**
         @brief Debug flag: if false, skip source term evaluation.
      */
      bool m_computeS;

      /**
         @brief If true, recompute the electric field consistently at each SDC node.
      */
      bool m_consistentE;

      /**
         @brief If true, recompute the RTE solution consistently at each SDC node.
      */
      bool m_consistentRTE;

      /**
         @brief True once a time step error estimate has been computed.
      */
      bool m_haveDtErr;

      /**
         @brief If true, write per-step SDC profile data to a file.
      */
      bool m_profileSteps;

      /**
         @brief Per-CDR-species error norm from the most recent correction comparison.
      */
      Vector<Real> m_cdrError;

      /**
         @brief Surface charge error norm from the most recent correction comparison.
      */
      Real m_sigmaError;

      /**
         @brief Maximum error across all species and surface charge.
      */
      Real m_maxError;

      /**
         @brief Error from the previous step, used by the PI controller.
      */
      Real m_preError;

      /**
         @brief SDC integration matrix Q_{mj}. Entry [m][j] is the weight of node j
         in the quadrature approximation of the integral from node 0 to node m.
      */
      Vector<Vector<Real>> m_qmj;

      /**
         @brief Vandermonde matrix used to compute the integration weights from the node positions.
      */
      Vector<Vector<Real>> m_vandermonde;

      /**
         @brief Quadrature node positions in [0, 1].
      */
      Vector<Real> m_nodes;

      /**
         @brief Absolute time at each SDC node: m_tm[m] = t^n + m_nodes[m] * dt.
      */
      Vector<Real> m_tm;

      /**
         @brief Sub-interval widths: m_dtm[m] = m_tm[m+1] - m_tm[m].
      */
      Vector<Real> m_dtm;

      /**
         @brief Node type selector: "uniform", "lobatto", or "chebyshev".
      */
      std::string m_whichNodes;

      /**
         @brief Return the SDC scratch storage associated with a given CDR solver.
         @param[in] a_solverIt CDR solver iterator pointing to the desired solver.
         @return Reference to the corresponding CdrStorage entry in m_cdrScratch.
      */
      RefCountedPtr<CdrStorage>&
      getCdrStorage(const CdrIterator<CdrSolver>& a_solverIt);

      /**
         @brief Return the SDC scratch storage associated with a given RTE solver.
         @param[in] a_solverIt RTE solver iterator pointing to the desired solver.
         @return Reference to the corresponding RtStorage entry in m_rteScratch.
      */
      RefCountedPtr<RtStorage>&
      getRtStorage(const RtIterator<RtSolver>& a_solverIt);

      /**
         @brief Allocate scratch storage for all CDR solvers.
      */
      static void
      allocateCdrStorage();

      /**
         @brief Allocate scratch storage for the field solver.
      */
      static void
      allocateFieldStorage();

      /**
         @brief Allocate scratch storage for all RTE solvers.
      */
      static void
      allocateRtStorage();

      /**
         @brief Allocate scratch storage for the surface charge solver.
      */
      static void
      allocateSigmaStorage();

      /**
         @brief Return the CDR phi^k values at SDC node m.
         @param[in] a_m SDC node index in [0, p].
         @return Vector of AMR data pointers, one per CDR species.
      */
      Vector<EBAMRCellData*>
      getCdrSolversPhiK(const int a_m);

      /**
         @brief Return the surface charge sigma^k at SDC node m.
         @param[in] a_m SDC node index in [0, p].
         @return Reference to the AMR EB data at node m.
      */
      EBAMRIVData&
      getSigmaSolverK(const int a_m);

      /**
         @brief Build quadrature nodes of type m_whichNodes with p subintervals.
         @param[in] a_p Number of subintervals (p+1 nodes).
      */
      void
      setupQuadratureNodes(const int a_p);

      /**
         @brief Build p+1 uniformly spaced nodes on [0, 1].
         @param[in] a_p Number of subintervals.
      */
      void
      setupUniformNodes(const int a_p);

      /**
         @brief Build p+1 Gauss-Lobatto nodes on [0, 1].
         @param[in] a_p Number of subintervals.
      */
      void
      setupLobattoNodes(const int a_p);

      /**
         @brief Build p+1 Chebyshev nodes on [0, 1].
         @param[in] a_p Number of subintervals.
      */
      void
      setupChebyshevNodes(const int a_p);

      /**
         @brief Compute the SDC integration matrix Q_{mj} from the current node positions.
         @param[in] a_p Number of subintervals.
      */
      void
      setupQmj(const int a_p);

      /**
         @brief Map the unit-interval node positions to [a_time, a_time + a_dt] and compute sub-interval widths.
         @param[in] a_time Start of the time step.
         @param[in] a_dt   Time step size.
      */
      void
      setupSubintervals(const Real a_time, const Real a_dt);

      /**
         @brief Compute the quadrature sum of a_integrand from node 0 to node a_m.
         @param[out] a_quad      Quadrature result (cell-centered).
         @param[in]  a_integrand Integrand at each SDC node.
         @param[in]  a_m         Upper integration limit (node index).
      */
      void
      quad(EBAMRCellData& a_quad, const Vector<EBAMRCellData>& a_integrand, const int a_m);

      /**
         @brief Compute the quadrature sum of a_integrand from node 0 to node a_m (EB-centered version).
         @param[out] a_quad      Quadrature result (EB-centered).
         @param[in]  a_integrand Integrand at each SDC node.
         @param[in]  a_m         Upper integration limit (node index).
      */
      void
      quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m);

      /**
         @brief Copy current CDR solver states into the SDC phi^k storage at node 0.
      */
      void
      copyCdrToPhiM0();

      /**
         @brief Copy current surface charge into the SDC sigma^k storage at node 0.
      */
      void
      copySigmaToM0();

      /**
         @brief Copy the SDC phi^k storage at the final node p back into the CDR solvers.
      */
      void
      copyPhiPToCdr();

      /**
         @brief Copy the SDC sigma^k storage at the final node p back into the surface charge solver.
      */
      void
      copySigmaPToSigma();

      /**
         @brief Advance RTE solvers transiently over a_dt.
         @param[in] a_dt Time step.
      */
      void
      integrateRtTransient(const Real a_dt);

      /**
         @brief Advance RTE solvers to the stationary solution.
      */
      void
      integrateRtStationary();

      /**
         @brief Perform one complete SDC sweep over all subintervals.
         @param[in] a_dt           Time step.
         @param[in] a_time         Start of the time step.
         @param[in] a_lagged_terms If true, use stored lagged operator values (correction sweep).
         If false, compute fresh values (predictor sweep).
      */
      void
      integrate(const Real a_dt, const Real a_time, const bool a_lagged_terms);

      /**
         @brief Advance CDR species from node m to m+1 using explicit advection+reaction and implicit diffusion.
         @param[in] a_dt           Full time step size (used to compute sub-interval widths from m_dtm).
         @param[in] a_m            Current SDC node index.
         @param[in] a_lagged_terms If true, include the lagged correction terms from the previous sweep.
      */
      void
      integrateAdvectionReaction(const Real a_dt, const int a_m, const bool a_lagged_terms);

      /**
         @brief Advance CDR species from node m to m+1 using explicit advection only (no diffusion).
         @param[in] a_dt           Full time step size.
         @param[in] a_m            Current SDC node index.
         @param[in] a_lagged_terms If true, include lagged correction terms.
      */
      void
      integrateAdvection(const Real a_dt, const int a_m, const bool a_lagged_terms);

      /**
         @brief Advance CDR species from node m to m+1 using implicit diffusion only.
         @param[in] a_dt           Full time step size.
         @param[in] a_m            Current SDC node index.
         @param[in] a_lagged_terms If true, include lagged correction terms.
      */
      void
      integrateDiffusion(const Real a_dt, const int a_m, const bool a_lagged_terms);

      /**
         @brief Return AMR pointers to the CDR error fields for all species.
         @return Vector of pointers to CdrStorage::m_error for each CDR species.
      */
      Vector<EBAMRCellData*>
      getCdrErrors();

      /**
         @brief Compute the combined operator F = F_AR + F_D at all SDC nodes.
         @details Called after each correction sweep to prepare the integrands for the next quadrature step.
      */
      void
      reconcileIntegrands();

      /**
         @brief Compute the diffusion operator F_D at the first SDC node (t^n) from current solver states.
      */
      void
      computeFD0();

      /**
         @brief Initialize error storage to zero before beginning correction sweeps.
      */
      void
      initializeErrors();

      /**
         @brief Accumulate the maximum correction norm into m_cdrError and m_sigmaError.
      */
      void
      finalizeErrors();

      /**
         @brief Evaluate the PI step-size controller and decide whether to accept or reject the step.
         @param[out] a_accept_step     Set to true if the step is accepted.
         @param[in]  a_dt              Time step used in this advance.
         @param[in]  a_num_corrections Number of corrections that were completed.
      */
      void
      computeNewDt(bool& a_accept_step, const Real a_dt, const int a_num_corrections);

      /**
         @brief Print a diagnostic summary of the adaptive stepping outcome.
         @param[in] a_first_dt    Time step proposed at the start of the step.
         @param[in] a_dt          Time step actually used.
         @param[in] a_new_dt      Proposed time step for the next advance.
         @param[in] a_corr        Number of corrections completed.
         @param[in] a_rej         Number of step rejections.
         @param[in] a_max_err     Maximum relative error across all species.
      */
      void
      adaptiveReport(const Real a_first_dt,
                     const Real a_dt,
                     const Real a_new_dt,
                     const int  a_corr,
                     const int  a_rej,
                     const Real a_max_err);

      /**
         @brief Solve the Poisson equation using the current CDR solver densities and surface charge.
      */
      void
      updateField();

      /**
         @brief Solve the Poisson equation from explicitly supplied densities and surface charge.
         @param[in] a_densities CDR species densities to use as space charge source.
         @param[in] a_sigma     Surface charge density.
      */
      void
      updateField(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma);

      /**
         @brief Recompute diffusion coefficients from the current electric field in scratch storage.
      */
      void
      updateDiffusionCoefficients();

      /**
         @brief Save the current CDR, RTE, field, and surface charge solver states to scratch.
      */
      void
      storeSolvers();

      /**
         @brief Restore CDR, RTE, field, and surface charge solver states from scratch.
      */
      void
      restoreSolvers();

      /**
         @brief Evaluate the reaction network source terms at SDC node m and store the result.
         @param[in] a_m    SDC node index.
         @param[in] a_time Current simulation time.
         @param[in] a_dt   Time step (passed to CdrPlasmaPhysics::advanceReactionNetwork).
      */
      void
      computeReactionNetwork(const int a_m, const Real a_time, const Real a_dt);

      /**
         @brief Compute the electric field from the current field solver state and store in m_fieldScratch.
      */
      void
      computeElectricFieldIntoScratch();

      /**
         @brief Extrapolate current CDR solver states to EB centroids.
      */
      void
      computeCdrEbStates();

      /**
         @brief Extrapolate explicitly supplied CDR states to EB centroids.
         @param[in] a_phis CDR cell-centered densities.
      */
      void
      computeCdrEbStates(const Vector<EBAMRCellData*>& a_phis);

      /**
         @brief Extrapolate current CDR solver states to domain faces.
      */
      void
      computeCdrDomainStates();

      /**
         @brief Extrapolate explicitly supplied CDR states to domain faces.
         @param[in] a_phis CDR cell-centered densities.
      */
      void
      computeCdrDomainStates(const Vector<EBAMRCellData*>& a_phis);

      /**
         @brief Compute cell-centered gradients of current CDR solver densities.
      */
      void
      computeCdrGradients();

      /**
         @brief Compute cell-centered gradients of explicitly supplied CDR densities.
         @param[in] a_phis CDR cell-centered densities.
      */
      void
      computeCdrGradients(const Vector<EBAMRCellData*>& a_phis);

      /**
         @brief Compute CDR boundary fluxes on EB faces from current solver states.
         @param[in] a_time Current simulation time.
      */
      void
      computeCdrFluxes(const Real a_time);

      /**
         @brief Compute CDR boundary fluxes on EB faces from explicitly supplied states.
         @param[in] a_phis CDR cell-centered densities.
         @param[in] a_time Current simulation time.
      */
      void
      computeCdrFluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time);

      /**
         @brief Compute CDR boundary fluxes on domain faces from current solver states.
         @param[in] a_time Current simulation time.
      */
      void
      computeCdrDomainFluxes(const Real a_time);

      /**
         @brief Compute CDR boundary fluxes on domain faces from explicitly supplied states.
         @param[in] a_phis CDR cell-centered densities.
         @param[in] a_time Current simulation time.
      */
      void
      computeCdrDomainFluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time);

      /**
         @brief Compute CDR drift velocities from the current field in scratch storage.
         @param[in] a_time Current simulation time.
      */
      void
      computeCdrVelo(const Real a_time);

      /**
         @brief Compute CDR drift velocities from explicitly supplied densities.
         @param[in] a_phis CDR cell-centered densities.
         @param[in] a_time Current simulation time.
      */
      void
      computeCdrVelo(const Vector<EBAMRCellData*>& a_phis, const Real a_time);

      /**
         @brief Compute the next time step based on CFL and adaptive-stepping constraints.
         @return Proposed time step.
      */
      Real
      computeDt() override;

      /**
         @brief Compute the surface charge boundary flux from current solver states.
      */
      void
      computeSigmaFlux();

      /**
         @brief Return the maximum CDR error norm across all species.
         @return Maximum of m_cdrError.
      */
      Real
      getMaxError();

      /**
         @brief Return the largest sub-interval width in the current SDC discretization.
         @return Maximum of m_dtm.
      */
      Real
      getMaxNodeDistance();

      /**
         @brief Write per-step SDC performance data to a profile file.
         @param[in] a_dt          Time step used.
         @param[in] a_error       Maximum relative error.
         @param[in] a_substeps    Number of sub-intervals (m_p).
         @param[in] a_corrections Number of correction sweeps completed.
         @param[in] a_rejections  Number of step rejections.
      */
      void
      writeStepProfile(const Real a_dt,
                       const Real a_error,
                       const int  a_substeps,
                       const int  a_corrections,
                       const int  a_rejections);

      /**
         @brief Parse the quadrature node type (m_whichNodes) and number of subintervals (m_p).
      */
      void
      parseNodes();

      /**
         @brief Parse the implicit diffusion coupling method (TGA or Crank-Nicholson).
      */
      void
      parseDiffusionCoupling();

      /**
         @brief Parse adaptive time-stepping parameters (thresholds, safety factors, growth limits).
      */
      void
      parseAdaptiveOptions();

      /**
         @brief Parse debug flags that toggle individual physics operators.
      */
      void
      parseDebugOptions();

      /**
         @brief Parse advection-related options (half-time extrapolation).
      */
      void
      parseAdvectionOptions();
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
