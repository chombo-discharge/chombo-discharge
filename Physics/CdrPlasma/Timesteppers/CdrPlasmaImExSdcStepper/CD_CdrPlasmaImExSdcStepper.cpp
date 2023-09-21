/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaImExSdcStepper.cpp
  @brief  Implementation of CD_CdrPlasmaImExSdcStepper.H
  @author Robert Marskar
*/

// Std includes
#include <fstream>
#include <iostream>
#include <iomanip>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CdrPlasmaImExSdcStepper.H>
#include <CD_CdrPlasmaImExSdcStorage.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_Timer.H>
#include <CD_CdrGodunov.H>
#include "CD_LaPackUtils.H"
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

typedef CdrPlasmaImExSdcStepper::CdrStorage   CdrStorage;
typedef CdrPlasmaImExSdcStepper::FieldStorage FieldStorage;
typedef CdrPlasmaImExSdcStepper::RtStorage    RtStorage;
typedef CdrPlasmaImExSdcStepper::SigmaStorage SigmaStorage;

CdrPlasmaImExSdcStepper::CdrPlasmaImExSdcStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics)
{
  CH_TIME("CdrPlasmaImExStepper::CdrPlasmaImExStepper()");

  m_className = "CdrPlasmaImExSdcStepper";
  m_physics   = a_physics;
}

CdrPlasmaImExSdcStepper::~CdrPlasmaImExSdcStepper()
{
  CH_TIME("CdrPlasmaImExSdcStepper::~CdrPlasmaImExSdcStepper()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::~CdrPlasmaImExSdcStepper()" << endl;
  }
}

void
CdrPlasmaImExSdcStepper::parseOptions()
{
  CH_TIME("CdrPlasmaImExSdcStepper::parseOptions()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::parseOptions()" << endl;
  }

  // Regular stuff from CdrPlasmaStepper that we almost always need
  this->parseVerbosity();
  this->parseSolverVerbosity();
  this->parseCFL();
  this->parseRelaxationTime();
  this->parseFastRadiativeTransfer();
  this->parseFastPoisson();
  this->parseMinDt();
  this->parseMaxDt();
  this->parseSourceComputation();

  // Specific to this class
  this->parseNodes();
  this->parseDiffusionCoupling();
  this->parseAdaptiveOptions();
  this->parseDebugOptions();
  this->parseAdvectionOptions();

  // Set up the quadrature nodes.
  this->setupQuadratureNodes(m_p);
  this->setupQmj(m_p);
}

void
CdrPlasmaImExSdcStepper::parseRuntimeOptions()
{
  CH_TIME("CdrPlasmaImExSdcStepper::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::parseRuntimeOptions()" << endl;
  }

  // Regular stuff from CdrPlasmaStepper that we almost always need
  this->parseVerbosity();
  this->parseSolverVerbosity();
  this->parseCFL();
  this->parseRelaxationTime();
  this->parseFastRadiativeTransfer();
  this->parseFastPoisson();
  this->parseMinDt();
  this->parseMaxDt();
  this->parseSourceComputation();

  // Parse nodes, diffusion couplnig and various other debugging flags.
  this->parseNodes();
  this->parseDiffusionCoupling();
  this->parseAdaptiveOptions();
  this->parseDebugOptions();
  this->parseAdvectionOptions();

  // Set up the new quadrature nodes.
  this->setupQuadratureNodes(m_p);
  this->setupQmj(m_p);

  // Solvers also parse run-time options.
  m_cdr->parseRuntimeOptions();
  m_rte->parseRuntimeOptions();
  m_fieldSolver->parseRuntimeOptions();
}

void
CdrPlasmaImExSdcStepper::parseNodes()
{
  CH_TIME("CdrPlasmaImExSdcStepper::parseNodes");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::parseNodes" << endl;
  }

  ParmParse   pp(m_className.c_str());
  std::string str;

  pp.get("quad_nodes", str);
  if (str == "lobatto") {
    m_whichNodes = "lobatto";
  }
  else if (str == "uniform") {
    m_whichNodes = "uniform";
  }
  else if (str == "chebyshev") {
    m_whichNodes = "chebyshev";
  }
  else {
    MayDay::Abort("CdrPlasmaImExSdcStepper::parseNodes - unknown node type requested");
  }

  pp.get("subintervals", m_p);
  pp.get("corr_iter", m_k);

  if (m_p < 1) {
    MayDay::Abort("CdrPlasmaImExSdcStepper::parseNodes - CdrPlasmaImExSdcStepper.subintervals cannot be < 1");
  }
  if (m_k < 0) {
    MayDay::Abort("CdrPlasmaImExSdcStepper::parseNodes - CdrPlasmaImExSdcStepper.corr_iter cannot be < 0");
  }
}

void
CdrPlasmaImExSdcStepper::parseDiffusionCoupling()
{
  CH_TIME("CdrPlasmaImExSdcStepper::parseDiffusionCoupling");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::parseDiffusionCoupling" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("use_tga", m_useTGA);
}

void
CdrPlasmaImExSdcStepper::parseAdaptiveOptions()
{
  CH_TIME("CdrPlasmaImExSdcStepper::parseAdaptiveOptions");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::parseAdaptiveOptions" << endl;
  }

  ParmParse   pp(m_className.c_str());
  std::string str;

  pp.get("error_norm", m_errorNorm);
  pp.get("min_corr", m_minCorr);
  pp.get("max_retries", m_maxRetries);
  pp.get("max_growth", m_maxDtGrowth);
  pp.get("decrease_safety", m_decreaseSafety);
  pp.get("min_cfl", m_minCFL);
  pp.get("max_cfl", m_maxCFL);
  pp.get("max_error", m_errThresh);
  pp.get("error_index", m_errorIdx);
  pp.get("safety", m_safety);
  pp.get("print_report", m_printReport);
  pp.get("adaptive_dt", m_adaptiveDt);

  m_minCorr = (!m_adaptiveDt) ? 0 : m_minCorr;
}

void
CdrPlasmaImExSdcStepper::parseDebugOptions()
{
  CH_TIME("CdrPlasmaImExSdcStepper::parseDebugOptions");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::parseDebugOptions" << endl;
  }

  ParmParse   pp(m_className.c_str());
  std::string str;

  pp.get("compute_D", str);
  m_computeD = (str == "true") ? true : false;
  pp.get("compute_v", str);
  m_computeV = (str == "true") ? true : false;
  pp.get("compute_S", str);
  m_computeS = (str == "true") ? true : false;

  pp.get("consistent_E", str);
  m_consistentE = (str == "true") ? true : false;
  pp.get("consistent_rte", str);
  m_consistentRTE = (str == "true") ? true : false;

  pp.get("do_advec_src", str);
  m_doAdvectionSource = (str == "true") ? true : false;
  pp.get("do_diffusion", str);
  m_doDiffusion = (str == "true") ? true : false;
  pp.get("do_rte", str);
  m_doRTE = (str == "true") ? true : false;
  pp.get("do_poisson", str);
  m_doPoisson = (str == "true") ? true : false;

  pp.get("profile_steps", str);
  m_profileSteps = (str == "true") ? true : false;
}

void
CdrPlasmaImExSdcStepper::parseAdvectionOptions()
{
  CH_TIME("CdrPlasmaImExSdcStepper::parseAdvectionOptions");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::parseAdvectionOptions" << endl;
  }

  ParmParse   pp(m_className.c_str());
  std::string str;

  m_extrapDt =
    0.5; // Relic of an ancient past. I don't see any reason why extrapolating to anything but the half interval
  // would make sense.

  pp.get("extrap_advect", str);
  m_extrapAdvect = (str == "true") ? true : false;
}

RefCountedPtr<CdrStorage>&
CdrPlasmaImExSdcStepper::getCdrStorage(const CdrIterator<CdrSolver>& a_solverit)
{
  return m_cdrScratch[a_solverit.index()];
}

RefCountedPtr<RtStorage>&
CdrPlasmaImExSdcStepper::getRtStorage(const RtIterator<RtSolver>& a_solverit)
{
  return m_rteScratch[a_solverit.index()];
}

Real
CdrPlasmaImExSdcStepper::getMaxNodeDistance()
{
  CH_TIME("CdrPlasmaImExSdcStepper::getMaxNodeDistance");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::getMaxNodeDistance" << endl;
  }

  Real max_dist = 0.0;
  for (int m = 0; m < m_p; m++) {
    max_dist = Max(max_dist, m_nodes[m + 1] - m_nodes[m]);
  }

  return max_dist;
}

void
CdrPlasmaImExSdcStepper::setupQuadratureNodes(const int a_p)
{
  CH_TIME("CdrPlasmaImExSdcStepper::setupQuadratureNodes");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::setupQuadratureNodes" << endl;
  }

  if (m_whichNodes == "uniform") {
    CdrPlasmaImExSdcStepper::setupUniformNodes(a_p);
  }
  else if (m_whichNodes == "lobatto") {
    CdrPlasmaImExSdcStepper::setupLobattoNodes(a_p);
  }
  else if (m_whichNodes == "chebyshev") {
    CdrPlasmaImExSdcStepper::setupChebyshevNodes(a_p);
  }
  else {
    MayDay::Abort("CdrPlasmaImExSdcStepper::setupQuadratureNodes - unknown nodes requested");
  }
}

void
CdrPlasmaImExSdcStepper::setupUniformNodes(const int a_p)
{
  CH_TIME("CdrPlasmaImExSdcStepper::setupUniformNodes");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::setupUniformNodes" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials.
  m_nodes.resize(1 + a_p);

  const Real delta = 2. / a_p;
  for (int m = 0; m <= a_p; m++) {
    m_nodes[m] = m * delta;
  }
}

void
CdrPlasmaImExSdcStepper::setupLobattoNodes(const int a_p)
{
  CH_TIME("CdrPlasmaImExSdcStepper::setupLobattoNodes");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::setupLobattoNodes" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials.
  m_nodes.resize(1 + a_p);

  if (a_p == 1) {
    m_nodes[0] = -1.0;
    m_nodes[1] = 1.0;
  }
  else if (a_p == 2) {
    m_nodes[0] = -1.0;
    m_nodes[1] = 0.0;
    m_nodes[2] = 1.0;
  }
  else if (a_p == 3) {
    m_nodes[0] = -1.0;
    m_nodes[1] = -1. / sqrt(5.);
    m_nodes[2] = 1. / sqrt(5.);
    m_nodes[3] = 1.0;
  }
  else if (a_p == 4) {
    m_nodes[0] = -1.0;
    m_nodes[1] = -sqrt(3. / 7);
    m_nodes[2] = 0.0;
    m_nodes[3] = sqrt(3. / 7);
    m_nodes[4] = 1.0;
  }
  else if (a_p == 5) {
    m_nodes[0] = -1.0;
    m_nodes[1] = -0.76595532;
    m_nodes[2] = -0.28532152;
    m_nodes[3] = 0.28532152;
    m_nodes[4] = 0.76595532;
    m_nodes[5] = 1.0;
  }
  else if (a_p == 6) {
    m_nodes[0] = -1.0;
    m_nodes[1] = -0.83022390;
    m_nodes[2] = -0.46884879;
    m_nodes[3] = 0.0;
    m_nodes[4] = 0.46884879;
    m_nodes[5] = 0.83022390;
    m_nodes[6] = 1.0;
  }
  else {
    MayDay::Abort(
      "CdrPlasmaImExSdcStepper::setupLobattoNodes - requested order exceeds 7. Compute your own damn nodes!");
  }
}

void
CdrPlasmaImExSdcStepper::setupChebyshevNodes(const int a_p)
{
  CH_TIME("CdrPlasmaImExSdcStepper::setupChebyshevNodes");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::setupChebyshevNodes" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials.
  m_nodes.resize(1 + a_p);
  m_nodes[0] = -1.0;
  for (int m = 1; m < a_p; m++) {
    m_nodes[m] = -cos((2 * m - 1) * Units::pi / (2 * (a_p - 1)));
  }
  m_nodes[a_p] = 1.0;
}

void
CdrPlasmaImExSdcStepper::setupQmj(const int a_p)
{
  CH_TIME("CdrPlasmaImExSdcStepper::setupQmj");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::setupQmj" << endl;
  }

  const int nnodes = 1 + a_p;

  // Resize the integration matrix, it should be p x (p+1)
  m_qmj.resize(a_p);
  for (int m = 0; m < a_p; m++) {
    m_qmj[m].resize(nnodes, 0.0);
  }

  // Generate the qmj matrix
  for (int j = 0; j < nnodes; j++) {

    // Set up the Vandermonde matrix (in Fortran order since we will call LaPack)
    double V[nnodes * nnodes];
    for (int j = 0; j < nnodes; j++) {
      for (int i = 0; i < nnodes; i++) {
        const int k = j * nnodes + i;
        V[k]        = pow(m_nodes[i], j);
      }
    }

    // Setup f = delta_kj. When we solve, this becomes the solution vector
    double cj[nnodes];
    for (int k = 0; k < nnodes; k++) {
      cj[k] = (k == j) ? 1.0 : 0.0;
    }

    // Solve V*c = f. This calls LAPACK
    int N    = nnodes;
    int NRHS = 1;
    int LDA  = nnodes;
    int IPIV[nnodes];
    int LDB  = nnodes;
    int INFO = 10;
    dgesv_(&N, &NRHS, V, &LDA, IPIV, cj, &LDB, &INFO);
    if (INFO != 0)
      MayDay::Abort("CdrPlasmaImExSdcStepper::setupQmj - could not compute weights");

    // Now construct qmj
    for (int m = 0; m < a_p; m++) {
      m_qmj[m][j] = 0.0;
      for (int k = 0; k < nnodes; k++) {
        m_qmj[m][j] += cj[k] * (pow(m_nodes[m + 1], k + 1) - pow(m_nodes[m], k + 1)) / (k + 1);
      }
    }
  }
}

void
CdrPlasmaImExSdcStepper::setupSubintervals(const Real a_time, const Real a_dt)
{
  CH_TIME("CdrPlasmaImExSdcStepper::setupSubintervals");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::setupSubintervals" << endl;
  }

  // m_nodes are Gauss-Lobatto nodes on [-1,1]. These must
  // be shifted to [t_n,t_n + a_dt]
  m_tm.resize(m_nodes.size());
  Vector<Real> shifted_nodes = m_nodes;
  for (int m = 0; m < shifted_nodes.size(); m++) {
    shifted_nodes[m] += 1.0;    // [0,2]
    shifted_nodes[m] *= 0.5;    // [0,1]
    shifted_nodes[m] *= a_dt;   // [0, a_dt]
    shifted_nodes[m] += a_time; // [a_time, a_time + a_dt]

    m_tm[m] = shifted_nodes[m];
  }

  // dtm = t_{m+1} - t_m. Order 1 is special since we only use the IMEX_SDC predictor from a second order formulation
  m_dtm.resize(m_tm.size() - 1);
  for (int m = 0; m < m_tm.size() - 1; m++) {
    m_dtm[m] = m_tm[m + 1] - m_tm[m];
  }
}

void
CdrPlasmaImExSdcStepper::quad(EBAMRCellData& a_quad, const Vector<EBAMRCellData>& a_integrand, const int a_m)
{
  CH_TIME("CdrPlasmaImExSdcStepper::quad");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::quad" << endl;
  }

  if (a_m < 0)
    MayDay::Abort("CdrPlasmaImExSdcStepper::quad - bad index a_m < 0");
  if (a_m >= m_p)
    MayDay::Abort("CdrPlasmaImExSdcStepper::quad - bad index a_m >= m_p");

  DataOps::setValue(a_quad, 0.0);
  for (int j = 0; j <= m_p; j++) {
    DataOps::incr(a_quad, a_integrand[j], m_qmj[a_m][j]);
  }
}

void
CdrPlasmaImExSdcStepper::quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m)
{
  CH_TIME("CdrPlasmaImExSdcStepper::quad");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::quad" << endl;
  }

  if (a_m < 0)
    MayDay::Abort("CdrPlasmaImExSdcStepper::quad - bad index a_m < 0");
  if (a_m >= m_p)
    MayDay::Abort("CdrPlasmaImExSdcStepper::quad - bad index a_m >= m_p");

  DataOps::setValue(a_quad, 0.0);
  for (int j = 0; j <= m_p; j++) {
    DataOps::incr(a_quad, a_integrand[j], m_qmj[a_m][j]);
  }
}

void
CdrPlasmaImExSdcStepper::copyPhiPToCdr()
{
  CH_TIME("CdrPlasmaImExSdcStepper::copyPhiPToCdr");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::copyPhiPToCdr" << endl;
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<CdrStorage>& storage = getCdrStorage(solver_it);

    EBAMRCellData&       phi  = solver->getPhi();
    const EBAMRCellData& phip = storage->getPhi()[m_p];
    DataOps::copy(phi, phip);
  }
}

void
CdrPlasmaImExSdcStepper::copySigmaPToSigma()
{
  CH_TIME("CdrPlasmaImExSdcStepper::copySigmaPToSigma");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::copySigmaPToSigma" << endl;
  }

  EBAMRIVData&       sigma  = m_sigma->getPhi();
  const EBAMRIVData& sigmap = m_sigmaScratch->getSigmaSolver()[m_p];
  DataOps::copy(sigma, sigmap);
}

Real
CdrPlasmaImExSdcStepper::advance(const Real a_dt)
{
  CH_TIME("CdrPlasmaImExSdcStepper::advance(Real)");
  if (m_verbosity > 2) {
    pout() << "CdrPlasmaImExSdcStepper::advance(Real)" << endl;
  }

  // ---------------------------------------------------------------------------------------------------
  // TLDR:  When we enter this routine, solvers SHOULD have been filled with valid ready and be ready
  //        advancement. If you think that this may not be the case, activate the debugging below
  // ---------------------------------------------------------------------------------------------------
  // Initialize integrations. If we do corrections, we need FD(phi_0) since this is implicit. If we do adaptive_dt, we should
  // also a
  CdrPlasmaImExSdcStepper::copyCdrToPhiM0();
  CdrPlasmaImExSdcStepper::copySigmaToM0();
  CdrPlasmaImExSdcStepper::computeFD0();
  CdrPlasmaImExSdcStepper::storeSolvers();

  // IMEX_SDC advance
  Real first_dt        = a_dt;
  Real actual_dt       = a_dt;
  int  num_reject      = 0;
  int  num_corrections = 0;
  bool accept_step     = false;
  bool retry_step      = true;

  m_maxError = 0.1234E5;
  Real t     = 0.0;
  while (!accept_step && retry_step) {
    num_corrections = 0;
    CdrPlasmaImExSdcStepper::setupSubintervals(m_time, actual_dt);

    // First SDC sweep. No lagged slopes here.
    CdrPlasmaImExSdcStepper::integrate(actual_dt, m_time, false);

    // SDC correction sweeps. Need to take care of lagged terms.
    for (int icorr = 0; icorr < Max(m_k, m_minCorr); icorr++) {
      num_corrections++;

      // Initialize error and reconcile integrands (i.e. make them quadrature-ready)
      CdrPlasmaImExSdcStepper::initializeErrors();
      CdrPlasmaImExSdcStepper::reconcileIntegrands();

      // SDC correction along whole interval
      CdrPlasmaImExSdcStepper::integrate(actual_dt, m_time, true);

      // Compute error and check if we need to keep iterating
      CdrPlasmaImExSdcStepper::finalizeErrors();
      if (m_maxError < m_errThresh && m_adaptiveDt && icorr >= m_minCorr)
        break; // No need in going beyond
    }

    // Compute a new time step. If it is smaller than the minimum allowed CFL step, accept the step anyways
    if (m_adaptiveDt) {
      CdrPlasmaImExSdcStepper::computeNewDt(accept_step, actual_dt, num_corrections);

      if (!accept_step) { // Step rejection, use the new dt for next step.
        actual_dt = m_newDt;
        num_reject++;

        retry_step = num_reject <= m_maxRetries;

        if (retry_step) {
          CdrPlasmaImExSdcStepper::restoreSolvers();
          CdrPlasmaImExSdcStepper::computeElectricFieldIntoScratch();
          CdrPlasmaImExSdcStepper::computeCdrGradients();
          CdrPlasmaImExSdcStepper::computeCdrVelo(m_time);
          CdrPlasmaStepper::computeCdrDiffusion(m_fieldScratch->getElectricFieldCell(),
                                                m_fieldScratch->getElectricFieldEb());
        }
      }
    }
    else {
      m_newDt     = 1.234567E89;
      accept_step = true;
    }
  }

  // Copy results back to solvers
  CdrPlasmaImExSdcStepper::copyPhiPToCdr();
  CdrPlasmaImExSdcStepper::copySigmaPToSigma();

  // Always recompute velocities and diffusion coefficients before the next time step. The Poisson and RTE equations
  // have been updated when we come in here.
  CdrPlasmaImExSdcStepper::computeCdrVelo(m_time + actual_dt);
  CdrPlasmaStepper::computeCdrDiffusion(m_fieldScratch->getElectricFieldCell(), m_fieldScratch->getElectricFieldEb());

  // Profile step
  if (m_printReport)
    CdrPlasmaImExSdcStepper::adaptiveReport(first_dt, actual_dt, m_newDt, num_corrections, num_reject, m_maxError);
  if (m_profileSteps)
    CdrPlasmaImExSdcStepper::writeStepProfile(actual_dt, m_maxError, m_p, num_corrections, num_reject);

  // Store current error.
  m_haveError = true;
  m_preError  = m_maxError;

  return actual_dt;
}

void
CdrPlasmaImExSdcStepper::copyCdrToPhiM0()
{
  CH_TIME("CdrPlasmaImExSdcStepper::copyCdrToPhiM0");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::copyCdrToPhiM0" << endl;
  }

  // CDR solvers
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<CdrStorage>& storage = getCdrStorage(solver_it);

    EBAMRCellData&       phi0 = storage->getPhi()[0];
    const EBAMRCellData& phi  = solver->getPhi();
    DataOps::copy(phi0, phi);
  }
}

void
CdrPlasmaImExSdcStepper::copySigmaToM0()
{
  CH_TIME("CdrPlasmaImExSdcStepper::copy_sigma_sigma_m0");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::copySigmaToM0" << endl;
  }

  // Copy sigma to starting state
  EBAMRIVData&       sigma0 = m_sigmaScratch->getSigmaSolver()[0];
  const EBAMRIVData& sigma  = m_sigma->getPhi();
  DataOps::copy(sigma0, sigma);
}

void
CdrPlasmaImExSdcStepper::computeFD0()
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeFD0");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeFD0" << endl;
  }

  if (m_k > 0) { // We only need this if we're actually doing any corrections....
    for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
      RefCountedPtr<CdrSolver>&  solver  = solver_it();
      RefCountedPtr<CdrStorage>& storage = getCdrStorage(solver_it);

      EBAMRCellData& phi_0 = storage->getPhi()[0]; // phi_0
      EBAMRCellData& FD_0  = storage->getFD()[0];  // FD(phi_0)

      if (solver->isDiffusive()) {
        solver->computeDivD(FD_0, phi_0, false, false, false); // Domain fluxes always come in through advection terms.

        // Shouldn't be necesary
        m_amr->conservativeAverage(FD_0, m_fluidRealm, m_cdr->getPhase());
        m_amr->interpGhost(FD_0, m_fluidRealm, m_cdr->getPhase());
      }
      else {
        DataOps::setValue(FD_0, 0.0);
      }
    }
  }
}

void
CdrPlasmaImExSdcStepper::integrate(const Real a_dt, const Real a_time, const bool a_lagged_terms)
{
  CH_TIME("CdrPlasmaImExSdcStepper::integrate");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::integrate" << endl;
  }

  // 1. The first time we enter this routine, velocities were updated.
  // 2. For further calls, source terms and velocities have been overwritten, but since the explicit
  //    operator slopes do not change, this is perfectly fine. We just increment with the lagged terms.
  Real t0, t1;
  Real total_time     = 0.0;
  Real setup_time     = 0.0;
  Real advect_time    = 0.0;
  Real diffusive_time = 0.0;

  // We begin with phi[0] = phi(t_n). Then update phi[m+1].
  Real time = a_time;
  for (int m = 0; m < m_p; m++) {

    // Update source terms every time we go through this
    CdrPlasmaImExSdcStepper::computeElectricFieldIntoScratch();
    CdrPlasmaImExSdcStepper::
      computeReactionNetwork(m,
                             a_time,
                             m_dtm[m]); // Ppdate the CDR and RTE source terms using the correct step size

    // Always update boundary conditions on the way in. All of these calls use the stuff that reside in the solvers,
    // which is what we need to do at the start of the time step. In principle, these things do not change
    // and so we could probably store them somewhere for increased performance.
    if (m == 0 &&
        !a_lagged_terms) { // This updates the CDR boundary conditions; since these are used to compute the slopes,
      t0 =
        Timer::wallClock(); // and the m=0 hyperbolic slopes do not change, we only need to do this for the predictor.
      CdrPlasmaImExSdcStepper::computeCdrEbStates();
      CdrPlasmaImExSdcStepper::computeCdrFluxes(a_time);
      CdrPlasmaImExSdcStepper::computeCdrDomainStates();
      CdrPlasmaImExSdcStepper::computeCdrDomainFluxes(a_time);
      CdrPlasmaImExSdcStepper::computeSigmaFlux();
      t1 = Timer::wallClock();

      total_time = -t0;
      setup_time = t1 - t0;
    }

    // This does the transient rte advance. Source terms were uåpdated in the computeReactionNetwork routine above.
    t0 = Timer::wallClock();
    if (!(m_rte->isStationary()))
      CdrPlasmaImExSdcStepper::integrateRtTransient(a_dt);

    // This computes phi_(m+1) = phi_m + dtm*FAR_m(phi_m) + lagged quadrature and lagged advection-reaction
    t0 = Timer::wallClock();
    CdrPlasmaImExSdcStepper::integrateAdvectionReaction(a_dt, m, a_lagged_terms);
    t1 = Timer::wallClock();
    advect_time += t1 - t0;

    // This does the diffusion advance. It also adds in the remaining lagged diffusion terms before the implicit diffusion solve
    t0 = Timer::wallClock();
    CdrPlasmaImExSdcStepper::integrateDiffusion(a_dt, m, a_lagged_terms);
    t1 = Timer::wallClock();
    diffusive_time += t1 - t0;

    // After the diffusion step we update the Poisson and *stationary* RTE equations
    Vector<EBAMRCellData*> cdr_densities_mp1 = CdrPlasmaImExSdcStepper::getCdrSolversPhiK(m + 1);
    EBAMRIVData&           sigma_mp1         = CdrPlasmaImExSdcStepper::getSigmaSolverK(m + 1);
    const Real             t_mp1             = m_tm[m + 1];

    // Update electric field and stationary RTE equations
    if (m_consistentE)
      CdrPlasmaImExSdcStepper::updateField(cdr_densities_mp1, sigma_mp1);
    if (m_consistentRTE) {
      if (m_rte->isStationary()) {
        CdrPlasmaImExSdcStepper::computeReactionNetwork(m + 1, time + m_dtm[m], m_dtm[m]);
        CdrPlasmaImExSdcStepper::integrateRtStationary();
      }
    }

    // If we need another step, we should update boundary conditions agains. We DONT do this on the last step
    // because this was also done on the way INTO this routine. If we've updated m=m_p, we either recompute
    // boundary conditions in the next SDC sweep, or we allow the next time step to take care of this.
    const int last = m == m_p - 1;
    if (!last) {
      if (m_computeS)
        CdrPlasmaImExSdcStepper::computeCdrGradients(cdr_densities_mp1);
      if (m_computeV)
        CdrPlasmaImExSdcStepper::computeCdrVelo(cdr_densities_mp1, t_mp1);
      if (m_computeD)
        CdrPlasmaImExSdcStepper::updateDiffusionCoefficients();

      // Update boundary conditions for cdr and sigma equations. We need them for the next step
      CdrPlasmaImExSdcStepper::computeCdrEbStates(cdr_densities_mp1);
      CdrPlasmaImExSdcStepper::computeCdrFluxes(cdr_densities_mp1, t_mp1);
      CdrPlasmaImExSdcStepper::computeCdrDomainStates(cdr_densities_mp1);
      CdrPlasmaImExSdcStepper::computeCdrDomainFluxes(cdr_densities_mp1, t_mp1);
      CdrPlasmaImExSdcStepper::computeSigmaFlux();
    }

    time += m_dtm[m];
  }
  t1 = Timer::wallClock();

  total_time += t1;

#if 0
  pout() << endl
	 << "setup time = " << setup_time << endl
	 << "advect_time = " << advect_time << endl
	 << "diffusive_time = " << diffusive_time << endl
	 << "total time = " << total_time << endl
	 << endl;
#endif
}

void
CdrPlasmaImExSdcStepper::integrateAdvectionReaction(const Real a_dt, const int a_m, const bool a_lagged_terms)
{
  CH_TIME("CdrPlasmaImExSdcStepper::integrateAdvectionReaction");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::integrateAdvectionReaction" << endl;
  }

  // Advance phi_(m+1) = phi_m + dtm*F_A. These routines do nothing
  // with the operator slopes for phi, but they do adjust the slopes m_Fsig (but not m_Fsum) for sigma. Incidentally,
  // if m=0 and a_lagged_terms=true, we can increment directly with the precomputed advection-reaction. This means that
  // we can skip the advective advance. The sigma advance is accordingly also skipped.
  const bool skip = (a_m == 0 && a_lagged_terms);
  const Real t0   = Timer::wallClock();
  if (!skip) {
    CdrPlasmaImExSdcStepper::integrateAdvection(a_dt, a_m, a_lagged_terms);
  }
  const Real t1 = Timer::wallClock();

  // Add in the reaction term and then compute the new operator slopes.
  // If this is the corrector and m=0, we skipped the advection advance because we can use the precomputed
  // advection-reaction operator slope. In this case phi_(m+1) is bogus and we need to recompute it. Otherwise,
  // phi_(m+1) = phi_m + dtm*FA_m, and we just increment with the reaction operator.
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<CdrStorage>& storage = getCdrStorage(solver_it);

    // phi_(m+1) = phi_M
    EBAMRCellData&       phi_m1  = storage->getPhi()[a_m + 1];
    EBAMRCellData&       scratch = storage->getScratch();
    const EBAMRCellData& phi_m   = storage->getPhi()[a_m];

    // Increment with operator slopes. m=0 and corrector is a special case where we skipped the advective advance,
    // choosing instead to use the old slopes (which did not change)
    if (skip) {                                            // Can use the old slopes
      const EBAMRCellData& FAR_m = storage->getFAR()[a_m]; // Slope, doesn't require recomputation.
      DataOps::copy(phi_m1, phi_m);
      DataOps::incr(phi_m1, FAR_m, m_dtm[a_m]);
      if (a_lagged_terms) {
        DataOps::copy(scratch, FAR_m);
      }
    }
    else { // If we made it here, phi_(m+1) = phi_m + dtm*FA(phi_m) through the integrateAdvection routine
      EBAMRCellData& FAR_m = storage->getFAR()[a_m]; // Currently the old slope
      EBAMRCellData& src   = solver->getSource();    // Updated source

      // Increment swith source and then compute slope. This has already been done
      DataOps::incr(phi_m1, src, m_dtm[a_m]); // phi_(m+1) = phi_m + dtm*(FA_m + FR_m)

      // This shouldn't be necessary
      m_amr->conservativeAverage(phi_m1, m_fluidRealm, m_cdr->getPhase());
      m_amr->interpGhost(phi_m1, m_fluidRealm, m_cdr->getPhase());

      if (a_lagged_terms) { // Back up the old slope first, we will need it for the lagged term
        DataOps::copy(scratch, FAR_m);
      }

      // Re-compute the advection-reaction slope for node t_m
      DataOps::copy(FAR_m, phi_m1);           // FAR_m = (phi_(m+1) - phi_m)/dtm
      DataOps::incr(FAR_m, phi_m, -1.0);      // :
      DataOps::scale(FAR_m, 1. / m_dtm[a_m]); // :

      // Shouldn't be necessary
      m_amr->conservativeAverage(FAR_m, m_fluidRealm, m_cdr->getPhase());
      m_amr->interpGhost(FAR_m, m_fluidRealm, m_cdr->getPhase());
    }

    // Now add in the lagged advection-reaction and quadrature terms. This is a bit weird, but we did overwrite
    // FAR_m above after the advection-reaction advance, but we also backed up the old term into scratch.
    if (a_lagged_terms) {
      DataOps::incr(phi_m1, scratch, -m_dtm[a_m]); // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k)
      CdrPlasmaImExSdcStepper::quad(scratch,
                                    storage->getF(),
                                    a_m); // Does the quadrature of the lagged operator slopes.
      DataOps::incr(phi_m1,
                    scratch,
                    0.5 * a_dt); // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k) + I_m^(m+1)
    }
  }
  const Real t2 = Timer::wallClock();

  // Add in the lagged terms for sigma. As above, m=0 and corrector is a special case where we just use the old slopes.
  EBAMRIVData&       sigma_m1 = m_sigmaScratch->getSigmaSolver()[a_m + 1];
  const EBAMRIVData& sigma_m  = m_sigmaScratch->getSigmaSolver()[a_m];
  if (skip) {
    const EBAMRIVData& Fsig_m = m_sigmaScratch->getFold()[a_m]; // Here, we should be able to use either Fold or Fnew
    DataOps::copy(sigma_m1, sigma_m);                           // since Fsig_0 is only computed once.
    DataOps::incr(sigma_m1, Fsig_m, m_dtm[a_m]);
  }

  const Real t3 = Timer::wallClock();
  if (a_lagged_terms) { // Add in the lagged terms. When we make it here, sigma_(m+1) = sigma_m + dtm*Fsig_m.
    EBAMRIVData& Fsig_lag = m_sigmaScratch->getFold()[a_m];
    DataOps::incr(sigma_m1, Fsig_lag, -m_dtm[a_m]);

    // Add in the quadrature term
    EBAMRIVData& scratch = m_sigmaScratch->getScratch();
    CdrPlasmaImExSdcStepper::quad(scratch, m_sigmaScratch->getFold(), a_m);
    DataOps::incr(sigma_m1, scratch, 0.5 * a_dt); // Mult by 0.5*a_dt due to scaling on [-1,1] for quadrature
  }
  const Real t4 = Timer::wallClock();

#if 0
  pout() << "integrateAdvectionReaction::" << endl;
  pout() << "t1-t0 = " << t1-t0 << endl;
  pout() << "t2-t1 = " << t2-t2 << endl;
  pout() << "t3-t2 = " << t3-t2 << endl;
  pout() << "t4-t3 = " << t4-t3 << endl;
#endif
}

void
CdrPlasmaImExSdcStepper::integrateAdvection(const Real a_dt, const int a_m, const bool a_lagged_terms)
{
  CH_TIME("CdrPlasmaImExSdcStepper::integrateAdvection");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::integrateAdvection" << endl;
  }

  // TLDR; This routine should do phi_(m+1) = phi_m + dtm*FA_m, and sigma_(m+1) = sigma_m + dt*Fsig_m.
  //       It also computes the sigma slope.
  //
  //       The lagged terms are not a part of this routine.

  if (a_m == 0 && a_lagged_terms) {
    MayDay::Abort(
      "CdrPlasmaImExSdcStepper::integrateAdvection - (m==0 && corrector==true) logic bust which should never happen");
  }

  // Advance cdr equations
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<CdrStorage>& storage = getCdrStorage(solver_it);

    EBAMRCellData& phi_m1  = storage->getPhi()[a_m + 1];
    EBAMRCellData& scratch = storage->getScratch();
    EBAMRCellData& phi_m   = storage->getPhi()[a_m];

    if (solver->isMobile()) {
      const Real extrap_dt = m_extrapAdvect ? 2.0 * m_extrapDt * m_dtm[a_m] : 0.0; // Factor of 2 due to EBPatchAdvect
      solver->computeDivF(scratch, phi_m, extrap_dt, false, true, true);           // scratch =  Div(v_m*phi_m^(k+1))

      DataOps::copy(phi_m1, phi_m);
      DataOps::incr(phi_m1, scratch, -m_dtm[a_m]);
      DataOps::floor(phi_m1, 0.0);
      m_amr->conservativeAverage(phi_m1, m_fluidRealm, m_cdr->getPhase());
      m_amr->interpGhost(phi_m1, m_fluidRealm, m_cdr->getPhase());
    }
    else {
      DataOps::copy(phi_m1, phi_m);
    }
  }

  // Update sigma. Also compute the new slope.
  EBAMRIVData&       sigma_m1 = m_sigmaScratch->getSigmaSolver()[a_m + 1];
  EBAMRIVData&       Fsig_new = m_sigmaScratch->getFnew()[a_m];
  const EBAMRIVData& sigma_m  = m_sigmaScratch->getSigmaSolver()[a_m];
#if 1 // New code after SigmaSolver -> SurfaceODESolver
  DataOps::copy(Fsig_new, m_sigma->getRHS());
#else
  m_sigma->computeRHS(Fsig_new); // Fills Fsig_new with BCs from CDR solvers
#endif
  DataOps::copy(sigma_m1, sigma_m);
  DataOps::incr(sigma_m1, Fsig_new, m_dtm[a_m]);
}

void
CdrPlasmaImExSdcStepper::integrateDiffusion(const Real a_dt, const int a_m, const bool a_lagged_terms)
{
  CH_TIME("CdrPlasmaImExSdcStepper::integrateDiffusion");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::integrateDiffusion" << endl;
  }

  // TLDR: We're solving
  //
  // phi_(m+1)^(k+1) = phi_(m)^(k+1,\ast) + dtm*FD_(m+1)^(k+1) + sources.
  //
  // This routine does not modify FD_(m+1)^k. This is replaced by FD_(m+1)^(k+1) later on.
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaImExSdcStepper::getCdrStorage(solver_it);

    if (solver->isDiffusive()) {
      EBAMRCellData&       phi_m1 = storage->getPhi()[a_m + 1]; // Advected solution. Possibly with lagged terms.
      const EBAMRCellData& phi_m  = storage->getPhi()[a_m];

      // Build the diffusion source term
      EBAMRCellData& source    = storage->getScratch();
      EBAMRCellData& init_soln = storage->getScratch2();
      DataOps::setValue(source, 0.0); // No source term

      DataOps::copy(init_soln, phi_m1); // Copy initial solutions
      if (a_lagged_terms) {
        const EBAMRCellData& FD_m1k = storage->getFD()[a_m + 1]; // FD_(m+1)^k. Lagged term.
        DataOps::incr(init_soln, FD_m1k, -m_dtm[a_m]);
      }
      m_amr->conservativeAverage(init_soln, m_fluidRealm, m_cdr->getPhase());
      m_amr->interpGhost(init_soln, m_fluidRealm, m_cdr->getPhase());
      DataOps::copy(phi_m1, phi_m);

      // Solve
      if (m_useTGA) {
        MayDay::Error("CdrPlasmaImExSdcStepper::integrateDiffusion - logic bust");
      }
      else {
        solver->advanceEuler(phi_m1, init_soln, source, m_dtm[a_m]); // No source.
      }
      m_amr->conservativeAverage(phi_m1, m_fluidRealm, m_cdr->getPhase());
      m_amr->interpGhost(phi_m1, m_fluidRealm, m_cdr->getPhase());
      DataOps::floor(phi_m1, 0.0);

      // Update the operator slope
      EBAMRCellData& FD_m1k = storage->getFD()[a_m + 1];
      DataOps::setValue(FD_m1k, 0.0);
      DataOps::incr(FD_m1k, phi_m1, 1.0);
      DataOps::incr(FD_m1k, init_soln, -1.0);
      DataOps::scale(FD_m1k, 1. / m_dtm[a_m]);

      m_amr->conservativeAverage(FD_m1k, m_fluidRealm, m_cdr->getPhase());
      m_amr->interpGhost(FD_m1k, m_fluidRealm, m_cdr->getPhase());
    }
    else {
      EBAMRCellData& FD_m1k = storage->getFD()[a_m + 1];
      DataOps::setValue(FD_m1k, 0.0);
    }
  }
}

void
CdrPlasmaImExSdcStepper::reconcileIntegrands()
{
  CH_TIME("CdrPlasmaImExSdcStepper::reconcileIntegrands");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::reconcileIntegrands" << endl;
  }

  // TLDR: When we come in here, all solutions (Poisson, CDR, RTE, Sigma) are known at node m_p. But we
  //       do need the extra slopes for the explicit operators

  Vector<EBAMRCellData*> cdr_densities_p = CdrPlasmaImExSdcStepper::getCdrSolversPhiK(m_p);
  EBAMRIVData&           sigma_p         = CdrPlasmaImExSdcStepper::getSigmaSolverK(m_p);
  const Real             t_p             = m_tm[m_p];

  // Update boundary conditions for cdr and sigma equations before getting the slope at the final node
  CdrPlasmaImExSdcStepper::computeCdrEbStates(cdr_densities_p);
  CdrPlasmaImExSdcStepper::computeCdrFluxes(cdr_densities_p, t_p);
  CdrPlasmaImExSdcStepper::computeCdrDomainStates(cdr_densities_p);
  CdrPlasmaImExSdcStepper::computeCdrDomainFluxes(cdr_densities_p, t_p);
  CdrPlasmaImExSdcStepper::computeSigmaFlux();

  // Now compute FAR_p - that wasn't done when we integrated
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaImExSdcStepper::getCdrStorage(solver_it);
    const int                  idx     = solver_it.index();

    // This has not been computed yet. Do it.
    EBAMRCellData&       FAR_p = storage->getFAR()[m_p];
    EBAMRCellData&       phi_p = *cdr_densities_p[idx];
    const EBAMRCellData& src   = solver->getSource();

    if (solver->isMobile()) {
      const Real extrap_dt = m_extrapAdvect ? 2.0 * m_extrapDt * m_dtm[m_p - 1]
                                            : 0.0;                     // Factor of 2 because of EBPatchAdvect
      solver->computeDivF(FAR_p, phi_p, extrap_dt, false, true, true); // FAR_p =  Div(v_p*phi_p)
      DataOps::scale(FAR_p, -1.0);                                     // FAR_p = -Div(v_p*phi_p)
    }
    else {
      DataOps::setValue(FAR_p, 0.0);
    }
    DataOps::incr(FAR_p, src, 1.0); // RHS = -Div(v_m*phi_m) + S_m = FAR(phi_m)

    // Build the integrand
    for (int m = 0; m <= m_p; m++) {
      EBAMRCellData& F_m   = storage->getF()[m];
      EBAMRCellData& FD_m  = storage->getFD()[m];
      EBAMRCellData& FAR_m = storage->getFAR()[m];

      DataOps::copy(F_m, FAR_m);
      if (solver->isDiffusive()) {
        DataOps::incr(F_m, FD_m, 1.0);
      }

      // Shouldn't be necessary
      m_amr->conservativeAverage(F_m, m_fluidRealm, m_cdr->getPhase());
      m_amr->interpGhost(F_m, m_fluidRealm, m_cdr->getPhase());
    }
  }

  // Compute Fsig_p - that wasn't done either
  EBAMRIVData& Fnew_p = m_sigmaScratch->getFnew()[m_p];
#if 1 // New code after SigmaSolver -> SurfaceODESolver
  DataOps::copy(Fnew_p, m_sigma->getRHS());
#else
  m_sigma->computeRHS(Fnew_p);
#endif
  for (int m = 0; m <= m_p; m++) {
    EBAMRIVData& Fold_m = m_sigmaScratch->getFold()[m];
    EBAMRIVData& Fnew_m = m_sigmaScratch->getFnew()[m];
    DataOps::copy(Fold_m, Fnew_m);
  }
}

void
CdrPlasmaImExSdcStepper::initializeErrors()
{
  CH_TIME("CdrPlasmaImExSdcStepper::corrector_initializeErrors");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::corrector_initializeErrors" << endl;
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaImExSdcStepper::getCdrStorage(solver_it);
    const int                  idx     = solver_it.index();

    // These should be zero
    if (idx == m_errorIdx || m_errorIdx < 0) {
      EBAMRCellData&       error     = storage->getError();
      const EBAMRCellData& phi_final = storage->getPhi()[m_p];

      DataOps::setValue(error, 0.0);
      DataOps::incr(error, phi_final, -1.0);
    }
  }

  EBAMRIVData&       error       = m_sigmaScratch->getError();
  const EBAMRIVData& sigma_final = m_sigmaScratch->getSigmaSolver()[m_p];
  DataOps::setValue(error, 0.0);
  DataOps::incr(error, sigma_final, -1.0);
}

void
CdrPlasmaImExSdcStepper::finalizeErrors()
{
  CH_TIME("CdrPlasmaImExSdcStepper::corrector_finalizeErrors");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::corrector_finalizeErrors" << endl;
  }

  const Real safety = 1.E-20;

  m_maxError = 0.0;
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaImExSdcStepper::getCdrStorage(solver_it);
    const int                  idx     = solver_it.index();

    // Compute error
    if (idx == m_errorIdx || m_errorIdx < 0) {
      EBAMRCellData&       error = storage->getError();
      const EBAMRCellData& phi_p = storage->getPhi()[m_p];
      DataOps::incr(error, phi_p, 1.0);

      // Compute norms. Only coarsest level
      Real      Lerr, Lphi;
      const int lvl = 0;
      Lerr          = DataOps::norm(*error[lvl], m_errorNorm);
      Lphi          = DataOps::norm(*phi_p[lvl], m_errorNorm);

      if (Lphi > 0.0) {
        m_cdrError[idx] = Lerr / Lphi;

        m_maxError = Max(m_cdrError[idx], m_maxError);
      }
#if 0 // Debug
      if(procID() == 0){
	std::cout << "Lerr = " << Lerr << "\t Lphi = " << Lphi << "\t Lerr/Lphi = " << Lerr/Lphi << std::endl;
      }
#endif
    }
  }

  // Override if
  if (m_errorIdx >= 0) {
    m_maxError = m_cdrError[m_errorIdx];
  }

  // Compute the surface charge conservation error
  EBAMRIVData&       error       = m_sigmaScratch->getError();
  const EBAMRIVData& sigma_final = m_sigmaScratch->getSigmaSolver()[m_p];
  DataOps::incr(error, sigma_final, 1.0);
  m_sigmaError = 0.0; // I don't think this is ever used...
}

void
CdrPlasmaImExSdcStepper::computeNewDt(bool& a_accept_step, const Real a_dt, const int a_num_corrections)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeNewDt");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeNewDt" << endl;
  }

  // If a_dt was the smallest possible CFL or hardcap time step, we just have to accept it
  const Real max_gl_dist = CdrPlasmaImExSdcStepper::getMaxNodeDistance();
  Real       dt_cfl      = 2.0 * m_dtCFL / max_gl_dist; // This is the smallest time step ON THE FINEST LEVEL

  // Try time step
  const Real rel_err    = (m_safety * m_errThresh) / m_maxError;
  const Real dt_adapt   = (m_maxError > 0.0) ? a_dt * pow(rel_err, 1.0 / (a_num_corrections + 1)) : m_maxDt;
  const Real min_dt_cfl = dt_cfl * m_minCFL;
  const Real max_dt_cfl = dt_cfl * m_maxCFL;

  if (m_maxError <= m_errThresh) { // Always accept, and compute new step
    a_accept_step = true;

    // Do not grow step too fast
    if (rel_err > 1.0) { // rel_err > 1 => dt_adapt > a_dt
      m_newDt = Min(m_maxDtGrowth * a_dt, dt_adapt);
    }
    else { // rel_err > 1 => dt_adapt < a_dt. This shrinks the error down to the safety factor.
      m_newDt = dt_adapt;
    }
    m_newDt = Max(m_newDt, m_minDt); // Don't go below hardcap
    m_newDt = Min(m_newDt, m_maxDt); // Don't go above other hardcap

    m_newDt = Max(m_newDt, dt_cfl * m_minCFL); // Don't drop below minimum CFL
    m_newDt = Min(m_newDt, dt_cfl * m_maxCFL); // Don't go above maximum CFL
  }
  else {
    a_accept_step = false;

    m_newDt = m_decreaseSafety * dt_adapt;      // Decrease time step a little bit extra to avoid another rejection
    if (a_dt <= min_dt_cfl || a_dt < m_minDt) { // Step already at minimum. Accept it anyways.
      a_accept_step = true;
    }

    m_newDt = Max(m_newDt, dt_cfl * m_minCFL); // Don't drop below minimum CFL
    m_newDt = Min(m_newDt, dt_cfl * m_maxCFL); // Don't go above maximum CFL

    m_newDt = Max(m_newDt, m_minDt); // Don't go below hardcap
    m_newDt = Min(m_newDt, m_maxDt); // Don't go above other hardcap
  }

#if 0 // Debug
  if(procID() == 0) std::cout << "accept = " << a_accept_step
			      << " dt = " << a_dt
			      << " new_dt = " << m_newDt
			      << " fraction = " << m_newDt/a_dt
			      << std::endl;
#endif

  m_haveDtErr = true;
}

void
CdrPlasmaImExSdcStepper::adaptiveReport(const Real a_first_dt,
                                        const Real a_dt,
                                        const Real a_new_dt,
                                        const int  a_corr,
                                        const int  a_rej,
                                        const Real a_max_err)
{
  CH_TIME("CdrPlasmaImExSdcStepper::adaptiveReport");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::adaptiveReport" << endl;
  }

  pout() << "\n";
  pout() << "CdrPlasmaImExSdcStepper::adaptiveReport breakdown" << endl;
  pout() << "--------------------------------\n";
  pout() << "\t Try dt       = " << a_first_dt << endl;
  pout() << "\t Advanced dt  = " << a_dt << endl;
  pout() << "\t New dt       = " << a_new_dt << endl;
  pout() << "\t Subintervals = " << m_p << endl;
  pout() << "\t Corrections  = " << a_corr << endl;
  pout() << "\t Rejections   = " << a_rej << endl;
  pout() << "\t Max error    = " << a_max_err << endl;
  pout() << "\n";
}

Real
CdrPlasmaImExSdcStepper::computeDt()
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeDt" << endl;
  }

  Real dt   = std::numeric_limits<Real>::max();
  Real a_dt = std::numeric_limits<Real>::max();

  int Nref = 1;
  for (int lvl = 0; lvl < m_amr->getFinestLevel(); lvl++) {
    Nref = Nref * m_amr->getRefinementRatios()[lvl];
  }
  const Real max_gl_dist = CdrPlasmaImExSdcStepper::getMaxNodeDistance();
  m_dtCFL                = m_cdr->computeAdvectionDt();

  Real dt_cfl = 2.0 * m_dtCFL / max_gl_dist;

  // Time step selection for non-adaptive stepping
  if (!m_adaptiveDt) {
    if (dt_cfl < dt) {
      dt         = m_cfl * dt_cfl;
      m_timeCode = TimeCode::Advection;
    }
  }
  else {
    Real new_dt;

    // Step should not exceed m_newDt. Also, it shoul
    if (m_haveDtErr) {
      new_dt = m_newDt;

      new_dt = Max(new_dt, dt_cfl * m_minCFL);
      new_dt = Min(new_dt, dt_cfl * m_maxCFL);
    }
    else {
      new_dt = m_minCFL * dt_cfl;
    }

    if (new_dt < dt) {
      dt         = new_dt;
      m_timeCode = TimeCode::Error;
    }
  }

  // EVERYTHING BELOW HERE IS "STANDARD"
  // -----------------------------------
  const Real dt_relax = m_relaxTime * this->computeRelaxationTime();
  if (dt_relax < dt) {
    dt         = dt_relax;
    m_timeCode = TimeCode::RelaxationTime;
  }

  if (dt < m_minDt) {
    dt         = m_minDt;
    m_timeCode = TimeCode::Hardcap;
  }

  if (dt > m_maxDt) {
    dt         = m_maxDt;
    m_timeCode = TimeCode::Hardcap;
  }

  a_dt = dt;

  // Copy the time code, it is needed for diagnostics
  m_timeCode = m_timeCode;

#if 0 // Debug
  if(procID() == 0){
    std::cout << "computeDt = " << a_dt << "\t m_newDt = " << m_newDt << std::endl; 
  }
#endif

  return a_dt;
}

void
CdrPlasmaImExSdcStepper::regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("CdrPlasmaImExSdcStepper::regridInternals(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::regridInternals(int, int, int)" << endl;
  }
}

void
CdrPlasmaImExSdcStepper::allocateInternals()
{
  CH_TIME("CdrPlasmaImExSdcStepper::allocateInternals");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::allocateInternals" << endl;
  }

  m_cdrError.resize(m_physics->getNumCdrSpecies());

  CdrPlasmaImExSdcStepper::allocateCdrStorage();
  CdrPlasmaImExSdcStepper::allocateFieldStorage();
  CdrPlasmaImExSdcStepper::allocateRtStorage();
  CdrPlasmaImExSdcStepper::allocateSigmaStorage();

  CdrPlasmaImExSdcStepper::setupQuadratureNodes(m_p);
  CdrPlasmaImExSdcStepper::setupQmj(m_p);
}

void
CdrPlasmaImExSdcStepper::allocateCdrStorage()
{
  const int ncomp       = 1;
  const int num_species = m_physics->getNumCdrSpecies();

  m_cdrScratch.resize(num_species);

  for (CdrIterator<CdrSolver> solver_it(*m_cdr); solver_it.ok(); ++solver_it) {
    const int idx     = solver_it.index();
    m_cdrScratch[idx] = RefCountedPtr<CdrStorage>(new CdrStorage(m_amr, m_fluidRealm, m_cdr->getPhase(), ncomp));
    m_cdrScratch[idx]->allocateStorage(m_p);
  }
}

void
CdrPlasmaImExSdcStepper::allocateFieldStorage()
{
  const int ncomp = 1;
  m_fieldScratch  = RefCountedPtr<FieldStorage>(new FieldStorage(m_amr, m_fluidRealm, m_cdr->getPhase(), ncomp));
  m_fieldScratch->allocateStorage(m_p);
}

void
CdrPlasmaImExSdcStepper::allocateRtStorage()
{
  const int ncomp       = 1;
  const int num_Photons = m_physics->getNumRtSpecies();
  m_rteScratch.resize(num_Photons);

  for (RtIterator<RtSolver> solver_it(*m_rte); solver_it.ok(); ++solver_it) {
    const int idx     = solver_it.index();
    m_rteScratch[idx] = RefCountedPtr<RtStorage>(new RtStorage(m_amr, m_fluidRealm, m_rte->getPhase(), ncomp));
    m_rteScratch[idx]->allocateStorage(m_p);
  }
}

void
CdrPlasmaImExSdcStepper::allocateSigmaStorage()
{
  const int ncomp = 1;
  m_sigmaScratch  = RefCountedPtr<SigmaStorage>(new SigmaStorage(m_amr, m_fluidRealm, m_cdr->getPhase(), ncomp));
  m_sigmaScratch->allocateStorage(m_p);
}

void
CdrPlasmaImExSdcStepper::deallocateInternals()
{
  CH_TIME("CdrPlasmaImExSdcStepper::deallocateInternals");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::deallocateInternals" << endl;
  }

  for (CdrIterator<CdrSolver> solver_it(*m_cdr); solver_it.ok(); ++solver_it) {
    const int idx = solver_it.index();
    m_cdrScratch[idx]->deallocateStorage();
    m_cdrScratch[idx] = RefCountedPtr<CdrStorage>(0);
  }

  for (RtIterator<RtSolver> solver_it(*m_rte); solver_it.ok(); ++solver_it) {
    const int idx = solver_it.index();
    m_rteScratch[idx]->deallocateStorage();
    m_rteScratch[idx] = RefCountedPtr<RtStorage>(0);
  }

  m_cdrScratch.resize(0);
  m_rteScratch.resize(0);

  m_fieldScratch->deallocateStorage();
  m_fieldScratch = RefCountedPtr<FieldStorage>(0);

  m_sigmaScratch->deallocateStorage();
  m_sigmaScratch = RefCountedPtr<SigmaStorage>(0);
}

void
CdrPlasmaImExSdcStepper::computeElectricFieldIntoScratch()
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeElectricFieldIntoScratch");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeElectricFieldIntoScratch" << endl;
  }

  EBAMRCellData& E_cell = m_fieldScratch->getElectricFieldCell();
  EBAMRFluxData& E_face = m_fieldScratch->getElectricFieldFace();
  EBAMRIVData&   E_eb   = m_fieldScratch->getElectricFieldEb();
  EBAMRIFData&   E_dom  = m_fieldScratch->getElectricFieldDomain();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();

  CdrPlasmaImExSdcStepper::computeElectricField(E_cell, m_cdr->getPhase(), phi);    // Compute cell-centered field
  CdrPlasmaImExSdcStepper::computeElectricField(E_face, m_cdr->getPhase(), E_cell); // Compute face-centered field
  CdrPlasmaImExSdcStepper::computeElectricField(E_eb, m_cdr->getPhase(), E_cell);   // EB-centered field

  CdrPlasmaStepper::extrapolateToDomainFaces(E_dom, m_cdr->getPhase(), E_cell);
}

void
CdrPlasmaImExSdcStepper::computeCdrGradients()
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrGradients");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrGradients" << endl;
  }

  CdrPlasmaImExSdcStepper::computeCdrGradients(m_cdr->getPhis());
}

void
CdrPlasmaImExSdcStepper::computeCdrGradients(const Vector<EBAMRCellData*>& a_phis)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrGradients");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrGradients" << endl;
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    const int                  idx     = solver_it.index();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaImExSdcStepper::getCdrStorage(solver_it);
    EBAMRCellData&             grad    = storage->getGradient();
    m_amr->computeGradient(grad, *a_phis[idx], m_fluidRealm, m_cdr->getPhase());
    //    m_amr->conservativeAverage(grad, m_fluidRealm, m_cdr->getPhase());
    m_amr->interpGhost(grad, m_fluidRealm, m_cdr->getPhase());
  }
}

void
CdrPlasmaImExSdcStepper::computeCdrVelo(const Real a_time)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrVelo");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrVelo" << endl;
  }

  CdrPlasmaImExSdcStepper::computeCdrVelo(m_cdr->getPhis(), a_time);
}

void
CdrPlasmaImExSdcStepper::computeCdrVelo(const Vector<EBAMRCellData*>& a_phis, const Real a_time)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrVelo(Vector<EBAMRCellData*>, Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrVelo(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->getVelocities();
  CdrPlasmaImExSdcStepper::computeCdrDriftVelocities(velocities,
                                                     a_phis,
                                                     m_fieldScratch->getElectricFieldCell(),
                                                     a_time);
}

void
CdrPlasmaImExSdcStepper::computeCdrEbStates()
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrEbStates");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrEbStates" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<CdrSolver>& solver  = solver_it();
    RefCountedPtr<CdrStorage>&      storage = CdrPlasmaImExSdcStepper::getCdrStorage(solver_it);

    cdr_states.push_back(&(solver->getPhi()));
    eb_states.push_back(&(storage->getEbState()));
    eb_gradients.push_back(&(storage->getEbGrad()));
    cdr_gradients.push_back(&(storage->getGradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  CdrPlasmaImExSdcStepper::extrapolateToEb(eb_states, m_cdr->getPhase(), cdr_states);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    const int idx = solver_it.index();
    DataOps::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux.
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_fluidRealm, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++) {
    CdrPlasmaImExSdcStepper::extrapolateToEb(eb_gradient, m_cdr->getPhase(), *cdr_gradients[i]);
    CdrPlasmaImExSdcStepper::projectFlux(*eb_gradients[i], eb_gradient);
  }
}

void
CdrPlasmaImExSdcStepper::computeCdrEbStates(const Vector<EBAMRCellData*>& a_phis)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrEbStates(vec)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrEbStates(vec)" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_gradients;

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaImExSdcStepper::getCdrStorage(solver_it);

    eb_states.push_back(&(storage->getEbState()));
    eb_gradients.push_back(&(storage->getEbGrad()));
    cdr_gradients.push_back(&(storage->getGradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  CdrPlasmaImExSdcStepper::extrapolateToEb(eb_states, m_cdr->getPhase(), a_phis);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    const int idx = solver_it.index();
    DataOps::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux.
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_fluidRealm, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < a_phis.size(); i++) {
    CdrPlasmaImExSdcStepper::extrapolateToEb(eb_gradient, m_cdr->getPhase(), *cdr_gradients[i]);
    CdrPlasmaImExSdcStepper::projectFlux(*eb_gradients[i], eb_gradient);
  }
}

void
CdrPlasmaImExSdcStepper::computeCdrDomainStates()
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrDomainStates");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrDomainStates" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<CdrSolver>& solver  = solver_it();
    RefCountedPtr<CdrStorage>&      storage = CdrPlasmaImExSdcStepper::getCdrStorage(solver_it);

    cdr_states.push_back(&(solver->getPhi()));
    domain_states.push_back(&(storage->getDomainState()));
    domain_gradients.push_back(&(storage->getDomainGrad()));
    cdr_gradients.push_back(&(storage->getGradient())); // Should already be computed
  }

  // Extrapolate states to the domain faces
  CdrPlasmaImExSdcStepper::extrapolateToDomainFaces(domain_states, m_cdr->getPhase(), cdr_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux.
  EBAMRIFData grad;
  m_amr->allocate(grad, m_fluidRealm, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++) {
    CdrPlasmaImExSdcStepper::extrapolateToDomainFaces(grad, m_cdr->getPhase(), *cdr_gradients[i]);
    CdrPlasmaImExSdcStepper::projectDomain(*domain_gradients[i], grad);
  }
}

void
CdrPlasmaImExSdcStepper::computeCdrDomainStates(const Vector<EBAMRCellData*>& a_phis)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrDomainStates");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrDomainStates" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_gradients;

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<CdrSolver>& solver  = solver_it();
    RefCountedPtr<CdrStorage>&      storage = this->getCdrStorage(solver_it);

    domain_states.push_back(&(storage->getDomainState()));
    domain_gradients.push_back(&(storage->getDomainGrad()));
    cdr_gradients.push_back(&(storage->getGradient()));
  }

  // Extrapolate states to the domain faces
  this->extrapolateToDomainFaces(domain_states, m_cdr->getPhase(), a_phis);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux.
  EBAMRIFData grad;
  m_amr->allocate(grad, m_fluidRealm, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < a_phis.size(); i++) {
    this->extrapolateToDomainFaces(grad, m_cdr->getPhase(), *cdr_gradients[i]);
    this->projectDomain(*domain_gradients[i], grad);
  }
}

void
CdrPlasmaImExSdcStepper::computeCdrFluxes(const Real a_time)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrFluxes");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrFluxes" << endl;
  }

  this->computeCdrFluxes(m_cdr->getPhis(), a_time);
}

void
CdrPlasmaImExSdcStepper::computeCdrFluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrFluxes(Vector<EBAMRCellData*>, Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrFluxes(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRIVData*> cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_densities;
  Vector<EBAMRIVData*> extrap_cdr_velocities;
  Vector<EBAMRIVData*> extrap_cdr_gradients;
  Vector<EBAMRIVData*> extrap_rte_fluxes;

  cdr_fluxes = m_cdr->getEbFlux();

  for (CdrIterator<CdrSolver> solver_it(*m_cdr); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrStorage>& storage = this->getCdrStorage(solver_it);

    EBAMRIVData& dens_eb = storage->getEbState();
    EBAMRIVData& velo_eb = storage->getEbVelo();
    EBAMRIVData& flux_eb = storage->getEbFlux();
    EBAMRIVData& grad_eb = storage->getEbGrad();

    extrap_cdr_densities.push_back(&dens_eb);  // Computed in computeCdrEbStates
    extrap_cdr_velocities.push_back(&velo_eb); // Not yet computed
    extrap_cdr_fluxes.push_back(&flux_eb);     // Not yet computed
    extrap_cdr_gradients.push_back(&grad_eb);  // Computed in computeCdrEbStates
  }

  // Extrapolate densities, velocities, and fluxes
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->getVelocities();
  CdrPlasmaStepper::computeExtrapolatedFluxes(extrap_cdr_fluxes, a_phis, cdr_velocities, m_cdr->getPhase());
  CdrPlasmaStepper::computeExtrapolatedVelocities(extrap_cdr_velocities, cdr_velocities, m_cdr->getPhase());

  // Compute RTE flux on the boundary
  for (RtIterator<RtSolver> solver_it(*m_rte); solver_it.ok(); ++solver_it) {
    RefCountedPtr<RtSolver>&  solver  = solver_it();
    RefCountedPtr<RtStorage>& storage = this->getRtStorage(solver_it);

    EBAMRIVData& flux_eb = storage->getEbFlux();
    solver->computeBoundaryFlux(flux_eb, solver->getPhi());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_fieldScratch->getElectricFieldEb();
  CdrPlasmaStepper::computeCdrFluxes(cdr_fluxes,
                                     extrap_cdr_fluxes,
                                     extrap_cdr_densities,
                                     extrap_cdr_velocities,
                                     extrap_cdr_gradients,
                                     extrap_rte_fluxes,
                                     E,
                                     a_time);
}

void
CdrPlasmaImExSdcStepper::computeCdrDomainFluxes(const Real a_time)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrDomainFluxes");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrDomainFluxes" << endl;
  }

  this->computeCdrDomainFluxes(m_cdr->getPhis(), a_time);
}

void
CdrPlasmaImExSdcStepper::computeCdrDomainFluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeCdrDomainFluxes(Vector<EBAMRCellData*>, Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeCdrDomainFluxes(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRIFData*> cdr_fluxes;
  Vector<EBAMRIFData*> extrap_cdr_fluxes;
  Vector<EBAMRIFData*> extrap_cdr_densities;
  Vector<EBAMRIFData*> extrap_cdr_velocities;
  Vector<EBAMRIFData*> extrap_cdr_gradients;
  Vector<EBAMRIFData*> extrap_rte_fluxes;

  Vector<EBAMRCellData*> cdr_velocities;
  Vector<EBAMRCellData*> cdr_gradients;

  cdr_fluxes     = m_cdr->getDomainFlux();
  cdr_velocities = m_cdr->getVelocities();
  for (CdrIterator<CdrSolver> solver_it(*m_cdr); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrStorage>& storage = this->getCdrStorage(solver_it);

    EBAMRIFData&   dens_domain = storage->getDomainState();
    EBAMRIFData&   velo_domain = storage->getDomainVelo();
    EBAMRIFData&   flux_domain = storage->getDomainFlux();
    EBAMRIFData&   grad_domain = storage->getDomainGrad();
    EBAMRCellData& gradient    = storage->getGradient();

    extrap_cdr_densities.push_back(&dens_domain);  // Has not been computed
    extrap_cdr_velocities.push_back(&velo_domain); // Has not been computed
    extrap_cdr_fluxes.push_back(&flux_domain);     // Has not been computed
    extrap_cdr_gradients.push_back(&grad_domain);  // Has not been computed
    cdr_gradients.push_back(&gradient);
  }

  // Compute extrapolated velocities and fluxes at the domain faces
  this->extrapolateToDomainFaces(extrap_cdr_densities, m_cdr->getPhase(), a_phis);
  this->extrapolateVelocitiesToDomainFaces(extrap_cdr_velocities, m_cdr->getPhase(), cdr_velocities);
  this->computeExtrapolatedDomainFluxes(extrap_cdr_fluxes, a_phis, cdr_velocities, m_cdr->getPhase());
  this->extrapolateVectorToDomainFaces(extrap_cdr_gradients, m_cdr->getPhase(), cdr_gradients);

  // Compute RTE flux on domain faces
  for (RtIterator<RtSolver> solver_it(*m_rte); solver_it.ok(); ++solver_it) {
    RefCountedPtr<RtSolver>&  solver  = solver_it();
    RefCountedPtr<RtStorage>& storage = this->getRtStorage(solver_it);

    EBAMRIFData& domain_flux = storage->getDomainFlux();
    solver->computeDomainFlux(domain_flux, solver->getPhi());
    extrap_rte_fluxes.push_back(&domain_flux);
  }

  const EBAMRIFData& E = m_fieldScratch->getElectricFieldDomain();

  // This fills the solvers' domain fluxes
  CdrPlasmaStepper::computeCdrDomainFluxes(cdr_fluxes,
                                           extrap_cdr_fluxes,
                                           extrap_cdr_densities,
                                           extrap_cdr_velocities,
                                           extrap_cdr_gradients,
                                           extrap_rte_fluxes,
                                           E,
                                           a_time);
}

void
CdrPlasmaImExSdcStepper::computeSigmaFlux()
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeSigmaFlux");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeSigmaFlux" << endl;
  }

  EBAMRIVData& flux = m_sigma->getRHS();
  DataOps::setValue(flux, 0.0);

  for (CdrIterator<CdrSolver> solver_it(*m_cdr); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<CdrSolver>&  solver      = solver_it();
    const RefCountedPtr<CdrSpecies>& spec        = solver_it.getSpecies();
    const EBAMRIVData&               solver_flux = solver->getEbFlux();

    DataOps::incr(flux, solver_flux, spec->getChargeNumber() * Units::Qe);
  }

  m_sigma->resetElectrodes(flux, 0.0);
}

void
CdrPlasmaImExSdcStepper::computeReactionNetwork(const int a_m, const Real a_time, const Real a_dt)
{
  CH_TIME("CdrPlasmaImExSdcStepper::computeReactionNetwork");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::computeReactionNetwork";
  }

  Vector<EBAMRCellData*> cdr_sources = m_cdr->getSources();
  Vector<EBAMRCellData*> rte_sources = m_rte->getSources();

  const Vector<EBAMRCellData*> cdr_densities = getCdrSolversPhiK(a_m);
  const Vector<EBAMRCellData*> rte_densities = m_rte->getPhis();
  const EBAMRCellData&         E             = m_fieldScratch->getElectricFieldCell();

  CdrPlasmaStepper::advanceReactionNetwork(cdr_sources, rte_sources, cdr_densities, rte_densities, E, a_time, a_dt);
}

void
CdrPlasmaImExSdcStepper::updateField()
{
  CH_TIME("CdrPlasmaImExSdcStepper::updateField(solver)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::updateField(solver)" << endl;
  }

  if (m_doPoisson) { // Solve Poisson equation
    if ((m_timeStep + 1) % m_fastPoisson == 0) {
      CdrPlasmaStepper::solvePoisson();
      this->computeElectricFieldIntoScratch();
    }
  }
}

void
CdrPlasmaImExSdcStepper::updateField(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma)
{
  CH_TIME("CdrPlasmaImExSdcStepper::updateField(full)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::updateField(full)" << endl;
  }

  if (m_doPoisson) { // Solve Poisson equation
    if ((m_timeStep + 1) % m_fastPoisson == 0) {
      CdrPlasmaStepper::solvePoisson(m_fieldSolver->getPotential(), m_fieldSolver->getRho(), a_densities, a_sigma);
      this->computeElectricFieldIntoScratch();
    }
  }
}

void
CdrPlasmaImExSdcStepper::integrateRtTransient(const Real a_dt)
{
  CH_TIME("CdrPlasmaImExSdcStepper::integrateRtTransient");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::integrateRtTransient" << endl;
  }

  if (m_doRTE) {
    if ((m_timeStep + 1) % m_fastRTE == 0) {
      if (!(m_rte->isStationary())) {
        for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
          RefCountedPtr<RtSolver>& solver = solver_it();
          solver->advance(a_dt);
        }
      }
    }
  }
}

void
CdrPlasmaImExSdcStepper::integrateRtStationary()
{
  CH_TIME("CdrPlasmaImExSdcStepper::integrateRtTransient");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::integrateRtTransient" << endl;
  }

  if (m_doRTE) {
    if ((m_timeStep + 1) % m_fastRTE == 0) {
      if ((m_rte->isStationary())) {
        for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
          RefCountedPtr<RtSolver>& solver = solver_it();
          solver->advance(0.0);
        }
      }
    }
  }
}

void
CdrPlasmaImExSdcStepper::updateDiffusionCoefficients()
{
  CH_TIME("CdrPlasmaImExSdcStepper::updateDiffusionCoefficients");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::updateDiffusionCoefficients" << endl;
  }
  CdrPlasmaStepper::computeCdrDiffusion(m_fieldScratch->getElectricFieldCell(), m_fieldScratch->getElectricFieldEb());
}

Vector<EBAMRCellData*>
CdrPlasmaImExSdcStepper::getCdrErrors()
{
  CH_TIME("CdrPlasmaImExSdcStepper::getCdrErrors");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::getCdrErrors" << endl;
  }

  Vector<EBAMRCellData*> ret;
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaImExSdcStepper::getCdrStorage(solver_it);
    ret.push_back(&(storage->getError()));
  }

  return ret;
}

Vector<EBAMRCellData*>
CdrPlasmaImExSdcStepper::getCdrSolversPhiK(const int a_m)
{
  CH_TIME("CdrPlasmaImExSdcStepper::getCdrSolversPhiK");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::getCdrSolversPhiK" << endl;
  }

  Vector<EBAMRCellData*> ret;
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaImExSdcStepper::getCdrStorage(solver_it);
    ret.push_back(&(storage->getPhi()[a_m]));
  }

  return ret;
}

EBAMRIVData&
CdrPlasmaImExSdcStepper::getSigmaSolverK(const int a_m)
{
  CH_TIME("CdrPlasmaImExSdcStepper::getSigmaSolverK");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::getSigmaSolverK)" << endl;
  }
  return m_sigmaScratch->getSigmaSolver()[a_m];
}

void
CdrPlasmaImExSdcStepper::writeStepProfile(const Real a_dt,
                                          const Real a_error,
                                          const int  a_substeps,
                                          const int  a_corrections,
                                          const int  a_rejections)
{
  CH_TIME("sissdc::writeStepProfile");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::writeStepProfile" << endl;
  }

  if (procID() == 0) {

    const std::string fname("CdrPlasmaImExSdcStepper_step_profile.txt");

    bool write_header;
    { // Write header if we must
      std::ifstream infile(fname);
      write_header = infile.peek() == std::ifstream::traits_type::eof() ? true : false;
    }

    // Write output
    std::ofstream f;
    f.open(fname, std::ios_base::app);
    const int width = 12;

    if (write_header) {
      f << std::left << std::setw(width) << "# Step"
        << "\t" << std::left << std::setw(width) << "Time"
        << "\t" << std::left << std::setw(width) << "dt"
        << "\t" << std::left << std::setw(width) << "Substeps"
        << "\t" << std::left << std::setw(width) << "Corrections"
        << "\t" << std::left << std::setw(width) << "Rejections"
        << "\t" << std::left << std::setw(width) << "Error"
        << "\t" << std::left << std::setw(width) << "CFL"
        << "\t" << endl;
    }

    f << std::left << std::setw(width) << m_timeStep << "\t" << std::left << std::setw(width) << m_time << "\t"
      << std::left << std::setw(width) << a_dt << "\t" << std::left << std::setw(width) << a_substeps << "\t"
      << std::left << std::setw(width) << a_corrections << "\t" << std::left << std::setw(width) << a_rejections << "\t"
      << std::left << std::setw(width) << m_maxError << "\t" << std::left << std::setw(width) << m_dtCFL << "\t"
      << endl;
  }
}

void
CdrPlasmaImExSdcStepper::storeSolvers()
{
  CH_TIME("CdrPlasmaImExSdcStepper::storeSolvers");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::storeSolvers" << endl;
  }

  if (m_k > 0 && m_adaptiveDt) {
    // IMEX_SDC does not manipulate cdr and sigma solvers until the end of the time step. Only need to do
    // Poisson and RTE here.

    // Poisson
    MFAMRCellData&       previous = m_fieldScratch->getPrevious();
    const MFAMRCellData& state    = m_fieldSolver->getPotential();
    DataOps::copy(previous, state);

    // RTE
    for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
      RefCountedPtr<RtStorage>&      storage = CdrPlasmaImExSdcStepper::getRtStorage(solver_it);
      const RefCountedPtr<RtSolver>& solver  = solver_it();

      EBAMRCellData&       previous = storage->getPrevious();
      const EBAMRCellData& state    = solver->getPhi();

      DataOps::copy(previous, state);
    }
  }
}

void
CdrPlasmaImExSdcStepper::restoreSolvers()
{
  CH_TIME("CdrPlasmaImExSdcStepper::restoreSolvers");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaImExSdcStepper::restoreSolvers" << endl;
  }

  // IMEX_SDC does not manipulate cdr and sigma solvers until the end of the time step. Only need to do
  // Poisson and RTE here.

  // Poisson
  MFAMRCellData&       state    = m_fieldSolver->getPotential();
  const MFAMRCellData& previous = m_fieldScratch->getPrevious();

  DataOps::copy(state, previous);

  // RTE
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<RtStorage>& storage = CdrPlasmaImExSdcStepper::getRtStorage(solver_it);
    RefCountedPtr<RtSolver>&  solver  = solver_it();

    EBAMRCellData& previous = storage->getPrevious();
    EBAMRCellData& state    = solver->getPhi();

    DataOps::copy(state, previous);
  }
}

#include <CD_NamespaceFooter.H>
