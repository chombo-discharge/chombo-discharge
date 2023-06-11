/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaGodunovStepper.cpp
  @brief  Implementation of CD_CdrPlasmaGodunovStepper.H
  @author Robert Marskar
*/

// Std includes
#include <limits>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaGodunovStorage.H>
#include <CD_DischargeIO.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

typedef CdrPlasmaGodunovStepper::CdrStorage   CdrStorage;
typedef CdrPlasmaGodunovStepper::FieldStorage FieldStorage;
typedef CdrPlasmaGodunovStepper::RtStorage    RtStorage;
typedef CdrPlasmaGodunovStepper::SigmaStorage SigmaStorage;

CdrPlasmaGodunovStepper::CdrPlasmaGodunovStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics)
{
  CH_TIME("CdrPlasmaGodunovStepper::CdrPlasmaGodunovStepper()");

  // Default settings
  m_className    = "CdrPlasmaGodunovStepper";
  m_physics      = a_physics;
  m_extrapAdvect = true;
  m_regridSlopes = true;
}

CdrPlasmaGodunovStepper::~CdrPlasmaGodunovStepper()
{
  CH_TIME("CdrPlasmaGodunovStepper::~CdrPlasmaGodunovStepper()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::~CdrPlasmaGodunovStepper" << endl;
  }
}

void
CdrPlasmaGodunovStepper::parseOptions()
{
  CH_TIME("CdrPlasmaGodunovStepper::parseOptions()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::parseOptions()" << endl;
  }

  // Parse various options.
  this->parseVerbosity();
  this->parseSolverVerbosity();
  this->parseFastPoisson();
  this->parseFastRadiativeTransfer();
  this->parseCFL();
  this->parseRelaxationTime();
  this->parseMinDt();
  this->parseMaxDt();
  this->parseSourceComputation();
  this->parseField();
  this->parseAdvection();
  this->parseDiffusion();
  this->parseFloor();
  this->parseDebug();
  this->parseProfile();
  this->parseFHD();
  this->parseRegridSlopes();
}

void
CdrPlasmaGodunovStepper::parseRuntimeOptions()
{
  CH_TIME("CdrPlasmaGodunovStepper::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::parseRuntimeOptions()" << endl;
  }

  this->parseVerbosity();
  this->parseSolverVerbosity();
  this->parseFastPoisson();
  this->parseFastRadiativeTransfer();
  this->parseCFL();
  this->parseRelaxationTime();
  this->parseMinDt();
  this->parseMaxDt();
  this->parseSourceComputation();
  this->parseField();
  this->parseAdvection();
  this->parseDiffusion();
  this->parseFloor();
  this->parseDebug();
  this->parseProfile();
  this->parseFHD();
  this->parseRegridSlopes();

  // Solvers also parse their runtime options.
  m_cdr->parseRuntimeOptions();
  m_rte->parseRuntimeOptions();
  m_fieldSolver->parseRuntimeOptions();
  m_sigma->parseRuntimeOptions();

  // Physics also parses run-time options
  m_physics->parseRuntimeOptions();
}

void
CdrPlasmaGodunovStepper::parseDiffusion()
{
  CH_TIME("CdrPlasmaGodunovStepper::parseDiffusion()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::parseDiffusion()" << endl;
  }

  ParmParse pp(m_className.c_str());

  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  m_useImplicitDiffusion.resize(numCdrSpecies, false);

  std::string str;
  pp.get("diffusion", str);
  if (str == "explicit") {
    m_diffusionAlgorithm = DiffusionAlgorithm::Explicit;

    // Set explicit diffusion for all species
    m_useImplicitDiffusion.resize(numCdrSpecies, false);
  }
  else if (str == "implicit") {
    m_diffusionAlgorithm = DiffusionAlgorithm::Implicit;

    // Set implicit diffusion for all species
    m_useImplicitDiffusion.resize(numCdrSpecies, true);
  }
  else if (str == "auto") {
    m_diffusionAlgorithm = DiffusionAlgorithm::Automatic;

    // Diffusion will be determined during time step computation, but for now we initialize it to explicit.
    m_useImplicitDiffusion.resize(numCdrSpecies, false);
  }
  else {
    MayDay::Error("CdrPlasmaGodunovStepper::parseDiffusion - unknown diffusion type requested");
  }

  m_diffusionOrder = -1;
  pp.get("diffusion_order", m_diffusionOrder);
  if (!(m_diffusionOrder == 1 || m_diffusionOrder == 2)) {
    MayDay::Error("CdrPlasmaGodunovStepper::parseDiffusion - unknown diffusion order requested. Must be 1 or 2. ");
  }
}

void
CdrPlasmaGodunovStepper::parseField()
{
  CH_TIME("CdrPlasmaGodunovStepper::parseField()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::parseField()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("field_coupling", str);
  if (str == "explicit") {
    m_fieldCoupling = FieldCoupling::Explicit;
  }
  else if (str == "semi_implicit") {
    m_fieldCoupling = FieldCoupling::SemiImplicit;
  }
  else {
    MayDay::Error("CdrPlasmaGodunovStepper::parseField - unknown transport algorithm requested");
  }
}

void
CdrPlasmaGodunovStepper::parseAdvection()
{
  CH_TIME("CdrPlasmaGodunovStepper::parseAdvection()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::parseAdvection()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("advection", str);

  if (str == "euler") {
    m_advectionSolver = AdvectionSolver::Euler;
  }
  else if (str == "rk2") {
    m_advectionSolver = AdvectionSolver::RK2;
  }
  else if (str == "muscl") {
    m_advectionSolver = AdvectionSolver::MUSCL;
  }
  else {
    MayDay::Error("CdrPlasmaGodunovStepper::parseAdvection - unknown argument");
  }
}

void
CdrPlasmaGodunovStepper::parseFloor()
{
  CH_TIME("CdrPlasmaGodunovStepper::parseFloor()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::parseFloor()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("floor_cdr", m_floor);
}

void
CdrPlasmaGodunovStepper::parseDebug()
{
  CH_TIME("CdrPlasmaGodunovStepper::parseDebug()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::parseDebug()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("debug", m_debug);
}

void
CdrPlasmaGodunovStepper::parseProfile()
{
  CH_TIME("CdrPlasmaGodunovStepper::parseProfile()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::parseProfile()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("profile", m_profile);
}

void
CdrPlasmaGodunovStepper::parseFHD()
{
  CH_TIME("CdrPlasmaGodunovStepper::parseFHD()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::parseFHD()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("fhd", m_fhd);
}

void
CdrPlasmaGodunovStepper::parseRegridSlopes()
{
  CH_TIME("CdrPlasmaGodunovStepper::parseRegridSlopes()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::parseRegridSlopes()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("use_regrid_slopes", m_regridSlopes);
}

RefCountedPtr<CdrStorage>&
CdrPlasmaGodunovStepper::getCdrStorage(const CdrIterator<CdrSolver>& a_solverIt)
{
  return m_cdrScratch[a_solverIt.index()];
}

RefCountedPtr<RtStorage>&
CdrPlasmaGodunovStepper::getRtStorage(const RtIterator<RtSolver>& a_solverIt)
{
  return m_rteScratch[a_solverIt.index()];
}

Real
CdrPlasmaGodunovStepper::advance(const Real a_dt)
{
  CH_TIME("CdrPlasmaGodunovStepper::advance(Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::advance(Real)" << endl;
  }

  // INFO: When we enter here, the CdrSolvers should have been filled with velocities and diffusion coefficients.
  //
  // The next steps are:
  //
  //    1. Advance transport
  //    2. Solve the Poisson equation (special option for semi-implicit transport)
  //    3. Solve the reactive problem.
  //    4. Solve the radiative transfer problem
  //    5. Put data back into solvers to prepare for the next time step.

  this->allocateScratch();

  m_timer = std::make_unique<Timer>("CdrPlasmaGodunovStepper::advance");

  Timer timer("CdrPlasmaGodunovStepper::advance");

  // 1. Solve the transport problem. Note that we call advanceTransport which holds the implementation. This differs for explicit and semi-implicit formulations.
  CdrPlasmaGodunovStepper::advanceTransport(a_dt);

  // 2. Solve the Poisson equation and compute the electric field. If we did a semi-implicit solve then the field has already been computed.
  if (m_fieldCoupling != FieldCoupling::SemiImplicit) {
    m_timer->startEvent("Poisson");
    if ((m_timeStep + 1) % m_fastPoisson == 0) {
      CdrPlasmaStepper::solvePoisson();
    }
    CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch();
    m_timer->stopEvent("Poisson");
  }

  // 3. Solve the reactive problem.
  m_timer->startEvent("Reactions");
  CdrPlasmaGodunovStepper::computeCdrGradients();
  CdrPlasmaGodunovStepper::computeSourceTerms(a_dt);
  m_timer->stopEvent("Reactions");

  // 3. Advance CDR equations with reactive terms.
  m_timer->startEvent("CDR reactions");
  CdrPlasmaGodunovStepper::advanceCdrReactions(a_dt);
  m_timer->stopEvent("CDR reactions");

  // 4. Solve the radiative transfer problem.
  m_timer->startEvent("Radiation");
  if ((m_timeStep + 1) % m_fastRTE == 0) {
    CdrPlasmaGodunovStepper::advanceRadiativeTransfer(a_dt);
  }
  m_timer->stopEvent("Radiation");

  // Do post step operations.
  m_timer->startEvent("Post-step");
  CdrPlasmaGodunovStepper::postStep();
  m_timer->stopEvent("Post-step");

  // 5. Update velocities and diffusion coefficients in order to prepare for the next time step.
  m_timer->startEvent("Velocities");
  CdrPlasmaGodunovStepper::computeCdrDriftVelocities(m_time + a_dt);
  m_timer->stopEvent("Velocities");

  m_timer->startEvent("Diffu-coeffs");
  CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients(m_time + a_dt);
  m_timer->stopEvent("Diffu-coeffs");

  m_timer->startEvent("Compute J");
  this->computeJ(m_currentDensity);
  m_timer->stopEvent("Compute J");

  if (m_profile) {
    m_timer->eventReport(pout(), false);
  }

  this->deallocateScratch();

  return a_dt;
}

void
CdrPlasmaGodunovStepper::preRegrid(const int a_lbase, const int a_finestLevel)
{
  CH_TIME("CdrPlasmaGodunovStepper::preRegrid(int, int)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::preRegrid(int, int)" << endl;
  }

  // Call the parent method -- this will put all solvers in pre-regrid mode.
  CdrPlasmaStepper::preRegrid(a_lbase, a_finestLevel);

  // When we regrid we must compute the field on the new grid using the semi-implicit update. However, we only have
  // the conductivities and space charge on the old grids, so we must store them and interpolate them to the new grids.
  if (m_fieldCoupling == FieldCoupling::SemiImplicit) {
    m_amr->allocate(m_scratchConductivity, m_realm, phase::gas, 1);
    m_amr->allocate(m_scratchSemiImplicitRho, m_realm, phase::gas, 1);

    DataOps::copy(m_scratchConductivity, m_conductivityFactorCell);
    DataOps::copy(m_scratchSemiImplicitRho, m_semiImplicitRho);
  }

  m_semiImplicitRho.clear();
  m_conductivityFactorCell.clear();
  m_conductivityFactorFace.clear();
  m_conductivityFactorEB.clear();
}

void
CdrPlasmaGodunovStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("CdrPlasmaGodunovStepper::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::regrid(int, int, int)" << endl;
  }

  // TLDR: If we are not using a semi-implicit scheme then we can just call the parent method. Modifications are needed for the
  //       semi-implicit scheme because we solve for the field using div((eps + dt*sigma^k/eps0)E^(k+1)) = -rho^(k+1)/eps0 but
  //       this means we need the conductivity and space charge at the previous time step for restoring the field on the new mesh.

  // Just use regular regrid when we start the simulation.
  if (m_fieldCoupling != FieldCoupling::SemiImplicit || m_timeStep == 0) {
    CdrPlasmaStepper::regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
  else {
    const Interval interv(0, 0);

    // Call parent method for setting up storage and regridding solvers.
    this->allocateInternals();
    this->regridSolvers(a_lmin, a_oldFinestLevel, a_newFinestLevel);
    this->regridInternals(a_lmin, a_oldFinestLevel, a_newFinestLevel);

    const EBCoarseToFineInterp::Type interpType = m_regridSlopes ? EBCoarseToFineInterp::Type::ConservativeMinMod
                                                                 : EBCoarseToFineInterp::Type::ConservativePWC;

    m_amr->interpToNewGrids(m_conductivityFactorCell,
                            m_scratchConductivity,
                            phase::gas,
                            a_lmin,
                            a_oldFinestLevel,
                            a_newFinestLevel,
                            interpType);

    m_amr->interpToNewGrids(m_semiImplicitRho,
                            m_scratchSemiImplicitRho,
                            phase::gas,
                            a_lmin,
                            a_oldFinestLevel,
                            a_newFinestLevel,
                            interpType);

    // Coarsen the conductivity and space charge from the last step and update ghost cells.
    m_amr->arithmeticAverage(m_conductivityFactorCell, m_realm, m_phase);
    m_amr->interpGhostPwl(m_conductivityFactorCell, m_realm, m_phase);

    m_amr->arithmeticAverage(m_semiImplicitRho, m_realm, m_phase);
    m_amr->interpGhostPwl(m_semiImplicitRho, m_realm, m_phase);

    // Set up the semi-implicit Poisson equation and solve it.
    m_fieldSolver->setupSolver();
    this->computeFaceConductivity(m_conductivityFactorFace, m_conductivityFactorEB, m_conductivityFactorCell);
    this->setupSemiImplicitPoisson(m_conductivityFactorFace, m_conductivityFactorEB, 1.0);

    const bool converged = this->solveSemiImplicitPoisson();

    // Debug if we don't converge.
    if (!converged) {
      pout() << "CdrPlasmaGodunovStepper::regrid - Poisson solver failed to converge during semi-implicit regrid."
             << endl;

#if 1 // For now, add a debug file.
#ifdef CH_USE_HDF5
      pout() << "CdrPlasmaGodunovStepper::regrid - I'm adding a debug file in 'semi_implicit_debug.hdf5'" << endl;
      EBAMRCellData data;
      m_amr->allocate(data, m_realm, m_phase, 2);

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        m_conductivityFactorCell[lvl]->copyTo(Interval(0, 0), *data[lvl], Interval(0, 0));
        m_semiImplicitRho[lvl]->copyTo(Interval(0, 0), *data[lvl], Interval(1, 1));
      }

      DischargeIO::writeEBHDF5(data, "semi_implicit_debug.hdf5");
#endif
#endif
    }

    // Now compute drift velocities and diffusion -- using the electric field from last time step on the new mesh.
    CdrPlasmaStepper::computeCdrDriftVelocities();
    CdrPlasmaStepper::computeCdrDiffusion();

    // If we're doing a stationary RTE solve, we also have to recompute source terms
    if (this->stationaryRTE()) { // Solve RTE equations by using data that exists inside solvers
      const Real dummyDt = 0.0;

      // Need new source terms for RTE equations
      this->advanceReactionNetwork(m_time, dummyDt);
      this->solveRadiativeTransfer(dummyDt); // Argument does not matter, it's a stationary solver.
    }
  }
}

void
CdrPlasmaGodunovStepper::postRegrid()
{
  CH_TIME("CdrPlasmaGodunovStepper::postRegrid()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::postRegrid()" << endl;
  }

  CdrPlasmaStepper::postRegrid();

  // Release memory.
  if (m_fieldCoupling == FieldCoupling::SemiImplicit) {
    m_amr->deallocate(m_scratchConductivity);
    m_amr->deallocate(m_scratchSemiImplicitRho);
  }
}

void
CdrPlasmaGodunovStepper::postCheckpointSetup()
{
  CH_TIME("CdrPlasmaGodunovStepper::postCheckpointSetup()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::postCheckpointSetup()" << endl;
  }

  // TLDR: Only the semi-implicit part overrides the parent method, and it is done because the field needs to be computed from
  //       a different equation.

  if (m_fieldCoupling == FieldCoupling::SemiImplicit) {

    // When we enter this routine we will already have called read the checkpoint data into the conductivityFactor and semiimplicit space charge. We need
    // to set up the field solver with those quantities rather than the regular space charge.
    m_fieldSolver->setupSolver();

    m_amr->arithmeticAverage(m_conductivityFactorCell, m_realm, m_phase);
    m_amr->interpGhostPwl(m_conductivityFactorCell, m_realm, m_phase);

    m_amr->arithmeticAverage(m_semiImplicitRho, m_realm, m_phase);
    m_amr->interpGhostPwl(m_semiImplicitRho, m_realm, m_phase);

    this->computeFaceConductivity(m_conductivityFactorFace, m_conductivityFactorEB, m_conductivityFactorCell);
    this->setupSemiImplicitPoisson(m_conductivityFactorFace, m_conductivityFactorEB, 1.0);

    // Now solve for the field on the new grids.
    const bool converged = this->solveSemiImplicitPoisson();
    if (!converged) {
      pout() << "CdrPlasmaGodunovStepper::postCheckpointSetup - Poisson solver failed to converge during restart."
             << endl;
    }

    // Now compute the drift and diffusion velocities and update the source terms. This will be OK because the
    // CdrPlasmaStepper functions fetch the electric field from the solver, but it already has the correct potential.
    CdrPlasmaStepper::computeCdrDriftVelocities();
    CdrPlasmaStepper::computeCdrDiffusion();

    // If we are doing a stationary RTE then we probably have to update the elliptic equations.
    if (this->stationaryRTE()) {
      const Real dummyDt = 0.0;

      // We also need new source terms for RTE equations.
      this->advanceReactionNetwork(m_time, dummyDt);
      this->solveRadiativeTransfer(dummyDt);
    }
  }
  else {
    CdrPlasmaStepper::postCheckpointSetup();
  }
}

void
CdrPlasmaGodunovStepper::regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("CdrPlasmaGodunovStepper::regridInternals(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::regridInternals(int, int, int)" << endl;
  }
}

bool
CdrPlasmaGodunovStepper::solveSemiImplicitPoisson()
{
  CH_TIME("CdrPlasmaGodunovStepper::solveSemiImplicitPoisson");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::solveSemiImplicitPoisson" << endl;
  }

  // Now set the semi-implicit space charge.
  MFAMRCellData& rho      = m_fieldSolver->getRho();
  EBAMRCellData  rhoPhase = m_amr->alias(m_phase, rho);

  DataOps::setValue(rho, 0.0);
  DataOps::copy(rhoPhase, m_semiImplicitRho);

  m_amr->arithmeticAverage(rho, m_realm);
  m_amr->interpGhostPwl(rho, m_realm);

  m_amr->interpToCentroids(rhoPhase, m_realm, m_phase);

  const bool converged = m_fieldSolver->solve(m_fieldSolver->getPotential(), rho, m_sigma->getPhi(), false);

  return converged;
}

void
CdrPlasmaGodunovStepper::allocateInternals()
{
  CH_TIME("CdrPlasmaGodunovStepper::allocateInternals()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::allocateInternals()" << endl;
  }

  constexpr int nComp = 1;

  CdrPlasmaStepper::allocateInternals();

  // Set up storage for semi-implicit field solves.
  m_amr->allocate(m_semiImplicitRho, m_realm, m_phase, nComp);
  m_amr->allocate(m_conductivityFactorCell, m_realm, m_phase, nComp);
  m_amr->allocate(m_conductivityFactorFace, m_realm, m_phase, nComp);
  m_amr->allocate(m_conductivityFactorEB, m_realm, m_phase, nComp);

  // Initialize values.
  DataOps::setValue(m_semiImplicitRho, 0.0);
  DataOps::setValue(m_conductivityFactorCell, 0.0);
  DataOps::setValue(m_conductivityFactorFace, 0.0);
  DataOps::setValue(m_conductivityFactorEB, 0.0);
}

void
CdrPlasmaGodunovStepper::allocateScratch()
{
  CH_TIME("CdrPlasmaGodunovStepper::allocateScratch()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::allocateScratch()" << endl;
  }

  // Number of CDR solvers and RTE solvers that we have.
  const int numCdrSpecies = m_physics->getNumCdrSpecies();
  const int numRteSpecies = m_physics->getNumRtSpecies();

  // Allocate storage for the CDR solvers.
  m_cdrScratch.resize(numCdrSpecies);
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    m_cdrScratch[idx] = RefCountedPtr<CdrStorage>(new CdrStorage(m_amr, m_realm, m_cdr->getPhase()));
    m_cdrScratch[idx]->allocateStorage();
  }

  // Allocate RTE storage.
  m_rteScratch.resize(numRteSpecies);
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    m_rteScratch[idx] = RefCountedPtr<RtStorage>(new RtStorage(m_amr, m_realm, m_rte->getPhase()));
    m_rteScratch[idx]->allocateStorage();
  }

  // Allocate storage for field solver.
  m_fieldScratch = RefCountedPtr<FieldStorage>(new FieldStorage(m_amr, m_realm, m_cdr->getPhase()));
  m_fieldScratch->allocateStorage();

  // Allocate storage for surface charge solver.
  m_sigmaScratch = RefCountedPtr<SigmaStorage>(new SigmaStorage(m_amr, m_realm, m_cdr->getPhase()));
  m_sigmaScratch->allocateStorage();
}

void
CdrPlasmaGodunovStepper::deallocateInternals()
{
  CH_TIME("CdrPlasmaGodunovStepper::deallocateInternals()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::deallocateInternals()" << endl;
  }

  // TLDR: This routine simply deallocates the transient memory used by CdrPlasmaGodunovStepper.

  // Run through CDR solvers and deallocate the transient memory assocaited with them.
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    m_cdrScratch[idx]->deallocateStorage();
    m_cdrScratch[idx] = RefCountedPtr<CdrStorage>(0);
  }

  // Run through RTE solvers and deallocate the transient memory assocaited with them.
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    m_rteScratch[idx]->deallocateStorage();
    m_rteScratch[idx] = RefCountedPtr<RtStorage>(0);
  }

  m_cdrScratch.resize(0);
  m_rteScratch.resize(0);

  m_fieldScratch->deallocateStorage();
  m_fieldScratch = RefCountedPtr<FieldStorage>(0);

  m_sigmaScratch->deallocateStorage();
  m_sigmaScratch = RefCountedPtr<SigmaStorage>(0);

  // Deallocate storage for the semi-implicit solve.
  m_amr->deallocate(m_semiImplicitRho);
  m_amr->deallocate(m_conductivityFactorCell);
  m_amr->deallocate(m_conductivityFactorFace);
  m_amr->deallocate(m_conductivityFactorEB);
}

void
CdrPlasmaGodunovStepper::deallocateScratch()
{
  CH_TIME("CdrPlasmaGodunovStepper::deallocateScratch()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::deallocateScratch()" << endl;
  }

  // TLDR: This routine simply deallocates the transient memory used by CdrPlasmaGodunovStepper.

  // Run through CDR solvers and deallocate the transient memory assocaited with them.
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    m_cdrScratch[idx]->deallocateStorage();
    m_cdrScratch[idx] = RefCountedPtr<CdrStorage>(0);
  }

  // Run through RTE solvers and deallocate the transient memory assocaited with them.
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

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
CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch()
{
  CH_TIME("CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch()" << endl;
  }

  // TLDR: This computes the electric field on the cell center, EB, and domain faces.

  EBAMRCellData& electricFieldCell   = m_fieldScratch->getElectricFieldCell();
  EBAMRIVData&   electricFieldEB     = m_fieldScratch->getElectricFieldEB();
  EBAMRIFData&   electricFieldDomain = m_fieldScratch->getElectricFieldDomain();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();

  CdrPlasmaStepper::computeElectricField(electricFieldCell, m_cdr->getPhase(), phi); // Compute cell-centered field
  CdrPlasmaStepper::computeElectricField(electricFieldEB, m_cdr->getPhase(), electricFieldCell); // EB-centered field
  CdrPlasmaStepper::extrapolateToDomainFaces(electricFieldDomain,
                                             m_cdr->getPhase(),
                                             electricFieldCell); // Domain centered field
}

void
CdrPlasmaGodunovStepper::computeCdrGradients()
{
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrGradients()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::computeCdrGradients()" << endl;
  }

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    // Fetch solver and associated storage.
    RefCountedPtr<CdrSolver>&  solver  = solverIt();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

    EBAMRCellData& scratch = storage->getScratch();
    EBAMRCellData& grad    = storage->getGradient();

    // Update the ghost cells so we can compute the gradient.
    m_amr->copyData(scratch,solver->getPhi());

    m_amr->arithmeticAverage(scratch, m_realm, m_phase);
    m_amr->interpGhostPwl(scratch, m_realm, m_phase);

    // Compute the gradient, coarsen it, and update the ghost cells.
    m_amr->computeGradient(grad, scratch, m_realm, phase::gas);

    m_amr->arithmeticAverage(grad, m_realm, m_cdr->getPhase());
    m_amr->interpGhost(grad, m_realm, m_cdr->getPhase());
  }
}

void
CdrPlasmaGodunovStepper::extrapolateWithSourceTerm(const Real a_dt)
{
  CH_TIME("CdrPlasmaGodunovStepper::extrapolateWithSourceTerm(Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::extrapolateWithSourceTerm(Real)" << endl;
  }

  // TLDR: If we extrapolate the advective derivative and include the source term in the extrapolation,
  //       the boundary conditions should be computed from phi + 0.5*source*dt rather than just phi. This probably
  //       does not matter, but for consistency we do it anyways.

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrSolver>&  solver  = solverIt();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

    const EBAMRCellData& state = solver->getPhi();
    //    const EBAMRCellData& source = solver->getSource();

    EBAMRCellData& extrap = storage->getExtrap();

    // Let n = n(t) + 0.5*a_dt * S(t) if we use the source term for extrapolation.
    DataOps::copy(extrap, state);
    if (m_extrapAdvect) {
      //      DataOps::incr(extrap, source, 0.5*a_dt);
    }
  }
}

void
CdrPlasmaGodunovStepper::extrapolateCdrToEB()
{
  CH_TIME("CdrPlasmaGodunovStepper::extrapolateCdrToEB()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::extrapolateCdrToEB()" << endl;
  }

  // TLDR: This routine is reponsible for computing the cell-centered states and gradients at the EB. This is necessary because
  //       the boundary condition routines require these things to be known at the EB. This is the routine that computes them. We
  //       will later fetch these quantities and pass them into our boundary condition routines.

  Vector<EBAMRCellData*>
                         cdrDensities; // Cell-centered densities used when extrapolating to EBs and domains when parsing boundary conditions.
  Vector<EBAMRCellData*> cdrGradients;   // Gradient of cell-centered densities
  Vector<EBAMRIVData*>   cdrDensitiesEB; // Extrapolation of cdrDensities to the EB
  Vector<EBAMRIVData*>   cdrGradientsEB; // Extrapolation of cdrGradients to the EB

  // Scratch storage for holding a gradient at the EB
  EBAMRIVData gradientEB;
  m_amr->allocate(gradientEB, m_realm, m_cdr->getPhase(), SpaceDim);

  // Run through the solvers and fetch the various data used for the extrapolation.
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

    // Note: For the CDR densities we use the extrap data holder in the scratch storage. The routine extrapolateWithSourceTerm will have
    //       been called prior to this routine. If we center advective discretizations at the half time step, we need to increment the
    //       edge centered states by 0.5*S*dt.

    // Populate the data.
    cdrDensities.push_back(&(storage->getExtrap()));
    cdrGradients.push_back(&(storage->getGradient()));
    cdrDensitiesEB.push_back(&(storage->getEbState()));
    cdrGradientsEB.push_back(&(storage->getEbGrad()));
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because no mass has been injected. That comes later.
  CdrPlasmaStepper::extrapolateToEb(cdrDensitiesEB, m_cdr->getPhase(), cdrDensities);
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    DataOps::floor(*cdrDensitiesEB[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux.
  for (int i = 0; i < cdrDensities.size(); i++) {

    // Extrapolate the EB to the gradient
    CdrPlasmaStepper::extrapolateToEb(gradientEB, m_cdr->getPhase(), *cdrGradients[i]);

    // And project it along the EB normal.
    CdrPlasmaStepper::projectFlux(*cdrGradientsEB[i], gradientEB);
  }
}

void
CdrPlasmaGodunovStepper::computeCdrFluxesEB()
{
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrFluxesEB()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::computeCdrFluxesEB()";
  }

  // TLDR: This is the main routine for computing the CDR fluxes on the EB, i.e. the boundary conditions. When we enter this routine
  //       we will have populated the states we use for extrapolation, and the gradients. The velocities on the EB and the extrapoalted
  //       fluxes on the EB will not have been computed, so we do those here.

  // Holds the CDR densities.
  Vector<EBAMRCellData*> cdrDensities;

  // These are things that must be populated and passed into our
  // nifty boundary condition framework.
  Vector<EBAMRIVData*> extrapCdrFluxesEB;
  Vector<EBAMRIVData*> extrapCdrDensitiesEB;
  Vector<EBAMRIVData*> extrapCdrVelocitiesEB;
  Vector<EBAMRIVData*> extrapCdrGradientsEB;
  Vector<EBAMRIVData*> extrapRteFluxesEB;

  // Run through CDR solvers and populate Vectors.
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrStorage>& storage = this->getCdrStorage(solverIt);

    EBAMRCellData& densityCell = storage->getExtrap();
    EBAMRIVData&   densityEB   = storage->getEbState();
    EBAMRIVData&   velocityEB  = storage->getEbVelo();
    EBAMRIVData&   fluxEB      = storage->getEbFlux();
    EBAMRIVData&   gradientEB  = storage->getEbGrad();

    cdrDensities.push_back(&densityCell);         // Will not touch this, but use it for extrapolating the fluxes.
    extrapCdrDensitiesEB.push_back(&densityEB);   // Computed in extrapolateCdrToEB
    extrapCdrVelocitiesEB.push_back(&velocityEB); // Not yet computed
    extrapCdrFluxesEB.push_back(&fluxEB);         // Not yet computed
    extrapCdrGradientsEB.push_back(&gradientEB);  // Computed in extrapolateCdrToEB
  }

  // Extrapolate the CDR fluxes and velocities to the EB. After this, we have all
  // the pertinent CDR quantities we need in our BC framework.
  Vector<EBAMRCellData*> cdrVelocities = m_cdr->getVelocities();
  CdrPlasmaStepper::computeExtrapolatedFluxes(extrapCdrFluxesEB, cdrDensities, cdrVelocities, m_cdr->getPhase());
  CdrPlasmaStepper::computeExtrapolatedVelocities(extrapCdrVelocitiesEB, cdrVelocities, m_cdr->getPhase());

  // Compute RTE flux on the boundary
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<RtSolver>&  solver  = solverIt();
    RefCountedPtr<RtStorage>& storage = this->getRtStorage(solverIt);

    EBAMRIVData& fluxEB = storage->getEbFlux();

    // Let the solver compute the boundary flux.
    solver->computeBoundaryFlux(fluxEB, solver->getPhi());

    // Add the extrapolated EB flux to the data holder.
    extrapRteFluxesEB.push_back(&fluxEB);
  }

  // This is where we put the result -- directly in the solvers.
  Vector<EBAMRIVData*> cdrFluxesEB = m_cdr->getEbFlux();

  // Electric field which has been extrapolated to the EB (in computeElectricFieldIntoScratch).
  const EBAMRIVData& electricFieldEB = m_fieldScratch->getElectricFieldEB();

  // Now call the parent method which does all the BC computations.
  CdrPlasmaStepper::computeCdrFluxes(cdrFluxesEB,
                                     extrapCdrFluxesEB,
                                     extrapCdrDensitiesEB,
                                     extrapCdrVelocitiesEB,
                                     extrapCdrGradientsEB,
                                     extrapRteFluxesEB,
                                     electricFieldEB,
                                     m_time);
}

void
CdrPlasmaGodunovStepper::extrapolateCdrToDomain()
{
  CH_TIME("CdrPlasmaGodunovStepper::extrapolateCdrToDomain()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::extrapolateCdrToDomain()" << endl;
  }

  // TLDR: This routine is reponsible for computing the cell-centered states and gradients at domainfaces. This is necessary because
  //       the boundary condition routines require these things to be known. This is the routine that computes them. We
  //       will later fetch these quantities and pass them into our boundary condition routines.

  Vector<EBAMRCellData*> cdrDensities;       // CDR densities on the cell center
  Vector<EBAMRCellData*> cdrGradients;       // CDR gradients on the cell center
  Vector<EBAMRIFData*>   cdrDensitiesDomain; // CDR densities on the domain faces
  Vector<EBAMRIFData*>   cdrGradientsDomain; // CDR gradients on the domain faces

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux.
  EBAMRIFData grad;
  m_amr->allocate(grad, m_realm, m_cdr->getPhase(), SpaceDim);

  // Run through the CDR solvers and populate the vectors.
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

    cdrDensities.push_back(&(storage->getExtrap()));   // Already known and computed in extrapolateWithSourceTerm
    cdrGradients.push_back(&(storage->getGradient())); // Should already be computed in computeGradients
    cdrDensitiesDomain.push_back(&(storage->getDomainState()));
    cdrGradientsDomain.push_back(&(storage->getDomainGrad()));
  }

  // Extrapolate the cell-centered states to the domain faces
  CdrPlasmaStepper::extrapolateToDomainFaces(cdrDensitiesDomain, m_cdr->getPhase(), cdrDensities);

  // Run through the solvers and extrapolate the gradients to the domain faces.
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    // Extrapolate gradient.
    CdrPlasmaStepper::extrapolateToDomainFaces(grad, m_cdr->getPhase(), *cdrGradients[idx]);

    // Project it onto the domain face.
    CdrPlasmaStepper::projectDomain(*cdrGradientsDomain[idx], grad);
  }
}

void
CdrPlasmaGodunovStepper::computeCdrDomainFluxes()
{
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrDomainFluxes()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::computeCdrDomainFluxes()" << endl;
  }

  // TLDR: This is the main routine for computing the CDR fluxes on the domain faces, i.e. the boundary conditions. When we enter this routine
  //       we will have populated the states we use for extrapolation, and the gradients. The velocities on the domain faces and the extrapolated
  //       fluxes on the domain faces will not have been computed, so we do those here.

  Vector<EBAMRCellData*> cdrDensities; // For holding the cell-centered states used for extrapolation.
  Vector<EBAMRCellData*> cdrGradients; // For holding the cell-centered gradients.

  Vector<EBAMRIFData*> extrapCdrFluxesDomain;     // Extrapolated fluxes to domain faces
  Vector<EBAMRIFData*> extrapCdrDensitiesDomain;  // Extrapolated densities to domain faces
  Vector<EBAMRIFData*> extrapCdrVelocitiesDomain; // Extrapolated velocities to domain faces
  Vector<EBAMRIFData*> extrapCdrGradientsDomain;  // Extrapolated gradients to domain faces
  Vector<EBAMRIFData*> extrapRteFluxesDomain;     // Extrapolated RTE fluxes to domain faces

  // CDR velocities -- these are known.
  Vector<EBAMRCellData*> cdrVelocities = m_cdr->getVelocities();

  // Run through the CDR solvers and populate the relevant vectors.
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrStorage>& storage = this->getCdrStorage(solverIt);

    EBAMRIFData&   densityDomain  = storage->getDomainState(); // Not yet computed.
    EBAMRIFData&   velocityDomain = storage->getDomainVelo();  // Not yet computed.
    EBAMRIFData&   fluxDomain     = storage->getDomainFlux();  // Not yet computed.
    EBAMRIFData&   gradDomain     = storage->getDomainGrad();  // Already computed.
    EBAMRCellData& gradCell       = storage->getGradient();    // Already computed.
    EBAMRCellData& densityCell    = storage->getExtrap();      // Already computed.

    extrapCdrDensitiesDomain.push_back(&densityDomain);   // Not yet computed.
    extrapCdrVelocitiesDomain.push_back(&velocityDomain); // Not yet computed.
    extrapCdrFluxesDomain.push_back(&fluxDomain);         // Not yet computed.
    extrapCdrGradientsDomain.push_back(&gradDomain);      // Already computed.
    cdrGradients.push_back(&gradCell);                    // Already computed.
    cdrDensities.push_back(&densityCell);                 // Already computed.
  }

  // The API says we must have the extrapolated fluxes, densities, velocities, and gradients on the domain faces. We've already
  // done the densities and gradients in extrapolateCdrToDomain, but we have not yet done the velocities and fluxes.
  // this->extrapolateToDomainFaces(extrapCdrDensities,         m_cdr->getPhase(), states);
  this->extrapolateVelocitiesToDomainFaces(extrapCdrVelocitiesDomain, m_cdr->getPhase(), cdrVelocities);
  this->computeExtrapolatedDomainFluxes(extrapCdrFluxesDomain, cdrDensities, cdrVelocities, m_cdr->getPhase());
  // this->extrapolateVectorToDomainFaces(extrapCdrGradients,  m_cdr->getPhase(), cdrGradients);

  // Now compute the RTE flux on domain faces as well.
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<RtSolver>&  solver  = solverIt();
    RefCountedPtr<RtStorage>& storage = this->getRtStorage(solverIt);

    EBAMRIFData& rteFluxDomain = storage->getDomainFlux();

    // Solver computes the domain fluxes and puts it in the data holder.
    solver->computeDomainFlux(rteFluxDomain, solver->getPhi());

    extrapRteFluxesDomain.push_back(&rteFluxDomain);
  }

  // This where we put the fluxes that we compute. They go directly in the solvers so that the solvers can compute divergences.
  Vector<EBAMRIFData*> cdrFluxesDomain = m_cdr->getDomainFlux();

  // Electric field on the domain edge/face.
  const EBAMRIFData& electricFieldDomain = m_fieldScratch->getElectricFieldDomain();

  // We have prepared everything that the API needs. Now call the parent method which fills the solvers' domain fluxes
  CdrPlasmaStepper::computeCdrDomainFluxes(cdrFluxesDomain,
                                           extrapCdrFluxesDomain,
                                           extrapCdrDensitiesDomain,
                                           extrapCdrVelocitiesDomain,
                                           extrapCdrGradientsDomain,
                                           extrapRteFluxesDomain,
                                           electricFieldDomain,
                                           m_time);
}

void
CdrPlasmaGodunovStepper::computeSigmaFlux()
{
  CH_TIME("CdrPlasmaGodunovStepper::computeSigmaFlux()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::computeSigmaFlux()" << endl;
  }

  // Reset the data holder that holds the solver charge flux.
  EBAMRIVData& flux = m_sigma->getRHS();
  DataOps::setValue(flux, 0.0);

  // Run through the CDR solvers and increment by their fluxes.
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<CdrSolver>&  solver  = solverIt();
    const RefCountedPtr<CdrSpecies>& species = solverIt.getSpecies();

    const int Z = species->getChargeNumber();

    if (Z != 0) {
      const EBAMRIVData& cdrSolverFluxEB = solver->getEbFlux();

      DataOps::incr(flux, cdrSolverFluxEB, Z * Units::Qe);
    }
  }

  // Reset the flux on electrode interface cells.
  m_sigma->resetElectrodes(flux, 0.0);
}

void
CdrPlasmaGodunovStepper::advanceTransport(const Real a_dt)
{
  CH_TIME("CdrPlasmaGodunovStepper::advanceTransport(Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::advanceTransport(Real)" << endl;
  }

  switch (m_fieldCoupling) {
  case FieldCoupling::Explicit: {
    this->advanceTransportExplicitField(a_dt);

    break;
  }
  case FieldCoupling::SemiImplicit: {
    this->advanceTransportSemiImplicit(a_dt);

    break;
  }
  default: {
    MayDay::Error("CdrPlasmaGodunovStepper::advanceTransport - logic bust. Must have ");

    break;
  }
  }
}

void
CdrPlasmaGodunovStepper::advanceTransportExplicitField(const Real a_dt)
{
  CH_TIME("CdrPlasmaGodunovStepper::advanceTransportExplicitField(Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::advanceTransportExplicitField(Real)" << endl;
  }

  // TLDR: This advances the CDR equations using an Euler rule. The right-hand side of the CDR equations can be explictly discretized, or
  //       with implicit diffusion. If we use implicit diffusion we first advance the advective problem to the end state and use that as
  //       an initial condition in the diffusion equation.

  // First, update everything we need for consistently computing boundary conditions on the EBs and
  // domain faces.
  m_timer->startEvent("Gradients and BCs");
  CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch(); // Compute the electric field
  CdrPlasmaGodunovStepper::computeCdrGradients();             // Compute CDR gradients
  CdrPlasmaGodunovStepper::extrapolateWithSourceTerm(a_dt);   // If we used advective extrapolation, BCs are more work.
  CdrPlasmaGodunovStepper::extrapolateCdrToEB();              // Extrapolate cell-centered stuff to EB centroids
  CdrPlasmaGodunovStepper::extrapolateCdrToDomain();          // Extrapolate cell-centered states to domain edges
  CdrPlasmaGodunovStepper::computeCdrFluxesEB();              // Extrapolate cell-centered fluxes to EB centroids
  CdrPlasmaGodunovStepper::computeCdrDomainFluxes();          // Extrapolate cell-centered fluxes to domain edges
  CdrPlasmaGodunovStepper::computeSigmaFlux();                // Update charge flux for sigma solver
  m_timer->stopEvent("Gradients and BCs");

  // Run through the CDR solvers and update them as
  //
  //    phi^(k+1) = phi^k - dt*div(J).
  //
  // If we use implicit diffusion, then we actually solve
  //
  // phi^(k+1) = phi^k - dt*div(F) + dt*div(

  m_timer->startEvent("Transport advance");
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    RefCountedPtr<CdrSolver>&  solver  = solverIt();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

    // Get CDR density and some scratch storage that we can use.
    EBAMRCellData& phi      = solver->getPhi();
    EBAMRCellData& scratch  = storage->getScratch();
    EBAMRCellData& scratch2 = storage->getScratch2();

    // Do the advective advance.
    if (solver->isMobile()) {
      switch (m_advectionSolver) {
      case AdvectionSolver::Euler: {
        // Advance is just phi^(k+1) = phi^k - dt*div(F)
        solver->computeDivF(scratch, phi, 0.0, false, true, true);

        DataOps::incr(phi, scratch, -a_dt);

        break;
      }
      case AdvectionSolver::RK2: {
        // Use Heun's method.
        solver->computeDivF(scratch, phi, 0.0, false, true, true);

        // Make phi = phi^k  - dt * div(F)
        DataOps::incr(phi, scratch, -a_dt);

        // Compute divF(phi)
        solver->computeDivF(scratch2, phi, 0.0, false, true, true);

        // Make phi = phi^k + 0.5*dt * (scratch1 + scratch2)
        DataOps::incr(phi, scratch, 0.5 * a_dt);
        DataOps::incr(phi, scratch2, -0.5 * a_dt);

        break;
      }
      case AdvectionSolver::MUSCL: {
        // Advance is phi^(k+1) = phi^k - dt*div(F). Almost like the Euler method except for the transverse slopes.
        solver->computeDivF(scratch, phi, a_dt, false, true, true);

        DataOps::incr(phi, scratch, -a_dt);

        break;
      }
      default: {
        MayDay::Error("CdrPlasmaGodunovStepper::advanceTransportExplicitField -- logic bust");
      }
      }

      // Floor mass or not?
      if (m_floor) {
        this->floorMass(phi, "CdrPlasmaGodunovStepper::advanceTransportExplicitField", solver);
      }

      // Coarsen the solution and update ghost cells.
      m_amr->arithmeticAverage(phi, m_realm, m_cdr->getPhase());
      m_amr->interpGhost(phi, m_realm, m_cdr->getPhase());
    }

    // Do the diffusion advance. This can be explicit or implicit.
    if (solver->isDiffusive()) {
      if (m_useImplicitDiffusion[idx]) {
        DataOps::copy(scratch, phi);
        DataOps::setValue(scratch2, 0.0);

        if (m_diffusionOrder == 1) {
          solver->advanceEuler(phi, scratch, scratch2, a_dt);
        }
        else if (m_diffusionOrder == 2) {
          solver->advanceCrankNicholson(phi, scratch, scratch2, a_dt);
        }
      }
      else {
        if (m_diffusionOrder == 1) {
          solver->computeDivD(scratch, phi, false, false, false);

          DataOps::incr(phi, scratch, a_dt);
        }
        else if (m_diffusionOrder == 2) {
          solver->computeDivD(scratch, phi, false, false, false);

          DataOps::incr(phi, scratch, a_dt);

          solver->computeDivD(scratch2, phi, false, false, false);

          DataOps::incr(phi, scratch, -0.5 * a_dt);
          DataOps::incr(phi, scratch2, 0.5 * a_dt);
        }
      }
    }

    // Floor mass or not?
    if (m_floor) {
      this->floorMass(phi, "CdrPlasmaGodunovStepper::advanceTransportExplicitField", solver);
    }

    // Coarsen the solution and update ghost cells.
    m_amr->arithmeticAverage(phi, m_realm, m_cdr->getPhase());
    m_amr->interpGhost(phi, m_realm, m_cdr->getPhase());
  }
  m_timer->stopEvent("Transport advance");

  // Advance the sigma equation. This may seem weird but we've kept the flux through the EB constant during the transport step, so it
  // doesn't matter if we did an Euler or Heun advance in the advective step.
  EBAMRIVData&       sigma = m_sigma->getPhi();
  const EBAMRIVData& rhs   = m_sigma->getRHS();

  DataOps::incr(sigma, rhs, a_dt);
}

void
CdrPlasmaGodunovStepper::advanceTransportSemiImplicit(const Real a_dt)
{
  CH_TIME("CdrPlasmaGodunovStepper::advanceTransportSemiImplicit(Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::advanceTransportSemiImplicit(Real)" << endl;
  }

  // Compute the conductivity first. We store it as sigma^k*a_dt/eps0
  m_timer->startEvent("Compute conductivity");
  this->computeCellConductivity(m_conductivityFactorCell);
  DataOps::scale(m_conductivityFactorCell, a_dt / Units::eps0);
  DataOps::floor(m_conductivityFactorCell, 0.0);

  m_amr->arithmeticAverage(m_conductivityFactorCell, m_realm, m_phase);
  m_amr->interpGhostPwl(m_conductivityFactorCell, m_realm, m_phase);

  // Average conductivity to faces and set up the semi-implicit poisson equation.
  this->computeFaceConductivity(m_conductivityFactorFace, m_conductivityFactorEB, m_conductivityFactorCell);
  m_timer->stopEvent("Compute conductivity");

  m_timer->startEvent("Setup Poisson");
  this->setupSemiImplicitPoisson(m_conductivityFactorFace, m_conductivityFactorEB, 1.0);
  m_timer->stopEvent("Setup Poisson");

  // Compute the modified right-hand side. We store this as rho^k - dt*e * sum(Z * div(D*grad(phi))).
  m_timer->startEvent("Compute rho");
  DataOps::setValue(m_semiImplicitRho, 0.0);
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<CdrSolver>&  solver  = solverIt();
    const RefCountedPtr<CdrSpecies>& species = solverIt.getSpecies();

    const int Z = species->getChargeNumber();

    if (Z != 0) {
      EBAMRCellData& phi = solver->getPhi();

      DataOps::incr(m_semiImplicitRho, phi, Z * Units::Qe);

      // If the solver is diffusive we must compute the diffusion term as well, and then increment by it.
      if (solver->isDiffusive()) {
#if 0
	RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);
												 
	EBAMRCellData& divDgradPhi = storage->getScratch();

	solver->computeDivD(divDgradPhi, phi, false, false, false);

	DataOps::incr(m_semiImplicitRho, divDgradPhi, Z*a_dt*Units::Qe);
#endif
      }
    }
  }
  m_timer->stopEvent("Compute rho");

  // Now solve the semi-implicit Poisson. Issue a warning if we didn't converge.
  m_timer->startEvent("Poisson solve");
  const bool converged = this->solveSemiImplicitPoisson();
  if (!converged) {
    pout() << "CdrPlasmaGodunovStepper::advanceTransportSemiImplicit --  Poisson solve did not converge!" << endl;
  }
  m_timer->stopEvent("Poisson solve");

  // Compute the electric field into scratch storage.
  m_timer->startEvent("Field and velocities");
  CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch();

  // Recompute velocities and diffusion coefficient using the electric field after the semi-implicit field solve.
  CdrPlasmaGodunovStepper::computeCdrDriftVelocities(m_time + a_dt);
  CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients(m_time + a_dt);
  m_timer->stopEvent("Field and velocities");

  // Now call the Euler transport method -- it will know what to do.
  this->advanceTransportExplicitField(a_dt);
}

void
CdrPlasmaGodunovStepper::advanceCdrReactions(const Real a_dt)
{
  CH_TIME("CdrPlasmaGodunovStepper::advanceCdrReactions(Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::advanceCdrReactions(Real)" << endl;
  }

  // After calling advanceReactionNetwork the CDR solvers have been filled with appropriate source terms. We now advance the
  // states over the time step a_dt using those source terms.
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrSolver>& solver = solverIt();

    EBAMRCellData&       phi = solver->getPhi();
    const EBAMRCellData& src = solver->getSource();

    DataOps::incr(phi, src, a_dt);

    // Floor mass if asked for it. If running in debug mode we compute the mass before and after flooring it.
    if (m_floor) {
      if (m_debug) {
        const Real massBefore = solver->computeMass();

        DataOps::floor(phi, 0.0);

        const Real massAfter   = solver->computeMass();
        const Real relMassDiff = (massAfter - massBefore) / massBefore;

        pout() << "CdrPlasmaGodunovStepper::advanceCdrReactions - injecting relative " << solver->getName()
               << " mass = " << relMassDiff << endl;
      }
      else {
        DataOps::floor(phi, 0.0);
      }
    }
  }
}

void
CdrPlasmaGodunovStepper::advanceRadiativeTransfer(const Real a_dt)
{
  CH_TIME("CdrPlasmaGodunovStepper::advanceRadiativeTransfer(Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::advanceRadiativeTransfer(Real)" << endl;
  }

  // Source terms should already be in place so solvers can just run their advance method.
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->advance(a_dt);
  }
}

void
CdrPlasmaGodunovStepper::computeCdrDriftVelocities(const Real a_time)
{
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrDriftVelocities(Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::computeCdrDriftVelocities(Real)" << endl;
  }

  // TLDR: We call the parent method, but using the scratch storage that holds the electric field.

  Vector<EBAMRCellData*> velocities = m_cdr->getVelocities();
  CdrPlasmaStepper::computeCdrDriftVelocities(velocities,
                                              m_cdr->getPhis(),
                                              m_fieldScratch->getElectricFieldCell(),
                                              a_time);
}

void
CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients(const Real a_time)
{
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients(Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients(Real)" << endl;
  }

  // TLDR: We call the parent method, but using the scratch storage that holds the electric field.

  CdrPlasmaStepper::computeCdrDiffusion(m_fieldScratch->getElectricFieldCell(), m_fieldScratch->getElectricFieldEB());
}

void
CdrPlasmaGodunovStepper::computeSourceTerms(const Real a_dt)
{
  CH_TIME("CdrPlasmaGodunovStepper::computeSourceTerms(Real)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::computeSourceTerms(Real)" << endl;
  }

  // We have already computed E and the gradients of the CDR equations. The API says we also need the gradients of the CDR equations, but we've already
  // computed those in computeCdrGradients. So we just collect everything and then call the parent method which fills the source terms.

  Vector<EBAMRCellData*> cdrSources   = m_cdr->getSources();
  Vector<EBAMRCellData*> rteSources   = m_rte->getSources();
  Vector<EBAMRCellData*> cdrDensities = m_cdr->getPhis();
  Vector<EBAMRCellData*> rteDensities = m_rte->getPhis();

  // Fill the gradient data holders. These have already been computed.
  Vector<EBAMRCellData*> cdrGradients;
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrStorage>& storage = getCdrStorage(solverIt);

    EBAMRCellData& gradient = storage->getGradient();
    cdrGradients.push_back(&gradient);
  }

  // Get the electric field -- we use the one in the scratch storage.
  const EBAMRCellData& electricField = m_fieldScratch->getElectricFieldCell();

  // Compute all source terms for both CDR and RTE equations.
  CdrPlasmaStepper::advanceReactionNetwork(cdrSources,
                                           rteSources,
                                           cdrDensities,
                                           cdrGradients,
                                           rteDensities,
                                           electricField,
                                           m_time,
                                           a_dt);
}

Real
CdrPlasmaGodunovStepper::computeDt()
{
  CH_TIME("CdrPlasmaGodunovStepper::computeDt()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::computeDt()" << endl;
  }

  // TLDR: This routine really depends on what algorithms we use:
  //
  //       Explicit or semi-implicit -> restrict by advection, diffusion, and relaxation time
  //       Partially implicit        -> restrict by advection and relaxation time
  //
  // Note that the semi-implicit scheme does not require restriction by the relaxation time, but users will take
  // care of that through the input script.

  Real dt = std::numeric_limits<Real>::max();

  // First, figure out what the transport time step must be for explicit and explicit-implicit methods.
  if (m_diffusionAlgorithm == DiffusionAlgorithm::Explicit) {
    const Real advectionDt = m_cdr->computeAdvectionDt();
    const Real diffusionDt = m_cdr->computeDiffusionDt();

    m_dtCFL = std::min(advectionDt, diffusionDt);
    dt      = m_cfl * m_dtCFL;

    if (advectionDt < diffusionDt) {
      m_timeCode = TimeCode::Advection;
    }
    else {
      m_timeCode = TimeCode::Diffusion;
    }

    // Turn off implicit diffusion for all species.
    for (int i = 0; i < m_useImplicitDiffusion.size(); i++) {
      m_useImplicitDiffusion[i] = false;
    }
  }
  else if (m_diffusionAlgorithm == DiffusionAlgorithm::Implicit) {
    m_dtCFL    = m_cdr->computeAdvectionDt();
    m_timeCode = TimeCode::Advection;

    dt = m_cfl * m_dtCFL;

    // Turn on implicit diffusion for all species.
    for (int i = 0; i < m_useImplicitDiffusion.size(); i++) {
      m_useImplicitDiffusion[i] = true;
    }
  }
  else if (m_diffusionAlgorithm == DiffusionAlgorithm::Automatic) {

    // When we run with auto-diffusion, we check which species can be done using explicit diffusion and which ones
    // that should use implicit diffusion (based on a user threshold).

    // First. Store the various time step restrictions for the CDR solvers
    std::vector<Real> solverDt;
    std::vector<Real> advectionDt;
    std::vector<Real> diffusionDt;
    std::vector<Real> advectionDiffusionDt;

    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
      solverDt.emplace_back(std::numeric_limits<Real>::max());

      advectionDt.emplace_back(solverIt()->computeAdvectionDt());
      diffusionDt.emplace_back(solverIt()->computeDiffusionDt());
      advectionDiffusionDt.emplace_back(solverIt()->computeAdvectionDiffusionDt());
    }

    // Next, run through the CDR solvers and switch to implicit diffusion for the solvers that satisfy the threshold.
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

      const Real dtA  = advectionDt[idx];
      const Real dtAD = advectionDiffusionDt[idx];

      // Check if this solver should use implicit or explicit diffusion.
      if (dtA / dtAD > m_implicitDiffusionThreshold) {
        solverDt[idx] = dtA;

        m_useImplicitDiffusion[idx] = true;
      }
      else {
        solverDt[idx] = dtAD;

        m_useImplicitDiffusion[idx] = false;
      }
    }

    // Figure out the smallest time step among the solvers
    Real minDt = std::numeric_limits<Real>::max();
    for (const auto& dt : solverDt) {
      minDt = std::min(dt, minDt);
    }

    // In this sweep, go through the solvers again and switch back to explicit diffusion if there's no point in using
    // implicit diffusion.
    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

      const Real dtAD = advectionDiffusionDt[idx];

      // Switch to explicit diffusion if we can.
      if (dtAD > minDt) {
        solverDt[idx]               = dtAD;
        m_useImplicitDiffusion[idx] = false;
      }
    }

    // Finally, we will have found the smallest time step and also figured out which species that implicit/explicit diffusion.
    dt         = minDt;
    m_timeCode = TimeCode::AdvectionDiffusion;

    for (const auto& implicit : m_useImplicitDiffusion) {
      if (implicit) {
        m_timeCode = TimeCode::Advection;
      }
    }

    m_dtCFL = dt / m_cfl;
  }

  // Next, limit by the relaxation time.
  const Real dtRelax = m_relaxTime * this->computeRelaxationTime();
  if (dtRelax < dt) {
    dt         = dtRelax;
    m_timeCode = TimeCode::RelaxationTime;
  }

  // Limit by lower hardcap.
  if (dt < m_minDt) {
    dt         = m_minDt;
    m_timeCode = TimeCode::Hardcap;
  }

  // Limit by upper hardcap.
  if (dt > m_maxDt) {
    dt         = m_maxDt;
    m_timeCode = TimeCode::Hardcap;
  }

  return dt;
}

void
CdrPlasmaGodunovStepper::floorMass(EBAMRCellData&                  a_data,
                                   const std::string               a_message,
                                   const RefCountedPtr<CdrSolver>& a_solver) const
{
  CH_TIME("CdrPlasmaGodunovStepper::floorMass(EBAMRCellData, std::string, RefCountedPtr<CdrSolver>)");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::floorMass(EBAMRCellData, std::string, RefCountedPtr<CdrSolver>)" << endl;
  }

  if (m_debug) {

    // Compute the mass before flooring it.
    const Real massBefore = a_solver->computeMass(a_data);

    DataOps::floor(a_data, 0.0);

    // Compute the mass and relative mass increase after flooring.
    const Real massAfter   = a_solver->computeMass(a_data);
    const Real relMassDiff = (massAfter - massBefore) / massBefore;

    // Print an error message
    pout() << a_message + " - injecting relative " << a_solver->getName() << " mass = " << relMassDiff << endl;
  }
  else {
    DataOps::floor(a_data, 0.0);
  }
}

void
CdrPlasmaGodunovStepper::postStep()
{
  CH_TIME("CdrPlasmaGodunovStepper::postStep()");
  if (m_verbosity > 5) {
    pout() << "CdrPlasmaGodunovStepper::postStep()" << endl;
  }

  // Nothing to see here.
}

#ifdef CH_USE_HDF5
void
CdrPlasmaGodunovStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("CdrPlasmaGodunovStepper::writeCheckpointData(HDF5Handle, int)");
  if (m_verbosity > 3) {
    pout() << "CdrPlasmaGodunovStepper::writeCheckpointData(HDF5Handle, int)" << endl;
  }

  // In addition to the checkpointed stuff from the parent class, the semi-implicit scheme also requires us to checkpoint
  // the factors used in the semi-implicit field solve. We write them here.

  CdrPlasmaStepper::writeCheckpointData(a_handle, a_lvl);

  write(a_handle, *m_semiImplicitRho[a_lvl], "semiImplicitRho");
  write(a_handle, *m_conductivityFactorCell[a_lvl], "conductivityFactor");
}
#endif

#ifdef CH_USE_HDF5
void
CdrPlasmaGodunovStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("CdrPlasmaGodunovStepper::readCheckpointData(HDF5Handle, int)");
  if (m_verbosity > 3) {
    pout() << "CdrPlasmaGodunovStepper::readCheckpointData(HDF5Handle, int)" << endl;
  }

  // In addition to the checkpointed stuff from the parent class, the semi-implicit scheme also requires us to checkpoint
  // the factors used in the semi-implicit field solve. We read them here.

  const Interval interv(0, 0);

  CdrPlasmaStepper::readCheckpointData(a_handle, a_lvl);

  read(a_handle, *m_semiImplicitRho[a_lvl], "semiImplicitRho", m_amr->getGrids(m_realm)[a_lvl], interv, false);
  read(a_handle,
       *m_conductivityFactorCell[a_lvl],
       "conductivityFactor",
       m_amr->getGrids(m_realm)[a_lvl],
       interv,
       false);
}
#endif

#include <CD_NamespaceFooter.H>
