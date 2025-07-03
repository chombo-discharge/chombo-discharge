/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AdvectionDiffusionStepper.cpp
  @brief  Implementation of CD_AdvectionDiffusionStepper.H
  @author Robert Marskar
  @data   March 2020
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_AdvectionDiffusionStepper.H>
#include <CD_AdvectionDiffusionSpecies.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::AdvectionDiffusion;

AdvectionDiffusionStepper::AdvectionDiffusionStepper()
{
  CH_TIME("AdvectionDiffusionStepper::AdvectionDiffusionStepper");

  ParmParse pp("AdvectionDiffusion");

  m_realm = Realm::Primal;
  m_phase = phase::gas;
  m_debug = false;

  pp.query("debug", m_debug);
  pp.get("verbosity", m_verbosity);
  pp.get("cfl", m_cfl);
  pp.get("advection", m_mobile);
  pp.get("diffusion", m_diffusive);

  m_minDt    = 0.0;
  m_maxDt    = std::numeric_limits<Real>::max();
  m_forceCFL = -1.0;

  // Parse the default velocity and diffusion coefficients
  Real         diffCo        = 0.0;
  Real         omega         = 0.0;
  Real         blobAmplitude = 0.0;
  Real         blobRadius    = 0.0;
  RealVect     blobCenter    = RealVect::Zero;
  Vector<Real> v             = Vector<Real>(SpaceDim, 0.0);

  pp.get("diffco", diffCo);
  pp.get("omega", omega);
  pp.get("blob_amplitude", blobAmplitude);
  pp.get("blob_radius", blobRadius);
  pp.getarr("blob_center", v, 0, SpaceDim);

  blobCenter = RealVect(D_DECL(v[0], v[1], v[2]));

  // Set the default velocity and diffusion fields.
  m_velocity = [omega](const RealVect& pos) -> RealVect {
    const Real r     = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
    const Real theta = atan2(pos[1], pos[0]);

    return RealVect(D_DECL(-r * omega * sin(theta), r * omega * cos(theta), 0.));
  };

  m_diffCo = [diffCo](const RealVect& pos) -> Real {
    return diffCo;
  };

  m_initialData = [r = blobRadius, a = blobAmplitude, c = blobCenter](const RealVect& x) -> Real {
    const Real d = (x - c).dotProduct(x - c);
    return a * exp(-d * d / (2 * r * r * r * r));
  };

  this->parseIntegrator();
}

AdvectionDiffusionStepper::AdvectionDiffusionStepper(RefCountedPtr<CdrSolver>& a_solver) : AdvectionDiffusionStepper()
{
  CH_TIME("AdvectionDiffusionStepper::AdvectionDiffusionStepper(full)");

  m_solver = a_solver;
}

AdvectionDiffusionStepper::~AdvectionDiffusionStepper()
{
  CH_TIME("AdvectionDiffusionStepper::~AdvectionDiffusionStepper");
}

void
AdvectionDiffusionStepper::parseRuntimeOptions()
{
  CH_TIME("AdvectionDiffusionStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::parseRuntimeOptions" << endl;
  }

  ParmParse pp("AdvectionDiffusion");

  pp.get("verbosity", m_verbosity);
  pp.get("min_dt", m_minDt);
  pp.get("max_dt", m_maxDt);
  pp.get("cfl", m_cfl);

  this->parseIntegrator();

  m_solver->parseRuntimeOptions();
}

void
AdvectionDiffusionStepper::parseIntegrator()
{
  CH_TIME("AdvectionDiffusionStepper::parseIntegrator");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::parseIntegrator" << endl;
  }

  ParmParse pp("AdvectionDiffusion");

  std::string str;

  pp.get("integrator", str);

  if (str == "heun") {
    m_integrator = Integrator::Heun;
  }
  else if (str == "imex") {
    m_integrator = Integrator::IMEX;
  }
  else {
    MayDay::Error("AdvectionDiffusionStepper::parseIntegrator -- logic bust");
  }
}

void
AdvectionDiffusionStepper::setupSolvers()
{
  CH_TIME("AdvectionDiffusionStepper::setupSolvers");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::setupSolvers" << endl;
  }

  CH_assert(!m_solver.isNull());

  // Instantiate the species.
  m_species = RefCountedPtr<AdvectionDiffusionSpecies>(
    new AdvectionDiffusionSpecies(m_initialData, m_mobile, m_diffusive));

  // Prep the solver.
  m_solver->setVerbosity(m_verbosity);
  m_solver->setSpecies(m_species);
  m_solver->parseOptions();
  m_solver->setPhase(m_phase);
  m_solver->setAmr(m_amr);
  m_solver->setComputationalGeometry(m_computationalGeometry);
  m_solver->setRealm(m_realm);

  if (!m_solver->isMobile() && !m_solver->isDiffusive()) {
    MayDay::Error("AdvectionDiffusionStepper::setupSolvers - can't turn off both advection AND diffusion");
  }
}

void
AdvectionDiffusionStepper::registerRealms()
{
  CH_TIME("AdvectionDiffusionStepper::registerRealms");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_realm);
}

void
AdvectionDiffusionStepper::registerOperators()
{
  CH_TIME("AdvectionDiffusionStepper::registerOperators");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::registerOperators" << endl;
  }

  // Let the solver do this -- it knows what it needs.
  m_solver->registerOperators();
}

void
AdvectionDiffusionStepper::allocate()
{
  CH_TIME("AdvectionDiffusionStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::allocate" << endl;
  }

  m_solver->allocate();
}

void
AdvectionDiffusionStepper::initialData()
{
  CH_TIME("AdvectionDiffusionStepper::initialData");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::initialData" << endl;
  }

  // Fill the solver with initial data from the species.
  m_solver->initialData();

  // Set velocity, diffusion coefficient, and boundary conditions.
  m_solver->setSource(0.0);
  m_solver->setEbFlux(0.0);
  if (m_solver->isDiffusive()) {
    m_solver->setDiffusionCoefficient(m_diffCo);
  }
  if (m_solver->isMobile()) {
    m_solver->setVelocity(m_velocity);
  }

  // Set flux functions
  auto fluxFunc = [](const RealVect a_pos, const Real a_time) {
    return 0.0;
  };

  //  m_solver->setDomainFlux(fluxFunc);
}

#ifdef CH_USE_HDF5
void
AdvectionDiffusionStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("AdvectionDiffusionStepper::writeCheckpointData");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::writeCheckpointData" << endl;
  }

  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
void
AdvectionDiffusionStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("AdvectionDiffusionStepper::readCheckpointData");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::readCheckpointData" << endl;
  }

  m_solver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

void
AdvectionDiffusionStepper::postCheckpointSetup()
{
  CH_TIME("AdvectionDiffusionStepper::postCheckpointSetup");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::postCheckpointSetup" << endl;
  }

  // Set velocity, diffusion coefficient, and boundary conditions.
  m_solver->setSource(0.0);
  m_solver->setEbFlux(0.0);
  if (m_solver->isDiffusive()) {
    m_solver->setDiffusionCoefficient(m_diffCo);
  }
  if (m_solver->isMobile()) {
    m_solver->setVelocity(m_velocity);
  }

  // Set flux functions
  auto fluxFunc = [](const RealVect a_pos, const Real a_time) {
    return 0.0;
  };

  //  m_solver->setDomainFlux(fluxFunc);
}

int
AdvectionDiffusionStepper::getNumberOfPlotVariables() const
{
  CH_TIME("AdvectionDiffusionStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::getNumberOfPlotVariables" << endl;
  }

  // Not plotting anything of our own, so return whatever the solver wants to plot.
  return m_solver->getNumberOfPlotVariables();
}

Vector<std::string>
AdvectionDiffusionStepper::getPlotVariableNames() const
{
  CH_TIME("AdvectionDiffusionStepper::getPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::getPlotVariableNames" << endl;
  }

  return m_solver->getPlotVariableNames();
}

void
AdvectionDiffusionStepper::writePlotData(LevelData<EBCellFAB>& a_output,
                                         int&                  a_icomp,
                                         const std::string     a_outputRealm,
                                         const int             a_level) const
{
  CH_TIME("AdvectionDiffusionStepper::writePlotData");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::writePlotData" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  m_solver->writePlotData(a_output, a_icomp, a_outputRealm, a_level);
}

Real
AdvectionDiffusionStepper::computeDt()
{
  CH_TIME("AdvectionDiffusionStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::computeDt" << endl;
  }

  // TLDR: If we run explicit advection but implicit diffusion then we are only limited by the advective CFL. Otherwise,
  //       if diffusion is also explicit we need the advection-diffusion limited time step.

  // A weird thing, but sometimes we want to be able to force the CFL so that
  // we override run-time configurations of the CFL number. This code does that.
  Real cfl = 0.0;
  if (m_forceCFL > 0.0) {
    cfl = m_forceCFL;
  }
  else {
    cfl = m_cfl;
  }

  Real dt = std::numeric_limits<Real>::max();

  switch (m_integrator) {
  case Integrator::Heun: {
    dt = cfl * m_solver->computeAdvectionDiffusionDt();

    break;
  }
  case Integrator::IMEX: {
    dt = cfl * m_solver->computeAdvectionDt();

    break;
  }
  default: {
    MayDay::Error("AdvectionDiffusionStepper::computeDt - logic bust");

    break;
  }
  }

  dt = std::max(dt, m_minDt);
  dt = std::min(dt, m_maxDt);

  return dt;
}

Real
AdvectionDiffusionStepper::advance(const Real a_dt)
{
  CH_TIME("AdvectionDiffusionStepper::advance");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::advance" << endl;
  }

  // State to be advanced.
  EBAMRCellData& state = m_solver->getPhi();

  // Compute mass before advance.
  const Real initialMass = m_solver->computeMass();

  switch (m_integrator) {
  case Integrator::Heun: {
    const bool conservativeOnly = false;
    const bool addEbFlux        = true;
    const bool addDomainFlux    = true;

    // Transient storage
    EBAMRCellData yp;
    EBAMRCellData k1;
    EBAMRCellData k2;

    m_amr->allocate(yp, m_realm, m_phase, 1);
    m_amr->allocate(k1, m_realm, m_phase, 1);
    m_amr->allocate(k2, m_realm, m_phase, 1);

    // Compute k1 coefficient
    m_solver->computeDivJ(k1, state, 0.0, conservativeOnly, addEbFlux, addDomainFlux);
    DataOps::copy(yp, state);
    DataOps::incr(yp, k1, -a_dt);

    // Compute k2 coefficient and final state
    m_solver->computeDivJ(k2, yp, 0.0, conservativeOnly, addEbFlux, addDomainFlux);
    DataOps::incr(state, k1, -0.5 * a_dt);
    DataOps::incr(state, k2, -0.5 * a_dt);

    break;
  }
  case Integrator::IMEX: {
    const bool addEbFlux     = true;
    const bool addDomainFlux = true;

    // Transient storage
    EBAMRCellData yp;
    EBAMRCellData k1;
    EBAMRCellData k2;

    m_amr->allocate(k1, m_realm, m_phase, 1);
    m_amr->allocate(k2, m_realm, m_phase, 1);

    if (m_solver->isDiffusive()) {

      // Compute the finite volume approximation to kappa*div(F). The second "hook" is a debugging hook that includes redistribution when computing kappa*div(F). It
      // exists only for debugging/assurance reasons.
      if (false) {
        const bool conservativeOnly = true;

        m_solver->computeDivF(k1, state, a_dt, conservativeOnly, addEbFlux, addDomainFlux);
      }
      else {
        const bool conservativeOnly = false;

        m_solver->computeDivF(k1, state, a_dt, conservativeOnly, addEbFlux, addDomainFlux);
      }

      DataOps::kappaScale(k1);
      DataOps::scale(k1, -1.0);

      // Use k1 as the old solution.
      DataOps::copy(k2, state);

      // Do the Euler solve.
      m_solver->advanceCrankNicholson(state, k2, k1, a_dt);
    }
    else { // Purely inviscid advance.
      m_solver->computeDivF(k1, state, a_dt, false, addEbFlux, addDomainFlux);

      DataOps::incr(state, k1, -a_dt);
    }

    break;
  }
  default:
    MayDay::Error("AdvectionDiffusionStepper - unknown integrator requested");
    break;
  }

  m_amr->conservativeAverage(state, m_realm, m_phase);
  m_amr->interpGhost(state, m_realm, m_phase);

  // Compute final mass.
  const Real finalMass = m_solver->computeMass();

  if (procID() == 0 && m_debug) {
    std::cout << "step = " << m_timeStep + 1 << "\t\t\t"
              << "mass conservation % = " << 100. * (finalMass - initialMass) / initialMass << std::endl;
  }

  return a_dt;
}

void
AdvectionDiffusionStepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("AdvectionDiffusionStepper::synchronizeSolverTimes");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_solver->setTime(a_step, a_time, a_dt);
}

void
AdvectionDiffusionStepper::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("AdvectionDiffusionStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::preRegrid" << endl;
  }

  m_solver->preRegrid(a_lbase, a_oldFinestLevel);
}

void
AdvectionDiffusionStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("AdvectionDiffusionStepper::regrid");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::regrid" << endl;
  }

  // Regrid CDR solver and set up the flow fields
  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_solver->setSource(0.0);
  m_solver->setEbFlux(0.0);

  if (m_solver->isDiffusive()) {
    m_solver->setDiffusionCoefficient(m_diffCo);
  }
  if (m_solver->isMobile()) {
    m_solver->setVelocity(m_velocity);
  }
}

void
AdvectionDiffusionStepper::setCFL(const Real a_cfl)
{
  CH_TIME("AdvectionDiffusionStepper::setCFL");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::setCFL" << endl;
  }

  m_forceCFL = a_cfl;
}

void
AdvectionDiffusionStepper::setInitialData(const std::function<Real(const RealVect& a_position)>& a_initData) noexcept
{
  CH_TIME("AdvectionDiffusionStepper::setInitialData");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::setInitialData" << endl;
  }

  m_initialData = a_initData;
}

void
AdvectionDiffusionStepper::setVelocity(const std::function<RealVect(const RealVect& a_position)>& a_velocity) noexcept
{
  CH_TIME("AdvectionDiffusionStepper::setVelocity");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::setVelocity" << endl;
  }

  m_velocity = a_velocity;
}

void
AdvectionDiffusionStepper::setDiffusionCoefficient(
  const std::function<Real(const RealVect& a_position)>& a_diffusion) noexcept
{
  CH_TIME("AdvectionDiffusionStepper::setDiffusionCoefficient");
  if (m_verbosity > 5) {
    pout() << "AdvectionDiffusionStepper::setDiffusionCoefficient" << endl;
  }

  m_diffCo = a_diffusion;
}

#include <CD_NamespaceFooter.H>
