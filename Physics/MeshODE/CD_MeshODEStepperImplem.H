/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MeshODEStepperImplem.H
  @brief  Implementation of CD_MeshODEStepper.H
  @author Robert Marskar
*/

#ifndef CD_MeshODEStepperImplem_H
#define CD_MeshODEStepperImplem_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_MeshODESolver.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::MeshODE;

template <size_t N>
MeshODEStepper<N>::MeshODEStepper()
{
  CH_TIME("MeshODEStepper::MeshODEStepper");

  m_verbosity = 10;
  m_realm     = Realm::Primal;
  m_phase     = phase::gas;
}

template <size_t N>
MeshODEStepper<N>::~MeshODEStepper()
{
  CH_TIME("MeshODEStepper::~MeshODEStepper");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::~MeshODEStepper" << endl;
  }
}

template <size_t N>
void
MeshODEStepper<N>::setupSolvers()
{
  CH_TIME("MeshODEStepper::setupSolvers");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::setupSolvers" << endl;
  }

  m_solver = RefCountedPtr<MeshODESolver<N>>(new MeshODESolver<N>(m_amr));

  m_solver->setPhase(m_phase);
  m_solver->setRealm(m_realm);
  m_solver->parseOptions();

  // Initialize functions.
  m_initialData = [](const RealVect& a_pos) -> std::array<Real, N> {
    std::array<Real, N> Y;
    for (auto& y : Y) {
      y = a_pos[0];
    }

    return Y;
  };

  m_rhsFunction = [](const std::array<Real, N>, const Real) -> std::array<Real, N> { return std::array<Real, N>{}; };
}

template <size_t N>
void
MeshODEStepper<N>::allocate()
{
  CH_TIME("MeshODEStepper::allocate()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::allocate()" << endl;
  }

  m_solver->allocate();
}

template <size_t N>
void
MeshODEStepper<N>::initialData()
{
  CH_TIME("MeshODEStepper::initialData()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::initialData()" << endl;
  }

  m_solver->setPhi(m_initialData);
}

template <size_t N>
void
MeshODEStepper<N>::postInitialize()
{
  CH_TIME("MeshODEStepper::postInitialize()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::postInitialize()" << endl;
  }

  m_solver->computeRHS(m_rhsFunction);
}

template <size_t N>
void
MeshODEStepper<N>::registerRealms()
{
  CH_TIME("MeshODEStepper::registerRealms()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::registerRealms()" << endl;
  }

  m_amr->registerRealm(Realm::Primal);
}

template <size_t N>
void
MeshODEStepper<N>::registerOperators()
{
  CH_TIME("MeshODEStepper::registerOperators()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::registerOperators()" << endl;
  }

  m_solver->registerOperators();
}

template <size_t N>
void
MeshODEStepper<N>::parseOptions()
{
  CH_TIME("MeshODEStepper::parseOptions()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::parseOptions()" << endl;
  }

  m_solver->parseOptions();

  this->parseIntegrator();
  this->parseVerbosity();
}

template <size_t N>
void
MeshODEStepper<N>::parseRuntimeOptions()
{
  CH_TIME("MeshODEStepper::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::parseRuntimeOptions()" << endl;
  }

  m_solver->parseRuntimeOptions();

  this->parseIntegrator();
  this->parseVerbosity();
}

#ifdef CH_USE_HDF5
template <size_t N>
void
MeshODEStepper<N>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("MeshODEStepper::writeCheckpointData(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::writeCheckpointData(HDF5Handle, int)" << endl;
  }

  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
template <size_t N>
void
MeshODEStepper<N>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("MeshODEStepper::readCheckpointData(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::readCheckpointData(HDF5Handle, int)" << endl;
  }

  m_solver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

template <size_t N>
void
MeshODEStepper<N>::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const
{
  CH_TIME("MeshODEStepper::writePlotData(EBAMRCellData, Vector<std::string>, int)");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::writePlotData(EBAMRCellData, Vector<std::string>, int)" << endl;
  }

  m_solver->writePlotData(a_output, a_icomp);
  a_plotVariableNames.append(m_solver->getPlotVariableNames());
}

template <size_t N>
int
MeshODEStepper<N>::getNumberOfPlotVariables() const
{
  CH_TIME("MeshODEStepper::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::getNumberOfPlotVariables()" << endl;
  }

  return m_solver->getNumberOfPlotVariables();
}

template <size_t N>
void
MeshODEStepper<N>::computeDt(Real& a_dt, TimeCode& a_timeCode)
{
  CH_TIME("MeshODEStepper::computeDt()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::computeDt()" << endl;
  }

  a_dt = std::numeric_limits<Real>::max();
}

template <size_t N>
Real
MeshODEStepper<N>::advance(const Real a_dt)
{
  CH_TIME("MeshODEStepper::avdvance(Real)");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::advance(Real)" << endl;
  }

  return a_dt;
}

template <size_t N>
void
MeshODEStepper<N>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("MeshODEStepper::synchronizeSolverTimes(int, Real, Real)");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::synchronizeSolverTimes(int, Real, Real)" << endl;
  }

  m_solver->setTime(a_step, a_time, a_dt);
}

template <size_t N>
void
MeshODEStepper<N>::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("MeshODEStepper::preRegrid(int, int)");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::preRegrid(int, int)" << endl;
  }

  m_solver->preRegrid(a_lmin, a_oldFinestLevel);
}

template <size_t N>
void
MeshODEStepper<N>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("MeshODEStepper::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::regrid(int, int, int)" << endl;
  }

  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
}

template <size_t N>
void
MeshODEStepper<N>::postRegrid()
{
  CH_TIME("MeshODEStepper::postRegrid()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::postRegrid()" << endl;
  }

  // Compute new right-hand side.
  m_solver->computeRHS(m_rhsFunction);
}

template <size_t N>
void
MeshODEStepper<N>::parseIntegrator()
{
  CH_TIME("MeshODEStepper::parseIntegrator()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::parseIntegrator()" << endl;
  }

  ParmParse pp("MeshODEStepper");

  std::string str;

  pp.get("integration", str);
  if (str == "euler") {
    m_algorithm = IntegrationAlgorithm::Euler;
  }
  else if (str == "rk2") {
    m_algorithm = IntegrationAlgorithm::RK2;
  }
  else if (str == "rk4") {
    m_algorithm = IntegrationAlgorithm::RK4;
  }
  else {
    MayDay::Error("MeshODEStepper::parseIntegrator -- logic bust");
  }
}

template <size_t N>
void
MeshODEStepper<N>::parseVerbosity()
{
  CH_TIME("MeshODEStepper::parseVerbosity()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::parseVerbosity()" << endl;
  }

  ParmParse pp("MeshODEStepper");

  pp.get("verbosity", m_verbosity);
}

template <size_t N>
void
MeshODEStepper<N>::advanceEuler(const Real a_dt)
{
  CH_TIME("MeshODEStepper::advanceEuler()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::advanceEuler()" << endl;
  }

  MayDay::Error("MeshODEStepper::advanceEuler -- not implemented");
}

template <size_t N>
void
MeshODEStepper<N>::advanceRK2(const Real a_dt)
{
  CH_TIME("MeshODEStepper::advanceRK2()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::advanceRK2()" << endl;
  }

  MayDay::Error("MeshODEStepper::advanceRK2 -- not implemented");
}

template <size_t N>
void
MeshODEStepper<N>::advanceRK4(const Real a_dt)
{
  CH_TIME("MeshODEStepper::advanceRK4()");
  if (m_verbosity > 5) {
    pout() << "MeshODEStepper::advanceRK4()" << endl;
  }

  MayDay::Error("MeshODEStepper::advanceRK4 -- not implemented");
}

#include <CD_NamespaceFooter.H>

#endif
