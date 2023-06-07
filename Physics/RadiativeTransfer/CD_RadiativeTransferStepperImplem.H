/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RadiativeTransferStepperImplem.H
  @brief  Implementation of CD_RadiativeTransferStepper.H
  @author Robert Marskar
*/

// Std includes
#include <math.h>

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>
#include <PolyGeom.H>

// Our includes
#include <CD_RadiativeTransferStepper.H>
#include <CD_RadiativeTransferSpecies.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::RadiativeTransfer;

template <typename T>
RadiativeTransferStepper<T>::RadiativeTransferStepper()
{
  CH_TIME("RadiativeTransferStepper::RadiativeTransferStepper");

  // Get shit for setting the source term
  ParmParse    pp("RadiativeTransferStepper");
  Vector<Real> v;

  pp.get("verbosity", m_verbosity);
  pp.get("realm", m_realm);
  pp.get("dt", m_dt);

  m_forceDt = -1.0;
}

template <typename T>
RadiativeTransferStepper<T>::~RadiativeTransferStepper()
{
  CH_TIME("RadiativeTransferStepper::~RadiativeTransferStepper");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::~RadiativeTransferStepper" << endl;
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::setupSolvers()
{
  CH_TIME("RadiativeTransferStepper::setupSolvers");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::setupSolvers" << endl;
  }

  m_species = RefCountedPtr<RtSpecies>(new RadiativeTransferSpecies());

  // Solver setup
  m_solver = RefCountedPtr<RtSolver>(new T());
  m_solver->setVerbosity(m_verbosity);
  m_solver->setRtSpecies(m_species);
  m_solver->parseOptions();
  m_solver->setPhase(phase::gas);
  m_solver->setAmr(m_amr);
  m_solver->setComputationalGeometry(m_computationalGeometry);
  m_solver->setRealm(m_realm);
}

template <typename T>
void
RadiativeTransferStepper<T>::parseRuntimeOptions()
{
  CH_TIME("RadiativeTransferStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::parseRuntimeOptions" << endl;
  }

  ParmParse pp("RadiativeTransferStepper");

  pp.get("verbosity", m_verbosity);
  pp.get("dt", m_dt);

  m_solver->parseRuntimeOptions();
}

template <typename T>
void
RadiativeTransferStepper<T>::allocate()
{
  CH_TIME("RadiativeTransferStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::allocate" << endl;
  }

  m_solver->allocate();
}

template <typename T>
void
RadiativeTransferStepper<T>::initialData()
{
  CH_TIME("RadiativeTransferStepper::initialData");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::initialData" << endl;
  }

  m_solver->initialData();
  this->setGaussianSource();

  if (m_solver->isStationary()) {
    m_solver->advance(0.0, false);
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::postInitialize()
{
  CH_TIME("RadiativeTransferStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::postInitialize" << endl;
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::forceDt(const Real a_dt)
{
  CH_TIME("RadiativeTransferStepper::forceDt");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::forceDt" << endl;
  }

  m_forceDt = a_dt;
}

#ifdef CH_USE_HDF5
template <typename T>
void
RadiativeTransferStepper<T>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("RadiativeTransferStepper::writeCheckpointData");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::writeCheckpointData" << endl;
  }

  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
template <typename T>
void
RadiativeTransferStepper<T>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("RadiativeTransferStepper::readCheckpointData");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::readCheckpointData" << endl;
  }

  m_solver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

template <typename T>
void
RadiativeTransferStepper<T>::postCheckpointSetup()
{
  CH_TIME("RadiativeTransferStepper::postCheckpointSetup");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::postCheckpointSetup" << endl;
  }

  this->setGaussianSource();
}

template <typename T>
int
RadiativeTransferStepper<T>::getNumberOfPlotVariables() const
{
  CH_TIME("RadiativeTransferStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::getNumberOfPlotVariables" << endl;
  }

  return m_solver->getNumberOfPlotVariables();
}

template <typename T>
Vector<std::string>
RadiativeTransferStepper<T>::getPlotVariableNames() const
{
  CH_TIME("RadiativeTransferStepper::getPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::getPlotVariableNames" << endl;
  }

  return m_solver->getPlotVariableNames();
}

template <typename T>
void
RadiativeTransferStepper<T>::writePlotData(LevelData<EBCellFAB>& a_output, int& a_icomp, const int a_level) const
{
  CH_TIME("RadiativeTransferStepper::writePlotData");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::writePlotData" << endl;
  }

  m_solver->writePlotData(a_output, a_icomp, a_level);
}

template <typename T>
Real
RadiativeTransferStepper<T>::computeDt()
{
  CH_TIME("RadiativeTransferStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::computeDt" << endl;
  }

  Real dt = m_dt;

  if (m_forceDt > 0.0) {
    dt = m_forceDt;
  }

  return dt;
}

template <typename T>
Real
RadiativeTransferStepper<T>::advance(const Real a_dt)
{
  CH_TIME("RadiativeTransferStepper::advance");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::advance" << endl;
  }

  m_solver->advance(a_dt);

  return a_dt;
}

template <typename T>
void
RadiativeTransferStepper<T>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("RadiativeTransferStepper::synchronizeSolverTimes");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_solver->setTime(m_timeStep, m_time, m_dt);
}

template <typename T>
void
RadiativeTransferStepper<T>::printStepReport()
{
  CH_TIME("RadiativeTransferStepper::printStepReport");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::printStepReport" << endl;
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::registerRealms()
{
  CH_TIME("RadiativeTransferStepper::registerRealms");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_realm);
}

template <typename T>
void
RadiativeTransferStepper<T>::registerOperators()
{
  CH_TIME("RadiativeTransferStepper::registerOperators");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::registerOperators" << endl;
  }

  m_solver->registerOperators();
}

template <typename T>
void
RadiativeTransferStepper<T>::preRegrid(const int a_base, const int a_oldFinestLevel)
{
  CH_TIME("RadiativeTransferStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::preRegrid" << endl;
  }

  m_solver->preRegrid(a_base, a_oldFinestLevel);
}

template <typename T>
void
RadiativeTransferStepper<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("RadiativeTransferStepper::regrid");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::regrid" << endl;
  }

  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
}

template <typename T>
void
RadiativeTransferStepper<T>::postRegrid()
{
  CH_TIME("RadiativeTransferStepper::postRegrid");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::postRegrid" << endl;
  }

  this->setGaussianSource();
}

template <typename T>
const EBAMRCellData&
RadiativeTransferStepper<T>::getPhi() const
{
  CH_TIME("RadiativeTransferStepper::getPhi");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::getPhi" << endl;
  }

  return m_solver->getPhi();
}

template <typename T>
void
RadiativeTransferStepper<T>::setGaussianSource()
{
  CH_TIME("RadiativeTransferStepper::setGaussianSource");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::setGaussianSource" << endl;
  }

  ParmParse pp("RadiativeTransferStepper");

  Vector<Real> v(SpaceDim);
  RealVect     blobCenter;
  Real         blobAmp;
  Real         blobRad;

  pp.get("blob_amplitude", blobAmp);
  pp.get("blob_radius", blobRad);
  pp.getarr("blob_center", v, 0, SpaceDim);
  blobCenter = RealVect(D_DECL(v[0], v[1], v[2]));

  auto sourceFunction = [blobAmp, blobRad, blobCenter](const RealVect a_position) -> Real {
    const RealVect dist  = a_position - blobCenter;
    const Real     dist2 = dist.dotProduct(dist);
    return blobAmp * exp(-dist2 / (2 * blobRad * blobRad));
  };

  m_solver->setSource(sourceFunction);
}

#include <CD_NamespaceFooter.H>
