/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrLayoutImplem.H
  @brief  Implementation of CD_CdrLayout.H
  @author Robert Marskar
*/

#ifndef CD_CdrLayoutImplem_H
#define CD_CdrLayoutImplem_H

// Our includes
#include <CD_CdrLayout.H>
#include <CD_CdrIterator.H>
#include <CD_Units.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

template <class T>
CdrLayout<T>::CdrLayout(const Vector<RefCountedPtr<CdrSpecies>>& a_species)
{
  CH_TIME("CdrLayout<T>::CdrLayout()");

  // Default settings.
  m_verbosity = -1;
  m_species   = a_species;
  m_realm     = Realm::Primal;

  m_solvers.resize(0);
}

template <class T>
CdrLayout<T>::~CdrLayout()
{
  CH_TIME("CdrLayout<T>::~CdrLayout()");
}

template <class T>
phase::which_phase
CdrLayout<T>::getPhase() const
{
  CH_TIME("CdrLayout<T>::getPhase()");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::getPhase()" << endl;
  }

  return m_phase;
}

template <class T>
CdrIterator<T>
CdrLayout<T>::iterator()
{
  CH_TIME("CdrLayout<T>::iterator()");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::iterator()" << endl;
  }

  return CdrIterator<T>(*this);
}

template <class T>
std::string
CdrLayout<T>::getRealm() const
{
  CH_TIME("CdrLayout<T>::getRealm()");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::getRealm()" << endl;
  }

  return m_realm;
}

template <class T>
void
CdrLayout<T>::setRealm(const std::string a_realm)
{
  CH_TIME("CdrLayout<T>::setRealm()");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::setRealm()" << endl;
  }

  m_realm = a_realm;

  for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->setRealm(a_realm);
  }
}

template <class T>
void
CdrLayout<T>::addSolver(RefCountedPtr<T>& a_solver)
{
  CH_TIME("CdrLayout<T>::addSolver(RefCountedPtr<T>");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::addSolver()" << endl;
  }

  m_solvers.push_back(a_solver);
}

template <class T>
void
CdrLayout<T>::parseOptions()
{
  CH_TIME("CdrLayout<T>::parseOptions()");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::parseOptions()" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->parseOptions();
  }
}

template <class T>
void
CdrLayout<T>::parseRuntimeOptions()
{
  CH_TIME("CdrLayout<T>::parseRuntimeOptions()");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::parseRuntimeOptions()" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->parseRuntimeOptions();
  }
}

template <class T>
void
CdrLayout<T>::allocate()
{
  CH_TIME("CdrLayout<T>::allocate()");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::allocate()" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->allocate();
  }
}

template <class T>
void
CdrLayout<T>::deallocate()
{
  CH_TIME("CdrLayout<T>::deallocate()");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::deallocate()" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->deallocate();
  }
}

template <class T>
void
CdrLayout<T>::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("CdrLayout<T>::preRegrid(int, int)");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::preRegrid(int, int)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->preRegrid(a_lbase, a_oldFinestLevel);
  }
}

template <class T>
void
CdrLayout<T>::initialData()
{
  CH_TIME("CdrLayout<T>::initialData()");
  if (m_verbosity > 6) {
    pout() << "CdrLayout<T>::initialData()" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->initialData();
  }
}

template <class T>
void
CdrLayout<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("CdrLayout<T>::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::regrid(int, int, int)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
}

template <class T>
void
CdrLayout<T>::registerOperators()
{
  CH_TIME("CdrLayout<T>::registerOperators()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::registerOperators()" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->registerOperators();
  }
}

template <class T>
void
CdrLayout<T>::setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh)
{
  CH_TIME("CdrLayout<T>::setAmr(RefCountedPtr<AmrMesh>)");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::setAmr(RefCountedPtr<AmrMesh>)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->setAmr(a_amrMesh);
  }
}

template <class T>
void
CdrLayout<T>::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("CdrLayout<T>::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->setComputationalGeometry(a_computationalGeometry);
  }
}

template <class T>
void
CdrLayout<T>::setPhase(const phase::which_phase a_phase)
{
  CH_TIME("CdrLayout<T>::setPhase(phase::which_phase)");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::setPhase(phase::which_phase)" << endl;
  }

  m_phase = a_phase;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->setPhase(m_phase);
  }
}

template <class T>
void
CdrLayout<T>::setVerbosity(const int a_verbosity)
{
  CH_TIME("CdrLayout<T>::setVerbosity(int)");
  m_verbosity = a_verbosity;
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::setVerbosity(int)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->setVerbosity(a_verbosity);
  }
}

template <class T>
void
CdrLayout<T>::setTime(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("CdrLayout<T>::setTime(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::setTime(int, int, int)" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->setTime(a_step, a_time, a_dt);
  }
}

template <class T>
void
CdrLayout<T>::writePlotFile()
{
  CH_TIME("CdrLayout<T>::writePlotFile()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::writePlotFile()" << endl;
  }

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->writePlotFile();
  }
}

template <class T>
Real
CdrLayout<T>::computeAdvectionDt()
{
  CH_TIME("CdrLayout<T>::computeAdvectionDt()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::computeAdvectionDt()" << endl;
  }

  Real dt = std::numeric_limits<Real>::max();

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    const Real curDt = solver_it()->computeAdvectionDt();

    dt = std::min(dt, curDt);
  }

  return dt;
}

template <class T>
Real
CdrLayout<T>::computeDiffusionDt()
{
  CH_TIME("CdrLayout<T>::computeDiffusionDt()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::computeDiffusionDt()" << endl;
  }

  Real dt = std::numeric_limits<Real>::max();

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    const Real curDt = solver_it()->computeDiffusionDt();

    dt = std::min(dt, curDt);
  }

  return dt;
}

template <class T>
Real
CdrLayout<T>::computeAdvectionDiffusionDt()
{
  CH_TIME("CdrLayout<T>::computeAdvectionDiffusionDt()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::computeAdvectionDiffusionDt()" << endl;
  }

  Real dt = std::numeric_limits<Real>::max();

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    const Real curDt = solver_it()->computeAdvectionDiffusionDt();

    dt = std::min(dt, curDt);
  }

  return dt;
}

template <class T>
Vector<RefCountedPtr<T>>&
CdrLayout<T>::getSolvers()
{
  CH_TIME("CdrLayout<T>::getSolvers()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::getSolvers()" << endl;
  }

  return (m_solvers);
}

template <class T>
Vector<RefCountedPtr<CdrSpecies>>&
CdrLayout<T>::getSpecies()
{
  CH_TIME("CdrLayout<T>::getSpecies()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::getSpecies()" << endl;
  }

  return (m_species);
}

template <class T>
Vector<EBAMRCellData*>
CdrLayout<T>::getPhis()
{
  CH_TIME("CdrLayout<T>::getPhis()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::getPhis()" << endl;
  }

  Vector<EBAMRCellData*> states;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    states.push_back(&(solver_it()->getPhi()));
  }

  return states;
}

template <class T>
Vector<EBAMRCellData*>
CdrLayout<T>::getSources()
{
  CH_TIME("CdrLayout<T>::getSources()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::getSources()" << endl;
  }

  Vector<EBAMRCellData*> sources;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    sources.push_back(&(solver_it()->getSource()));
  }

  return sources;
}

template <class T>
Vector<EBAMRCellData*>
CdrLayout<T>::getVelocities()
{
  CH_TIME("CdrLayout<T>::getVelocities()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::getVelocities()" << endl;
  }

  Vector<EBAMRCellData*> velocities;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    velocities.push_back(&(solver_it()->getCellCenteredVelocity()));
  }

  return velocities;
}

template <class T>
Vector<EBAMRFluxData*>
CdrLayout<T>::getFaceCenteredDiffusionCoefficient()
{
  CH_TIME("CdrLayout<T>::getFaceCenteredDiffusionCoefficient()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::getFaceCenteredDiffusionCoefficient()" << endl;
  }

  Vector<EBAMRFluxData*> diffCo;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    diffCo.push_back(&(solver_it()->getFaceCenteredDiffusionCoefficient()));
  }

  return diffCo;
}

template <class T>
Vector<EBAMRIVData*>
CdrLayout<T>::getEbCenteredDiffusionCoefficient()
{
  CH_TIME("CdrLayout<T>::getEbCenteredDiffusionCoefficient");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::getEbCenteredDiffusionCoefficient" << endl;
  }

  Vector<EBAMRIVData*> diffCo;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    diffCo.push_back(&(solver_it()->getEbCenteredDiffusionCoefficient()));
  }

  return diffCo;
}

template <class T>
Vector<EBAMRIVData*>
CdrLayout<T>::getEbFlux()
{
  CH_TIME("CdrLayout<T>::getEbFlux()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::getEbFlux()" << endl;
  }

  Vector<EBAMRIVData*> ebFlux;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    ebFlux.push_back(&(solver_it()->getEbFlux()));
  }

  return ebFlux;
}

template <class T>
Vector<EBAMRIFData*>
CdrLayout<T>::getDomainFlux()
{
  CH_TIME("CdrLayout<T>::getDomainFlux()");
  if (m_verbosity > 5) {
    pout() << "CdrLayout<T>::getDomainFlux()" << endl;
  }

  Vector<EBAMRIFData*> domainFlux;

  for (CdrIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    domainFlux.push_back(&(solver_it()->getDomainFlux()));
  }

  return domainFlux;
}

template <class T, class S>
CdrFactory<T, S>::CdrFactory()
{
  CH_TIME("CdrFactory<T, S>::CdrFactory()");
}

template <class T, class S>
CdrFactory<T, S>::~CdrFactory()
{
  CH_TIME("CdrFactory<T, S>::~CdrFactory()");
}

template <class T, class S>
RefCountedPtr<CdrLayout<T>>
CdrFactory<T, S>::newLayout(const Vector<RefCountedPtr<CdrSpecies>>& a_species) const
{
  CH_TIME("CdrFactory<T, S>::newLayout(Vector<RefCountedPtr<CdrSpecies> >)");

  // Instatiate the return CdrLayout object first.
  auto cdr = RefCountedPtr<CdrLayout<T>>(new CdrLayout<T>(a_species));

  // Cast solvers and instantiate them
  for (int i = 0; i < a_species.size(); i++) {

    auto solver = RefCountedPtr<T>(static_cast<T*>(new S()));

    solver->setSpecies(a_species[i]);
    solver->setPhase(phase::gas);
    solver->setVerbosity(-1);

    cdr->addSolver(solver);
  }

  return cdr;
}

#include <CD_NamespaceFooter.H>

#endif
