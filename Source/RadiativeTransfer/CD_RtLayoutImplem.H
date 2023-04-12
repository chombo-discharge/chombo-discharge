/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtLayoutImplem.H
  @brief  Implementation of CD_RtLayout.H
  @author Robert Marskar
*/

#ifndef CD_RtLayoutImplem_H
#define CD_RtLayoutImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_RtIterator.H>
#include <CD_RtLayout.H>
#include <CD_NamespaceHeader.H>

template <class T>
RtLayout<T>::RtLayout(const Vector<RefCountedPtr<RtSpecies>>& a_RtSpecies)
{
  CH_TIME("RtLayout<T>::RtLayout(Vector)");

  // Default settings.
  m_verbosity = -1;
  m_species   = a_RtSpecies;
  m_solvers.resize(0);
}

template <class T>
RtLayout<T>::~RtLayout()
{
  CH_TIME("RtLayout<T>::~RtLayout");
}

template <class T>
RtIterator<T>
RtLayout<T>::iterator()
{
  return RtIterator<T>(*this);
}

template <class T>
const std::string
RtLayout<T>::getRealm() const
{
  return m_realm;
}

template <class T>
void
RtLayout<T>::setRealm(const std::string a_realm)
{
  m_realm = a_realm;

  for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->setRealm(m_realm);
  }
}

template <class T>
void
RtLayout<T>::addSolver(RefCountedPtr<T> a_solver)
{
  m_solvers.push_back(a_solver);
}

template <class T>
void
RtLayout<T>::parseOptions()
{
  CH_TIME("RtLayout<T>::parseOptions");
  if (m_verbosity > 6) {
    pout() << "RtLayout<T>::parseOptions" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->parseOptions();
  }
}

template <class T>
void
RtLayout<T>::parseRuntimeOptions()
{
  CH_TIME("RtLayout<T>::parseRuntimeOptions");
  if (m_verbosity > 6) {
    pout() << "RtLayout<T>::parseRuntimeOptions" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->parseRuntimeOptions();
  }
}

template <class T>
void
RtLayout<T>::allocate()
{
  CH_TIME("RtLayout<T>::allocate");
  if (m_verbosity > 6) {
    pout() << "RtLayout<T>::allocate" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->allocate();
  }
}

template <class T>
void
RtLayout<T>::preRegrid(const int a_base, const int a_finestLevel)
{
  CH_TIME("RtLayout<T>::preRegrid");
  if (m_verbosity > 6) {
    pout() << "RtLayout<T>::preRegrid" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->preRegrid(a_base, a_finestLevel);
  }
}

template <class T>
void
RtLayout<T>::deallocate()
{
  CH_TIME("RtLayout<T>::deallocate");
  if (m_verbosity > 6) {
    pout() << "RtLayout<T>::deallocate" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->deallocate();
  }
}

template <class T>
void
RtLayout<T>::setAmr(const RefCountedPtr<AmrMesh>& a_amr)
{
  CH_TIME("RtLayout<T>::setAmr");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::setAmr" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->setAmr(a_amr);
  }
}

template <class T>
void
RtLayout<T>::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("RtLayout<T>::setComputationalGeometry");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::setComputationalGeometry" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->setComputationalGeometry(a_computationalGeometry);
  }
}

template <class T>
void
RtLayout<T>::setPhase(phase::which_phase a_phase)
{
  CH_TIME("RtLayout<T>::setPhase");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::setPhase" << endl;
  }

  m_phase = a_phase;

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->setPhase(a_phase);
  }
}

template <class T>
void
RtLayout<T>::setVerbosity(const int a_verbosity)
{
  CH_TIME("RtLayout<T>::setVerbosity");

  m_verbosity = a_verbosity;
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::setVerbosity" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->setVerbosity(a_verbosity);
  }
}

template <class T>
void
RtLayout<T>::sanityCheck()
{
  CH_TIME("RtLayout<T>::sanityCheck");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::sanityCheck" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->sanityCheck();
  }
}

template <class T>
void
RtLayout<T>::setTime(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("RtLayout<T>::setTime");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::setTime" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->setTime(a_step, a_time, a_dt);
  }
}

template <class T>
void
RtLayout<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("RtLayout<T>::regrid");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::regrid" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
}

template <class T>
void
RtLayout<T>::registerOperators()
{
  CH_TIME("RtLayout<T>::registerOperators");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::registerOperators" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->registerOperators();
  }
}

template <class T>
void
RtLayout<T>::setSource(const EBAMRCellData& a_source)
{
  CH_TIME("RtLayout<T>::setSource(ebamrcell)");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::setSource(ebamrcell)" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->setSource(a_source);
  }
}

template <class T>
void
RtLayout<T>::setSource(const Real a_source)
{
  CH_TIME("RtLayout<T>::setSource(constant)");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::setSource(constant)" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->setSource(a_source);
  }
}

template <class T>
void
RtLayout<T>::setStationary(const bool a_stationary)
{
  CH_TIME("RtLayout<T>::setStationary");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::setStationary" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->setStationary(a_stationary);
  }
}

template <class T>
void
RtLayout<T>::writePlotFile()
{
  CH_TIME("RtLayout<T>::writePlotFile");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::writePlotFile" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->writePlotFile();
  }
}

template <class T>
void
RtLayout<T>::advance(const Real a_dt)
{
  CH_TIME("RtLayout<T>::advance");
  if (m_verbosity > 6) {
    pout() << "RtLayout<T>::advance" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->advance(a_dt, solver->getPhi(), solver->getSource());
  }
}

template <class T>
void
RtLayout<T>::initialData()
{
  CH_TIME("RtLayout<T>::initialData");
  if (m_verbosity > 6) {
    pout() << "RtLayout<T>::initialData" << endl;
  }

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    solver->initialData();
  }
}

template <class T>
bool
RtLayout<T>::isStationary()
{
  CH_TIME("RtLayout<T>::isStationary");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::isStationary" << endl;
  }

  bool stationary = true;

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();

    if (solver->isStationary() == false) {
      stationary = false;
    }
  }

  return stationary;
}

template <class T>
phase::which_phase
RtLayout<T>::getPhase()
{
  CH_TIME("RtLayout<T>::getPhase");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::getPhase" << endl;
  }

  return m_phase;
}

template <class T>
Vector<RefCountedPtr<T>>&
RtLayout<T>::getSolvers()
{
  return m_solvers;
}

template <class T>
Vector<RefCountedPtr<RtSpecies>>&
RtLayout<T>::getSpecies()
{
  return m_species;
}

template <class T>
Vector<EBAMRCellData*>
RtLayout<T>::getSources()
{
  CH_TIME("RtLayout<T>::getSources");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::getSources" << endl;
  }

  Vector<EBAMRCellData*> sources;

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    sources.push_back(&(solver->getSource()));
  }

  return sources;
}

template <class T>
Vector<EBAMRCellData*>
RtLayout<T>::getPhis()
{
  CH_TIME("RtLayout<T>::getPhis");
  if (m_verbosity > 5) {
    pout() << "RtLayout<T>::getPhis" << endl;
  }

  Vector<EBAMRCellData*> states;

  for (RtIterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<T>& solver = solver_it();
    states.push_back(&(solver->getPhi()));
  }

  return states;
}

template <class T, class S>
RtFactory<T, S>::RtFactory()
{}

template <class T, class S>
RtFactory<T, S>::~RtFactory()
{}

template <class T, class S>
RefCountedPtr<RtLayout<T>>
RtFactory<T, S>::newLayout(const Vector<RefCountedPtr<RtSpecies>>& a_species) const
{
  // Build rte layout
  auto rte = RefCountedPtr<RtLayout<T>>(new RtLayout<T>(a_species));
  auto spe = a_species;

  // Cast solvers and instantiate them
  for (int i = 0; i < a_species.size(); i++) {
    auto solver = RefCountedPtr<T>(static_cast<T*>(new S()));
    solver->setRtSpecies(spe[i]);
    solver->setPhase(phase::gas);
    solver->setVerbosity(-1);
    rte->addSolver(solver);
  }

  return rte;
}

#include <CD_NamespaceFooter.H>

#endif
