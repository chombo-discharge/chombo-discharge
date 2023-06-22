/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoLayoutImplem.H
  @brief  Implementation of CD_ItoLayout.H
  @author Robert Marskar
*/

#ifndef CD_ItoLayoutImplem_H
#define CD_ItoLayoutImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_ItoLayout.H>
#include <CD_ItoIterator.H>
#include <CD_NamespaceHeader.H>

template <class T>
ItoLayout<T>::ItoLayout(const Vector<RefCountedPtr<ItoSpecies>>& a_species)
{
  this->define(a_species);
  m_solvers.resize(0);
}

template <class T>
ItoLayout<T>::~ItoLayout()
{}

template <class T>
ItoIterator<T>
ItoLayout<T>::iterator()
{
  return ItoIterator<T>(*this);
}

template <class T>
const std::string
ItoLayout<T>::getRealm() const
{
  return m_realm;
}

template <class T>
void
ItoLayout<T>::setRealm(const std::string a_realm)
{
  m_realm = a_realm;

  for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->setRealm(m_realm);
  }
}

template <class T>
void
ItoLayout<T>::define(const Vector<RefCountedPtr<ItoSpecies>>& a_species)
{
  m_isDefined = true;
}

template <class T>
void
ItoLayout<T>::parseOptions()
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->parseOptions();
  }
}

template <class T>
void
ItoLayout<T>::parseRuntimeOptions()
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->parseRuntimeOptions();
  }
}

template <class T>
void
ItoLayout<T>::allocate()
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->allocate();
  }
}

template <class T>
void
ItoLayout<T>::addSolver(RefCountedPtr<T>& a_solver)
{
  m_solvers.push_back(a_solver);
}

template <class T>
void
ItoLayout<T>::interpolateVelocities()
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->interpolateVelocities();
  }
}

template <class T>
void
ItoLayout<T>::preRegrid(const int a_lbase, const int a_finestLevel)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->preRegrid(a_lbase, a_finestLevel);
  }
}

template <class T>
void
ItoLayout<T>::initialData()
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->initialData();
  }
}

template <class T>
void
ItoLayout<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
}

template <class T>
void
ItoLayout<T>::registerOperators()
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->registerOperators();
  }
}

template <class T>
void
ItoLayout<T>::setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->setAmr(a_amrMesh);
  }
}

template <class T>
void
ItoLayout<T>::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->setComputationalGeometry(a_computationalGeometry);
  }
}

template <class T>
void
ItoLayout<T>::setPhase(phase::which_phase a_phase)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->setPhase(a_phase);
  }
}

template <class T>
void
ItoLayout<T>::setVerbosity(const int a_verbosity)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->setVerbosity(a_verbosity);
  }
}

template <class T>
void
ItoLayout<T>::setTime(const int a_step, const Real a_time, const Real a_dt)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->setTime(a_step, a_time, a_dt);
  }
}

template <class T>
void
ItoLayout<T>::organizeParticlesByCell(const ItoSolver::WhichContainer a_whichContainer)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->organizeParticlesByCell(a_whichContainer);
  }
}

template <class T>
void
ItoLayout<T>::organizeParticlesByPatch(const ItoSolver::WhichContainer a_whichContainer)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->organizeParticlesByPatch(a_whichContainer);
  }
}

template <class T>
void
ItoLayout<T>::makeSuperparticles(const ItoSolver::WhichContainer a_whichContainer, const int a_ppc)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->makeSuperparticles(a_whichContainer, a_ppc);
  }
}

template <class T>
void
ItoLayout<T>::makeSuperparticles(const ItoSolver::WhichContainer a_whichContainer, const Vector<int> a_ppc)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->makeSuperparticles(a_whichContainer, a_ppc);
  }
}

template <class T>
void
ItoLayout<T>::depositParticles()
{
  this->depositParticles(ItoSolver::WhichContainer::Bulk);
}

template <class T>
void
ItoLayout<T>::depositParticles(const ItoSolver::WhichContainer a_whichContainer)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->depositParticles(a_whichContainer);
  }
}

template <class T>
void
ItoLayout<T>::remap()
{
  this->remap(ItoSolver::WhichContainer::Bulk);
}

template <class T>
void
ItoLayout<T>::remap(const ItoSolver::WhichContainer a_whichContainer)
{
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    iter()->remap(a_whichContainer);
  }
}

template <class T>
size_t
ItoLayout<T>::getNumParticles(const ItoSolver::WhichContainer a_whichContainer, const bool a_localOnly)
{
  size_t ret = 0;
  for (auto iter = this->iterator(); iter.ok(); ++iter) {
    ret += iter()->getNumParticles(a_whichContainer, a_localOnly);
  }

  return ret;
}

template <class T>
Real
ItoLayout<T>::computeDt()
{
  Real minDt = std::numeric_limits<Real>::max();

  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    const Real thisDt = iter()->computeDt();
    minDt             = std::min(minDt, thisDt);
  }

  return minDt;
}

template <class T>
Real
ItoLayout<T>::computeAdvectiveDt()
{
  Real minDt = std::numeric_limits<Real>::max();

  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    const Real thisDt = iter()->computeAdvectiveDt();
    minDt             = std::min(minDt, thisDt);
  }

  return minDt;
}

template <class T>
Real
ItoLayout<T>::computeHopDt(const Real a_maxCellsToMove)
{
  Real minDt = std::numeric_limits<Real>::max();

  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    const Real thisDt = iter()->computeHopDt(a_maxCellsToMove);
    minDt             = std::min(minDt, thisDt);
  }

  return minDt;
}

template <class T>
Real
ItoLayout<T>::computeDiffusiveDt()
{
  Real minDt = std::numeric_limits<Real>::max();

  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    const Real thisDt = iter()->computeDiffusiveDt();
    minDt             = std::min(minDt, thisDt);
  }

  return minDt;
}

template <class T>
Vector<RefCountedPtr<T>>&
ItoLayout<T>::getSolvers()
{
  return m_solvers;
}

template <class T>
Vector<RefCountedPtr<ItoSpecies>>&
ItoLayout<T>::getSpecies()
{
  return m_species;
}

template <class T>
Vector<EBAMRCellData*>
ItoLayout<T>::getVelocityFunctions()
{
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    ret[iter.index()] = &(iter()->getVelocityFunction());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*>
ItoLayout<T>::getDensities()
{
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    ret[iter.index()] = &(iter()->getPhi());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*>
ItoLayout<T>::getDiffusionFunctions()
{
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    ret[iter.index()] = &(iter()->getDiffusionFunction());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*>
ItoLayout<T>::getMobilityFunctions()
{
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    ret[iter.index()] = &(iter()->getMobilityFunction());
  }

  return ret;
}

template <class T>
Vector<ParticleContainer<ItoParticle>*>
ItoLayout<T>::getParticles(const ItoSolver::WhichContainer a_container)
{
  Vector<ParticleContainer<ItoParticle>*> ret(m_solvers.size(), nullptr);
  for (ItoIterator<T> iter = this->iterator(); iter.ok(); ++iter) {
    ret[iter.index()] = &(iter()->getParticles(a_container));
  }

  return ret;
}

template <class T>
phase::which_phase
ItoLayout<T>::getPhase() const
{
  return m_phase;
}

template <class T, class S>
ItoFactory<T, S>::ItoFactory()
{}

template <class T, class S>
ItoFactory<T, S>::~ItoFactory()
{}

template <class T, class S>
RefCountedPtr<ItoLayout<T>>
ItoFactory<T, S>::newLayout(const Vector<RefCountedPtr<ItoSpecies>>& a_species) const
{

  auto itoLayout = RefCountedPtr<ItoLayout<T>>(new ItoLayout<T>(a_species));

  // Cast solvers and instantiate them
  for (int i = 0; i < a_species.size(); i++) {
    RefCountedPtr<T> solver = RefCountedPtr<T>(static_cast<T*>(new S()));
    solver->setSpecies(a_species[i]);
    solver->setVerbosity(-1);
    itoLayout->addSolver(solver);
  }

  return itoLayout;
}

#include <CD_NamespaceFooter.H>

#endif
