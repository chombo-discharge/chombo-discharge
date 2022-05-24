/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MeshODESolverImplem.H
  @brief  Implementation of CD_MeshODESolver.H
  @author Robert Marskar
*/

#ifndef CD_MeshODESolverImplem_H
#define CD_MeshODESolverImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MeshODESolver.H>

template <size_t N>
MeshODESolver<N>::MeshODESolver()
{
  CH_TIME("MeshODESolver::MeshODESolver");

  // Default settings. 
  m_verbosity = -1;
  m_name      = "MeshODESolver";
  m_className = "MeshODESolver";
  m_phase     = phase::gas;
  m_realm     = Realm::Primal;

  this->parseOptions();
}

template <size_t N>
MeshODESolver<N>::~MeshODESolver()
{
  CH_TIME("MeshODESolver::~MeshODESolver");
}

template <size_t N>
void
MeshODESolver<N>::parseOptions() noexcept
{
  CH_TIME("MeshODESolver::parseOptions");

  ParmParse pp(m_className.c_str());
  
  pp.get("verbosity", m_verbosity);
  pp.get("use_regrid_slopes", m_regridSlopes);
  
  if(m_verbosity > 5) {
    pout() << m_name + "::parseOptions()" << endl;
  }
}

template <size_t N>
void
MeshODESolver<N>::parseRuntimeOptions() noexcept
{
  CH_TIME("MeshODESolver::parseRuntimeOptions()");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions()" << endl;
  }
  
  ParmParse pp(m_className.c_str());
  
  pp.get("verbosity", m_verbosity);
  pp.get("use_regrid_slopes", m_regridSlopes);  
}

template <size_t N>
EBAMRCellData&
MeshODESolver<N>::getPhi() noexcept
{
  CH_TIME("MeshODESolver::getPhi()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getPhi()" << endl;
  }
  
  return m_phi;
}

template <size_t N>
const EBAMRCellData&
MeshODESolver<N>::getPhi() const noexcept
{
  CH_TIME("MeshODESolver::getPhi()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getPhi()" << endl;
  }
  
  return m_phi;
}

template <size_t N>
EBAMRCellData&
MeshODESolver<N>::getRHS() noexcept
{
  CH_TIME("MeshODESolver::getRHS()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getRHS()" << endl;
  }
  
  return m_RHS;
}

template <size_t N>
const EBAMRCellData&
MeshODESolver<N>::getRHS() const noexcept
{
  CH_TIME("MeshODESolver::getRHS()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getRHS()" << endl;
  }
  
  return m_RHS;
}

template <size_t N>
void
MeshODESolver<N>::setPhi(const std::function<Real(const RealVect& a_pos)>& a_phiFunc, const size_t a_comp) noexcept
{
  CH_TIME("MeshODESolver::setPhi(std::function, size_t)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setPhi(std::function, size_t)" << endl;
  }

  DataOps::setValue(m_phi, a_phiFunc, m_amr->getProbLo(), m_amr->getDx(), a_comp);
}

template <size_t N>
void
MeshODESolver<N>::setRHS(const std::function<Real(const RealVect& a_pos)>& a_srcFunc, const size_t a_comp) noexcept
{
  CH_TIME("MeshODESolver::setRHS(std::function, size_t)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setRHS(std::function, size_t)" << endl;
  }

  DataOps::setValue(m_rhs, a_phiFunc, m_amr->getProbLo(), m_amr->getDx(), a_comp);  
}

template <size_t N>
void
MeshODESolver<N>::preRegrid(const int a_lbase, const int a_oldFinestLevel) noexcept {
  CH_TIME("MeshODESolver::preRegrid(int, int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::preRegrid(int, int)" << endl;
  }

  m_amr->allocate(m_cache, m_realm, m_phase, N);

  m_cache.copy(m_phi);
}

template <size_t N>
void
MeshODESolver<N>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept
{
  CH_TIME("MeshODESolver::regrid(int, int, int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::regrid(int, int, int)" << endl;
  }

  this->allocateInternals();

  m_amr->interpToNewGrids(m_phi, m_cache, m_phase, a_lmin, a_oldFinestLevel, a_newFinestLevel, m_regridSlopes);

  m_amr->averageDown(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);  

  m_amr->deallocate(m_cache);
}

template <size_t N>
std::string
MeshODESolver<N>::getRealm() const noexcept
{
  CH_TIME("MeshODESolver::getRealm()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getRealm()" << endl;
  }

  return m_realm;
}

template <size_t N>
std::string
MeshODESolver<N>::getName() const noexcept
{
  CH_TIME("MeshODESolver::getName()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getName()" << endl;
  }

  return m_name;
}

template <size_t N>
Vector<std::string>
MeshODESolver<N>::getPlotVariableNames() const noexcept
{
  CH_TIME("MeshODESolver::getPlotVariableNames()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getPlotVariableNames()" << endl;
  }

  Vector<std::string> plotVarNames(0);

  if(m_plotPhi) {
    for (size_t i = 0; i < N; i++) {
      plotVarNames.push_back(m_name + " phi-" + std::to_string(i));
    }
  }

  if(m_plotSource) {
    for (size_t i = 0; i < N; i++) {
      plotVarNames.push_back(m_name + " source-" + std::to_string(i));
    }
  }

  return plotVarNames;
}

template <size_t N>
int
MeshODESolver<N>::getNumberOfPlotVariables() const noexcept
{
  CH_TIME("MeshODESolver::getNumberOfPlotVariables()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables()" << endl;
  }
  
  int numPlotVars = 0;

  if(m_plotPhi) {
    numPlotVars += int(N);
  }
  
  if(m_plotSource) {
    numPlotVars += int(N);
  }

  return numPlotVars;
}

template <size_t N>
void
MeshODESolver<N>::registerOperators() const noexcept
{
  CH_TIME("MeshODESolver::registerOperators()");
  if(m_verbosity > 5) {
    pout() << m_name + "::registerOperators()" << endl;
  }
  
  CH_assert(!m_amr.isNull());

  m_amr->registerOperator(s_eb_coar_ave, m_realm, m_phase);
  m_amr->registerOperator(s_eb_fill_patch, m_realm, m_phase);
  m_amr->registerOperator(s_eb_fine_interp, m_realm, m_phase);    
}

template <size_t N>
void
MeshODESolver<N>::setRealm(const std::string a_realm) noexcept
{
  CH_TIME("MeshODESolver::setRealm(std::string)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setRealm(std::string)" << endl;
  }

  m_realm.assign(a_realm);  
}

template <size_t N>
void
MeshODESolver<N>::setPhase(phase::which_phase a_phase) noexcept
{
  CH_TIME("MeshODESolver::setPhase(phase)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setPhase(phase)" << endl;
  }

  m_phase = a_phase;  
}

template <size_t N>
void
MeshODESolver<N>::setVerbosity(const int a_verbosity) noexcept
{
  CH_TIME("MeshODESolver::setVerbosity(int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setVerbosity(int)" << endl;
  }

  m_verbosity = a_verbosity;
}

template <size_t N>
void
MeshODESolver<N>::setTime(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("MeshODESolver::setTime(int, Real, Real)");
  if(m_verbosity > 5) {
    pout() << m_name + "::setTime(int, Real, Real)" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;  
}

template <size_t N>
void
MeshODESolver<N>::writePlotFile() noexcept
{
  CH_TIME("MeshODESolver::writePlotFile()");
  if(m_verbosity > 5) {
    pout() << m_name + "::writePlotFile()" << endl;
  }

  // Number of output components and their names
  const int                 numPlotVars  = this->getNumberOfPlotVariables();
  const Vector<std::string> plotVarNames = this->getPlotVariableNames();

  // Allocate storage
  EBAMRCellData output;
  m_amr->allocate(output, m_realm, m_phase, numPlotVars);
  DataOps::setValue(output, 0.0);

  int icomp = 0;
  this->writePlotData(output, icomp);

  char filename[100];
  sprintf(filename, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_timeStep, SpaceDim);
  std::string fname(filename);

  // Alias, because Chombo's EBAMRIO wants raw pointers (but we stick to smart pointers, like God intended).
  Vector<LevelData<EBCellFAB>*> outputPtr;
  m_amr->alias(outputPtr, output);

#ifdef CH_USE_HDF5
  constexpr int numPlotGhost = 0;

  DischargeIO::writeEBHDF5(fname,
                           plotVarNames,
                           m_amr->getGrids(m_realm),
                           outputPtr,
                           m_amr->getDomains(),
                           m_amr->getDx(),
                           m_amr->getRefinementRatios(),
                           m_dt,
                           m_time,
                           m_amr->getProbLo(),
                           m_amr->getFinestLevel() + 1,
                           numPlotGhost);
#endif  
}

template <size_t N>
void
MeshODESolver<N>::writePlotData(EBAMRCellData& a_output, int& a_icomp) noexcept
{
  CH_TIME("MeshODESolver::writePlotData(EBAMRCellData, int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::writePlotData(EBAMRCellData, int)" << endl;
  }

  if(m_plotPhi) {
    this->writeData(a_output, a_icomp, m_phi, true);
  }

  if(m_plotSource) {
    this->writeData(a_output, a_icomp, m_source, true);
  }
}

#ifdef CH_USE_HDF5
template <size_t N>
void
MeshODESolver<N>::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const noexcept
{
  CH_TIME("MeshODESolver::writeCheckpointLevel(HDF5Handle, int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel(HDF5Handle, int)" << endl;
  }

  write(a_handle, *m_phi[a_level], m_name);
  write(a_handle, *m_source[a_level], m_name + "_src");
}
#endif

#ifdef CH_USE_HDF5
template <size_t N>
void
MeshODESolver<N>::readCheckpointLevel(HDF5Handle& a_handle, const int a_level) noexcept
{
  CH_TIME("MeshODESolver::writeCheckpointLevel(HDF5Handle, int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel(HDF5Handle, int)" << endl;
  }

  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name, m_amr->getGrids(m_realm)[a_level], interv, false);
  read<EBCellFAB>(a_handle, *m_source[a_level], m_name + "_src", m_amr->getGrids(m_realm)[a_level], interv, false);
}
#endif

template <size_t N>
void
MeshODESolver::writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp)
{
  CH_TIME("MeshODESolver::writeData(EBAMRCellData, int, EBAMRCellData, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeData(EBAMRCellData, int, EBAMRCellData, bool)" << endl;
  }

  // Number of components we are working with.
  const int ncomp = a_data[0]->nComp();

  // Component ranges.
  const Interval srcInterv(0, ncomp - 1);
  const Interval dstInterv(a_comp, a_comp + ncomp - 1);

  // Copy data to scratch and interpolate scratch to cell centroids if we are asked to.
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, ncomp);
  DataOps::copy(scratch, a_data);

  if (a_interp) {
    m_amr->interpToCentroids(scratch, m_realm, phase::gas);
  }

  m_amr->averageDown(scratch, m_realm, m_phase);
  m_amr->interpGhost(scratch, m_realm, m_phase);

  // Need a more general copy method because we can't call DataOps::copy (because realms might not be the same) and
  // we can't call EBAMRData<T>::copy either (because components don't align). So -- general type of copy here.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    if (m_realm == a_output.getRealm()) {
      scratch[lvl]->localCopyTo(srcInterv, *a_output[lvl], dstInterv);
    }
    else {
      scratch[lvl]->copyTo(srcInterv, *a_output[lvl], dstInterv);
    }
  }

  DataOps::setCoveredValue(a_output, a_comp, 0.0);

  a_comp += ncomp;
}

#include <CD_NamespaceFooter.H>

#endif
