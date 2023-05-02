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
#include <ParmParse.H>

// Our includes
#include <CD_MeshODESolver.H>
#include <CD_DischargeIO.H>
#include <CD_BoxLoops.H>
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

template <size_t N>
MeshODESolver<N>::MeshODESolver()
{
  CH_TIME("MeshODESolver::MeshODESolver()");

  // Default settings.
  m_verbosity = -1;
  m_name      = "MeshODESolver";
  m_className = "MeshODESolver";
  m_phase     = phase::gas;
  m_realm     = Realm::Primal;

  this->parseOptions();
}

template <size_t N>
MeshODESolver<N>::MeshODESolver(const RefCountedPtr<AmrMesh>& a_amr) noexcept : MeshODESolver<N>()
{
  CH_TIME("MeshODESolver::MeshODESolver(AMR)");

  m_amr = a_amr;
}

template <size_t N>
MeshODESolver<N>::~MeshODESolver()
{
  CH_TIME("MeshODESolver::~MeshODESolver");
}

template <size_t N>
void
MeshODESolver<N>::setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh) noexcept
{
  CH_TIME("MeshODESolver::setAmr");
  if (m_verbosity > 5) {
    pout() << m_name + "::setAmr" << endl;
  }

  m_amr = a_amrMesh;
}

template <size_t N>
void
MeshODESolver<N>::parseOptions() noexcept
{
  CH_TIME("MeshODESolver::parseOptions");

  ParmParse pp(m_className.c_str());

  pp.get("verbosity", m_verbosity);
  pp.get("use_regrid_slopes", m_regridSlopes);

  this->parsePlotVariables();

  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions()" << endl;
  }
}

template <size_t N>
void
MeshODESolver<N>::parseRuntimeOptions() noexcept
{
  CH_TIME("MeshODESolver::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("verbosity", m_verbosity);
  pp.get("use_regrid_slopes", m_regridSlopes);

  this->parsePlotVariables();
}

template <size_t N>
void
MeshODESolver<N>::parsePlotVariables() noexcept
{
  CH_TIME("MeshODESolver::parsePlotVariables()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePlotVariables()" << endl;
  }

  m_plotPhi = false;
  m_plotRHS = false;

  ParmParse pp(m_className.c_str());
  const int num = pp.countval("plt_vars");

  if (num > 0) {
    Vector<std::string> str(num);
    pp.getarr("plt_vars", str, 0, num);

    for (int i = 0; i < num; i++) {
      if (str[i] == "phi") {
        m_plotPhi = true;
      }
      else if (str[i] == "rhs") {
        m_plotRHS = true;
      }
    }
  }
}

template <size_t N>
EBAMRCellData&
MeshODESolver<N>::getPhi() noexcept
{
  CH_TIME("MeshODESolver::getPhi()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPhi()" << endl;
  }

  return m_phi;
}

template <size_t N>
const EBAMRCellData&
MeshODESolver<N>::getPhi() const noexcept
{
  CH_TIME("MeshODESolver::getPhi()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPhi()" << endl;
  }

  return m_phi;
}

template <size_t N>
EBAMRCellData&
MeshODESolver<N>::getRHS() noexcept
{
  CH_TIME("MeshODESolver::getRHS()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getRHS()" << endl;
  }

  return m_rhs;
}

template <size_t N>
const EBAMRCellData&
MeshODESolver<N>::getRHS() const noexcept
{
  CH_TIME("MeshODESolver::getRHS()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getRHS()" << endl;
  }

  return m_rhs;
}

template <size_t N>
void
MeshODESolver<N>::setPhi(const std::function<Real(const RealVect& a_pos)>& a_phiFunc, const size_t a_comp) noexcept
{
  CH_TIME("MeshODESolver::setPhi(std::function, size_t)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhi(std::function, size_t)" << endl;
  }

  DataOps::setValue(m_phi, a_phiFunc, m_amr->getProbLo(), m_amr->getDx(), a_comp);
}

template <size_t N>
void
MeshODESolver<N>::setPhi(const std::function<std::array<Real, N>(const RealVect& a_pos)>& a_phiFunc) noexcept
{
  CH_TIME("MeshODESolver::setPhi(std::function)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhi(std::function)" << endl;
  }

  constexpr int comp = 0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const Real&              dx     = m_amr->getDx()[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const RealVect&          probLo = m_amr->getProbLo();

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      EBCellFAB&           phi        = (*m_phi[lvl])[dit()];
      const EBISBox&       ebisbox    = ebisl[dit()];
      FArrayBox&           phiFAB     = phi.getFArrayBox();
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];

      // Regular kernel. Only do cells that are not covered by a finer grid.
      auto regularKernel = [&](const IntVect& iv) -> void {
        std::array<Real, N> y{};

        if (validCells(iv, comp) && ebisbox.isRegular(iv)) {
          const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;

          y = a_phiFunc(pos);
        }

        // Put right-hand side in mesh data. This is set to zero if the cell is covered by a finer grid cell.
        for (size_t i = 0; i < N; i++) {
          phiFAB(iv, i) = y[i];
        }
      };

      // Irregular kernel. Only do cells that are not covered by a finer grid.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        std::array<Real, N> y{};
        const IntVect       iv = vof.gridIndex();

        if (validCells(iv, comp)) {
          const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);

          y = a_phiFunc(pos);
        }

        // Put right-hand side in mesh data. This is set to zero if the cell is covered by a finer grid cell.
        for (size_t i = 0; i < N; i++) {
          phi(vof, i) = y[i];
        }
      };

      // Kernel regions.
      const Box    cellBox = dbl[dit()];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      // Run the kernels.
      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

template <size_t N>
void
MeshODESolver<N>::setRHS(const std::function<Real(const RealVect& a_pos)>& a_srcFunc, const size_t a_comp) noexcept
{
  CH_TIME("MeshODESolver::setRHS(std::function, size_t)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setRHS(std::function, size_t)" << endl;
  }

  DataOps::setValue(m_rhs, a_srcFunc, m_amr->getProbLo(), m_amr->getDx(), a_comp);
}

template <size_t N>
void
MeshODESolver<N>::computeRHS(const RHSFunction& a_rhsFunction) noexcept
{
  CH_TIME("MeshODESolver::computeRHS(std::function<std::array<Real, N>(...)>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeRHS(std::function<std::array<Real, N>(...)>)" << endl;
  }

  this->computeRHS(m_rhs, a_rhsFunction);
}

template <size_t N>
void
MeshODESolver<N>::computeRHS(EBAMRCellData& a_rhs, const RHSFunction& a_rhsFunction) const noexcept
{
  CH_TIME("MeshODESolver::computeRHS(EBAMRCellData, std::function<std::array<Real, N>(...)>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeRHS(EBAMRCellData, std::function<std::array<Real, N>(...)>)" << endl;
  }

  constexpr int comp = 0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      EBCellFAB&       rhs     = (*a_rhs[lvl])[dit()];
      const EBCellFAB& phi     = (*m_phi[lvl])[dit()];
      const EBISBox&   ebisbox = ebisl[dit()];

      FArrayBox&       rhsFAB = rhs.getFArrayBox();
      const FArrayBox& phiFAB = phi.getFArrayBox();

      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];

      // Regular kernel. Only do cells that are not covered by a finer grid.
      auto regularKernel = [&](const IntVect& iv) -> void {
        std::array<Real, N> y{};
        std::array<Real, N> f{};

        if (validCells(iv, comp) && ebisbox.isRegular(iv)) {

          // Make y from the mesh data.
          for (size_t i = 0; i < N; i++) {
            y[i] = phiFAB(iv, i);
          }

          // Compute the right-hand side.
          f = a_rhsFunction(y, m_time);
        }

        // Put right-hand side in mesh data. This is set to zero if the cell is covered by a finer grid cell.
        for (size_t i = 0; i < N; i++) {
          rhsFAB(iv, i) = f[i];
        }
      };

      // Irregular kernel. Only do cells that are not covered by a finer grid.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        std::array<Real, N> y{};
        std::array<Real, N> f{};

        const IntVect iv = vof.gridIndex();

        if (validCells(iv, comp)) {

          // Make y from the mesh data.
          for (size_t i = 0; i < N; i++) {
            y[i] = phi(vof, i);
          }

          // Compute the right-hand side.
          f = a_rhsFunction(y, m_time);
        }

        // Put right-hand side in mesh data. This is set to zero if the cell is covered by a finer grid cell.
        for (size_t i = 0; i < N; i++) {
          rhs(vof, i) = f[i];
        }
      };

      // Kernel regions.
      const Box    cellBox = dbl[dit()];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      // Run the kernels.
      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

template <size_t N>
void
MeshODESolver<N>::allocate() noexcept
{
  CH_TIME("MeshODESolver::allocate()");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate()" << endl;
  }

  m_amr->allocate(m_phi, m_realm, m_phase, N);
  m_amr->allocate(m_rhs, m_realm, m_phase, N);
}

template <size_t N>
void
MeshODESolver<N>::preRegrid(const int a_lbase, const int a_oldFinestLevel) noexcept
{
  CH_TIME("MeshODESolver::preRegrid(int, int)");
  if (m_verbosity > 5) {
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
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid(int, int, int)" << endl;
  }

  this->allocate();

  const EBCoarseToFineInterp::Type interpType = m_regridSlopes ? EBCoarseToFineInterp::Type::ConservativeMinMod
                                                               : EBCoarseToFineInterp::Type::ConservativePWC;

  m_amr->interpToNewGrids(m_phi, m_cache, m_phase, a_lmin, a_oldFinestLevel, a_newFinestLevel, interpType);

  m_amr->conservativeAverage(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);

  m_amr->deallocate(m_cache);
}

template <size_t N>
std::string
MeshODESolver<N>::getRealm() const noexcept
{
  CH_TIME("MeshODESolver::getRealm()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getRealm()" << endl;
  }

  return m_realm;
}

template <size_t N>
std::string
MeshODESolver<N>::getName() const noexcept
{
  CH_TIME("MeshODESolver::getName()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getName()" << endl;
  }

  return m_name;
}

template <size_t N>
void
MeshODESolver<N>::setName(const std::string& a_name) noexcept
{
  CH_TIME("MeshODESolver::setName()");
  if (m_verbosity > 5) {
    pout() << m_name + "::setName()" << endl;
  }

  m_name = a_name;
}

template <size_t N>
void
MeshODESolver<N>::registerOperators() const noexcept
{
  CH_TIME("MeshODESolver::registerOperators()");
  if (m_verbosity > 5) {
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
  if (m_verbosity > 5) {
    pout() << m_name + "::setRealm(std::string)" << endl;
  }

  m_realm.assign(a_realm);
}

template <size_t N>
void
MeshODESolver<N>::setPhase(phase::which_phase a_phase) noexcept
{
  CH_TIME("MeshODESolver::setPhase(phase)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhase(phase)" << endl;
  }

  m_phase = a_phase;
}

template <size_t N>
void
MeshODESolver<N>::setVerbosity(const int a_verbosity) noexcept
{
  CH_TIME("MeshODESolver::setVerbosity(int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setVerbosity(int)" << endl;
  }

  m_verbosity = a_verbosity;
}

template <size_t N>
void
MeshODESolver<N>::setTime(const int a_step, const Real a_time, const Real a_dt) noexcept
{
  CH_TIME("MeshODESolver::setTime(int, Real, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setTime(int, Real, Real)" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;
}

template <size_t N>
void
MeshODESolver<N>::writePlotFile() const noexcept
{
  CH_TIME("MeshODESolver::writePlotFile()");
  if (m_verbosity > 5) {
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
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->writePlotData(*output[lvl], icomp, lvl);
  }

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
int
MeshODESolver<N>::getNumberOfPlotVariables() const noexcept
{
  CH_TIME("MeshODESolver::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables()" << endl;
  }

  int numPlotVars = 0;

  if (m_plotPhi) {
    numPlotVars += int(N);
  }

  if (m_plotRHS) {
    numPlotVars += int(N);
  }

  return numPlotVars;
}

template <size_t N>
Vector<std::string>
MeshODESolver<N>::getPlotVariableNames() const noexcept
{
  CH_TIME("MeshODESolver::getPlotVariableNames()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPlotVariableNames()" << endl;
  }

  Vector<std::string> plotVarNames(0);

  if (m_plotPhi) {
    for (size_t i = 0; i < N; i++) {
      plotVarNames.push_back(m_name + " phi-" + std::to_string(i));
    }
  }

  if (m_plotRHS) {
    for (size_t i = 0; i < N; i++) {
      plotVarNames.push_back(m_name + " source-" + std::to_string(i));
    }
  }

  return plotVarNames;
}

template <size_t N>
void
MeshODESolver<N>::writePlotData(LevelData<EBCellFAB>& a_output, int& a_icomp, const int a_level) const noexcept
{
  CH_TIME("MeshODESolver::writePlotData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  if (m_plotPhi) {
    this->writeData(a_output, a_icomp, m_phi, a_level, false, true);
  }

  if (m_plotRHS) {
    this->writeData(a_output, a_icomp, m_rhs, a_level, false, true);
  }
}

template <size_t N>
void
MeshODESolver<N>::writeData(LevelData<EBCellFAB>& a_output,
                            int&                  a_comp,
                            const EBAMRCellData&  a_data,
                            const int             a_level,
                            const bool            a_interpToCentroids,
                            const bool            a_interpGhost) const noexcept

{
  CH_TIMERS("MeshODESolver::writeData");
  CH_TIMER("MeshODESolver::writeData::allocate", t1);
  CH_TIMER("MeshODESolver::writeData::local_copy", t2);
  CH_TIMER("MeshODESolver::writeData::interp_ghost", t3);
  CH_TIMER("MeshODESolver::writeData::interp_centroid", t4);
  CH_TIMER("MeshODESolver::writeData::define_copier", t5);
  CH_TIMER("MeshODESolver::writeData::final_copy", t6);
  if (m_verbosity > 5) {
    pout() << m_name + "::writeData" << endl;
  }

  // Number of components we are working with.
  const int numComp = a_data[a_level]->nComp();

  // Component ranges that we copy to/from.
  const Interval srcInterv(0, numComp - 1);
  const Interval dstInterv(a_comp, a_comp + numComp - 1);

  CH_START(t1);
  LevelData<EBCellFAB> scratch;
  m_amr->allocate(scratch, m_realm, m_phase, a_level, numComp);
  CH_STOP(t1);

  CH_START(t2);
  a_data[a_level]->localCopyTo(scratch);
  scratch.exchange();
  CH_START(t2);

  // Interpolate ghost cells
  CH_START(t3);
  if (a_level > 0 && a_interpGhost) {
    for (int icomp = 0; icomp < numComp; icomp++) {
      LevelData<EBCellFAB> coarData;

      aliasLevelData<EBCellFAB>(coarData, &(*a_data[a_level - 1]), Interval(icomp, icomp));

      m_amr->interpGhost(scratch, coarData, a_level, m_realm, m_phase);
    }
  }
  CH_STOP(t3);

  CH_START(t4);
  if (a_interpToCentroids) {
    m_amr->interpToCentroids(scratch, m_realm, m_phase, a_level);
  }
  CH_STOP(t4);

  // Need a more general copy method because we can't call DataOps::copy (because realms might not be the same) and
  // we can't call EBAMRData<T>::copy either (because components don't align). So -- general type of copy here.
  CH_START(t5);
  Copier copier;
  copier.ghostDefine(scratch.disjointBoxLayout(),
                     a_output.disjointBoxLayout(),
                     m_amr->getDomains()[a_level],
                     scratch.ghostVect(),
                     a_output.ghostVect());
  CH_STOP(t5);

  DataOps::setCoveredValue(scratch, 0.0);

  CH_START(t6);
  scratch.copyTo(srcInterv, a_output, dstInterv, copier);
  a_output.exchange();
  CH_STOP(t6);

  a_comp += numComp;
}

#ifdef CH_USE_HDF5
template <size_t N>
void
MeshODESolver<N>::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const noexcept
{
  CH_TIME("MeshODESolver::writeCheckpointLevel(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel(HDF5Handle, int)" << endl;
  }

  write(a_handle, *m_phi[a_level], m_name);
  write(a_handle, *m_rhs[a_level], m_name + "_src");
}
#endif

#ifdef CH_USE_HDF5
template <size_t N>
void
MeshODESolver<N>::readCheckpointLevel(HDF5Handle& a_handle, const int a_level) noexcept
{
  CH_TIME("MeshODESolver::writeCheckpointLevel(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel(HDF5Handle, int)" << endl;
  }

  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name, m_amr->getGrids(m_realm)[a_level], Interval(0, N - 1), false);
  read<EBCellFAB>(a_handle,
                  *m_rhs[a_level],
                  m_name + "_src",
                  m_amr->getGrids(m_realm)[a_level],
                  Interval(0, N - 1),
                  false);
}
#endif

#include <CD_NamespaceFooter.H>

#endif
