/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SurfaceODESolverImplem.H
  @brief  Implementation of CD_SurfaceODESolver.H
  @author Robert Marskar
*/

#ifndef CD_SurfaceODESolverImplem_H
#define CD_SurfaceODESolverImplem_H

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_SurfaceODESolver.H>
#include <CD_DischargeIO.H>
#include <CD_NamespaceHeader.H>

template <int N>
SurfaceODESolver<N>::SurfaceODESolver()
{
  CH_TIME("SurfaceODESolver::SurfaceODESolver");

  m_className = "SurfaceODESolver";

  m_verbosity = -1;

  this->setVerbosity(-1);
  this->setRealm(Realm::Primal);
  this->setPhase(phase::gas);
  this->setName(m_className);
  this->setTime(0, 0.0, 0.0);
}

template <int N>
SurfaceODESolver<N>::SurfaceODESolver(const RefCountedPtr<AmrMesh>& a_amr) : SurfaceODESolver<N>()
{
  CH_TIME("SurfaceODESolver::SurfaceODESolver");
  if (m_verbosity > 5) {
    pout() << m_name + "::SurfaceODESolver()" << endl;
  }

  m_amr = a_amr;
}

template <int N>
SurfaceODESolver<N>::~SurfaceODESolver()
{
  CH_TIME("SurfaceODESolver::~SurfaceODESolver");
  if (m_verbosity > 5) {
    pout() << m_name + "::~SurfaceODESolver()" << endl;
  }
}

template <int N>
void
SurfaceODESolver<N>::parseOptions() noexcept
{
  CH_TIME("SurfaceODESolver::parseOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions" << endl;
  }

  this->parseVerbosity();
  this->parseRegrid();
  this->parsePlotVariables();
}

template <int N>
void
SurfaceODESolver<N>::parseRuntimeOptions() noexcept
{
  CH_TIME("SurfaceODESolver::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  this->parseVerbosity();
  this->parseRegrid();
  this->parsePlotVariables();
}

template <int N>
void
SurfaceODESolver<N>::parseVerbosity() noexcept
{
  CH_TIME("SurfaceODESolver::parseVerbosity");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseVerbosity" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("verbosity", m_verbosity);
}

template <int N>
void
SurfaceODESolver<N>::parseRegrid() noexcept
{
  CH_TIME("SurfaceODESolver::parseRegrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRegrid" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("regrid", str);

  if (str == "conservative") {
    m_conservativeRegrid = true;
  }
  else if (str == "arithmetic") {
    m_conservativeRegrid = false;
  }
  else {
    const std::string err = "SurfaceODESolver<N>::parseRegrid() - argument '" + str + "' not recognized";
    MayDay::Error(err.c_str());
  }
}

template <int N>
void
SurfaceODESolver<N>::parsePlotVariables() noexcept
{
  CH_TIME("SurfaceODESolver::parsePlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePlotVariables" << endl;
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

template <int N>
void
SurfaceODESolver<N>::setAmr(const RefCountedPtr<AmrMesh>& a_amr) noexcept
{
  CH_TIME("SurfaceODESolver::setAmr");
  if (m_verbosity > 5) {
    pout() << m_name + "::setAmr" << endl;
  }

  CH_assert(!(a_amr.isNull()));

  m_amr = a_amr;
}

template <int N>
void
SurfaceODESolver<N>::setRealm(const std::string a_realm) noexcept
{
  CH_TIME("SurfaceODESolver::setRealm");
  if (m_verbosity > 5) {
    pout() << m_name + "::setRealm" << endl;
  }

  m_realm = a_realm;
}

template <int N>
void
SurfaceODESolver<N>::setName(const std::string a_name) noexcept
{
  CH_TIME("SurfaceODESolver::setName");
  if (m_verbosity > 5) {
    pout() << "SurfaceODESolver::setName" << endl;
  }

  m_name = a_name;
}

template <int N>
std::string
SurfaceODESolver<N>::getRealm() const noexcept
{
  CH_TIME("SurfaceODESolver::getRealm");
  if (m_verbosity > 5) {
    pout() << "SurfaceODESolver::getRealm" << endl;
  }

  return m_realm;
}

template <int N>
void
SurfaceODESolver<N>::setPhase(const phase::which_phase a_phase) noexcept
{
  CH_TIME("SurfaceODESolver::setPhase");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhase" << endl;
  }

  m_phase = a_phase;
}

template <int N>
phase::which_phase
SurfaceODESolver<N>::getPhase() const noexcept
{
  CH_TIME("SurfaceODESolver::getPhase");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPhase" << endl;
  }

  return m_phase;
}

template <int N>
void
SurfaceODESolver<N>::setVerbosity(const int a_verbosity) noexcept
{
  CH_TIME("SurfaceODESolver::setVerbosity");
  if (m_verbosity > 5) {
    pout() << m_name + "::setVerbosity" << endl;
  }

  m_verbosity = a_verbosity;
}

template <int N>
int
SurfaceODESolver<N>::getVerbosity() const noexcept
{
  CH_TIME("SurfaceODESolver::getVerbosity");
  if (m_verbosity > 5) {
    pout() << m_name + "::getVerbosity" << endl;
  }

  return m_verbosity;
}

template <int N>
void
SurfaceODESolver<N>::setTime(const int a_step, const Real a_time, const Real a_dt) noexcept
{
  CH_TIME("SurfaceODESolver::setTime");
  if (m_verbosity > 5) {
    pout() << m_name + "::setTime" << endl;
  }

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

template <int N>
Real
SurfaceODESolver<N>::computeMass(const int a_comp) const noexcept
{
  CH_TIME("SurfaceODESolver::computeMass(int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeMass(int)" << endl;
  }

  return this->computeMass(m_phi, a_comp);
}

template <int N>
Real
SurfaceODESolver<N>::computeMass(const EBAMRIVData& a_data, const int a_comp) const noexcept
{
  CH_TIME("SurfaceODESolver::computeMass(EBAMRIVData, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeMass(EBAMRIVData, int)" << endl;
  }

  Real dataSum = 0.0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx     = m_amr->getDx()[lvl];
    const Real               dxArea = std::pow(dx, SpaceDim - 1);

    CH_assert(a_data[lvl]->nComp() > a_comp);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox&         ebisbox    = ebisl[dit()];
      const BaseFab<bool>&   validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];
      const BaseIVFAB<Real>& data       = (*a_data[lvl])[dit()];

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();

        if (validCells(iv)) {
          dataSum += data(vof, a_comp) * ebisbox.bndryArea(vof) * dxArea;
        }
      };

      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  return ParallelOps::sum(dataSum);
}

template <int N>
int
SurfaceODESolver<N>::getStep() const noexcept
{
  CH_TIME("SurfaceODESolver::getStep");
  if (m_verbosity > 5) {
    pout() << m_name + "::getStep" << endl;
  }

  return m_step;
}

template <int N>
Real
SurfaceODESolver<N>::getTime() const noexcept
{
  CH_TIME("SurfaceODESolver::getTime");
  if (m_verbosity > 5) {
    pout() << m_name + "::getTime" << endl;
  }

  return m_time;
}

template <int N>
Real
SurfaceODESolver<N>::getDt() const noexcept
{
  CH_TIME("SurfaceODESolver::getDt");
  if (m_verbosity > 5) {
    pout() << m_name + "::getDt" << endl;
  }

  return m_dt;
}

template <int N>
void
SurfaceODESolver<N>::setPhi(const Real a_phi)
{
  CH_TIME("SurfaceODESolver::setPhi(Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhi(Real)" << endl;
  }

  DataOps::setValue(m_phi, a_phi);
}

template <int N>
void
SurfaceODESolver<N>::setPhi(const std::array<Real, N>& a_phi)
{
  CH_TIME("SurfaceODESolver::setPhi(std::array<Real, N>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhi(std::array<Real, N>)" << endl;
  }

  for (int i = 0; i < N; i++) {
    DataOps::setValue(m_phi, a_phi[i]);
  }
}

template <int N>
void
SurfaceODESolver<N>::setPhi(const EBAMRIVData& a_phi)
{
  CH_TIME("SurfaceODESolver::setPhi(EBAMRIVData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhi(EBAMRIVData)" << endl;
  }

  DataOps::copy(m_phi, a_phi);
}

template <int N>
EBAMRIVData&
SurfaceODESolver<N>::getPhi() noexcept
{
  CH_TIME("SurfaceODESolver::getPhi");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPhi" << endl;
  }

  return m_phi;
}

template <int N>
const EBAMRIVData&
SurfaceODESolver<N>::getPhi() const noexcept
{
  CH_TIME("SurfaceODESolver::getPhi");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPhi" << endl;
  }

  return m_phi;
}

template <int N>
void
SurfaceODESolver<N>::setRHS(const Real a_rhs)
{
  CH_TIME("SurfaceODESolver::setRHS(Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setRHS(Real)" << endl;
  }

  DataOps::setValue(m_rhs, a_rhs);
}

template <int N>
void
SurfaceODESolver<N>::setRHS(const std::array<Real, N>& a_rhs)
{
  CH_TIME("SurfaceODESolver::setRHS(std::array<Real, N>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setRHS(std::array<Real, N>)" << endl;
  }

  for (int i = 0; i < N; i++) {
    DataOps::setValue(m_rhs, a_rhs[i]);
  }
}

template <int N>
void
SurfaceODESolver<N>::setRHS(const EBAMRIVData& a_rhs)
{
  CH_TIME("SurfaceODESolver::setRHS(EBAMRIVData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setRHS(EBAMRIVData)" << endl;
  }

  DataOps::copy(m_rhs, a_rhs);
}

template <int N>
EBAMRIVData&
SurfaceODESolver<N>::getRHS()
{
  CH_TIME("SurfaceODESolver::getRHS");
  if (m_verbosity > 5) {
    pout() << m_name + "::getRHS" << endl;
  }

  return m_rhs;
}

template <int N>
const EBAMRIVData&
SurfaceODESolver<N>::getRHS() const noexcept
{
  CH_TIME("SurfaceODESolver::getRHS");
  if (m_verbosity > 5) {
    pout() << m_name + "::getRHS" << endl;
  }

  return m_rhs;
}

template <int N>
void
SurfaceODESolver<N>::allocate() noexcept
{
  CH_TIME("SurfaceODESolver::allocate");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate" << endl;
  }

  CH_assert(!(m_amr.isNull()));

  m_amr->allocate(m_phi, m_realm, m_phase, N);
  m_amr->allocate(m_rhs, m_realm, m_phase, N);

  this->defineVoFIterators();
}

template <int N>
void
SurfaceODESolver<N>::deallocate() noexcept
{
  CH_TIME("SurfaceODESolver::deallocate");
  if (m_verbosity > 5) {
    pout() << m_name + "::deallocate" << endl;
  }

  m_phi.clear();
  m_rhs.clear();
}

template <int N>
void
SurfaceODESolver<N>::defineVoFIterators() noexcept
{
  CH_TIME("SurfaceODESolver::defineVoFIterators");
  if (m_verbosity > 5) {
    pout() << m_name + "::defineVoFIterators" << endl;
  }

  using LDVoFs = LayoutData<VoFIterator>;

  CH_assert(!(m_amr.isNull()));

  const int finestLevel = m_amr->getFinestLevel();

  m_dielectricVoFs.resize(1 + finestLevel);
  m_electrodeVoFs.resize(1 + finestLevel);

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const MFLevelGrid&       mflg  = *m_amr->getMFLevelGrid(m_realm)[lvl];

    m_dielectricVoFs[lvl] = RefCountedPtr<LDVoFs>(new LDVoFs(dbl));
    m_electrodeVoFs[lvl]  = RefCountedPtr<LDVoFs>(new LDVoFs(dbl));

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box&     cellBox = dbl[dit()];
      const EBISBox& ebisBox = ebisl[dit()];
      const EBGraph& ebGraph = ebisBox.getEBGraph();

      VoFIterator& dielectricVoFIt = (*m_dielectricVoFs[lvl])[dit()];
      VoFIterator& electrodeVoFIt  = (*m_electrodeVoFs[lvl])[dit()];

      IntVectSet dielectricIVS = mflg.interfaceRegion(cellBox, dit());
      IntVectSet electrodeIVS  = ebisBox.getIrregIVS(cellBox) - dielectricIVS;

      dielectricVoFIt.define(dielectricIVS, ebGraph);
      electrodeVoFIt.define(electrodeIVS, ebGraph);
    }
  }
}

template <int N>
void
SurfaceODESolver<N>::registerOperators() noexcept
{
  CH_TIME("SurfaceODESolver::registerOperators");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerOperators" << endl;
  }

  CH_assert(!(m_amr.isNull()));

  m_amr->registerOperator(s_eb_coar_ave, m_realm, m_phase);
}

template <int N>
void
SurfaceODESolver<N>::preRegrid(const int a_lbase, const int a_oldFinestLevel) noexcept
{
  CH_TIME("SurfaceODESolver::registerOperators");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerOperators" << endl;
  }

  CH_assert(!(m_amr.isNull()));

  m_amr->allocate(m_cache, m_realm, m_phase, N);

  DataOps::copy(m_cache, m_phi);

  m_phi.clear();
}

template <int N>
void
SurfaceODESolver<N>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept
{
  CH_TIME("SurfaceODESolver::regrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerOperators" << endl;
  }

  this->allocate();

  const EBCoarseToFineInterp::Type interpType = m_conservativeRegrid ? EBCoarseToFineInterp::Type::ConservativePWC
                                                                     : EBCoarseToFineInterp::Type::PWC;

  m_amr->interpToNewGrids(m_phi, m_cache, m_phase, a_lmin, a_oldFinestLevel, a_newFinestLevel, interpType);

  m_cache.clear();
}

template <int N>
void
SurfaceODESolver<N>::resetElectrodes(const Real a_value) noexcept
{
  CH_TIME("SurfaceODESolver::resetElectrodes");
  if (m_verbosity > 5) {
    pout() << m_name + "::resetElectrodes" << endl;
  }

  this->resetElectrodes(m_phi, a_value);
}

template <int N>
void
SurfaceODESolver<N>::resetElectrodes(EBAMRIVData& a_phi, const Real a_value) const noexcept
{
  CH_TIME("SurfaceODESolver::resetElectrodes(EBAMRIVData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::resetElectrodes(EBAMRIVData)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      BaseIVFAB<Real>& phi   = (*a_phi[lvl])[dit()];
      VoFIterator&     vofit = (*m_electrodeVoFs[lvl])[dit()];

      auto kernel = [&](const VolIndex& vof) -> void {
        for (int comp = 0; comp < N; comp++) {
          phi(vof, comp) = a_value;
        }
      };

      BoxLoops::loop(vofit, kernel);
    }
  }
}

template <int N>
void
SurfaceODESolver<N>::resetDielectrics(const Real a_value) noexcept
{
  CH_TIME("SurfaceODESolver::resetDielectrics");
  if (m_verbosity > 5) {
    pout() << m_name + "::resetDielectrics" << endl;
  }

  this->resetDielectrics(m_phi, a_value);
}

template <int N>
void
SurfaceODESolver<N>::resetDielectrics(EBAMRIVData& a_phi, const Real a_value) const noexcept
{
  CH_TIME("SurfaceODESolver::resetDielectrics(EBAMRIVData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::resetDielectrics(EBAMRIVData)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      BaseIVFAB<Real>& phi   = (*a_phi[lvl])[dit()];
      VoFIterator&     vofit = (*m_dielectricVoFs[lvl])[dit()];

      auto kernel = [&](const VolIndex& vof) -> void {
        for (int comp = 0; comp < N; comp++) {
          phi(vof, comp) = a_value;
        }
      };

      BoxLoops::loop(vofit, kernel);
    }
  }
}

#ifdef CH_USE_HDF5
template <int N>
void
SurfaceODESolver<N>::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const noexcept
{
  CH_TIME("SurfaceODESolver::writeCheckpointLevel");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel" << endl;
  }

  const std::string id = "SurfaceODESolver::" + m_name + "::phi";

  CH_assert(!(m_amr.isNull()));

  // HDF5 doesn't like irregular data so we copy the data onto an EBCellFAB, which also has storage for multi-cells
  EBAMRCellData scratch;

  m_amr->allocate(scratch, m_realm, m_phase, N);

  DataOps::setValue(*scratch[a_level], 0.0);
  DataOps::incr(*scratch[a_level], *m_phi[a_level], 1.0);

  write(a_handle, *scratch[a_level], id);
}
#endif

#ifdef CH_USE_HDF5

template <int N>
void
SurfaceODESolver<N>::readCheckpointLevel(HDF5Handle& a_handle, const int a_level) noexcept
{
  CH_TIME("SurfaceODESolver::readCheckpointLevel");
  if (m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevel" << endl;
  }

  const std::string id = "SurfaceODESolver::" + m_name + "::phi";

  CH_assert(!(m_amr.isNull()));

  // HDF5 doesn't like irregular data so we wrote the data onto an EBCellFAB, which also has storage for multi-cells. We
  // also need to read from that...
  EBAMRCellData scratch;

  m_amr->allocate(scratch, m_realm, m_phase, N);

  read<EBCellFAB>(a_handle, *scratch[a_level], id, m_amr->getGrids(m_realm)[a_level], Interval(0, N - 1), false);

  DataOps::setValue(*m_phi[a_level], 0.0);
  DataOps::incr(*m_phi[a_level], *scratch[a_level], 1.0);
}
#endif

template <int N>
void
SurfaceODESolver<N>::writePlotFile() const noexcept
{
  CH_TIME("SurfaceODESolver::writePlotFile");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotFile" << endl;
  }

  const int numPlotVars = this->getNumberOfPlotVariables();

  if (numPlotVars > 0) {
    const Vector<std::string> plotVarNames = this->getPlotVariableNames();

    CH_assert(plotVarNames.size() == numPlotVars);
    CH_assert(!(m_amr.isNull()));

    EBAMRCellData output;
    m_amr->allocate(output, m_realm, m_phase, numPlotVars);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      int icomp = 0;
      this->writePlotData(*output[lvl], icomp, m_realm, lvl);
    }

    // Filename
    char filename[100];
    sprintf(filename, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_step, SpaceDim);
    std::string fname(filename);

    // Alias, becuase Chombo is a bit anal about smart pointers.
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
}

template <int N>
void
SurfaceODESolver<N>::writePlotData(LevelData<EBCellFAB>& a_output,
                                   int&                  a_comp,
                                   const std::string     a_outputRealm,
                                   const int             a_level) const noexcept
{
  CH_TIME("SurfaceODESolver::writePlotData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  // HDF5 doesn't do irregular data the way we want so we copy into EBCellFABs and
  // use that as output.
  LevelData<EBCellFAB> scratch;
  m_amr->allocate(scratch, m_realm, m_phase, a_level, N);

  // Copy m_phi to a_output
  if (m_plotPhi) {
    DataOps::setValue(scratch, 0.0);
    DataOps::incr(scratch, *m_phi[a_level], 1.0);

    const Interval srcInterv(0, N - 1);
    const Interval dstInterv(a_comp, a_comp + N - 1);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, dstInterv, srcInterv);

    a_output.exchange(dstInterv);

    a_comp += N;
  }

  // Copy m_rhs to a_output
  if (m_plotRHS) {
    DataOps::setValue(scratch, 0.0);
    DataOps::incr(scratch, *m_rhs[a_level], 1.0);

    const Interval srcInterv(0, N - 1);
    const Interval dstInterv(a_comp, a_comp + N - 1);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, dstInterv, srcInterv);

    a_output.exchange(dstInterv);

    a_comp += N;
  }
}

template <int N>
int
SurfaceODESolver<N>::getNumberOfPlotVariables() const noexcept
{
  CH_TIME("SurfaceODESolver::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int num = 0;

  if (m_plotPhi) {
    num += N;
  }

  if (m_plotRHS) {
    num += N;
  }

  return num;
}

template <int N>
Vector<std::string>
SurfaceODESolver<N>::getPlotVariableNames() const noexcept
{
  CH_TIME("SurfaceODESolver::getPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPlotVariablesNames" << endl;
  }

  Vector<std::string> plotVarNames;

  if (m_plotPhi) {
    for (int i = 0; i < N; i++) {
      if (N > 1) {
        plotVarNames.push_back(m_name + " phi-" + std::to_string(i));
      }
      else {
        plotVarNames.push_back(m_name + " phi");
      }
    }
  }

  if (m_plotRHS) {
    for (int i = 0; i < N; i++) {
      if (N > 1) {
        plotVarNames.push_back(m_name + " rhs-" + std::to_string(i));
      }
      else {
        plotVarNames.push_back(m_name + " rhs");
      }
    }
  }

  return plotVarNames;
}

#include <CD_NamespaceFooter.H>

#endif
