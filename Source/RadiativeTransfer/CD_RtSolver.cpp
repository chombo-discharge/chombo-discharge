/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtSolver.cpp
  @brief  Implementation of CD_RtSolver.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_RtSolver.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

constexpr int RtSolver::m_comp;
constexpr int RtSolver::m_nComp;

RtSolver::RtSolver()
{
  CH_TIME("RtSolver::RtSolver");

  // Default settings
  m_verbosity = -1;
  m_name      = "RtSolver";
  m_className = "RtSolver";
}

RtSolver::~RtSolver() { CH_TIME("RtSolver::~RtSolver"); }

std::string
RtSolver::getName()
{
  return m_name;
}

const std::string
RtSolver::getRealm() const
{
  return m_realm;
}

bool
RtSolver::isStationary()
{
  return m_stationary;
}

bool
RtSolver::advance(const Real a_dt, const bool a_zeroPhi)
{
  CH_TIME("RtSolver::advance(dt)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advance(dt)" << endl;
  }

  const bool converged = this->advance(a_dt, m_phi, a_zeroPhi);

  return converged;
}

bool
RtSolver::advance(const Real a_dt, EBAMRCellData& a_phi, const bool a_zeroPhi)
{
  CH_TIME("RtSolver::advance(dt, phi, bool )");
  if (m_verbosity > 5) {
    pout() << m_name + "::advance(dt, phi, bool)" << endl;
  }

  const bool converged = this->advance(a_dt, a_phi, m_source, a_zeroPhi);

  return converged;
}

void
RtSolver::setRealm(const std::string a_realm)
{
  CH_TIME("RtSolver::setRealm");
  if (m_verbosity > 5) {
    pout() << m_name + "::setRealm" << endl;
  }

  m_realm = a_realm;
}

void
RtSolver::setRtSpecies(const RefCountedPtr<RtSpecies>& a_RtSpecies)
{
  CH_TIME("RtSolver::setRtSpecies");
  if (m_verbosity > 5) {
    pout() << m_name + "::setRtSpecies" << endl;
  }

  m_rtSpecies = a_RtSpecies;
  m_name      = m_rtSpecies->getName();
}

void
RtSolver::setPhase(const phase::which_phase a_phase)
{
  CH_TIME("RtSolver::setPhase");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhase" << endl;
  }

  m_phase = a_phase;
}

void
RtSolver::sanityCheck()
{
  CH_TIME("RtSolver::sanityCheck");
  if (m_verbosity > 5) {
    pout() << m_name + "::sanityCheck" << endl;
  }

  CH_assert(!m_computationalGeometry.isNull());
  CH_assert(!m_amr.isNull());
  CH_assert(!m_rtSpecies.isNull());
  CH_assert(!m_ebis.isNull());
}

void
RtSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry> a_computationalGeometry)
{
  CH_TIME("RtSolver::setComputationalGeometry");
  if (m_verbosity > 5) {
    pout() << m_name + "::setComputationalGeometry" << endl;
  }

  m_computationalGeometry = a_computationalGeometry;

  const RefCountedPtr<MultiFluidIndexSpace> MultiFluidIndexSpace = m_computationalGeometry->getMfIndexSpace();

  this->setEbIndexSpace(MultiFluidIndexSpace->getEBIndexSpace(m_phase));
}

void
RtSolver::setEbIndexSpace(const RefCountedPtr<EBIndexSpace>& a_ebis)
{
  CH_TIME("RtSolver::setEbIndexSpace");
  if (m_verbosity > 5) {
    pout() << m_name + "::setEbIndexSpace" << endl;
  }

  m_ebis = a_ebis;
}

void
RtSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr)
{
  CH_TIME("RtSolver::setAmr");
  if (m_verbosity > 5) {
    pout() << m_name + "::setAmr" << endl;
  }

  m_amr = a_amr;
}

void
RtSolver::setTime(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("RtSolver::setTime");
  if (m_verbosity > 5) {
    pout() << m_name + "::setTime" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;
}

void
RtSolver::setStationary(const bool a_stationary)
{
  CH_TIME("RtSolver::setStationary");
  if (m_verbosity > 5) {
    pout() << m_name + "::set_stationry" << endl;
  }

  m_stationary = a_stationary;
}

void
RtSolver::setVerbosity(const int a_verbosity)
{
  CH_TIME("RtSolver::setVerbosity");

  m_verbosity = a_verbosity;
  if (m_verbosity > 5) {
    pout() << m_name + "::setVerbosity" << endl;
  }
}

void
RtSolver::setSource(const EBAMRCellData& a_source)
{
  CH_TIME("RtSolver::setSource(ebamrcell)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setSource(ebamrcell)" << endl;
  }

  m_amr->copyData(m_source, a_source);

  m_amr->conservativeAverage(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void
RtSolver::setSource(const std::function<Real(const RealVect a_pos)> a_source)
{
  CH_TIME("RtSolver::setSource(function)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setSource(function)" << endl;
  }

  constexpr int defaultComponent = 0;

  DataOps::setValue(m_source, a_source, m_amr->getProbLo(), m_amr->getDx(), defaultComponent);
}

void
RtSolver::setSource(const Real a_source)
{
  CH_TIME("RtSolver::setSource(constant)");
  if (m_verbosity > 5) {
    pout() << m_name + "::setSource(constant)" << endl;
  }

  DataOps::setValue(m_source, a_source);

  m_amr->conservativeAverage(m_source, m_realm, m_phase);
  m_amr->interpGhost(m_source, m_realm, m_phase);
}

void
RtSolver::initialData()
{
  CH_TIME("RtSolver::initialData");
  if (m_verbosity > 5) {
    pout() << m_name + "::initialData" << endl;
  }

  DataOps::setValue(m_phi, 0.0);
}

int
RtSolver::getNumberOfPlotVariables() const
{
  CH_TIME("RtSolver::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int numPlotVars = 0;

  if (m_plotPhi) {
    numPlotVars = numPlotVars + 1;
  }
  if (m_plotSource) {
    numPlotVars = numPlotVars + 1;
  }

  return numPlotVars;
}

Vector<std::string>
RtSolver::getPlotVariableNames() const
{
  CH_TIME("RtSolver::getPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPlotVariableNames" << endl;
  }

  Vector<std::string> plotVars(0);

  if (m_plotPhi) {
    plotVars.push_back(m_name + " phi");
  }
  if (m_plotSource) {
    plotVars.push_back(m_name + " source");
  }

  return plotVars;
}

void
RtSolver::writePlotData(LevelData<EBCellFAB>& a_output,
                        int&                  a_comp,
                        const std::string     a_outputRealm,
                        const int             a_level) const noexcept
{
  CH_TIME("RtSolver::writePlotData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  if (m_plotPhi) {
    this->writeData(a_output, a_comp, m_phi, a_outputRealm, a_level, true, true);
  }
  if (m_plotSource) {
    this->writeData(a_output, a_comp, m_source, a_outputRealm, a_level, false, true);
  }
}

void
RtSolver::writeData(LevelData<EBCellFAB>& a_output,
                    int&                  a_comp,
                    const EBAMRCellData&  a_data,
                    const std::string     a_outputRealm,
                    const int             a_level,
                    const bool            a_interpToCentroids,
                    const bool            a_interpGhost) const noexcept

{
  CH_TIMERS("RtSolver::writeData");
  CH_TIMER("RtSolver::writeData::allocate", t1);
  CH_TIMER("RtSolver::writeData::local_copy", t2);
  CH_TIMER("RtSolver::writeData::interp_ghost", t3);
  CH_TIMER("RtSolver::writeData::interp_centroid", t4);
  CH_TIMER("RtSolver::writeData::final_copy", t5);
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
  m_amr->copyData(scratch, *a_data[a_level], a_level, m_realm, m_realm);

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

  DataOps::setCoveredValue(scratch, 0.0);

  CH_START(t5);
  m_amr->copyData(a_output,
                  scratch,
                  a_level,
                  a_outputRealm,
                  m_realm,
                  dstInterv,
                  srcInterv,
                  CopyStrategy::ValidGhost,
                  CopyStrategy::ValidGhost);
  CH_STOP(t5);

  a_comp += numComp;
}

Real
RtSolver::getTime() const
{
  CH_TIME("RtSolver::getTime");
  if (m_verbosity > 5) {
    pout() << m_name + "::getTime" << endl;
  }

  return m_time;
}

phase::which_phase
RtSolver::getPhase()
{
  CH_TIME("RtSolver::getPhase");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPhase" << endl;
  }

  return m_phase;
}

EBAMRCellData&
RtSolver::getPhi()
{
  return m_phi;
}

EBAMRCellData&
RtSolver::getSource()
{
  return m_source;
}

EBAMRFluxData&
RtSolver::getKappa()
{
  return m_kappa;
}

EBAMRIVData&
RtSolver::getKappaEb()
{
  return m_kappaEB;
}

RefCountedPtr<RtSpecies>&
RtSolver::getSpecies()
{
  return m_rtSpecies;
}

void
RtSolver::parseVerbosity() noexcept
{
  CH_TIME("RtSolver::parseVerbosity");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseVerbosity" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("verbosity", m_verbosity);
}

void
RtSolver::computeLoads(Vector<long long>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level) const noexcept
{
  CH_TIME("RtSolver::computeLoads");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeLoads" << endl;
  }

  a_loads.resize(a_dbl.size(), 0LL);

  for (DataIterator dit(a_dbl); dit.ok(); ++dit) {
    a_loads[dit().intCode()] = a_dbl[dit()].numPts();
  }
}

#include <CD_NamespaceFooter.H>
