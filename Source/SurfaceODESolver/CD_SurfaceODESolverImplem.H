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
#include <CD_NamespaceHeader.H>

template <size_t N>
SurfaceODESolver<N>::SurfaceODESolver()
{
  CH_TIME("SurfaceODESolver::SurfaceODESolver");

  m_className = "SurfaceODESolver";

  this->setVerbosity(-1);
  this->setRealm(Realm::Primal);
  this->setPhase(phase::gas);
}

template <size_t N>
SurfaceODESolver<N>::SurfaceODESolver(const RefCountedPtr<AmrMesh>&               a_amr,
                                      const RefCountedPtr<ComputationalGeometry>& a_compGeom)
  : SurfaceODESolver<N>()
{
  CH_TIME("SurfaceODESolver::SurfaceODESolver");
  if (m_verbosity > 5) {
    pout() << m_name + "::SurfaceODESolver()" << endl;
  }

  m_amr      = a_amr;
  m_compGeom = a_compGeom;
}

template <size_t N>
SurfaceODESolver<N>::~SurfaceODESolver()
{
  CH_TIME("SurfaceODESolver::~SurfaceODESolver");
  if (m_verbosity > 5) {
    pout() << m_name + "::~SurfaceODESolver()" << endl;
  }
}

template <size_t N>
void
SurfaceODESolver<N>::parseOptions()
{
  CH_TIME("SurfaceODESolver::parseOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions" << endl;
  }

  this->parseVerbosity();
  this->parseRegrid();
  this->parsePlotVariables();
}

template <size_t N>
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

template <size_t N>
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

template <size_t N>
void
SurfaceODESolver<N>::parseRegrid() noexcept
{
  CH_TIME("SurfaceODESolver::parseRegrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRegrid" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("conservative_regrid", m_conservativeRegrid);
}

template <size_t N>
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

  if (pp.countval("plt_vars") > 0) {
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
void
SurfaceODESolver<N>::setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh) noexcept
{
  CH_TIME("SurfaceODESolver::setAmr");
  if (m_verbosity > 5) {
    pout() << m_name + "::setAmr" << endl;
  }

  CH_assert(!(a_amrMesh.isNull()));

  m_amr = a_amr;
}

template <size_t N>
void
SurfaceODESolver<N>::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_compGeom) noexcept
{
  CH_TIME("SurfaceODESolver::setComputationalGeometry");
  if (m_verbosity > 5) {
    pout() << m_name + "::setComputationalGeometry" << endl;
  }

  CH_assert(!(a_compGeom.isNull()));

  m_compGeom = a_compGeom;
  m_mfis     = m_comp->getMfIndexSpace();
}

template <size_t N>
void
SurfaceODESolver<N>::setRealm(const std::string a_realm) noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::setName(const std::string a_name) noexcept
{
  CH_TIME("SurfaceODESolver::");
}

const std::string
SurfaceODESolver<N>::getRealm() const noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::setPhase(phase::which_phase a_phase) noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::setVerbosity(const int a_verbosity) noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::getVerbosity() const noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::setTime(const int a_step, const Real a_time, const Real a_dt) noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
Real
SurfaceODESolver<N>::getStep() const noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
Real
SurfaceODESolver<N>::getTime() const noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
Real
SurfaceODESolver<N>::getTimeStep() const noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
phase::which_phase
SurfaceODESolver<N>::getPhase() const noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::setPhi(const Real a_phi)
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::setPhi(const std::array<Real, N>& a_phi)
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::setPhi(const EBAMRIVData& a_phi)
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
EBAMRIVData&
SurfaceODESolver<N>::getPhi() noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
const EBAMRIVData&
SurfaceODESolver<N>::getPhi() const noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::setRHS(const Real a_rhs)
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::setRHS(const std::array<Real, N>& a_rhs)
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::setRHS(const EBAMRIVData& a_rhs)
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
EBAMRIVData&
SurfaceODESolver<N>::getRHS()
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
const EBAMRIVData&
SurfaceODESolver<N>::getRHS() const noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::allocateInternals() noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::deallocateInternals() noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::registerOperators() noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::preRegrid(const int a_lbase, const int a_oldFinestLevel) noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::resetElectrodes(const Real a_value) noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::resetElectrodes(EBAMRIVData& a_phi, const Real a_value) const noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::resetDielectrics(const Real a_value) noexcept
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
void
SurfaceODESolver<N>::resetDielectrics(EBAMRIVData& a_phi, const Real a_value) const noexcept
{
  CH_TIME("SurfaceODESolver::");
}

#ifdef CH_USE_HDF5

template <size_t N>
void
SurfaceODESolver<N>::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const
{
  CH_TIME("SurfaceODESolver::");
}
#endif

#ifdef CH_USE_HDF5

template <size_t N>
void
SurfaceODESolver<N>::readCheckpointLevel(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("SurfaceODESolver::");
}
#endif

template <size_t N>
void
SurfaceODESolver<N>::writePlotData(EBAMRCellData& a_output, int& a_comp)
{
  CH_TIME("SurfaceODESolver::");
}

template <size_t N>
int
SurfaceODESolver<N>::getNumberOfPlotVariables()
{
  CH_TIME("SurfaceODESolver::getNumberOfPlotVariables");
}

template <size_t N>
Vector<std::string>
SurfaceODESolver<N>::getPlotVariableNames() const
{
  CH_TIME("SurfaceODESolver::getPlotVariableNames");
}

#include <CD_NamespaceFooter.H>

#endif
