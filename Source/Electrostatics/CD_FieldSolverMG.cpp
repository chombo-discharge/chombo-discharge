/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_FieldSolverMG.cpp
 * @brief  Implementation of CD_FieldSolverMG.H
 * @author Robert Marskar
 */

// Std includes
#include <algorithm>

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_FieldSolverMG.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

FieldSolverMG::FieldSolverMG()
  : m_isSolverSetup(false),
    m_bcOrder(1),
    m_bcWeight(0),
    m_hasSolverPermittivities(false)
#ifdef CH_USE_PETSC
    ,
    m_amgFallback(false),
    m_amgLatched(false),
    m_amgRetryInterval(1),
    m_solvesSinceProbe(0)
#endif
{
  CH_TIME("FieldSolverMG::FieldSolverMG()");

  m_className = "FieldSolverMG";

  m_gmg = RefCountedPtr<FieldSolverGMG>(new FieldSolverGMG());
}

FieldSolverMG::~FieldSolverMG()
{
  CH_TIME("FieldSolverMG::~FieldSolverMG()");
}

void
FieldSolverMG::parseOptions()
{
  CH_TIME("FieldSolverMG::parseOptions()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::parseOptions()" << endl;
  }

  this->parseVerbosity();
  this->parseDomainBc();
  this->parsePlotVariables();
  this->parseRegridSlopes();
  this->parseSolverSettings();

  // The interface condition is fixed after the first parse, as in FieldSolverGMG. FieldSolverAMG reads it in
  // parseSolverSettings.
  m_gmg->parseJumpBC();
}

void
FieldSolverMG::parseRuntimeOptions()
{
  CH_TIME("FieldSolverMG::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::parseRuntimeOptions()" << endl;
  }

  this->parseVerbosity();
  this->parsePlotVariables();
  this->parseRegridSlopes();
  this->parseSolverSettings();
}

void
FieldSolverMG::parseSolverSettings()
{
  CH_TIME("FieldSolverMG::parseSolverSettings()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::parseSolverSettings()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("bc_order", m_bcOrder);
  pp.get("bc_weight", m_bcWeight);

  if (m_bcOrder < 1 || m_bcWeight < 0) {
    MayDay::Abort("FieldSolverMG::parseSolverSettings -- 'bc_order' must be at least 1 and 'bc_weight' non-negative");
  }

  // The geometric multigrid solver reads its settings under this class's name. Its stencil settings come from
  // bc_order/bc_weight, and it never replaces the settings by its defaults, which would override the input.
  m_gmg->m_className = m_className;
  m_gmg->m_verbosity = m_verbosity;

  m_gmg->parseMultigridSettings();
  m_gmg->parseKrylovSettings("gmg_");
  m_gmg->parseKappaSource();

  m_gmg->m_multigridUseDefaultSettings = false;
  m_gmg->m_multigridBcOrder            = m_bcOrder;
  m_gmg->m_multigridBcWeight           = m_bcWeight;
  m_gmg->m_multigridJumpOrder          = m_bcOrder;
  m_gmg->m_multigridJumpWeight         = m_bcWeight;

  bool amgFallback      = false;
  int  amgRetryInterval = 1;

  pp.get("amg_fallback", amgFallback);
  pp.get("amg_retry_interval", amgRetryInterval);

#ifdef CH_USE_PETSC
  m_amgFallback      = amgFallback;
  m_amgRetryInterval = std::max(1, amgRetryInterval);

  if (m_amgFallback) {
    if (m_amg.isNull()) {
      m_amg = RefCountedPtr<FieldSolverAMG>(new FieldSolverAMG());
    }

    // The algebraic multigrid solver also reads its settings under this class's name, and uses the same EB stencils
    // as geometric multigrid.
    m_amg->m_className = m_className;
    m_amg->m_verbosity = m_verbosity;

    m_amg->parseSolverSettings();

    m_amg->m_ebOrder  = m_bcOrder;
    m_amg->m_ebWeight = m_bcWeight;
  }
  else {
    m_amg        = RefCountedPtr<FieldSolverAMG>();
    m_amgLatched = false;

    m_phiStart.clear();
  }
#else
  (void)amgRetryInterval;

  if (amgFallback) {
    pout() << "FieldSolverMG -- 'amg_fallback' is ignored because chombo-discharge was built without PETSc" << endl;
  }
#endif
}

void
FieldSolverMG::syncSolvers()
{
  CH_TIME("FieldSolverMG::syncSolvers()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::syncSolvers()" << endl;
  }

  m_gmg->aliasState(*this);

#ifdef CH_USE_PETSC
  if (!m_amg.isNull()) {
    m_amg->aliasState(*this);
  }
#endif
}

void
FieldSolverMG::registerOperators()
{
  CH_TIME("FieldSolverMG::registerOperators()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::registerOperators()" << endl;
  }

  CH_assert(!m_amr.isNull());

  // The operators of geometric multigrid include all those that algebraic multigrid needs.
  this->syncSolvers();

  m_gmg->registerOperators();
}

void
FieldSolverMG::allocate()
{
  CH_TIME("FieldSolverMG::allocate()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::allocate()" << endl;
  }

  FieldSolver::allocate();

#ifdef CH_USE_PETSC
  if (m_amgFallback) {
    m_amr->allocate(m_phiStart, m_realm, m_nComp);
  }
#endif

  this->syncSolvers();
}

void
FieldSolverMG::deallocate()
{
  CH_TIME("FieldSolverMG::deallocate()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::deallocate()" << endl;
  }

  FieldSolver::deallocate();

  m_hasSolverPermittivities = false;

  m_solverPermCell.clear();
  m_solverPermFace.clear();
  m_solverPermEB.clear();

#ifdef CH_USE_PETSC
  m_phiStart.clear();
#endif

  // The solvers must not keep the old data alive.
  this->syncSolvers();
}

void
FieldSolverMG::setupSolver()
{
  CH_TIME("FieldSolverMG::setupSolver()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::setupSolver()" << endl;
  }

  this->syncSolvers();

  m_gmg->setupSolver();

#ifdef CH_USE_PETSC
  if (!m_amg.isNull()) {
    m_amg->invalidateSolver();
  }
#endif

  // As in FieldSolverGMG and FieldSolverAMG, a setup restores the permittivities.
  m_hasSolverPermittivities = false;

  m_solverPermCell.clear();
  m_solverPermFace.clear();
  m_solverPermEB.clear();

  m_isSolverSetup = true;
}

void
FieldSolverMG::invalidateSolver()
{
  CH_TIME("FieldSolverMG::invalidateSolver()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::invalidateSolver()" << endl;
  }

  m_gmg->invalidateSolver();

#ifdef CH_USE_PETSC
  if (!m_amg.isNull()) {
    m_amg->invalidateSolver();
  }
#endif

  m_isSolverSetup = false;
}

void
FieldSolverMG::setSolverPermittivities(const MFAMRCellData& a_permittivityCell,
                                       const MFAMRFluxData& a_permittivityFace,
                                       const MFAMRIVData&   a_permittivityEB)
{
  CH_TIME("FieldSolverMG::setSolverPermittivities()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::setSolverPermittivities()" << endl;
  }

  if (!m_isSolverSetup) {
    MayDay::Error("FieldSolverMG::setSolverPermittivities -- must set up solver first!");
  }

  m_gmg->setSolverPermittivities(a_permittivityCell, a_permittivityFace, a_permittivityEB);

  // Keep the coefficients for algebraic multigrid, which may be set up later. The copies share the caller's storage.
  m_hasSolverPermittivities = true;

  m_solverPermCell = a_permittivityCell;
  m_solverPermFace = a_permittivityFace;
  m_solverPermEB   = a_permittivityEB;

#ifdef CH_USE_PETSC
  if (!m_amg.isNull() && m_amg->m_isSolverSetup) {
    m_amg->setSolverPermittivities(a_permittivityCell, a_permittivityFace, a_permittivityEB);
  }
#endif
}

void
FieldSolverMG::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("FieldSolverMG::preRegrid(int, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::preRegrid(int, int)" << endl;
  }

  FieldSolver::preRegrid(a_lbase, a_oldFinestLevel);

  this->invalidateSolver();
}

void
FieldSolverMG::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("FieldSolverMG::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::regrid(int, int, int)" << endl;
  }

  FieldSolver::regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  this->invalidateSolver();
}

bool
FieldSolverMG::solve(MFAMRCellData& a_phi, const MFAMRCellData& a_rho, const EBAMRIVData& a_sigma, const bool a_zeroPhi)
{
  CH_TIME("FieldSolverMG::solve(MFAMRCellData, MFAMRCellData, EBAMRIVData, bool)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::solve(MFAMRCellData, MFAMRCellData, EBAMRIVData, bool)" << endl;
  }

  CH_assert(m_isVoltageSet);

  if (!m_isSolverSetup) {
    this->setupSolver();
  }

  // Pass on the time, which the boundary conditions and the covered-cell potential depend on.
  this->syncSolvers();

#ifdef CH_USE_PETSC
  std::string reason;

  const bool useAMG = m_amgFallback && this->isAMGSupported(reason);

  if (useAMG) {
    // While latched onto algebraic multigrid, geometric multigrid is only re-tried every m_amgRetryInterval solves.
    bool tryGMG = true;

    if (m_amgLatched) {
      m_solvesSinceProbe++;

      if (m_solvesSinceProbe < m_amgRetryInterval) {
        tryGMG = false;
      }
      else {
        m_solvesSinceProbe = 0;
      }
    }

    if (a_zeroPhi) {
      DataOps::setValue(a_phi, 0.0);
    }

    if (tryGMG) {
      DataOps::copy(m_phiStart, a_phi);

      if (m_gmg->solve(a_phi, a_rho, a_sigma, false)) {
        if (m_amgLatched) {
          pout() << "FieldSolverMG::solve -- geometric multigrid converged again" << endl;
        }

        m_amgLatched       = false;
        m_solvesSinceProbe = 0;

        return true;
      }

      // Algebraic multigrid starts from the same initial guess as geometric multigrid did.
      DataOps::copy(a_phi, m_phiStart);

      if (!m_amgLatched) {
        pout() << "FieldSolverMG::solve -- geometric multigrid did not converge; falling back to algebraic multigrid"
               << endl;
      }
    }

    const bool converged = this->solveAMG(a_phi, a_rho, a_sigma);

    if (converged) {
      m_amgLatched = true;
    }
    else {
      pout() << "FieldSolverMG::solve -- algebraic multigrid did not converge either" << endl;
    }

    return converged;
  }

  const bool converged = m_gmg->solve(a_phi, a_rho, a_sigma, a_zeroPhi);

  if (!converged && m_amgFallback) {
    pout() << "FieldSolverMG::solve -- geometric multigrid did not converge, and the algebraic multigrid fallback is "
           << "not available: " << reason << endl;
  }

  return converged;
#else
  return m_gmg->solve(a_phi, a_rho, a_sigma, a_zeroPhi);
#endif
}

#ifdef CH_USE_PETSC
bool
FieldSolverMG::isAMGSupported(std::string& a_reason) const
{
  CH_TIME("FieldSolverMG::isAMGSupported(std::string)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::isAMGSupported(std::string)" << endl;
  }

  const ProblemDomain& domain = m_amr->getDomains()[0];

  for (int dir = 0; dir < SpaceDim; dir++) {
    if (domain.isPeriodic(dir)) {
      a_reason = "the domain is periodic";

      return false;
    }
  }

  // The EB stencils of algebraic multigrid reach as far as their order, which the coarse-fine interpolators must
  // cover.
  for (int lvl = 1; lvl <= m_amr->getFinestLevel(); lvl++) {
    const RefCountedPtr<EBMultigridInterpolator>& interp = m_amr->getMultigridInterpolator(m_realm, phase::gas)[lvl];

    if (!interp.isNull() && m_bcOrder > interp->getGhostCF()) {
      a_reason = "bc_order is larger than the coarse-fine interpolation radius (AmrMesh.mg_interp_radius)";

      return false;
    }
  }

  return true;
}

bool
FieldSolverMG::solveAMG(MFAMRCellData& a_phi, const MFAMRCellData& a_rho, const EBAMRIVData& a_sigma)
{
  CH_TIME("FieldSolverMG::solveAMG(MFAMRCellData, MFAMRCellData, EBAMRIVData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::solveAMG(MFAMRCellData, MFAMRCellData, EBAMRIVData)" << endl;
  }

  CH_assert(!m_amg.isNull());

  // Algebraic multigrid is set up on demand, with the permittivities that geometric multigrid currently uses.
  if (!m_amg->m_isSolverSetup) {
    m_amg->setupSolver();

    if (m_hasSolverPermittivities) {
      m_amg->setSolverPermittivities(m_solverPermCell, m_solverPermFace, m_solverPermEB);
    }
  }

  const bool converged = m_amg->solve(a_phi, a_rho, a_sigma, false);

  // As in FieldSolverGMG: optionally fill the covered cells with the electrode potential, for visualization.
  ParmParse pp(m_className.c_str());

  bool fillElectrodePotential = false;

  pp.query("fill_electrode_potential", fillElectrodePotential);

  if (fillElectrodePotential) {
    this->fillCoveredPotential(a_phi);
  }

  return converged;
}
#endif

void
FieldSolverMG::computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential) const
{
  CH_TIME("FieldSolverMG::computeElectricField(MFAMRCellData, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::computeElectricField(MFAMRCellData, MFAMRCellData)" << endl;
  }

  m_gmg->computeElectricField(a_E, a_potential);
}

void
FieldSolverMG::computeElectricField(MFAMRFluxData& a_E, const MFAMRCellData& a_potential) const
{
  CH_TIME("FieldSolverMG::computeElectricField(MFAMRFluxData, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::computeElectricField(MFAMRFluxData, MFAMRCellData)" << endl;
  }

  m_gmg->computeElectricField(a_E, a_potential);
}

void
FieldSolverMG::computeElectricField(EBAMRCellData&           a_E,
                                    const phase::which_phase a_phase,
                                    const MFAMRCellData&     a_potential) const
{
  CH_TIME("FieldSolverMG::computeElectricField(EBAMRCellData, phase, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::computeElectricField(EBAMRCellData, phase, MFAMRCellData)" << endl;
  }

  m_gmg->computeElectricField(a_E, a_phase, a_potential);
}

void
FieldSolverMG::computeElectricField(EBAMRFluxData&           a_E,
                                    const phase::which_phase a_phase,
                                    const MFAMRCellData&     a_potential) const
{
  CH_TIME("FieldSolverMG::computeElectricField(EBAMRFluxData, phase, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::computeElectricField(EBAMRFluxData, phase, MFAMRCellData)" << endl;
  }

  m_gmg->computeElectricField(a_E, a_phase, a_potential);
}

Vector<long long>
FieldSolverMG::computeLoads(const DisjointBoxLayout& a_dbl, const int a_level)
{
  CH_TIME("FieldSolverMG::computeLoads(DisjointBoxLayout, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMG::computeLoads(DisjointBoxLayout, int)" << endl;
  }

  this->syncSolvers();

  return m_gmg->computeLoads(a_dbl, a_level);
}

#include <CD_NamespaceFooter.H>
