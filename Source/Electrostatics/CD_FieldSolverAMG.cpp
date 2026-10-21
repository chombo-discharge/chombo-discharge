/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_FieldSolverAMG.cpp
 * @brief  Implementation of CD_FieldSolverAMG.H
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_FieldSolverAMG.H>
#include <CD_PetscHelmholtzElectrostaticDomainBC.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

constexpr Real FieldSolverAMG::m_alpha;
constexpr Real FieldSolverAMG::m_beta;

FieldSolverAMG::FieldSolverAMG() : m_isSolverSetup(false), m_ksp(nullptr)
{
  CH_TIME("FieldSolverAMG::FieldSolverAMG()");

  m_className = "FieldSolverAMG";
}

FieldSolverAMG::~FieldSolverAMG()
{
  CH_TIME("FieldSolverAMG::~FieldSolverAMG()");

  this->destroySolver();
}

void
FieldSolverAMG::parseOptions()
{
  CH_TIME("FieldSolverAMG::parseOptions()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::parseOptions()" << endl;
  }

  this->parseVerbosity();
  this->parseDomainBc();
  this->parsePlotVariables();
  this->parseSolverSettings();
  this->parseRegridSlopes();
}

void
FieldSolverAMG::parseRuntimeOptions()
{
  CH_TIME("FieldSolverAMG::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::parseRuntimeOptions()" << endl;
  }

  this->parseVerbosity();
  this->parsePlotVariables();
  this->parseRegridSlopes();
}

void
FieldSolverAMG::parseSolverSettings()
{
  CH_TIME("FieldSolverAMG::parseSolverSettings()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::parseSolverSettings()" << endl;
  }

  ParmParse pp(m_className.c_str());

  // The Krylov solver is given by our name for it, and stored as the PETSc type.
  std::string kspType;

  pp.get("ksp_type", kspType);

  if (kspType == "bicgstab") {
    m_kspType = KSPBCGS;
  }
  else if (kspType == "gmres") {
    m_kspType = KSPGMRES;
  }
  else {
    MayDay::Abort("FieldSolverAMG::parseSolverSettings -- 'ksp_type' must be 'bicgstab' or 'gmres'");
  }

  pp.get("pc_type", m_pcType);
  pp.get("ksp_rtol", m_kspRtol);
  pp.get("ksp_atol", m_kspAtol);
  pp.get("ksp_max_iter", m_kspMaxIter);
}

bool
FieldSolverAMG::solve(MFAMRCellData&       a_potential,
                      const MFAMRCellData& a_rho,
                      const EBAMRIVData& /*a_sigma*/,
                      const bool a_zeroPhi)
{
  CH_TIME("FieldSolverAMG::solve(MFAMRCellData, MFAMRCellData, EBAMRIVData, bool)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::solve(MFAMRCellData, MFAMRCellData, EBAMRIVData, bool)" << endl;
  }

  CH_assert(a_potential[0]->nComp() == 1);
  CH_assert(a_rho[0]->nComp() == 1);

  if (!m_isSolverSetup) {
    this->setupSolver();
  }

  const RefCountedPtr<PetscGrid>& petscGrid = m_amr->getPetscGrid(m_realm);

  // Right-hand side rho/eps0, with the boundary data part of the domain boundary fluxes moved into it. There are no
  // embedded boundaries yet, so the surface charge does not enter.
  MFAMRCellData rhoByEps0;
  m_amr->allocate(rhoByEps0, m_realm, m_nComp);
  DataOps::copy(rhoByEps0, a_rho);
  DataOps::scale(rhoByEps0, 1.0 / Units::eps0);

  Vec b;
  Vec x;
  Vec r;

  petscGrid->create(b);
  petscGrid->create(x);
  petscGrid->create(r);

  petscGrid->putChomboInPetsc(b, rhoByEps0);
  m_petscOp->addBoundaryConditionsToRHS(b);

  if (a_zeroPhi) {
    petscGrid->setValue(x, 0.0);
  }
  else {
    petscGrid->putChomboInPetsc(x, a_potential);
  }

  PetscCallAbort(Chombo_MPI::comm, KSPSetInitialGuessNonzero(m_ksp, a_zeroPhi ? PETSC_FALSE : PETSC_TRUE));
  PetscCallAbort(Chombo_MPI::comm, KSPSolve(m_ksp, b, x));

  KSPConvergedReason reason;
  PetscInt           numIter;

  PetscCallAbort(Chombo_MPI::comm, KSPGetConvergedReason(m_ksp, &reason));
  PetscCallAbort(Chombo_MPI::comm, KSPGetIterationNumber(m_ksp, &numIter));

  const bool converged = (reason > 0);

  // Residual r = b - A*x, in the max norm, which is what FieldSolverGMG reports.
  PetscReal residNorm;
  PetscReal rhsNorm;

  PetscCallAbort(Chombo_MPI::comm, MatResidual(m_petscOp->getMatrix(), b, x, r));
  PetscCallAbort(Chombo_MPI::comm, VecNorm(r, NORM_INFINITY, &residNorm));
  PetscCallAbort(Chombo_MPI::comm, VecNorm(b, NORM_INFINITY, &rhsNorm));

  if (m_verbosity > 1) {
    pout() << "FieldSolverAMG::solve -- " << (converged ? "converged" : "did not converge") << " in " << numIter
           << " iterations (reason = " << reason << "), max-norm residual = " << residNorm
           << ", max-norm right-hand side = " << rhsNorm << endl;
  }

  petscGrid->putPetscInChombo(a_potential, x);
  petscGrid->putPetscInChombo(m_residue, r);

  petscGrid->destroy(b);
  petscGrid->destroy(x);
  petscGrid->destroy(r);

  // Coarsen and update ghost cells before computing the field, as FieldSolverGMG does.
  m_amr->conservativeAverage(a_potential, m_realm);
  m_amr->interpGhostPwl(a_potential, m_realm);

  this->computeElectricField(m_electricField, a_potential);

  return converged;
}

void
FieldSolverAMG::setupSolver()
{
  CH_TIME("FieldSolverAMG::setupSolver()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::setupSolver()" << endl;
  }

  this->destroySolver();

  auto domainBc = RefCountedPtr<PetscHelmholtzDomainBC>(new PetscHelmholtzElectrostaticDomainBC(m_domainBc));

  m_petscOp = RefCountedPtr<PetscHelmholtzOp>(new PetscHelmholtzOp());
  m_petscOp->define(m_amr->getPetscGrid(m_realm),
                    domainBc,
                    m_permittivityCell,
                    m_permittivityFace,
                    m_amr->getProbLo(),
                    m_amr->getDx(),
                    m_alpha,
                    m_beta);

  // The settings from the input file are the defaults. They can be overridden with command-line options that carry
  // the prefix 'fieldsolveramg_', e.g. -fieldsolveramg_ksp_monitor.
  PC pc;

  PetscCallAbort(Chombo_MPI::comm, KSPCreate(Chombo_MPI::comm, &m_ksp));
  PetscCallAbort(Chombo_MPI::comm, KSPSetOperators(m_ksp, m_petscOp->getMatrix(), m_petscOp->getMatrix()));
  PetscCallAbort(Chombo_MPI::comm, KSPSetType(m_ksp, m_kspType.c_str()));
  PetscCallAbort(Chombo_MPI::comm, KSPGetPC(m_ksp, &pc));
  PetscCallAbort(Chombo_MPI::comm, PCSetType(pc, m_pcType.c_str()));
  // Keep PETSc's divergence tolerance. PETSC_DEFAULT is deprecated for this from PETSc 3.22, which has PETSC_CURRENT.
#if PETSC_VERSION_GE(3, 22, 0)
  constexpr PetscReal keepDivergenceTol = PETSC_CURRENT;
#else
  constexpr PetscReal keepDivergenceTol = PETSC_DEFAULT;
#endif

  PetscCallAbort(Chombo_MPI::comm, KSPSetTolerances(m_ksp, m_kspRtol, m_kspAtol, keepDivergenceTol, m_kspMaxIter));
  PetscCallAbort(Chombo_MPI::comm, KSPSetOptionsPrefix(m_ksp, "fieldsolveramg_"));
  PetscCallAbort(Chombo_MPI::comm, KSPSetFromOptions(m_ksp));

  m_isSolverSetup = true;
}

void
FieldSolverAMG::destroySolver()
{
  CH_TIME("FieldSolverAMG::destroySolver()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::destroySolver()" << endl;
  }

  if (m_ksp != nullptr) {
    PetscCallAbort(Chombo_MPI::comm, KSPDestroy(&m_ksp));
  }

  m_petscOp       = RefCountedPtr<PetscHelmholtzOp>();
  m_isSolverSetup = false;
}

void
FieldSolverAMG::setSolverPermittivities(const MFAMRCellData& /*a_permittivityCell*/,
                                        const MFAMRFluxData& /*a_permittivityFace*/,
                                        const MFAMRIVData& /*a_permittivityEB*/)
{
  CH_TIME("FieldSolverAMG::setSolverPermittivities()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::setSolverPermittivities()" << endl;
  }

  MayDay::Abort("FieldSolverAMG::setSolverPermittivities -- not supported yet");
}

void
FieldSolverAMG::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("FieldSolverAMG::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::regrid(int, int, int)" << endl;
  }

  FieldSolver::regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  this->destroySolver();
}

void
FieldSolverAMG::registerOperators()
{
  CH_TIME("FieldSolverAMG::registerOperators()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::registerOperators()" << endl;
  }

  CH_assert(!m_amr.isNull());

  m_amr->registerOperator(s_eb_gradient, m_realm, phase::gas);
  m_amr->registerOperator(s_eb_gradient, m_realm, phase::solid);
  m_amr->registerOperator(s_eb_coar_ave, m_realm, phase::gas);
  m_amr->registerOperator(s_eb_coar_ave, m_realm, phase::solid);
  m_amr->registerOperator(s_eb_fill_patch, m_realm, phase::gas);
  m_amr->registerOperator(s_eb_fill_patch, m_realm, phase::solid);
  m_amr->registerOperator(s_eb_fine_interp, m_realm, phase::gas);
  m_amr->registerOperator(s_eb_fine_interp, m_realm, phase::solid);
  m_amr->registerOperator(s_eb_irreg_interp, m_realm, phase::gas);
  m_amr->registerOperator(s_eb_irreg_interp, m_realm, phase::solid);
  m_amr->registerOperator(s_eb_multigrid, m_realm, phase::gas);
  m_amr->registerOperator(s_eb_multigrid, m_realm, phase::solid);
}

void
FieldSolverAMG::computeElectricField(MFAMRCellData& a_electricField, const MFAMRCellData& a_potential) const
{
  CH_TIME("FieldSolverAMG::computeElectricField(MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::computeElectricField(MFAMRCellData)" << endl;
  }

  CH_assert(a_electricField[0]->nComp() == SpaceDim);
  CH_assert(a_potential[0]->nComp() == 1);

  // As in FieldSolverGMG: fill ghost cells with the multigrid interpolator, which is consistent with the
  // discretization, then compute the gradient.
  MFAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_nComp);
  m_amr->copyData(scratch, a_potential);
  m_amr->interpGhostMG(scratch, m_realm);

  m_amr->computeGradient(a_electricField, scratch, m_realm);
  DataOps::scale(a_electricField, -1.0);

  m_amr->conservativeAverage(a_electricField, m_realm);
  m_amr->interpGhost(a_electricField, m_realm);
}

void
FieldSolverAMG::computeElectricField(MFAMRFluxData& a_electricField, const MFAMRCellData& a_potential) const
{
  CH_TIME("FieldSolverAMG::computeElectricField(MFAMRFluxData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::computeElectricField(MFAMRFluxData)" << endl;
  }

  CH_assert(a_electricField[0]->nComp() == SpaceDim);
  CH_assert(a_potential[0]->nComp() == 1);

  MFAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_nComp);
  m_amr->copyData(scratch, a_potential);
  m_amr->interpGhostMG(scratch, m_realm);

  m_amr->computeGradient(a_electricField, scratch, m_realm);
  DataOps::scale(a_electricField, -1.0);
}

void
FieldSolverAMG::computeElectricField(EBAMRCellData&           a_electricField,
                                     const phase::which_phase a_phase,
                                     const MFAMRCellData&     a_potential) const
{
  CH_TIME("FieldSolverAMG::computeElectricField(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::computeElectricField(EBAMRCellData)" << endl;
  }

  CH_assert(a_electricField[0]->nComp() == SpaceDim);
  CH_assert(a_potential[0]->nComp() == 1);

  EBAMRCellData scratch;
  EBAMRCellData potentialPhase;

  m_amr->allocatePointer(potentialPhase, m_realm);
  m_amr->allocate(scratch, m_realm, a_phase, m_nComp);
  m_amr->alias(potentialPhase, a_phase, a_potential);

  m_amr->copyData(scratch, potentialPhase);
  m_amr->interpGhostMG(scratch, m_realm, a_phase);

  m_amr->computeGradient(a_electricField, scratch, m_realm, a_phase);
  DataOps::scale(a_electricField, -1.0);

  m_amr->conservativeAverage(a_electricField, m_realm, a_phase);
  m_amr->interpGhost(a_electricField, m_realm, a_phase);
}

void
FieldSolverAMG::computeElectricField(EBAMRFluxData&           a_electricField,
                                     const phase::which_phase a_phase,
                                     const MFAMRCellData&     a_potential) const
{
  CH_TIME("FieldSolverAMG::computeElectricField(EBAMRFluxData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::computeElectricField(EBAMRFluxData)" << endl;
  }

  CH_assert(a_electricField[0]->nComp() == SpaceDim);
  CH_assert(a_potential[0]->nComp() == 1);

  EBAMRCellData scratch;
  EBAMRCellData potentialPhase;

  m_amr->allocatePointer(potentialPhase, m_realm);
  m_amr->allocate(scratch, m_realm, a_phase, m_nComp);
  m_amr->alias(potentialPhase, a_phase, a_potential);

  m_amr->copyData(scratch, potentialPhase);
  m_amr->interpGhostMG(scratch, m_realm, a_phase);

  m_amr->computeGradient(a_electricField, scratch, m_realm, a_phase);
  DataOps::scale(a_electricField, -1.0);
}

#include <CD_NamespaceFooter.H>

#endif
