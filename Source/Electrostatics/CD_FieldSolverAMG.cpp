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
#include <CD_PetscHelmholtzElectrostaticEBBC.H>
#include <CD_PetscHelmholtzNaturalJumpBC.H>
#include <CD_PetscHelmholtzSaturationChargeJumpBC.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

constexpr Real FieldSolverAMG::m_alpha;
constexpr Real FieldSolverAMG::m_beta;

FieldSolverAMG::FieldSolverAMG()
  : m_isSolverSetup(false),
    m_ksp(nullptr),
    m_petscVerbosity(-1),
    m_monitorZeroResidNorm(0.0),
    m_monitorFirstNorm(0.0),
    m_monitorPreviousNorm(0.0)
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

  // The outer Krylov solver is any PETSc KSP type (checked by PETSc in setupSolver), with 'bicgstab' accepted for
  // PETSc's 'bcgs'.
  pp.get("ksp_type", m_kspType);

  if (m_kspType == "bicgstab") {
    m_kspType = KSPBCGS;
  }

  pp.get("ksp_gmres_restart", m_kspGmresRestart);
  pp.get("pc_type", m_pcType);

  // The smoother on each level of a multigrid preconditioner (e.g. gamg): its KSP and PC types, and the number of
  // smoothing steps.
  pp.get("amg_smoother_ksp", m_amgSmootherKsp);
  pp.get("amg_smoother_pc", m_amgSmootherPc);
  pp.get("amg_smoother_its", m_amgSmootherIts);

  // GAMG's coarsening: the graph threshold below which weak connections are dropped, and the number of levels that
  // are coarsened aggressively. Other preconditioners ignore them.
  pp.get("gamg_threshold", m_gamgThreshold);
  pp.get("gamg_aggressive_coarsening", m_gamgAggressiveCoarsening);
  pp.get("ksp_rtol", m_kspRtol);
  pp.get("ksp_atol", m_kspAtol);
  pp.get("ksp_max_iter", m_kspMaxIter);
  pp.get("petsc_verbosity", m_petscVerbosity);
  pp.get("kappa_source", m_kappaSource);

  // EB face stencils, for the electrodes and the interface between the gas and the dielectrics.
  std::string jumpBc;

  pp.get("eb_order", m_ebOrder);
  pp.get("eb_weight", m_ebWeight);
  pp.get("jump_bc", jumpBc);
  pp.get("interface_points", m_interfacePoints);

  if (jumpBc == "natural") {
    m_jumpBcType = JumpBCType::Natural;
  }
  else if (jumpBc == "saturation_charge") {
    m_jumpBcType = JumpBCType::SaturationCharge;
  }
  else {
    MayDay::Abort("FieldSolverAMG::parseSolverSettings -- 'jump_bc' must be 'natural' or 'saturation_charge'");
  }
}

bool
FieldSolverAMG::solve(MFAMRCellData&       a_potential,
                      const MFAMRCellData& a_rho,
                      const EBAMRIVData&   a_sigma,
                      const bool           a_zeroPhi)
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

  // Right-hand side rho/eps0, with the boundary data parts of the domain, electrode and interface fluxes moved into
  // it. The surface charge enters through the interface fluxes.
  MFAMRCellData rhoByEps0;
  m_amr->allocate(rhoByEps0, m_realm, m_nComp);
  DataOps::copy(rhoByEps0, a_rho);
  DataOps::scale(rhoByEps0, 1.0 / Units::eps0);

  // The operator rows are kappa-weighted, so the right-hand side must be too. As in FieldSolverGMG, the space charge
  // is either kappa-weighted already or weighted here.
  if (m_kappaSource) {
    DataOps::kappaScale(rhoByEps0,
                        m_amr->getVofIterator(m_realm, phase::gas),
                        m_amr->getVofIterator(m_realm, phase::solid));
  }

  Vec b;
  Vec x;
  Vec r;

  petscGrid->create(b);
  petscGrid->create(x);
  petscGrid->create(r);

  // The jump condition is eps_0 * D_0 + eps_1 * D_1 = sigma/eps0, as in FieldSolverGMG.
  EBAMRIVData sigmaByEps0;
  m_amr->allocate(sigmaByEps0, m_realm, phase::gas, m_nComp);
  DataOps::copy(sigmaByEps0, a_sigma);
  DataOps::scale(sigmaByEps0, 1.0 / Units::eps0, m_amr->getVofIterator(m_realm, phase::gas));

  petscGrid->putChomboInPetsc(b, rhoByEps0);
  m_petscOp->addBoundaryConditionsToRHS(b, sigmaByEps0);

  // The matrix rows are weighted per level (PetscHelmholtzOp::getRowWeight), and the right-hand side must be too. The
  // unweighted norm of b is kept for reporting the residual as it was before the weighting.
  PetscReal zeroResidNormUnweighted;

  PetscCallAbort(Chombo_MPI::comm, VecNorm(b, NORM_2, &zeroResidNormUnweighted));

  m_petscOp->scaleRows(b, false);

  if (a_zeroPhi) {
    petscGrid->setValue(x, 0.0);
  }
  else {
    petscGrid->putChomboInPetsc(x, a_potential);
  }

  // Residual r = b - A*x. As in FieldSolverGMG's Krylov solvers, the residual is reported in the 2-norm and relative
  // to the zero-residual b - A*0 = b, which does not depend on the initial guess. This is the true residual of the
  // matrix, not the (possibly preconditioned) residual that the KSP tests convergence on.
  PetscReal residNorm;
  PetscReal zeroResidNorm;

  PetscCallAbort(Chombo_MPI::comm, VecNorm(b, NORM_2, &zeroResidNorm));

  m_monitorZeroResidNorm = zeroResidNorm;

  PetscCallAbort(Chombo_MPI::comm, KSPSetInitialGuessNonzero(m_ksp, a_zeroPhi ? PETSC_FALSE : PETSC_TRUE));
  PetscCallAbort(Chombo_MPI::comm, KSPSolve(m_ksp, b, x));

  KSPConvergedReason reason;
  PetscInt           numIter;

  PetscCallAbort(Chombo_MPI::comm, KSPGetConvergedReason(m_ksp, &reason));
  PetscCallAbort(Chombo_MPI::comm, KSPGetIterationNumber(m_ksp, &numIter));

  const bool converged = (reason > 0);

  PetscReal residNormUnweighted;

  PetscCallAbort(Chombo_MPI::comm, MatResidual(m_petscOp->getMatrix(), b, x, r));
  PetscCallAbort(Chombo_MPI::comm, VecNorm(r, NORM_2, &residNorm));

  // Remove the row weights, so that the residual is that of the unweighted divergence.
  m_petscOp->scaleRows(r, true);

  PetscCallAbort(Chombo_MPI::comm, VecNorm(r, NORM_2, &residNormUnweighted));

  if (m_verbosity > 1) {
    KSPType kspType;

    PetscCallAbort(Chombo_MPI::comm, KSPGetType(m_ksp, &kspType));

    const Real relResid           = (zeroResidNorm > 0.0) ? residNorm / zeroResidNorm : residNorm;
    const Real relResidUnweighted = (zeroResidNormUnweighted > 0.0) ? residNormUnweighted / zeroResidNormUnweighted
                                                                    : residNormUnweighted;

    pout() << "FieldSolverAMG - " << kspType << (converged ? " converged" : " did NOT converge")
           << " (rel.res = " << relResid << ", unweighted rel.res = " << relResidUnweighted << ", |r| = " << residNorm
           << ", |r_zero| = " << zeroResidNorm << ", iters = " << numIter
           << ", reason = " << KSPConvergedReasons[reason] << ")" << endl;
  }

  // The solver and the preconditioner, e.g. the levels of the algebraic multigrid hierarchy. PETSc writes this to
  // standard output.
  if (m_petscVerbosity >= 3) {
    PetscCallAbort(Chombo_MPI::comm, KSPView(m_ksp, PETSC_VIEWER_STDOUT_(Chombo_MPI::comm)));
  }

  // Closures whose jump data is an output (saturation charge) compute it from the solution and keep it.
  m_petscJumpBc->updateJumpData(x, m_permittivityEB);

  petscGrid->putPetscInChombo(a_potential, x);
  petscGrid->putPetscInChombo(m_residue, r);

  petscGrid->destroy(b);
  petscGrid->destroy(x);
  petscGrid->destroy(r);

  // Coarsen and update ghost cells before computing the field, as FieldSolverGMG does.
  m_amr->conservativeAverage(a_potential, m_realm);
  m_amr->interpGhostPwl(a_potential, m_realm);

  this->computeElectricField(m_electricField, a_potential);

  // With the saturation charge the surface charge density is an output of the solve, which FieldSolverGMG also copies
  // into m_sigma.
  if (m_jumpBcType == JumpBCType::SaturationCharge) {
    const auto& jumpBc = static_cast<const PetscHelmholtzSaturationChargeJumpBC&>(*m_petscJumpBc);

    DataOps::copy(m_sigma, jumpBc.getSigma());
    DataOps::scale(m_sigma, Units::eps0, m_amr->getVofIterator(m_realm, phase::gas));

    m_amr->conservativeAverage(m_sigma, m_realm, phase::gas);
  }

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

  // Electrode and interface conditions. One order and one weight apply to all EB face stencils.
  m_petscEbBc = RefCountedPtr<PetscHelmholtzEBBC>(new PetscHelmholtzElectrostaticEBBC(m_ebBc));

  switch (m_jumpBcType) {
  case JumpBCType::Natural: {
    m_petscJumpBc = RefCountedPtr<PetscHelmholtzJumpBC>(new PetscHelmholtzNaturalJumpBC());

    break;
  }
  case JumpBCType::SaturationCharge: {
    m_petscJumpBc = RefCountedPtr<PetscHelmholtzJumpBC>(new PetscHelmholtzSaturationChargeJumpBC(phase::gas));

    break;
  }
  default: {
    MayDay::Abort("FieldSolverAMG::setupSolver -- logic bust");

    break;
  }
  }

  m_petscEbBc->setOrder(m_ebOrder);
  m_petscEbBc->setWeight(m_ebWeight);
  m_petscJumpBc->setOrder(m_ebOrder);
  m_petscJumpBc->setWeight(m_ebWeight);

  // Impose the interface conditions on the dielectric surface rather than at the EB face centroids.
  if (m_interfacePoints) {
    m_petscJumpBc->setInterface(m_computationalGeometry->getSolidImplicitFunction());
  }

  // The coarse-fine interpolators of both phases on each level, in the phase order of the PetscGrid.
  const RefCountedPtr<EBIndexSpace>& ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  const int finestLevel = m_amr->getFinestLevel();
  const int numPhases   = m_multifluidIndexSpace->numPhases();

  Vector<MFMultigridInterpolator> interpolators(1 + finestLevel);

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    Vector<RefCountedPtr<EBMultigridInterpolator>> interpPhases(numPhases);

    if (!(ebisGas.isNull())) {
      interpPhases[phase::gas] = m_amr->getMultigridInterpolator(m_realm, phase::gas)[lvl];
    }
    if (!(ebisSol.isNull())) {
      interpPhases[phase::solid] = m_amr->getMultigridInterpolator(m_realm, phase::solid)[lvl];
    }

    interpolators[lvl].define(interpPhases);
  }

  m_petscOp = RefCountedPtr<PetscHelmholtzOp>(new PetscHelmholtzOp());
  m_petscOp->define(m_amr->getPetscGrid(m_realm),
                    interpolators,
                    domainBc,
                    m_petscEbBc,
                    m_petscJumpBc,
                    m_permittivityCell,
                    m_permittivityFace,
                    m_permittivityEB,
                    m_amr->getProbLo(),
                    m_amr->getDx(),
                    m_alpha,
                    m_beta);

  // Cut cells whose least-squares fits were all rejected have no flux through their EB face, and their potential
  // floats. The electrode and the interface conditions compute the stencils of the same cells, so one count covers
  // both.
  const long long numWithoutStencil = m_petscEbBc->getNumCellsWithoutStencil();

  if (numWithoutStencil > 0 && m_verbosity > 0) {
    pout() << "FieldSolverAMG::setupSolver -- " << numWithoutStencil
           << " cut cells (counted per phase) have no EB stencil, and no flux through their EB face" << endl;
  }

  // The settings from the input file are the defaults. They can be overridden with PETSc options that carry the prefix
  // 'fieldsolveramg_', given through the PETSC_OPTIONS environment variable (or, for options without a value, on the
  // command line), e.g. -fieldsolveramg_ksp_monitor.
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

  // No effect unless the KSP is GMRES or FGMRES.
  PetscCallAbort(Chombo_MPI::comm, KSPGMRESSetRestart(m_ksp, m_kspGmresRestart));

  // The level smoothers of a multigrid preconditioner are only created when the preconditioner is set up, so their
  // settings go into the options database, where they are picked up from. They are written on every setup, so that a
  // change of the input options takes effect, except where the user gave the option to PETSc (PETSC_OPTIONS or the
  // command line), which is recorded before this class first writes it. All FieldSolverAMG instances share the
  // 'fieldsolveramg_' prefix, and the smoothers read the options when the preconditioner is set up, at the first solve.
  // So with two instances whose amg_smoother_* or gamg_* inputs differ, the last setupSolver before a solve decides:
  // last writer wins. With identical inputs, the normal case, this does not matter. Preconditioners that are not
  // multigrid ignore these options.
  auto isUserOption = [](const char* a_name) -> bool {
    PetscBool isSet = PETSC_FALSE;

    PetscCallAbort(Chombo_MPI::comm, PetscOptionsHasName(nullptr, nullptr, a_name, &isSet));

    return isSet == PETSC_TRUE;
  };

  static const bool userSmootherKsp = isUserOption("-fieldsolveramg_mg_levels_ksp_type");
  static const bool userSmootherPc  = isUserOption("-fieldsolveramg_mg_levels_pc_type");
  static const bool userSmootherIts = isUserOption("-fieldsolveramg_mg_levels_ksp_max_it");
  static const bool userThreshold   = isUserOption("-fieldsolveramg_pc_gamg_threshold");
  static const bool userAggressive  = isUserOption("-fieldsolveramg_pc_gamg_aggressive_coarsening");

  if (!userSmootherKsp) {
    PetscCallAbort(Chombo_MPI::comm,
                   PetscOptionsSetValue(nullptr, "-fieldsolveramg_mg_levels_ksp_type", m_amgSmootherKsp.c_str()));
  }
  if (!userSmootherPc) {
    PetscCallAbort(Chombo_MPI::comm,
                   PetscOptionsSetValue(nullptr, "-fieldsolveramg_mg_levels_pc_type", m_amgSmootherPc.c_str()));
  }
  if (!userSmootherIts) {
    PetscCallAbort(
      Chombo_MPI::comm,
      PetscOptionsSetValue(nullptr, "-fieldsolveramg_mg_levels_ksp_max_it", std::to_string(m_amgSmootherIts).c_str()));
  }
  if (!userThreshold) {
    PetscCallAbort(
      Chombo_MPI::comm,
      PetscOptionsSetValue(nullptr, "-fieldsolveramg_pc_gamg_threshold", std::to_string(m_gamgThreshold).c_str()));
  }
  if (!userAggressive) {
    PetscCallAbort(Chombo_MPI::comm,
                   PetscOptionsSetValue(nullptr,
                                        "-fieldsolveramg_pc_gamg_aggressive_coarsening",
                                        std::to_string(m_gamgAggressiveCoarsening).c_str()));
  }

  PetscCallAbort(Chombo_MPI::comm, KSPSetOptionsPrefix(m_ksp, "fieldsolveramg_"));
  PetscCallAbort(Chombo_MPI::comm, KSPSetFromOptions(m_ksp));

  if (m_petscVerbosity >= 1) {
    PetscCallAbort(Chombo_MPI::comm, KSPMonitorSet(m_ksp, FieldSolverAMG::monitorResidual, this, nullptr));
  }

  m_isSolverSetup = true;
}

PetscErrorCode
FieldSolverAMG::monitorResidual(KSP a_ksp, PetscInt a_iteration, PetscReal a_residualNorm, void* a_context)
{
  CH_TIME("FieldSolverAMG::monitorResidual");

  FieldSolverAMG* solver = static_cast<FieldSolverAMG*>(a_context);

  if (a_iteration == 0) {
    solver->m_monitorFirstNorm    = a_residualNorm;
    solver->m_monitorPreviousNorm = a_residualNorm;
  }

  const Real scale = (solver->m_monitorFirstNorm > 0.0) ? solver->m_monitorFirstNorm : 1.0;

  pout() << "      FieldSolverAMG:: iteration = " << a_iteration << ", rel. residual = " << a_residualNorm / scale;

  if (a_iteration > 0) {
    pout() << ", rate = " << ((a_residualNorm > 0.0) ? solver->m_monitorPreviousNorm / a_residualNorm : 0.0);
  }

  // The true residual is built from the current solution, which for GMRES costs a solution update and a matrix-vector
  // product.
  if (solver->m_petscVerbosity >= 2) {
    Vec       residual = nullptr;
    PetscReal trueNorm;

    PetscCall(KSPBuildResidual(a_ksp, nullptr, nullptr, &residual));
    PetscCall(VecNorm(residual, NORM_2, &trueNorm));
    PetscCall(VecDestroy(&residual));

    const Real zeroResid = solver->m_monitorZeroResidNorm;

    pout() << ", true rel. residual = " << ((zeroResid > 0.0) ? trueNorm / zeroResid : trueNorm);
  }

  pout() << endl;

  solver->m_monitorPreviousNorm = a_residualNorm;

  return PETSC_SUCCESS;
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

  if (!m_petscOp.isNull()) {
    m_petscOp->clear();
  }

  m_petscOp       = RefCountedPtr<PetscHelmholtzOp>();
  m_petscEbBc     = RefCountedPtr<PetscHelmholtzEBBC>();
  m_petscJumpBc   = RefCountedPtr<PetscHelmholtzJumpBC>();
  m_isSolverSetup = false;
}

void
FieldSolverAMG::setSolverPermittivities(const MFAMRCellData& a_permittivityCell,
                                        const MFAMRFluxData& a_permittivityFace,
                                        const MFAMRIVData&   a_permittivityEB)
{
  CH_TIME("FieldSolverAMG::setSolverPermittivities()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::setSolverPermittivities()" << endl;
  }

  // As in FieldSolverGMG, the solver must be set up first; the coefficients then replace the ones it was set up with
  // until the next setup (e.g. after a regrid). The matrix is refilled with the same sparsity pattern, which depends on
  // the geometry only, and handing it to the KSP again makes PETSc rebuild the preconditioner at the next solve.
  if (!m_isSolverSetup) {
    MayDay::Error("FieldSolverAMG::setSolverPermittivities -- must set up solver first!");
  }

  m_petscOp->setCoefficients(a_permittivityCell, a_permittivityFace, a_permittivityEB);
  m_petscOp->fillMatrix();

  PetscCallAbort(Chombo_MPI::comm, KSPSetOperators(m_ksp, m_petscOp->getMatrix(), m_petscOp->getMatrix()));
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
