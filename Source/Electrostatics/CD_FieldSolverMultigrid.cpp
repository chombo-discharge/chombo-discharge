/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   FieldSolverMultigrid.cpp
  @brief  Implementation of FieldSolverMultigrid.H
  @author Robert Marskar
  @todo   Once the new operator is in, check the computeLoads routine. 
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_FieldSolverMultigrid.H>
#include <CD_DataOps.H>
#include <CD_MultifluidAlias.H>
#include <CD_MFMultigridInterpolator.H>
#include <CD_MFHelmholtzElectrostaticDomainBCFactory.H>
#include <CD_MFHelmholtzElectrostaticEBBCFactory.H>
#include <CD_MFHelmholtzJumpBCFactory.H>
#include <CD_MFHelmholtzSaturationChargeJumpBCFactory.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

constexpr Real FieldSolverMultigrid::m_alpha;
constexpr Real FieldSolverMultigrid::m_beta;

FieldSolverMultigrid::FieldSolverMultigrid() : FieldSolver()
{
  CH_TIME("FieldSolverMultigrid::FieldSolverMultigrid()");

  // Default settings
  m_isSolverSetup = false;
  m_className     = "FieldSolverMultigrid";
}

FieldSolverMultigrid::~FieldSolverMultigrid()
{
  CH_TIME("FieldSolverMultigrid::~FieldSolverMultigrid()");
}

void
FieldSolverMultigrid::parseOptions()
{
  CH_TIME("FieldSolverMultigrid::parseOptions()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::parseOptions()" << endl;
  }

  this->parseVerbosity();
  this->parseDomainBc();
  this->parsePlotVariables();
  this->parseMultigridSettings();
  this->parseKappaSource();
  this->parseJumpBC();
  this->parseRegridSlopes();
}

void
FieldSolverMultigrid::parseRuntimeOptions()
{
  CH_TIME("FieldSolverMultigrid::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::parseRuntimeOptions()" << endl;
  }

  this->parseVerbosity();
  this->parseMultigridSettings();
  this->parseKappaSource();
  this->parsePlotVariables();
  this->parseRegridSlopes();
}

void
FieldSolverMultigrid::parseKappaSource()
{
  CH_TIME("FieldSolverMultigrid::parseKappaSource()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::parseKappaSource()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("kappa_source", m_kappaSource);
}

void
FieldSolverMultigrid::parseMultigridSettings()
{
  CH_TIME("FieldSolverMultigrid::parseMultigridSettings()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::parseMultigridSettings()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("gmg_verbosity", m_multigridVerbosity);
  pp.get("gmg_pre_smooth", m_multigridPreSmooth);
  pp.get("gmg_post_smooth", m_multigridPostSmooth);
  pp.get("gmg_bott_smooth", m_multigridBottomSmooth);
  pp.get("gmg_max_iter", m_multigridMaxIterations);
  pp.get("gmg_min_iter", m_multigridMinIterations);
  pp.get("gmg_exit_tol", m_multigridExitTolerance);
  pp.get("gmg_exit_hang", m_multigridExitHang);
  pp.get("gmg_min_cells", m_minCellsBottom);
  pp.get("gmg_drop_order", m_domainDropOrder);
  pp.get("gmg_bc_order", m_multigridBcOrder);
  pp.get("gmg_bc_weight", m_multigridBcWeight);
  pp.get("gmg_jump_order", m_multigridJumpOrder);
  pp.get("gmg_jump_weight", m_multigridJumpWeight);

  // Fetch the desired bottom solver from the input script. We look for things like FieldSolverMultigrid.gmg_bottom_solver = bicgstab or '= simple <number>'
  // where <number> is the number of relaxation for the smoothing solver.
  const int num = pp.countval("gmg_bottom_solver");
  if (num == 1) {
    pp.get("gmg_bottom_solver", str);
    if (str == "bicgstab") {
      m_bottomSolverType = BottomSolverType::BiCGStab;
    }
    else if (str == "gmres") {
      m_bottomSolverType = BottomSolverType::GMRES;
    }
    else {
      MayDay::Error(
        "FieldSolverMultigrid::parseMultigridSettings() - logic bust, you've specified one parameter and I expected either 'bicgstab' or 'gmres'");
    }
  }
  else if (num == 2) {
    int numSmooth;
    pp.get("gmg_bottom_solver", str, 0);
    pp.get("gmg_bottom_solver", numSmooth, 1);
    if (str == "simple") {
      m_bottomSolverType = BottomSolverType::Simple;
      m_mfsolver.setNumSmooths(numSmooth);
    }
    else {
      MayDay::Error(
        "FieldSolverMultigrid::parseMultigridSettings() - logic bust, you've specified two parameters and I expected 'simple <number>'");
    }
  }
  else {
    MayDay::Error(
      "FieldSolverMultigrid::parseMultigridSettings() - logic bust in bottom solver. You must specify ' = bicgstab', ' = gmres', or ' = simple <number>'");
  }

  // Get a string for the multigrid smoother. This must either be "jacobi", "red_black", or "multi_color".
  pp.get("gmg_smoother", str);
  if (str == "jacobi") {
    m_multigridRelaxMethod = MFHelmholtzOp::Smoother::PointJacobi;
  }
  else if (str == "red_black") {
    m_multigridRelaxMethod = MFHelmholtzOp::Smoother::GauSaiRedBlack;
  }
  else if (str == "multi_color") {
    m_multigridRelaxMethod = MFHelmholtzOp::Smoother::GauSaiMultiColor;
  }
  else {
    MayDay::Error("FieldSolverMultigrid::parseMultigridSettings() - unsupported relaxation method requested");
  }

  // Cycle type
  pp.get("gmg_cycle", str);
  if (str == "vcycle") {
    m_multigridType = MultigridType::VCycle;
  }
  else {
    MayDay::Error(
      "FieldSolverMultigrid::parseMultigridSettings - unsupported multigrid cycle type requested. Only vcycle supported for now. ");
  }

  // No lower than 2.
  if (m_minCellsBottom < 2) {
    m_minCellsBottom = 2;
  }

  // Things won't run unless this is fulfilled.
  CH_assert(m_minCellsBottom >= 2);
  CH_assert(m_multigridPreSmooth >= 0);
  CH_assert(m_multigridPostSmooth >= 0);
  CH_assert(m_multigridBottomSmooth >= 0);
  CH_assert(m_multigridBcOrder > 0);
  CH_assert(m_multigridBcWeight >= 0);
  CH_assert(m_multigridJumpOrder > 0);
  CH_assert(m_multigridJumpWeight >= 0);
}

void
FieldSolverMultigrid::parseJumpBC()
{
  CH_TIME("FieldSolverMultigrid::parseJumpBC()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::parseJumpBC()" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("jump_bc", str);

  if (str == "natural") {
    m_jumpBcType = JumpBCType::Natural;
  }
  else if (str == "saturation_charge") {
    m_jumpBcType = JumpBCType::SaturationCharge;
  }
  else {
    const std::string errorString = "FieldSolverMultigrid::parseJumpBC -- error, argument '" + str + "' not recognized";
    MayDay::Error(errorString.c_str());
  }
}

bool
FieldSolverMultigrid::solve(MFAMRCellData&       a_phi,
                            const MFAMRCellData& a_rho,
                            const EBAMRIVData&   a_sigma,
                            const bool           a_zeroPhi)
{
  CH_TIMERS("FieldSolverMultigrid::solve");
  CH_TIMER("FieldSolverMultigrid::solve::alloc_temps", t1);
  CH_TIMER("FieldSolverMultigrid::solve::set_temps", t2);
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::solve(MFAMRCellData, MFAMRCellData, EBAMRIVData, bool)" << endl;
  }

  // TLDR: This is the main solve routine. The operator factory was set up as kappa*L(phi) = -kappa*div(eps*grad(phi))
  //       which we use to solve the Poisson equation div(eps*grad(phi)) = -rho/eps0.
  //
  //       So, we must scale rho (which is defined on centroids) by kappa and divide by eps0. The minus-sign provided in the operator
  //       through the b-coefficient. (So really, we're actually solving kappa*[-div(eps*grad(phi))] = kappa*(rho/eps0)

  CH_assert(m_isVoltageSet);
  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(a_rho[0]->nComp() == 1);
  CH_assert(a_sigma[0]->nComp() == 1);

  bool converged = false;

  // Set up multigrid solver if it is not already done.
  if (!m_isSolverSetup) {
    this->setupSolver();
  }

  // Define temporaries; the incoming data might need to be scaled but we don't want to
  // alter it directly.
  CH_START(t1);
  MFAMRCellData zero;
  MFAMRCellData kappaRhoByEps0;
  EBAMRIVData   sigmaByEps0;

  m_amr->allocate(zero, m_realm, m_nComp);
  m_amr->allocate(kappaRhoByEps0, m_realm, m_nComp);
  m_amr->allocate(sigmaByEps0, m_realm, phase::gas, m_nComp);
  CH_STOP(t1);

  // Scale data as appropriate.
  CH_START(t2);
  DataOps::setValue(zero, 0.0);

  // Do the scaled space charge density.
  DataOps::copy(kappaRhoByEps0, a_rho);
  DataOps::scale(kappaRhoByEps0, 1. / (Units::eps0));

  // Special flag for when a_rho is on the centroid but was not scaled by kappa on input. The multigrid operator solves
  // kappa*L(phi) = kappa*rho so the right-hand side must be kappa-weighted.
  if (m_kappaSource) {
    DataOps::kappaScale(kappaRhoByEps0);
  }

  m_amr->conservativeAverage(a_phi, m_realm);
  m_amr->interpGhost(a_phi, m_realm);
  m_amr->arithmeticAverage(kappaRhoByEps0, m_realm);
  m_amr->interpGhost(kappaRhoByEps0, m_realm);

  // Do the scaled surface charge
  DataOps::copy(sigmaByEps0, a_sigma);
  DataOps::scale(sigmaByEps0, 1. / (Units::eps0));
  CH_STOP(t2);

  // Factory needs knowledge of the new surface charge -- it passes this data by reference to the multigrid operators (MFHelmholtzOps).
  m_helmholtzOpFactory->setJump(sigmaByEps0, 1.0);

  // Aliasing, because Chombo is not too smart about smart pointers.
  Vector<LevelData<MFCellFAB>*> phi; // Raw pointers of a_phi
  Vector<LevelData<MFCellFAB>*> rhs; // Raw pointers of m_kappaRhoByEps0
  Vector<LevelData<MFCellFAB>*> res; // Raw pointers of m_residue
  Vector<LevelData<MFCellFAB>*> zer; // Raw pointers of m_zero

  m_amr->alias(phi, a_phi);
  m_amr->alias(rhs, kappaRhoByEps0);
  m_amr->alias(res, m_residue);
  m_amr->alias(zer, zero);

  // GMG solve. Use phi = zero as initial metric. We want to reduce this by m_multigridExitTolerance, in which case
  // multigrid has "converged".
  const int coarsestLevel = 0;
  const int finestLevel   = m_amr->getFinestLevel();

  // This is the residue rho - L(phi)
  const Real phiResid = m_multigridSolver->computeAMRResidual(phi, rhs, finestLevel, 0);

  // This is the residue rho - L(phi=0)
  const Real zeroResid = m_multigridSolver->computeAMRResidual(zer, rhs, finestLevel, 0);

  // Convergence criterion.
  const Real convergedResid = zeroResid * m_multigridExitTolerance;

  // If the residue rho - L(phi) is too large then we must get a new solution.
  if (phiResid > convergedResid) {
    m_multigridSolver->m_convergenceMetric = zeroResid;
    m_multigridSolver->solveNoInitResid(phi, res, rhs, finestLevel, coarsestLevel, a_zeroPhi);

    const int status = m_multigridSolver->m_exitStatus; // 1 => Initial norm sufficiently reduced
    if (status == 1 || status == 8) {                   // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else {
    converged = true;
  }

  m_multigridSolver->revert(phi, rhs, finestLevel, 0);

  // Coarsen/update ghosts before computing the field.
  m_amr->conservativeAverage(a_phi, m_realm);
  m_amr->interpGhostPwl(a_phi, m_realm);

  this->computeElectricField(m_electricField, a_phi);

  // If we are also solving for the saturation charge we get that solution from the factory (it can be a free parameter in the Helmholtz solve).
  if (m_jumpBcType == JumpBCType::SaturationCharge) {
    const EBAMRIVData& factorySigma = m_helmholtzOpFactory->getSigma();
    DataOps::copy(m_sigma, factorySigma);
    DataOps::scale(m_sigma, Units::eps0);

    m_amr->conservativeAverage(m_sigma, m_realm, phase::gas);
  }

  return converged;
}

void
FieldSolverMultigrid::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("FieldSolverMultigrid::preRegrid(int, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::preRegrid(int, int)" << endl;
  }

  FieldSolver::preRegrid(a_lbase, a_oldFinestLevel);

  m_multigridSolver.freeMem();
  m_helmholtzOpFactory.freeMem();
}

void
FieldSolverMultigrid::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("FieldSolverMultigrid::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::regrid(int, int, int)" << endl;
  }

  FieldSolver::regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  m_isSolverSetup = false;
}

void
FieldSolverMultigrid::registerOperators()
{
  CH_TIME("FieldSolverMultigrid::registerOperators()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::registerOperators()" << endl;
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
  m_amr->registerOperator(s_eb_flux_reg, m_realm, phase::gas);
  m_amr->registerOperator(s_eb_flux_reg, m_realm, phase::solid);
  m_amr->registerOperator(s_eb_multigrid, m_realm, phase::gas);
  m_amr->registerOperator(s_eb_multigrid, m_realm, phase::solid);
}

void
FieldSolverMultigrid::setupSolver()
{
  CH_TIME("FieldSolverMultigrid::setupSolver()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::setupSolver()" << endl;
  }

  this->setupHelmholtzFactory(); // Set up the operator factory
  this->setupMultigrid();        // Set up the AMR multigrid solver

  m_isSolverSetup = true;
}

void
FieldSolverMultigrid::setSolverPermittivities(const MFAMRCellData& a_permittivityCell,
                                              const MFAMRFluxData& a_permittivityFace,
                                              const MFAMRIVData&   a_permittivityEB)
{
  CH_TIME("FieldSolverMultigrid::setSolverPermittivities()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::setSolverPermittivities()" << endl;
  }

  if (!m_isSolverSetup) {
    MayDay::Error("FieldSolverMultigrid::setSolverPermittivities -- must set up solver first!");
  }

  // Get the AMR operators and update the coefficients.
  Vector<AMRLevelOp<LevelData<MFCellFAB>>*>& operatorsAMR = m_multigridSolver->getAMROperators();

  CH_assert(operatorsAMR.size() == 1 + m_amr->getFinestLevel());

  // Set coefficients for the AMR levels.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    CH_assert(!(operatorsAMR[lvl] == nullptr));

    MFHelmholtzOp& op = static_cast<MFHelmholtzOp&>(*operatorsAMR[lvl]);

    op.setAcoAndBco(a_permittivityCell[lvl], a_permittivityFace[lvl], a_permittivityEB[lvl]);
  }

  // Get the deeper multigrid levels and coarsen onto that data as well. Strictly speaking, we don't
  // have to do this but it facilitates multigrid convergence and is therefore good practice. The operator
  // factory has routines for the coefficients that belong to the multigrid levels. The factory does not
  // have access to the operator, so we fetch those using AMRMultiGrid and call setAcoAndBco from there.
  // access to the operator.
  m_helmholtzOpFactory->coarsenCoefficientsMG();
  Vector<Vector<MGLevelOp<LevelData<MFCellFAB>>*>> operatorsMG = m_multigridSolver->getOperatorsMG();

  for (int amrLevel = 0; amrLevel < operatorsMG.size(); amrLevel++) {
    for (int mgLevel = 0; mgLevel < operatorsMG[amrLevel].size(); mgLevel++) {
      CH_assert(!(operatorsMG[amrLevel][mgLevel] == nullptr));

      MFHelmholtzOp& op = static_cast<MFHelmholtzOp&>(*operatorsMG[amrLevel][mgLevel]);

      op.setAcoAndBco(op.getAcoef(), op.getBcoef(), op.getBcoefIrreg());
    }
  }
}

void
FieldSolverMultigrid::setPermittivities()
{
  CH_TIME("FieldSolverMultigrid::setPermittivities()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::setPermittivities()" << endl;
  }

  // Parent method fills permittivities over the "valid" region.
  FieldSolver::setPermittivities();

#if 1
  // With EBHelmholtzOp/MFHelmholtzOp, the stencils can reach out of grid
  // patches and into ghost faces when computing the centroid flux on a cut-cell
  // face. To be on the safe side, we fill ghost cells for the cell-centered
  // permittivity and average that onto the grid faces. Importantly, we ensure
  // that we also fill one layer of "ghost faces".

  const Average  average  = Average::Arithmetic;
  const int      tanGhost = 1;
  const Interval interv   = Interval(0, 0);

  EBAMRCellData permCellGas;
  EBAMRCellData permCellSol;
  EBAMRFluxData permFluxGas;
  EBAMRFluxData permFluxSol;
  EBAMRIVData   permEBGas;
  EBAMRIVData   permEBSol;

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Do the gas-phase.
  if (!(ebisGas.isNull())) {
    permCellGas = m_amr->alias(phase::gas, m_permittivityCell);
    permFluxGas = m_amr->alias(phase::gas, m_permittivityFace);
    permEBGas   = m_amr->alias(phase::gas, m_permittivityEB);

    // Coarsen cell and EB data.
    m_amr->average(permCellGas, m_realm, phase::gas, average);
    m_amr->average(permEBGas, m_realm, phase::gas, average);

    // Interpolate cell-centered permittivities to ghost cells; then average
    // cell-centered data to faces, including one ghost face.
    m_amr->interpGhost(permCellGas, m_realm, phase::gas);

    DataOps::averageCellToFace(permFluxGas, permCellGas, m_amr->getDomains(), tanGhost, interv, interv, average);

    m_amr->average(permFluxGas, m_realm, phase::gas, average);
  }

  if (!(ebisSol.isNull())) {
    permCellSol = m_amr->alias(phase::solid, m_permittivityCell);
    permFluxSol = m_amr->alias(phase::solid, m_permittivityFace);
    permEBSol   = m_amr->alias(phase::solid, m_permittivityEB);

    // Coarsen cell and EB data.
    m_amr->average(permCellSol, m_realm, phase::solid, average);
    m_amr->average(permEBSol, m_realm, phase::solid, average);

    // Interpolate cell-centered permittivities to ghost cells; then average
    // cell-centered data to faces, including one ghost face.
    m_amr->interpGhost(permCellSol, m_realm, phase::solid);

    DataOps::averageCellToFace(permFluxSol, permCellSol, m_amr->getDomains(), tanGhost, interv, interv, average);

    m_amr->average(permFluxSol, m_realm, phase::solid, average);
  }
#endif
}

void
FieldSolverMultigrid::setupHelmholtzFactory()
{
  CH_TIME("FieldSolverMultigrid::setupHelmholtzFactory()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::setupHelmholtzFactory()" << endl;
  }

  // TLDR: This routine sets up a Helmholtz factory for creating MFHelmholtzOps which are used by Chombo's AMRMultiGrid. We
  //       should already have registered the operators when we come to this routine.

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // First, set up multifluid wrappers. These wrappers are just code-saving constructs for multifluid problems -- they
  // are essentially just a size-two vector containing the operators for each phase.
  const int finestLevel = m_amr->getFinestLevel();
  const int numPhases   = m_multifluidIndexSpace->numPhases();

  Vector<MFLevelGrid>             mflg(1 + finestLevel);
  Vector<MFMultigridInterpolator> mfInterp(1 + finestLevel);
  Vector<MFReflux>                mfFluxReg(1 + finestLevel);
  Vector<MFCoarAve>               mfCoarAve(1 + finestLevel);

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    Vector<EBLevelGrid>                            eblgPhases(numPhases);
    Vector<RefCountedPtr<EBMultigridInterpolator>> interpPhases(numPhases);
    Vector<RefCountedPtr<EBReflux>>                fluxRegPhases(numPhases);
    Vector<RefCountedPtr<EBCoarAve>>               avePhases(numPhases);

    if (!ebisGas.isNull()) {
      eblgPhases[phase::gas] = *(m_amr->getEBLevelGrid(m_realm, phase::gas)[lvl]);
    }
    if (!ebisSol.isNull()) {
      eblgPhases[phase::solid] = *(m_amr->getEBLevelGrid(m_realm, phase::solid)[lvl]);
    }

    if (!ebisGas.isNull()) {
      interpPhases[phase::gas] = (m_amr->getMultigridInterpolator(m_realm, phase::gas)[lvl]);
    }
    if (!ebisSol.isNull()) {
      interpPhases[phase::solid] = (m_amr->getMultigridInterpolator(m_realm, phase::solid)[lvl]);
    }

    if (!ebisGas.isNull()) {
      fluxRegPhases[phase::gas] = (m_amr->getFluxRegister(m_realm, phase::gas)[lvl]);
    }
    if (!ebisSol.isNull()) {
      fluxRegPhases[phase::solid] = (m_amr->getFluxRegister(m_realm, phase::solid)[lvl]);
    }

    if (!ebisGas.isNull()) {
      avePhases[phase::gas] = (m_amr->getCoarseAverage(m_realm, phase::gas)[lvl]);
    }
    if (!ebisSol.isNull()) {
      avePhases[phase::solid] = (m_amr->getCoarseAverage(m_realm, phase::solid)[lvl]);
    }

    mflg[lvl].define(m_multifluidIndexSpace, eblgPhases);
    mfInterp[lvl].define(interpPhases);
    mfFluxReg[lvl].define(fluxRegPhases);
    mfCoarAve[lvl].define(avePhases);
  }

  // Find the coarsest domain used for multigrid. The user will have specified the minimum number of cells in any
  // coordinate direction (m_minCellsBottom), and we coarsen until we have a domain which satisfies that constraint.
  // MFHelmholtzOpFactory will not generate multigrid operators below that depth.
  ProblemDomain bottomDomain = m_amr->getDomains()[0];
  while (bottomDomain.domainBox().shortside() >= 2 * m_minCellsBottom) {
    bottomDomain.coarsen(2);
  }

  // Number of ghost cells in data holders. This is needed because EBHelmholtzOp uses restriction/prolongation objects from Chombo, and they
  // need to know explicit object size to optimize irregular access.
  const IntVect ghostPhi = m_amr->getNumberOfGhostCells() * IntVect::Unit;
  const IntVect ghostRhs = m_amr->getNumberOfGhostCells() * IntVect::Unit;

  // BC factories. Fortunately, MFHelmholtzOp/EBHelmholtzOp have sane interfaces which allowed us to code up
  // boundary condition objects and factories.
  auto ebbcFactory = RefCountedPtr<MFHelmholtzElectrostaticEBBCFactory>(
    new MFHelmholtzElectrostaticEBBCFactory(m_multigridBcOrder, m_multigridBcWeight, m_ebBc));

  auto domainBcFactory = RefCountedPtr<MFHelmholtzDomainBCFactory>(
    new MFHelmholtzElectrostaticDomainBCFactory(m_domainBc));

  // Set the BC jump factory. This is either the "natural" factory or the saturation charge BC.
  RefCountedPtr<MFHelmholtzJumpBCFactory> jumpBcFactory;
  switch (m_jumpBcType) {
  case JumpBCType::Natural: {
    jumpBcFactory = RefCountedPtr<MFHelmholtzJumpBCFactory>(new MFHelmholtzJumpBCFactory());

    break;
  }
  case JumpBCType::SaturationCharge: {
    jumpBcFactory = RefCountedPtr<MFHelmholtzJumpBCFactory>(new MFHelmholtzSaturationChargeJumpBCFactory(phase::gas));

    break;
  }
  }

  // Drop order stuff for EB stencils -- can facilitate safer GMG relaxations
  // on deeper multigrid levels.
  ebbcFactory->setDomainDropOrder(m_domainDropOrder);
  jumpBcFactory->setDomainDropOrder(m_domainDropOrder);

  // Create the factory. Note that we pass m_permittivityCell in through the a-coefficient, but we also set alpha to zero
  // so there is no diagonal term in the operator after all.
  m_helmholtzOpFactory = RefCountedPtr<MFHelmholtzOpFactory>(
    new MFHelmholtzOpFactory(m_multifluidIndexSpace,
                             m_dataLocation,
                             m_alpha,
                             m_beta,
                             m_amr->getProbLo(),
                             mflg,
                             mfInterp,
                             mfFluxReg,
                             mfCoarAve,
                             m_amr->getRefinementRatios(),
                             m_amr->getDx(),
                             m_permittivityCell.getData(), // Dummy argument (recall m_alpha = 0.0)
                             m_permittivityFace.getData(),
                             m_permittivityEB.getData(),
                             domainBcFactory,
                             ebbcFactory,
                             jumpBcFactory,
                             ghostPhi,
                             ghostRhs,
                             m_multigridRelaxMethod,
                             bottomDomain,
                             m_multigridJumpOrder,
                             m_multigridJumpWeight,
                             m_amr->getMaxBoxSize()));
}

void
FieldSolverMultigrid::setupMultigrid()
{
  CH_TIME("FieldSolverMultigrid::setupMultigrid()");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::setupMultigrid()" << endl;
  }

  // This routine sets up a standard Chombo multigrid solver using the AMRMultiGrid template. Fortunately, MFHelmholtzOp implements
  // the required functions and we have already set up the factory.

  // Select the bottom solver -- the user will have specified this when parseMultigridSettings() was called.
  LinearSolver<LevelData<MFCellFAB>>* bottomSolver = nullptr;
  switch (m_bottomSolverType) {
  case BottomSolverType::Simple: {
    bottomSolver = &m_mfsolver;

    break;
  }
  case BottomSolverType::BiCGStab: {
    bottomSolver = &m_bicgstab;

    break;
  }
  case BottomSolverType::GMRES: {
    bottomSolver = &m_gmres;

    break;
  }
  default: {
    MayDay::Error("FieldSolverMultigrid::setupMultigrid - logic bust in bottom solver");

    break;
  }
  }

  // Select the multigrid type. Again, the user will have specified this when parseMultigridSettings was called.
  int gmgType;
  switch (m_multigridType) {
  case MultigridType::VCycle: {
    gmgType = 1;

    break;
  }
  case MultigridType::WCycle: {
    gmgType = 2;

    break;
  }
  default: {
    MayDay::Error("FieldSolverMultigrid::setupMultigrid - logic bust in multigrid type selection");

    break;
  }
  }

  // AMRMultiGrid requires the finest level and the coarsest domain.
  const int           finestLevel    = m_amr->getFinestLevel();
  const ProblemDomain coarsestDomain = m_amr->getDomains()[0];

  // Define the Chombo multigrid solver
  m_multigridSolver = RefCountedPtr<AMRMultiGrid<LevelData<MFCellFAB>>>(new AMRMultiGrid<LevelData<MFCellFAB>>);
  m_multigridSolver->define(coarsestDomain, *m_helmholtzOpFactory, bottomSolver, 1 + finestLevel);
  m_multigridSolver->setSolverParameters(m_multigridPreSmooth,
                                         m_multigridPostSmooth,
                                         m_multigridBottomSmooth,
                                         gmgType,
                                         m_multigridMaxIterations,
                                         m_multigridExitTolerance,
                                         m_multigridExitHang,
                                         1.E-99);

  // Minimum number of iterations.
  m_multigridSolver->m_imin = m_multigridMinIterations;

  // Multigrid verbosity, in case user wants to see convergence rates etc.
  m_multigridSolver->m_verbosity = m_multigridVerbosity;

  // Create some dummy storage for multigrid initialization. This is needed because
  // AMRMultiGrid must allocate the operators.
  MFAMRCellData dummy1;
  MFAMRCellData dummy2;

  m_amr->allocate(dummy1, m_realm, m_nComp);
  m_amr->allocate(dummy2, m_realm, m_nComp);

  DataOps::setValue(dummy1, 0.0);
  DataOps::setValue(dummy2, 0.0);

  // Aliasing because AMRMultigrid@Chombo is not too clever when it comes to smart pointers.
  Vector<LevelData<MFCellFAB>*> phi;
  Vector<LevelData<MFCellFAB>*> rhs;

  m_amr->alias(phi, dummy1);
  m_amr->alias(rhs, dummy2);

  // Init the solver. This instantiates the all the operators in AMRMultiGrid so we can just call "solve"
  m_multigridSolver->init(phi, rhs, finestLevel, 0);
}

Vector<long long>
FieldSolverMultigrid::computeLoads(const DisjointBoxLayout& a_dbl, const int a_level)
{
  CH_TIME("FieldSolverMultigrid::computeLoads(DisjointBoxLayout, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::computeLoads(DisjointBoxLayout, int)" << endl;
  }

  // This routine estimates the computational loads on a grid level. It does that by creating an operator (which we later delete)
  // on the grid level corresponding to a_level. It then uses a TimedDataIterator in the input DisjointBoxLayout to estimate the
  // compute time spent when visiting each patch during a typical smoothing step. This step is defined in MFHelmholtzOp::computeOperatorLoads
  // and consists of coarse-fine interpolation, updating BCs, and calling the smoothing kernel (e.g. red-black Gauss-seidel). Fortunately
  // this will provide this class with the opportunity to use introspective load balancing.

  constexpr int numApply = 50;

  // Set up the multigrid solver if that has not already been done.
  if (!m_isSolverSetup) {
    this->setupSolver();
  }

  // Allocate a dummy storage that we can operate on.
  MFAMRCellData dummy;
  m_amr->allocate(dummy, m_realm, m_nComp);

  // Create the operator and run a couple of kernels that involve coarse-fine interpolation, BC updates, and smoothing kernels.
  auto oper = (MFHelmholtzOp*)m_helmholtzOpFactory->MGnewOp(m_amr->getDomains()[a_level], 0, false);

  const Vector<long long> loads = oper->computeOperatorLoads(*dummy[a_level], numApply);

  delete oper;

  return loads;
}

void
FieldSolverMultigrid::computeElectricField(MFAMRCellData& a_electricField, const MFAMRCellData& a_potential) const
{
  CH_TIME("FieldSolverMultigrid::computeElectricField(MFAMRCellData, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::computeElectricField(MFAMRCellData, MFAMRCellData)" << endl;
  }

  CH_assert(a_electricField[0]->nComp() == SpaceDim);
  CH_assert(a_potential[0]->nComp() == 1);

  // Update ghost cells. Use scratch storage for this. Also, we use the multigrid interpolator to do this
  // because that is what is consistent with the Helmholtz discretization. Hence the call to interpGhostMG rather
  // than interpGhostPwl or interpGhost here.
  MFAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_nComp);
  m_amr->copyData(scratch, a_potential);
  m_amr->interpGhostMG(scratch, m_realm);

  // Compute the cell-centered gradient everywhere.
  m_amr->computeGradient(a_electricField, scratch, m_realm);
  DataOps::scale(a_electricField, -1.0);

  // Coarsen solution and update ghost cells.
  m_amr->conservativeAverage(a_electricField, m_realm);
  m_amr->interpGhost(a_electricField, m_realm);
}

void
FieldSolverMultigrid::computeElectricField(MFAMRFluxData& a_electricField, const MFAMRCellData& a_potential) const
{
  CH_TIME("FieldSolverMultigrid::computeElectricField(MFAMRFluxData, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::computeElectricField(MFAMRFluxData, MFAMRCellData)" << endl;
  }

  CH_assert(a_electricField[0]->nComp() == SpaceDim);
  CH_assert(a_potential[0]->nComp() == 1);

  // Update ghost cells. Use scratch storage for this. Also, we use the multigrid interpolator to do this
  // because that is what is consistent with the Helmholtz discretization. Hence the call to interpGhostMG rather
  // than interpGhostPwl or interpGhost here.
  MFAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_nComp);
  m_amr->copyData(scratch, a_potential);
  m_amr->interpGhostMG(scratch, m_realm);

  // Compute the cell-centered gradient everywhere.
  m_amr->computeGradient(a_electricField, scratch, m_realm);
  DataOps::scale(a_electricField, -1.0);
}

void
FieldSolverMultigrid::computeElectricField(EBAMRCellData&           a_electricField,
                                           const phase::which_phase a_phase,
                                           const MFAMRCellData&     a_potential) const
{
  CH_TIME("FieldSolverMultigrid::computeElectricField(EBAMRCellData, phase, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::computeElectricField(EBAMRCellData, phase, MFAMRCellData)" << endl;
  }

  CH_assert(a_electricField[0]->nComp() == SpaceDim);
  CH_assert(a_potential[0]->nComp() == 1);

  // Update ghost cells. Use scratch storage for this. Also, we use the multigrid interpolator to do this
  // because that is what is consistent with the Helmholtz discretization. Hence the call to interpGhostMG rather
  // than interpGhostPwl or interpGhost here.
  EBAMRCellData scratch;
  EBAMRCellData potentialPhase;

  m_amr->allocatePointer(potentialPhase, m_realm);
  m_amr->allocate(scratch, m_realm, a_phase, m_nComp);
  m_amr->alias(potentialPhase, a_phase, a_potential);

  m_amr->copyData(scratch, potentialPhase);
  m_amr->interpGhostMG(scratch, m_realm, a_phase);

  // Use EBGradient for computing the gradient.
  m_amr->computeGradient(a_electricField, scratch, m_realm, a_phase);
  DataOps::scale(a_electricField, -1.0);

  // Coarsen solution and update ghost cells.
  m_amr->conservativeAverage(a_electricField, m_realm, a_phase);
  m_amr->interpGhost(a_electricField, m_realm, a_phase);
}

void
FieldSolverMultigrid::computeElectricField(EBAMRFluxData&           a_electricField,
                                           const phase::which_phase a_phase,
                                           const MFAMRCellData&     a_potential) const
{
  CH_TIME("FieldSolverMultigrid::computeElectricField(EBAMRFluxData, phase, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverMultigrid::computeElectricField(EBAMRFluxData, phase, MFAMRCellData)" << endl;
  }

  CH_assert(a_electricField[0]->nComp() == SpaceDim);
  CH_assert(a_potential[0]->nComp() == 1);

  // Update ghost cells. Use scratch storage for this. Also, we use the multigrid interpolator to do this
  // because that is what is consistent with the Helmholtz discretization. Hence the call to interpGhostMG rather
  // than interpGhostPwl or interpGhost here.
  EBAMRCellData scratch;
  EBAMRCellData potentialPhase;

  m_amr->allocatePointer(potentialPhase, m_realm);
  m_amr->allocate(scratch, m_realm, a_phase, m_nComp);
  m_amr->alias(potentialPhase, a_phase, a_potential);

  m_amr->copyData(scratch, potentialPhase);
  m_amr->interpGhostMG(scratch, m_realm, a_phase);

  // Use EBGradient for computing the gradient.
  m_amr->computeGradient(a_electricField, scratch, m_realm, a_phase);
  DataOps::scale(a_electricField, -1.0);
}

#include <CD_NamespaceFooter.H>
