/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EddingtonSP1.cpp
  @brief  Implementation of CD_EddingtonSP1.H
  @author Robert Marskar
  @todo   Stencil weights and order need to be input parameters!
*/

// Std includes
#include <chrono>
#include <time.h>

// Chombo includes
#include <ParmParse.H>
#include <EBAMRIO.H>

// Our includes
#include <CD_EddingtonSP1.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_DischargeIO.H>
#include <CD_EBHelmholtzDirichletEBBCFactory.H>
#include <CD_EBHelmholtzNeumannEBBCFactory.H>
#include <CD_EBHelmholtzLarsenEBBCFactory.H>
#include <CD_EBHelmholtzEddingtonSP1DomainBCFactory.H>
#include <CD_NamespaceHeader.H>

constexpr Real EddingtonSP1::m_alpha;
constexpr Real EddingtonSP1::m_beta;

Real
EddingtonSP1::s_defaultDomainBcFunction(const RealVect a_position, const Real a_time)
{
  return 1.0;
}

EddingtonSP1::EddingtonSP1() : RtSolver()
{

  // Default settings
  m_name      = "EddingtonSP1";
  m_className = "EddingtonSP1";

  m_verbosity     = -1;
  m_isSolverSetup = false;
  m_dataLocation  = Location::Cell::Center;
  m_regridSlopes  = true;

  // This fills m_domainBcFunctions with s_defaultDomainBcFunction on every domain side.
  this->setDefaultDomainBcFunctions();
}

EddingtonSP1::~EddingtonSP1()
{}

void
EddingtonSP1::parseOptions()
{
  CH_TIME("EddingtonSP1::parseOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions" << endl;
  }

  this->parseVerbosity();         // Parse verbosity
  this->parseDomainBC();          // Parses domain bc
  this->parseEBBC();              // Parse ebbc
  this->parseReflection();        // Parses "reflection coefficients"
  this->parseStationary();        // Parse stationary solver
  this->parsePlotVariables();     // Parses plot variables
  this->parseMultigridSettings(); // Parses solver parameters for geometric multigrid
  this->parseKappaScale();        // Parses kappa-scaling
  this->parseRegridSlopes();      // Slopes on/off when regridding
}

void
EddingtonSP1::parseRuntimeOptions()
{
  CH_TIME("EddingtonSP1::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  this->parseVerbosity();         // Parse verbosity
  this->parsePlotVariables();     // Parses plot variables
  this->parseMultigridSettings(); // Parses solver parameters for geometric multigrid
  this->parseKappaScale();        // Parses kappa-scaling
  this->parseRegridSlopes();      // Slopes on/off when regridding
}

void
EddingtonSP1::setDefaultDomainBcFunctions()
{
  CH_TIME("EddingtonSP1::setDefaultDomainBcFunctions");
  if (m_verbosity > 5) {
    pout() << m_name + "::setDefaultDomainBcFunctions" << endl;
  }

  m_domainBcFunctions.clear();
  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      m_domainBcFunctions.emplace(std::make_pair(dir, sit()), EddingtonSP1::s_defaultDomainBcFunction);
    }
  }
}

void
EddingtonSP1::setDomainSideBcFunction(const int                               a_dir,
                                      const Side::LoHiSide                    a_side,
                                      const EddingtonSP1DomainBc::BcFunction& a_function)
{
  CH_TIME("EddingtonSP1::setDomainSideBcFunction");
  if (m_verbosity > 4) {
    pout() << m_name + "::setDomainSideBcFunction" << endl;
  }

  const EddingtonSP1DomainBc::DomainSide domainSide = std::make_pair(a_dir, a_side);

  m_domainBcFunctions.at(domainSide) = a_function;
}

std::string
EddingtonSP1::makeBcString(const int a_dir, const Side::LoHiSide a_side) const
{
  CH_TIME("EddingtonSP1::makeBcString");
  if (m_verbosity > 5) {
    pout() << m_name + "::makeBcString" << endl;
  }

  // TLDR: Returns string in format "bc.x.lo"
  std::string strDir;
  std::string strSide;

  if (a_dir == 0) {
    strDir = "x";
  }
  else if (a_dir == 1) {
    strDir = "y";
  }
  else if (a_dir == 2) {
    strDir = "z";
  }

  if (a_side == Side::Lo) {
    strSide = "lo";
  }
  else if (a_side == Side::Hi) {
    strSide = "hi";
  }

  const std::string ret = std::string("bc.") + strDir + std::string(".") + strSide;

  return ret;
}

EddingtonSP1DomainBc::BcType
EddingtonSP1::parseBcString(const std::string a_str) const
{
  CH_TIME("EddingtonSP1::parseBcString");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseBcString" << endl;
  }

  EddingtonSP1DomainBc::BcType ret;

  if (a_str == "dirichlet") {
    ret = EddingtonSP1DomainBc::BcType::Dirichlet;
  }
  else if (a_str == "neumann") {
    ret = EddingtonSP1DomainBc::BcType::Neumann;
  }
  else if (a_str == "larsen") {
    ret = EddingtonSP1DomainBc::BcType::Larsen;
  }
  else {
    MayDay::Error("EddingtonSP1DomainBc::BcType - unknown BC type!");
  }

  return ret;
}

void
EddingtonSP1::parseDomainBC()
{
  CH_TIME("EddingtonSP1::parseDomainBC");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseDomainBC" << endl;
  }

  // TLDR: This routine might seem big and complicated. What we are doing is that we are creating one function object which returns some value
  //       anywhere in space and time on a domain edge (face). The code below simply creates those functions and associates them with an edge (face in 3D)
  //
  //       Because of the stationary interface of EBHelmholtzOp and the transient interface of the radiative transfer solvers, we capture the time by solver
  //       time (RtSolver::m_time) by reference in these functions. That allows us to pass the functions in as stationary-looking functions in EBHelmholtzOp,
  //       yet retain transient BCs.
  //
  //       For flexibility we want to be able to specify the BC functions in two forms, using a multiplier or not. The specification for this is e.g.
  //       "dirichlet <number>" in which case the bc function is multiplied by <number>. Use "dirichlet_custom" to use the bc functions directly.

  ParmParse pp(m_className.c_str());

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const EddingtonSP1DomainBc::DomainSide domainSide = std::make_pair(dir, sit());
      const std::string                      bcString   = this->makeBcString(dir, sit());
      const int                              num        = pp.countval(bcString.c_str());

      std::string str;

      EddingtonSP1DomainBc::BcType      bcType; // Will be set based on what comes in through the input script
      EddingtonSP1DomainBc::BcFunction& bcFunc = m_domainBcFunctions.at(
        domainSide); // Will be set based on what comes in through the input script

      std::function<Real(const RealVect, const Real)> curFunc;

      switch (num) {
      case 1: { // If only providing dirichlet_custom, neumann_custom, larsen_custom, the second value is overridden and only the function is used as BC.
        pp.get(bcString.c_str(), str, 0);

        // Set the function. Capture solver time by reference.
        curFunc = [&bcFunc, &time = this->m_time](const RealVect a_pos, const Real a_time) {
          return bcFunc(a_pos, time);
        };

        // Ensure that we actually have a valid BC.
        if (str == "dirichlet_custom") {
          bcType = EddingtonSP1DomainBc::BcType::Dirichlet;
        }
        else if (str == "neumann_custom") {
          bcType = EddingtonSP1DomainBc::BcType::Neumann;
        }
        else if (str == "larsen_custom") {
          bcType = EddingtonSP1DomainBc::BcType::Larsen;
        }
        else {
          MayDay::Error(
            "EddingtonSP1::parseDomainBC -- got only one argument but this argument was not in the form of e.g. 'neumann_custom'!");
        }

        break;
      }
      case 2: { // Had two arguments in input, e.g. neumann 0.0. In this case we set the BC to be the specified function times the value.
        Real val;

        pp.get(bcString.c_str(), str, 0);
        pp.get(bcString.c_str(), val, 1);

        bcType = this->parseBcString(str);

        curFunc = [&bcFunc, &time = this->m_time, val](const RealVect a_pos, const Real a_time) {
          return bcFunc(a_pos, time) * val;
        };

        break;
      }
      default:
        MayDay::Error("EddingtonSP1::parseDomainBC -- unknown argument to domain bc");
        break;
      }

      m_domainBc.setBc(domainSide, std::make_pair(bcType, curFunc));
    }
  }
}

void
EddingtonSP1::parseEBBC()
{
  CH_TIME("EddingtonSP1::parseEBBC");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseEBBC" << endl;
  }

  ParmParse pp(m_className.c_str());

  const std::string bcString = "ebbc";
  const int         num      = pp.countval(bcString.c_str());

  if (num == 2) {
    Real        val;
    std::string str;

    pp.get(bcString.c_str(), str, 0);
    pp.get(bcString.c_str(), val, 1);

    EBBCType bcType;

    if (str == "dirichlet") {
      bcType = EBBCType::Dirichlet;
    }
    else if (str == "neumann") {
      bcType = EBBCType::Neumann;
    }
    else if (str == "larsen") {
      bcType = EBBCType::Larsen;
    }
    else {
      MayDay::Error("EddingtonSP1::parseEBBC -- not dirichlet/neumann/larsen!");
    }

    m_ebbc = std::make_pair(bcType, val);
  }
  else {
    MayDay::Error(
      "EddingtonSP1::parseEBBC -- bad input argument. Should be in the form 'larsen 0.0', 'dirichlet 0.0', etc");
  }
}

void
EddingtonSP1::parseStationary()
{
  CH_TIME("EddingtonSP1::parseStationary");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseStationary" << endl;
  }

  ParmParse   pp(m_className.c_str());
  std::string str;

  pp.get("stationary", m_stationary);
}

void
EddingtonSP1::parsePlotVariables()
{
  CH_TIME("EddingtonSP1::parsePlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePlotVariables" << endl;
  }

  m_plotPhi    = false;
  m_plotSource = false;

  ParmParse           pp(m_className.c_str());
  const int           num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++) {
    if (str[i] == "phi")
      m_plotPhi = true;
    else if (str[i] == "src")
      m_plotSource = true;
  }
}

void
EddingtonSP1::parseKappaScale()
{
  ParmParse pp(m_className.c_str());

  pp.get("kappa_scale", m_kappaScale);
}

void
EddingtonSP1::parseMultigridSettings()
{
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
  pp.get("gmg_ebbc_order", m_multigridBcOrder);
  pp.get("gmg_ebbc_weight", m_multigridBcWeight);

  // Fetch the desired bottom solver from the input script. We look for things like EddingtonSP1.gmg_bottom_solver = bicgstab or '= simple <number>'
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
        "EddingtonSP1::parseMultigridSettings - logic bust, you've specified one parameter and I expected either 'bicgstab' or 'gmres'");
    }
  }
  else if (num == 2) {
    int numSmooth;
    pp.get("gmg_bottom_solver", str, 0);
    pp.get("gmg_bottom_solver", numSmooth, 1);
    if (str == "simple") {
      m_bottomSolverType = BottomSolverType::Simple;
      m_simpleSolver.setNumSmooths(numSmooth);
    }
    else {
      MayDay::Error(
        "EddingtonSP1::parseMultigridSettings - logic bust, you've specified two parameters and I expected 'simple <number>'");
    }
  }
  else {
    MayDay::Error(
      "EddingtonSP1::parseMultigridSettings - logic bust in bottom solver. You must specify ' = bicgstab', ' = gmres', or ' = simple <number>'");
  }

  // Relaxation type
  pp.get("gmg_smoother", str);
  if (str == "jacobi") {
    m_multigridRelaxMethod = EBHelmholtzOp::Smoother::PointJacobi;
  }
  else if (str == "red_black") {
    m_multigridRelaxMethod = EBHelmholtzOp::Smoother::GauSaiRedBlack;
  }
  else if (str == "multi_color") {
    m_multigridRelaxMethod = EBHelmholtzOp::Smoother::GauSaiMultiColor;
  }
  else {
    MayDay::Error("EddingtonSP1::parseMultigridSettings - unknown relaxation method requested");
  }

  // Cycle type
  pp.get("gmg_cycle", str);
  if (str == "vcycle") {
    m_multigridType = MultigridType::VCycle;
  }
  else {
    MayDay::Error("EddingtonSP1::parseMultigridSettings - unknown cycle type requested");
  }

  // No lower than 2.
  if (m_minCellsBottom < 2) {
    m_minCellsBottom = 2;
  }
}

void
EddingtonSP1::parseReflection()
{
  CH_TIME("EddingtonSP1::parse_reflectivity");
  if (m_verbosity > 5) {
    pout() << m_name + "::parse_reflectivity" << endl;
  }

  ParmParse pp(m_className.c_str());

  Real r;
  pp.get("reflectivity", r);

  m_reflectCoefOne = r / (2.0);
  m_reflectCoefTwo = r / (3.0);
}

void
EddingtonSP1::parseRegridSlopes()
{
  CH_TIME("EddingtonSP1::parseRegridSlopes()");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRegridSlopes()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("use_regrid_slopes", m_regridSlopes);
}

void
EddingtonSP1::preRegrid(const int a_base, const int a_oldFinestLevel)
{
  CH_TIME("EddingtonSP1::preRegrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::preRegrid" << endl;
  }

  m_amr->allocate(m_cachePhi, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_cacheSrc, m_realm, m_phase, m_nComp);

  m_amr->copyData(m_cachePhi, m_phi);
  m_amr->copyData(m_cacheSrc, m_source);

  this->deallocate();
}

void
EddingtonSP1::allocate()
{
  CH_TIME("EddingtonSP1::allocate");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate" << endl;
  }

  m_amr->allocate(m_helmAco, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_helmBco, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_helmBcoIrreg, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_phi, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_source, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_resid, m_realm, m_phase, m_nComp);

  DataOps::setValue(m_resid, 0.0);
  DataOps::setValue(m_phi, 0.0);
  DataOps::setValue(m_source, 0.0);

  this->setHelmholtzCoefficients();
}

void
EddingtonSP1::deallocate()
{
  CH_TIME("EddingtonSP1::deallocate");
  if (m_verbosity > 5) {
    pout() << m_name + "::deallocate" << endl;
  }

  m_helmAco.clear();
  m_helmBco.clear();
  m_helmBcoIrreg.clear();
  m_phi.clear();
  m_source.clear();
  m_resid.clear();
}

void
EddingtonSP1::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("EddingtonSP1::regrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid" << endl;
  }

  const Interval interv(m_comp, m_comp);

  // Allocate storage.
  this->allocate();

  const EBCoarseToFineInterp::Type interpType = m_regridSlopes ? EBCoarseToFineInterp::Type::ConservativeMinMod
                                                               : EBCoarseToFineInterp::Type::ConservativePWC;

  // Regrid phi and source
  m_amr->interpToNewGrids(m_phi, m_cachePhi, m_phase, a_lmin, a_oldFinestLevel, a_newFinestLevel, interpType);
  m_amr->interpToNewGrids(m_source, m_cacheSrc, m_phase, a_lmin, a_oldFinestLevel, a_newFinestLevel, interpType);

  // Coarsen and update ghost cells.
  m_amr->conservativeAverage(m_phi, m_realm, m_phase);
  m_amr->interpGhost(m_phi, m_realm, m_phase);

  m_isSolverSetup = false;

  // Deallocate the scratch data.
  m_cachePhi.clear();
  m_cacheSrc.clear();
}

void
EddingtonSP1::registerOperators()
{
  CH_TIME("EddingtonSP1::registerOperators");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerOperators" << endl;
  }

  if (m_amr.isNull()) {
    MayDay::Error("EddingtonSP1::registerOperators - need to set AmrMesh!");
  }
  else {
    m_amr->registerOperator(s_eb_coar_ave, m_realm, m_phase);
    m_amr->registerOperator(s_eb_fill_patch, m_realm, m_phase);
    m_amr->registerOperator(s_eb_flux_reg, m_realm, m_phase);
    m_amr->registerOperator(s_eb_gradient, m_realm, m_phase);
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, m_phase);
    m_amr->registerOperator(s_eb_fine_interp, m_realm, m_phase);
    m_amr->registerOperator(s_eb_multigrid, m_realm, m_phase);
  }
}

bool
EddingtonSP1::advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zeroPhi)
{
  CH_TIME("EddingtonSP1::advance(ebamrcell, ebamrcell)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advance(ebamrcell, ebamrcell)" << endl;
  }

  // TLDR: We are solving the equation
  //
  // (1/c) d(Phi)/dt + kappa*Phi - div(1/(3*kappa) grad(Phi)) = eta/c
  //
  // Our EBHelmholtzOp takes the A-coefficient as and the b-coefficient as 1/(3*kappa). Other than that,
  // the only weird thing about our discretization is that we solve for the "time" c*dt and scale the isotropic
  // term by c

  bool converged = false;

  if (!m_isSolverSetup) {
    this->setupSolver();
  }

  EBAMRCellData zero;
  EBAMRCellData scaledSource;

  m_amr->allocate(zero, m_realm, m_phase, 1);
  m_amr->allocate(scaledSource, m_realm, m_phase, 1);

  // Modify the source term.  Operator is scaled by kappa and source term might also have to be scaled by kappa
  // Source should be scaled by 1./c0 (due to the way we do the EBHelmholtzOp coefficients)
  DataOps::copy(scaledSource, a_source); // Copy source term
  DataOps::scale(scaledSource, 1. / Units::c);

  if (m_stationary) {

    // If we're doing a stationary solve, we must scale the source term by kappa (unless it's otherwise been done).
    if (m_kappaScale) {
      DataOps::kappaScale(scaledSource);
    }

    // Aliasing, because Chombo is not too smart.
    Vector<LevelData<EBCellFAB>*> phi;
    Vector<LevelData<EBCellFAB>*> rhs;
    Vector<LevelData<EBCellFAB>*> res;
    Vector<LevelData<EBCellFAB>*> zer;

    m_amr->alias(phi, a_phi);
    m_amr->alias(rhs, scaledSource);
    m_amr->alias(res, m_resid);
    m_amr->alias(zer, zero);

    const int coarsestLevel = 0;
    const int finestLevel   = m_amr->getFinestLevel();

    // Compute the residual and determine if we must enter multigrid.
    const Real phiResid  = m_multigridSolver->computeAMRResidual(phi, rhs, finestLevel, coarsestLevel);
    const Real zeroResid = m_multigridSolver->computeAMRResidual(zer, rhs, finestLevel, coarsestLevel);

    if (phiResid > zeroResid * m_multigridExitTolerance) {
      // Residual is too large, solve.
      m_multigridSolver->m_convergenceMetric = zeroResid;
      m_multigridSolver->solveNoInitResid(phi, res, rhs, finestLevel, coarsestLevel, a_zeroPhi);

      const int status = m_multigridSolver->m_exitStatus; // 1 => Initial norm sufficiently reduced
      if (status == 1 || status == 8 || status == 9) {    // 8 => Norm sufficiently small
        converged = true;
      }
    }
    else {
      // Solution is already good enough
      converged = true;
    }
    m_multigridSolver->revert(phi, rhs, finestLevel, 0);

    DataOps::setCoveredValue(a_phi, 0, 0.0);
  }
  else {
    this->advanceEuler(a_phi, scaledSource, Units::c * a_dt, a_zeroPhi);

    const int status = m_multigridSolver->m_exitStatus; // 1 => Initial norm sufficiently reduced
    if (status == 1 || status == 8 || status == 9) {    // 8 => Norm sufficiently small
      converged = true;
    }
  }

  DataOps::setCoveredValue(a_phi, 0.0);

  m_amr->conservativeAverage(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);

  return converged;
}

void
EddingtonSP1::advanceEuler(EBAMRCellData&       a_phi,
                           const EBAMRCellData& a_source,
                           const Real           a_dt,
                           const bool           a_zeroPhi) noexcept
{
  CH_TIME("EddingtonSP1::advanceEuler");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceEuler" << endl;
  }

  // Set up the multigrid solver
  if (!m_isSolverSetup) {
    this->setupSolver();
  }

  EBAMRCellData zero;
  EBAMRCellData scratch;

  m_amr->allocate(zero, m_realm, m_phase, 1);
  m_amr->allocate(scratch, m_realm, m_phase, 1);

  // Make the right-hand side for the Euler equation advance. First compute kappa*acoef*(phi^k) + dt*kappa*S, which
  // is the right-hand side of our system.
  DataOps::setValue(zero, 0.0);
  DataOps::copy(scratch, a_phi);

  // Compute the diagonal scaling of the operator.
  Vector<AMRLevelOp<LevelData<EBCellFAB>>*> amrOps = m_multigridSolver->getAMROperators();
  Vector<MGLevelOp<LevelData<EBCellFAB>>*>  mgOps  = m_multigridSolver->getAllOperators();
  for (int i = 0; i < amrOps.size(); i++) {
    TGAHelmOp<LevelData<EBCellFAB>>* helmholtzOperator = (TGAHelmOp<LevelData<EBCellFAB>>*)amrOps[i];

    helmholtzOperator->diagonalScale(*scratch[i], false);
  }

  // Set coefficients for multigrid solve.
  for (int i = 0; i < mgOps.size(); i++) {
    TGAHelmOp<LevelData<EBCellFAB>>* helmholtzOperator = (TGAHelmOp<LevelData<EBCellFAB>>*)mgOps[i];

    helmholtzOperator->setAlphaAndBeta(1.0, -a_dt);
  }

  DataOps::incr(scratch, a_source, a_dt);
  if (m_kappaScale) {
    DataOps::kappaScale(scratch);
  }

  Vector<LevelData<EBCellFAB>*> newPhi;
  Vector<LevelData<EBCellFAB>*> eulerRHS;
  Vector<LevelData<EBCellFAB>*> zer;

  m_amr->alias(newPhi, a_phi);
  m_amr->alias(eulerRHS, scratch);
  m_amr->alias(zer, zero);

  const int coarsestLevel = 0;
  const int finestLevel   = m_amr->getFinestLevel();

  // Figure out how far away we are from a "converged" solution and set the convergence metric. Then solve.
  const Real zeroResid = m_multigridSolver->computeAMRResidual(zer, eulerRHS, finestLevel, coarsestLevel);

  m_multigridSolver->m_convergenceMetric = zeroResid;
  m_multigridSolver->solve(newPhi, eulerRHS, finestLevel, coarsestLevel, a_zeroPhi);
}

void
EddingtonSP1::setupSolver()
{
  CH_TIME("EddingtonSP1::setupSolver");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupSolver" << endl;
  }

  this->setHelmholtzCoefficients(); // Set coefficients, kappa, aco, bco
  this->setupHelmholtzFactory();    // Set up the Helmholtz operator factory
  this->setupMultigrid();           // Set up the AMR multigrid solver

  m_isSolverSetup = true;
}

void
EddingtonSP1::setHelmholtzCoefficients()
{
  CH_TIME("EddingtonSP1::setHelmholtzCoefficients");
  if (m_verbosity > 5) {
    pout() << m_name + "::setHelmholtzCoefficients" << endl;
  }

  // This loop fills aco with kappa and bco_irreg with 1./kappa
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    LevelData<EBCellFAB>&       helmAco      = *m_helmAco[lvl];
    LevelData<EBFluxFAB>&       helmBco      = *m_helmBco[lvl];
    LevelData<BaseIVFAB<Real>>& helmBcoIrreg = *m_helmBcoIrreg[lvl];

    // Fill data in each grid patch.
    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      this->setHelmholtzCoefficientsBox(helmAco[din], helmBco[din], helmBcoIrreg[din], lvl, din);
    }
  }
}

void
EddingtonSP1::setHelmholtzCoefficientsBox(EBCellFAB&       a_helmAco,
                                          EBFluxFAB&       a_helmBco,
                                          BaseIVFAB<Real>& a_helmBcoIrreg,
                                          const int        a_lvl,
                                          const DataIndex& a_dit)
{
  CH_TIME("EddingtonSP1::setHelmholtzCoefficientsBox");
  if (m_verbosity > 10) {
    pout() << m_name + "::setHelmholtzCoefficientsBox" << endl;
  }

  CH_assert(a_helmAco.nComp() == 1);
  CH_assert(a_helmBco.nComp() == 1);
  CH_assert(a_helmBcoIrreg.nComp() == 1);
  CH_assert(a_helmAco.box() == a_helmBco.box());

  // Setting this to trigger debugging.
#ifndef NDEBUG
  a_helmAco.setVal(std::numeric_limits<Real>::max());
  a_helmBco.setVal(std::numeric_limits<Real>::max());
  a_helmBcoIrreg.setVal(std::numeric_limits<Real>::max());
#endif

  // TLDR: This routine is for setting coefficients in the Helmholtz operator. These coefficients are set as A = kappa, B = 1/(3*kappa). We happen to know that
  //       the interior face stencils are interpolated using the neigboring face, so we must fill the "ghost faces" around the B-coefficient grid patch.

  const RealVect probLo = m_amr->getProbLo();
  const Real     dx     = m_amr->getDx()[a_lvl];

  const EBISBox& ebisbox = a_helmAco.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();

  const Box cellBox    = m_amr->getGrids(m_realm)[a_lvl][a_dit];
  const Box helmAcoBox = a_helmAco.box() & m_amr->getDomains()[a_lvl];
  const Box helmBcoBox = a_helmBco.box() & m_amr->getDomains()[a_lvl];

  // Regular A-coefficient kernel. Recall that the A-coefficient only affects the diagonal part of the stencil so there's no need to fill anything
  // outside of the cell-centered grid patch.
  BaseFab<Real>& helmAcoReg       = a_helmAco.getSingleValuedFAB();
  auto           regularAcoKernel = [&](const IntVect& iv) -> void {
    const RealVect pos   = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
    const Real     kappa = m_rtSpecies->getAbsorptionCoefficient(pos);

    helmAcoReg(iv, m_comp) = kappa;
  };

  BoxLoops::loop(cellBox, regularAcoKernel); // Fill single-valued a-coefficient.

  // Regular B-coefficient. Recall that EBHelmholtzOp sets up the face centroid fluxes by interpolating with neighboring face-centered fluxes. The interpolating stencil
  // will have a radius of 1, so we need to fill one of the ghost faces outside of the grid patch. Only the ones that are "tangential" to the face direction
  // are necessary to fill.
  for (int dir = 0; dir < SpaceDim; dir++) {
    EBFaceFAB& helmBcoFace = a_helmBco[dir];

    // Kernel region -- make a cell-centered box, grown by one in every direction except 'dir'. This box will also contain
    // the "ghost faces".
    Box grownCellBox = cellBox;
    for (int otherDir = 0; otherDir < SpaceDim; otherDir++) {
      if (otherDir != dir) {
        grownCellBox.grow(otherDir, 1);
      }
    }
    grownCellBox &= m_amr->getDomains()[a_lvl];

    // Cut-cells in the grown box.
    const IntVectSet irregIVS = ebisbox.getIrregIVS(grownCellBox);

    // Actual kernel region. The face-centered box which also contains the "ghost faces" and a FaceIterator for
    // doing the irregular face stuff.
    const Box    faceBox = surroundingNodes(grownCellBox, dir);
    FaceIterator faceit(irregIVS, ebgraph, dir, FaceStop::SurroundingWithBoundary);

    // Regular kernel.
    BaseFab<Real>& helmBcoReg       = helmBcoFace.getSingleValuedFAB();
    auto           regularBcoKernel = [&](const IntVect& iv) -> void {
      const RealVect pos   = probLo + dx * ((RealVect(iv) + 0.5 * RealVect::Unit) - 0.5 * BASISREALV(dir));
      const Real     kappa = m_rtSpecies->getAbsorptionCoefficient(pos);

      helmBcoReg(iv, m_comp) = 1. / (3.0 * kappa);
    };

    // Irregular Bco-kernel. Does exactly the same.
    auto irregularBcoKernel = [&](const FaceIndex& face) -> void {
      const RealVect pos   = probLo + Location::position(Location::Face::Center, face, ebisbox, dx);
      const Real     kappa = m_rtSpecies->getAbsorptionCoefficient(pos);

      helmBcoFace(face, m_comp) = 1. / (3.0 * kappa);
    };

    // Run the kernels
    BoxLoops::loop(faceBox, regularBcoKernel);
    BoxLoops::loop(faceit, irregularBcoKernel);
  }

  // Fill A-coefficient in cut-cells and the EB-centered B-coefficient. Again, the A-part is diagonal so we don't need to fill anything outside the grid patch. For
  // the B-coefficient on the EB face then the stencil in the cut-cell will not reach into neighboring EB faces (really, it shouldn't!) so no need for filling
  // things outside the grid patch here, either.
  // Kernel region for cut-cells
  VoFIterator& vofit           = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
  auto         irregularKernel = [&](const VolIndex& vof) -> void {
    const RealVect pos   = probLo + Location::position(m_dataLocation, vof, ebisbox, dx);
    const Real     kappa = m_rtSpecies->getAbsorptionCoefficient(pos);

    a_helmAco(vof, m_comp)      = kappa;
    a_helmBcoIrreg(vof, m_comp) = 1. / (3.0 * kappa);
  };

  // Run the irregular kernel
  BoxLoops::loop(vofit, irregularKernel);
}

void
EddingtonSP1::setupHelmholtzFactory()
{
  CH_TIME("EddingtonSP1::setupHelmholtzFactory");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupHelmholtzFactory" << endl;
  }

  const Vector<RefCountedPtr<EBLevelGrid>>&             levelGrids   = m_amr->getEBLevelGrid(m_realm, m_phase);
  const Vector<RefCountedPtr<EBCoarAve>>&               coarAve      = m_amr->getCoarseAverage(m_realm, m_phase);
  const Vector<RefCountedPtr<EBReflux>>&                fluxReg      = m_amr->getFluxRegister(m_realm, m_phase);
  const Vector<RefCountedPtr<EBMultigridInterpolator>>& interpolator = m_amr->getMultigridInterpolator(m_realm,
                                                                                                       m_phase);

  // Coarsest domain used for multigrid. The user specifies the minimum number of cells in any
  // coordinate direction, and we coarsen until we have a domain which satisfies that constraint.
  ProblemDomain bottomDomain = m_amr->getDomains()[0];
  while (bottomDomain.domainBox().shortside() >= 2 * m_minCellsBottom) {
    bottomDomain.coarsen(2);
  }

  // Number of ghost cells in data holders
  const IntVect ghostPhi = m_amr->getNumberOfGhostCells() * IntVect::Unit;
  const IntVect ghostRhs = m_amr->getNumberOfGhostCells() * IntVect::Unit;

  // Handle to boundary condition factories for domain and EB in an EBHelmholtzOp context.
  RefCountedPtr<EBHelmholtzDomainBCFactory> domainBcFactory;
  RefCountedPtr<EBHelmholtzEBBCFactory>     ebbcFactory;

  // Domain bcs use a wrapper for fetching the bc function on a specified edge
  domainBcFactory = RefCountedPtr<EBHelmholtzDomainBCFactory>(
    new EBHelmholtzEddingtonSP1DomainBCFactory(m_domainBc, m_rtSpecies, m_reflectCoefOne, m_reflectCoefTwo));

  // For EBs we run with only one type of BC (for now). This code sets up the embedded boundary conditions. In parseEBBC we happened to
  // read a line from the input script in the format = <type> <value> which was put into a function. Associate that function here.
  switch (m_ebbc.first) {
  case EBBCType::Dirichlet:
    ebbcFactory = RefCountedPtr<EBHelmholtzDirichletEBBCFactory>(
      new EBHelmholtzDirichletEBBCFactory(m_multigridBcOrder, m_multigridBcWeight, m_ebbc.second));
    break;
  case EBBCType::Neumann:
    ebbcFactory = RefCountedPtr<EBHelmholtzEBBCFactory>(new EBHelmholtzNeumannEBBCFactory(m_ebbc.second));
    break;
  case EBBCType::Larsen: {
    constexpr int order  = 1; // Need to become input parameters
    constexpr int weight = 0;
    ebbcFactory          = RefCountedPtr<EBHelmholtzEBBCFactory>(
      new EBHelmholtzLarsenEBBCFactory(order, weight, m_rtSpecies, m_reflectCoefOne, m_reflectCoefTwo, m_ebbc.second));
    break;
  }
  default:
    MayDay::Error("EddingtonSP1::setupHelmholtzFactory - logic bust in EB BC factory");
    break;
  }

  // Set up the operator
  const Real relaxFactor = 1.0;

  m_helmholtzOpFactory = RefCountedPtr<EBHelmholtzOpFactory>(new EBHelmholtzOpFactory(m_dataLocation,
                                                                                      m_alpha,
                                                                                      m_beta,
                                                                                      m_amr->getProbLo(),
                                                                                      levelGrids,
                                                                                      m_amr->getValidCells(m_realm),
                                                                                      interpolator,
                                                                                      fluxReg,
                                                                                      coarAve,
                                                                                      m_amr->getRefinementRatios(),
                                                                                      m_amr->getDx(),
                                                                                      m_helmAco.getData(),
                                                                                      m_helmBco.getData(),
                                                                                      m_helmBcoIrreg.getData(),
                                                                                      domainBcFactory,
                                                                                      ebbcFactory,
                                                                                      ghostPhi,
                                                                                      ghostRhs,
                                                                                      m_multigridRelaxMethod,
                                                                                      relaxFactor,
                                                                                      bottomDomain,
                                                                                      m_amr->getMaxBoxSize()));
}

void
EddingtonSP1::setupMultigrid()
{
  CH_TIME("EddingtonSP1::setupMultigrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupMultigrid" << endl;
  }

  // Select the bottom solver
  LinearSolver<LevelData<EBCellFAB>>* botsolver = NULL;
  if (m_bottomSolverType == BottomSolverType::Simple) {
    botsolver = &m_simpleSolver;
  }
  else if (m_bottomSolverType == BottomSolverType::BiCGStab) {
    botsolver = &m_bicgstab;
  }
  else if (m_bottomSolverType == BottomSolverType::GMRES) {
    botsolver           = &m_gmres;
    m_gmres.m_verbosity = 0; // Shut up.
  }

  // Make m_multigridType into an int for multigrid
  int gmgType;
  switch (m_multigridType) {
  case MultigridType::VCycle:
    gmgType = 1;
    break;
  case MultigridType::WCycle:
    gmgType = 2;
    break;
  default:
    MayDay::Error("EddingtonSP1::setupMultigrid -- logic bust in multigrid type selection");
  }

  const int           finestLevel    = m_amr->getFinestLevel();
  const ProblemDomain coarsestDomain = m_amr->getDomains()[0];

  // Define AMRMultiGrid
  m_multigridSolver = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB>>>(new AMRMultiGrid<LevelData<EBCellFAB>>());

  m_multigridSolver->define(coarsestDomain, *m_helmholtzOpFactory, botsolver, 1 + finestLevel);
  m_multigridSolver->setSolverParameters(m_multigridPreSmooth,
                                         m_multigridPostSmooth,
                                         m_multigridBottomSmooth,
                                         gmgType,
                                         m_multigridMaxIterations,
                                         m_multigridExitTolerance,
                                         m_multigridExitHang,
                                         1.E-99); // Residue set through other means

  m_multigridSolver->m_imin      = m_multigridMinIterations;
  m_multigridSolver->m_verbosity = m_multigridVerbosity;

  // Dummies for init
  EBAMRCellData dummy1;
  EBAMRCellData dummy2;

  m_amr->allocate(dummy1, m_realm, m_phase, m_nComp);
  m_amr->allocate(dummy2, m_realm, m_phase, m_nComp);

  DataOps::setValue(dummy1, 0.0);
  DataOps::setValue(dummy2, 0.0);

  // Aliasing
  Vector<LevelData<EBCellFAB>*> phi;
  Vector<LevelData<EBCellFAB>*> rhs;

  m_amr->alias(phi, dummy1);
  m_amr->alias(rhs, dummy2);

  // Init solver. This instantiates all the operators in AMRMultiGrid so we can just call "solve"
  m_multigridSolver->init(phi, rhs, finestLevel, 0);
}

void
EddingtonSP1::computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi)
{
  CH_TIME("EddingtonSP1::computeBoundaryFlux");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeBoundaryFlux" << endl;
  }

  // TLDR: Equations say boundary flux => F_n = c*Psi/2, so we
  //       extrapolate the solution to the boundary and multiply by c/2.

  const int finestLevel = m_amr->getFinestLevel();

  m_amr->interpToEB(a_ebFlux, a_phi, m_realm, m_phase);
  m_amr->conservativeAverage(a_ebFlux, m_realm, m_phase);

  DataOps::scale(a_ebFlux, 0.5 * Units::c);
}

void
EddingtonSP1::computeDomainFlux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_data)
{
  CH_TIME("EddingtonSP1::computeDomainFlux");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDomainFlux" << endl;
  }

  // Equations say boundary flux => F_n = c*Psi/2.
  // This code simply extrapolates to the boundary, using linear interpolation if it can.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const EBCellFAB&     data     = (*a_data[lvl])[din];
      const EBISBox&       ebisbox  = ebisl[din];
      const BaseFab<Real>& data_fab = data.getSingleValuedFAB();

      for (int dir = 0; dir < SpaceDim; dir++) {
        for (SideIterator sit; sit.ok(); ++sit) {
          BaseIFFAB<Real>& extrap = (*a_domainflux[lvl])[din](dir, sit());

          const IntVectSet& ivs     = extrap.getIVS();
          const EBGraph&    ebgraph = extrap.getEBGraph();

          // Extrapolate to the boundary. Use face-centered stuff for all faces (also multivalued ones)
          const FaceStop::WhichFaces crit = FaceStop::AllBoundaryOnly;
          FaceIterator               faceit(ivs, ebgraph, dir, crit);

          auto kernel = [&](const FaceIndex& face) -> void {
            const int sgn = sign(sit()); // Lo = -1, Hi = 1

            const VolIndex& vof = face.getVoF(flip(sit()));
            const IntVect   iv0 = vof.gridIndex();
            const IntVect   iv1 = iv0 - sgn * BASISV(dir);

            if (ebisbox.isCovered(iv0)) { // Just provide some bogus data because the face
              extrap(face, m_comp) = 0.0;
            }
            else {
              if (!ebisbox.isCovered(iv1)) { // linear extrapolation
                extrap(face, m_comp) = 1.5 * data_fab(iv0, m_comp) - 0.5 * data_fab(iv1, m_comp); // Should be ok
              }
              else { // Not enough cells available, use cell-centered only
                extrap(face, m_comp) = data_fab(iv0, m_comp);
              }
            }

            // Necessary scaling
            extrap(face, m_comp) = 0.5 * Units::c * extrap(face, m_comp);
          };

          BoxLoops::loop(faceit, kernel);
        }
      }
    }
  }
}

void
EddingtonSP1::computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi)
{
  CH_TIME("EddingtonSP1::computeFlux");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeFlux" << endl;
  }

  // Equations say flux = -c/(3*kappa)*grad(Psi).
  //
  // We happen that our discretization uses Aco = kappa, so we just divide by it rather than computing
  // it on the mesh again.

  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, m_nComp);

  DataOps::copy(scratch, a_phi);
  m_amr->interpGhostMG(scratch, m_realm, m_phase);
  m_amr->computeGradient(a_flux, a_phi, m_realm, m_phase); // flux = grad(phi)

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    DataOps::divideByScalar(*a_flux[lvl], *m_helmAco[lvl]);   // flux = grad(phi)/(c*kappa)
    DataOps::scale(*a_flux[lvl], -Units::c * Units::c / 3.0); // flux = -c*grad(phi)/3.
  }

  m_amr->conservativeAverage(a_flux, m_realm, m_phase);
  m_amr->interpGhost(a_flux, m_realm, m_phase);
}

void
EddingtonSP1::computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi)
{
  CH_TIME("EddingtonSP1::computeDensity");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDensity" << endl;
  }

  // For EddingtonSP1, the solution vector is only the isotropic part of the RTE, so just copy over.
  const Interval interv(m_comp, m_comp);
  m_amr->copyData(a_isotropic, a_phi, interv, interv, CopyStrategy::ValidGhost, CopyStrategy::ValidGhost);
}

void
EddingtonSP1::writePlotFile()
{
  CH_TIME("EddingtonSP1::writePlotFile");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotFile" << endl;
  }

  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_timeStep, SpaceDim);

  const int      ncomps = 3 + SpaceDim;
  Vector<string> names(ncomps);
  names[0] = "density";
  names[1] = "x-flux";
  names[2] = "y-flux";
  if (SpaceDim == 3) {
    names[3] = "z-flux";
    names[4] = "isotropic source";
    names[5] = "residue";
  }
  else {
    names[3] = "isotropic source";
    names[4] = "residue";
  }

  // Compute the flux
  EBAMRCellData flux;
  m_amr->allocate(flux, m_realm, m_phase, SpaceDim);
  this->computeFlux(flux, m_phi);

  // Allocate output storage
  EBAMRCellData output;
  m_amr->allocate(output, m_realm, m_phase, ncomps, 1);

  // Copy data to put data holder
  m_amr->copyData(output, m_phi, Interval(0, 0), Interval(0, 0), CopyStrategy::ValidGhost, CopyStrategy::ValidGhost);
  m_amr->copyData(output,
                  flux,
                  Interval(0, SpaceDim - 1),
                  Interval(1, SpaceDim),
                  CopyStrategy::ValidGhost,
                  CopyStrategy::ValidGhost);
  m_amr->copyData(output,
                  m_source,
                  Interval(0, 0),
                  Interval(1 + SpaceDim, 1 + SpaceDim),
                  CopyStrategy::ValidGhost,
                  CopyStrategy::ValidGhost);
  m_amr->copyData(output,
                  m_resid,
                  Interval(0, 0),
                  Interval(2 + SpaceDim, 2 + SpaceDim),
                  CopyStrategy::ValidGhost,
                  CopyStrategy::ValidGhost);

  // Transform to centroid-centered
  m_amr->interpToCentroids(output, m_realm, phase::gas);

  // Alias this stuff
  Vector<LevelData<EBCellFAB>*> output_ptr;
  m_amr->alias(output_ptr, output);

  Vector<Real> covered_values(ncomps, 0.0);
  string       fname(file_char);

#ifdef CH_USE_HDF5
  constexpr int numPlotGhost = 0;

  DischargeIO::writeEBHDF5(fname,
                           names,
                           m_amr->getGrids(m_realm),
                           output_ptr,
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

#ifdef CH_USE_HDF5
void
EddingtonSP1::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const
{
  CH_TIME("EddingtonSP1::writeCheckpointLevel");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel" << endl;
  }

  // Write state vector
  write(a_handle, *m_phi[a_level], m_name + "_phi");
  write(a_handle, *m_source[a_level], m_name + "_src");
}
#endif

#ifdef CH_USE_HDF5
void
EddingtonSP1::readCheckpointLevel(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("EddingtonSP1::readCheckpointLevel");
  if (m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevel" << endl;
  }

  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name + "_phi", m_amr->getGrids(m_realm)[a_level], Interval(0, 0), false);
  read<EBCellFAB>(a_handle,
                  *m_source[a_level],
                  m_name + "_src",
                  m_amr->getGrids(m_realm)[a_level],
                  Interval(0, 0),
                  false);
}
#endif

#include <CD_NamespaceFooter.H>
