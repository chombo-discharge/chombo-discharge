/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*! 
  @file   CD_FieldSolver.cpp
  @brief  Implementation of CD_FieldSolver.H
  @author Robert Marskar
*/

// Std includes
#include <iostream>

// Chombo includes
#include <EBArith.H>
#include <ParmParse.H>
#include <MFAMRIO.H>
#include <EBAMRIO.H>

// Our includes
#include <CD_FieldSolver.H>
#include <CD_MultifluidAlias.H>
#include <CD_DataOps.H>
#include <CD_DischargeIO.H>
#include <CD_Units.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int FieldSolver::m_comp;
constexpr int FieldSolver::m_nComp;

FieldSolver::FieldSolver()
{
  CH_TIME("FieldSolver::FieldSolver()");

  // Default settings.
  m_className    = "FieldSolver";
  m_realm        = Realm::Primal;
  m_isVoltageSet = false;
  m_regridSlopes = true;
  m_verbosity    = -1;

  this->setDataLocation(Location::Cell::Center);
  this->setDefaultDomainBcFunctions();
}

FieldSolver::~FieldSolver() {}

void
FieldSolver::setDataLocation(const Location::Cell a_dataLocation)
{
  CH_TIME("FieldSolver::setDataLocation(Location::Cell)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setDataLocation(Location::Cell)" << endl;
  }

  // This sets the data location for the FieldSolver -- we issue a warning regarding centroid discretizations because
  // the rest of chombo-discharge has not caught up with the centroid formulation.

  switch (a_dataLocation) {
  case Location::Cell::Center: {
    m_dataLocation = a_dataLocation;
    m_faceLocation = Location::Face::Center;

    break;
  }
  case Location::Cell::Centroid: {
    MayDay::Error(
      "FieldSolver::setDataLocation - centroid discretization is not yet fully supported in chombo-discharge");

    m_dataLocation = a_dataLocation;
    m_faceLocation = Location::Face::Centroid;

    break;
  }
  default: {
    MayDay::Error("FieldSolver::setDataLocation - location must be either cell center or cell centroid");

    break;
  }
  }
}

bool
FieldSolver::solve(const bool a_zeroPhi)
{
  CH_TIME("FieldSolver::solve(bool)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::solve(bool)" << endl;
  }

  CH_assert(m_isVoltageSet);

  const bool converged = this->solve(m_potential, m_rho, a_zeroPhi);

  return converged;
}

bool
FieldSolver::solve(MFAMRCellData& a_potential, const bool a_zeroPhi)
{
  CH_TIME("FieldSolver::solve(MFAMRCellData, bool)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::solve(MFAMRCellData, bool)" << endl;
  }

  CH_assert(m_isVoltageSet);

  const bool converged = this->solve(a_potential, m_rho, a_zeroPhi);

  return converged;
}

void
FieldSolver::setSolverPermittivities(const MFAMRCellData& a_permittivityCell,
                                     const MFAMRFluxData& a_permittivityFace,
                                     const MFAMRIVData&   a_permittivityEB)
{
  CH_TIME("FieldSolver::setSolverPermittivities");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setSolverPermittivities" << endl;
  }
}

void
FieldSolver::computeElectricField()
{
  CH_TIME("FieldSolver::computeElectricField()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::computeElectricField()" << endl;
  }

  this->computeElectricField(m_electricField, m_potential);
}

void
FieldSolver::allocateInternals()
{
  CH_TIME("FieldSolver::allocateInternals()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::allocateInternals()" << endl;
  }

  m_amr->allocate(m_potential, m_realm, m_nComp);
  m_amr->allocate(m_rho, m_realm, m_nComp);
  m_amr->allocate(m_residue, m_realm, m_nComp);
  m_amr->allocate(m_permittivityCell, m_realm, m_nComp);
  m_amr->allocate(m_permittivityFace, m_realm, m_nComp);
  m_amr->allocate(m_permittivityEB, m_realm, m_nComp);
  m_amr->allocate(m_electricField, m_realm, SpaceDim);

  m_amr->allocate(m_sigma, m_realm, phase::gas, m_nComp);

  DataOps::setValue(m_potential, 0.0);
  DataOps::setValue(m_rho, 0.0);
  DataOps::setValue(m_sigma, 0.0);
  DataOps::setValue(m_residue, 0.0);
  DataOps::setValue(m_electricField, 0.0);

  // Set permittivities
  this->setPermittivities();
}

void
FieldSolver::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("FieldSolver::preRegrid(int, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::preRegrid(int, int)" << endl;
  }

  // TLDR: This version does a backup of m_potential on the old grid. This is later used for interpolating onto
  //       the new grid (in the regrid routine).

  m_amr->allocate(m_cache, m_realm, m_nComp);

  for (int lvl = 0; lvl <= a_oldFinestLevel; lvl++) {
    m_potential[lvl]->localCopyTo(*m_cache[lvl]);
  }
}

void
FieldSolver::computeDisplacementField(MFAMRCellData& a_displacementField, const MFAMRCellData& a_electricField)
{
  CH_TIME("FieldSolver::computeDisplacementField(MFAMRCellData, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::computeDisplacementField(MFAMRCellData, MFAMRCellData)" << endl;
  }

  CH_assert(a_displacementField[0]->nComp() == SpaceDim);
  CH_assert(a_electricField[0]->nComp() == SpaceDim);

  // TLDR: This computes the displacement field D = eps*E on both phases. The data is either cell-centered or centroid-centered, and the
  //       permittivity can be spatially varying (in the dielectric). So, for the gas phase we only need to compute eps0*E while for the
  //       dielectric phase we have to iterate through each cell and find the corresponding permittivity.

  const Vector<Dielectric>& dielectrics = m_computationalGeometry->getDielectrics();

  // Make D = eps0*E
  DataOps::copy(a_displacementField, a_electricField);
  DataOps::scale(a_displacementField, Units::eps0);

  if (m_multifluidIndexSpace->numPhases() > 1 && dielectrics.size() > 0) {

    // Lambda which returns the relative permittivity at some position
    auto permittivity = [&dielectrics](const RealVect& pos) -> Real {
      Real minDist = std::numeric_limits<Real>::infinity();
      int  closest = 0;

      for (int i = 0; i < dielectrics.size(); i++) {
        const RefCountedPtr<BaseIF> func = dielectrics[i].getImplicitFunction();

        const Real curDist = func->value(pos);

        if (std::abs(curDist) <= std::abs(minDist)) {
          minDist = curDist;
          closest = i;
        }
      }

      return dielectrics[closest].getPermittivity(pos);
    };

    // Iterate through data
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real               dx     = m_amr->getDx()[lvl];
      const RealVect           probLo = m_amr->getProbLo();
      const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, phase::solid)[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const Box        cellBox = dbl[dit()];
        const EBISBox&   ebisbox = ebisl[dit()];
        const EBGraph&   ebgraph = ebisbox.getEBGraph();
        const IntVectSet ivs     = ebisbox.getIrregIVS(cellBox);

        // Get handle to data on the solid phase
        MFCellFAB& D    = (*a_displacementField[lvl])[dit()];
        EBCellFAB& Dsol = D.getPhase(phase::solid);
        FArrayBox& Dreg = Dsol.getFArrayBox();

        // Iteration space for irregular cells.
        VoFIterator vofit(ivs, ebgraph);

        // Regular kernel
        auto regularKernel = [&](const IntVect& iv) -> void {
          if (ebisbox.isRegular(iv)) {
            const RealVect pos    = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
            const Real     epsRel = permittivity(pos);

            for (int comp = 0; comp < SpaceDim; comp++) {
              Dreg(iv, comp) *= epsRel;
            }
          }
        };

        // Irregular kernel
        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const RealVect pos    = probLo + Location::position(m_dataLocation, vof, ebisbox, dx);
          const Real     epsRel = permittivity(pos);

          for (int comp = 0; comp < SpaceDim; comp++) {
            Dsol(vof, comp) *= epsRel;
          }
        };

        // Launch kernels.
        BoxLoops::loop(cellBox, regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }
  }
}

Real
FieldSolver::computeEnergy(const MFAMRCellData& a_electricField)
{
  CH_TIME("FieldSolver::computeEnergy(MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::computeEnergy(MFAMRCellData)" << endl;
  }

  CH_assert(a_electricField[0]->nComp() == SpaceDim);

  // TLDR: This routine computes Int(E*D dV) over the entire domain. Since we use conservative averaging, we coarsen E*D onto
  //       the coarsest grid level and do the integratino there.

  const bool reallyMultiPhase = (m_multifluidIndexSpace->numPhases() > 1);

  // Allocate storage first, because we have no scratch storage for this.
  MFAMRCellData D;
  MFAMRCellData EdotD;

  m_amr->allocate(D, m_realm, SpaceDim);
  m_amr->allocate(EdotD, m_realm, 1);

  // Compute D = eps*E and then the dot product E*D.
  this->computeDisplacementField(D, a_electricField);
  DataOps::dotProduct(EdotD, D, a_electricField);
  m_amr->conservativeAverage(EdotD, m_realm);

  // This defines a lambda which computes the energy in a specified phase
  auto phaseEnergy = [&EdotD, &amr = this->m_amr](const phase::which_phase a_phase) -> Real {
    const EBAMRCellData phaseEdotD = amr->alias(a_phase, EdotD);

    const Real energy = DataOps::norm(*phaseEdotD[0], 1);

    return energy;
  };

  // Contributions from each phase.
  const Real Ug = phaseEnergy(phase::gas);
  const Real Us = reallyMultiPhase ? phaseEnergy(phase::solid) : 0.0;

  return 0.5 * (Ug + Us);
}

Real
FieldSolver::computeCapacitance()
{
  CH_TIME("FieldSolver::computeCapacitance()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::computeCapacitance()" << endl;
  }

  CH_assert(m_isVoltageSet);

  // TLDR: The energy density U = 0.5*C*V^2, or U = Int(E*D dV). We set the potential to one and solve the Poisson equation.
  //       We then compute U from E*D without sources and use C = 2*U/(V*V). The caveat to this approach is that the solver
  //       may have been set with a zero voltage, non-zero space charge and non-zero surface charge. Here, we do a "clean" solve
  //       with voltage set to 1, and without charges.

  MFAMRCellData phi;
  MFAMRCellData E;
  MFAMRCellData source;
  EBAMRIVData   sigma;

  m_amr->allocate(phi, m_realm, 1);
  m_amr->allocate(E, m_realm, SpaceDim);
  m_amr->allocate(source, m_realm, 1);
  m_amr->allocate(sigma, m_realm, phase::gas, 1);

  DataOps::setValue(phi, 0.0);
  DataOps::setValue(E, 0.0);
  DataOps::setValue(source, 0.0);
  DataOps::setValue(sigma, 0.0);

  // Do a backup of the voltage.
  auto voltageBackup = m_voltage;
  auto voltageOne    = [](const Real a_time) -> Real {
    return 1.0;
  };

  // Set the voltage to one and use that to compute the potential/E-field
  this->setVoltage(voltageOne);

  // Do a solve without a source term.
  this->solve(phi, source, sigma, true);
  this->computeElectricField(E, phi);

  // The energy is U = 0.5*CV^2 so we can just compute the electrostatic energy and revert that expression.
  const Real U = this->computeEnergy(E); // Electrostatic energy
  const Real C = 2.0 * U;                // C = 2*U/(V*V) but voltage is 1.

  // Set the voltage back to what it was.
  this->setVoltage(voltageBackup);

  return C;
}

void
FieldSolver::deallocateInternals()
{
  CH_TIME("FieldSolver::deallocateInternals()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::deallocateInternals()" << endl;
  }

  m_amr->deallocate(m_potential);
  m_amr->deallocate(m_rho);
  m_amr->deallocate(m_residue);
  m_amr->deallocate(m_sigma);
}

void
FieldSolver::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("FieldSolver::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::regrid(int, int, int)" << endl;
  }

  const Interval interv(m_comp, m_comp);

  // Reallocate internals
  this->allocateInternals();

  // Regrid potential data holder
  m_amr->interpToNewGrids(m_potential, m_cache, a_lmin, a_oldFinestLevel, a_newFinestLevel, m_regridSlopes);

  // Synchronize over levels.
  m_amr->conservativeAverage(m_potential, m_realm);
  m_amr->interpGhost(m_potential, m_realm);

  // Recompute E from the new potential.
  this->computeElectricField();

  // Set permittivities
  this->setPermittivities();

  // Deallocate the scratch storage.
  m_amr->deallocate(m_cache);
}

void
FieldSolver::setRho(const Real a_rho)
{
  CH_TIME("FieldSolver::setRho(Real");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setRho(Real)" << endl;
  }

  DataOps::setValue(m_rho, a_rho);
}

void
FieldSolver::setRho(const std::function<Real(const RealVect)>& a_rho)
{
  CH_TIME("FieldSolver::setRho(std::function<Real(const RealVect)>)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setRho(std::function<Real(const RealVect)>))" << endl;
  }

  DataOps::setValue(m_rho, a_rho, m_amr->getProbLo(), m_amr->getDx(), m_comp);
}

void
FieldSolver::setSigma(const Real a_sigma)
{
  CH_TIME("FieldSolver::setSigma(Real");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setSigma(Real)" << endl;
  }

  DataOps::setValue(m_sigma, a_sigma);
}

void
FieldSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("FieldSolver::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)" << endl;
  }

  CH_assert(!a_computationalGeometry.isNull());

  m_computationalGeometry = a_computationalGeometry;
  m_multifluidIndexSpace  = m_computationalGeometry->getMfIndexSpace();

  this->setDefaultEbBcFunctions();
}

void
FieldSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr)
{
  CH_TIME("FieldSolver::setAmr(RefCountedPtr<AmrMesh>)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setAmr(RefCountedPtr<AmrMesh>)" << endl;
  }

  CH_assert(!a_amr.isNull());

  m_amr = a_amr;
}

void
FieldSolver::setVoltage(std::function<Real(const Real a_time)> a_voltage)
{
  CH_TIME("FieldSolver::setVoltage(std::function<Real(const Real a_time)>)");
  if (m_verbosity > 4) {
    pout() << "FieldSolver::setVoltage(std::function<Real(const Real a_time)>)" << endl;
  }

  m_voltage      = a_voltage;
  m_isVoltageSet = true;
}

void
FieldSolver::setElectrodeDirichletFunction(const int a_electrode, const ElectrostaticEbBc::BcFunction& a_function)
{
  CH_TIME("FieldSolver::setElectrodeDirichletFunction(int, ElectrostaticEbBc::BcFunction)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setElectrodeDirichletFunction(int, ElectrostaticEbBc::BcFunction)" << endl;
  }

  m_ebBc.setEbBc(a_electrode, a_function);
}

void
FieldSolver::setVerbosity(const int a_verbosity)
{
  CH_TIME("FieldSolver::setVerbosity(int)");
  m_verbosity = a_verbosity;

  if (m_verbosity > 4) {
    pout() << "FieldSolver::setVerbosity(int)" << endl;
  }
}

void
FieldSolver::setTime(const int a_timeStep, const Real a_time, const Real a_dt)
{
  CH_TIME("FieldSolver::setTime(int, Real, Real)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setTime(int, Real, Real)" << endl;
  }

  m_timeStep = a_timeStep;
  m_time     = a_time;
  m_dt       = a_dt;
}

void
FieldSolver::setRealm(const std::string a_realm)
{
  CH_TIME("FieldSolver::setRealm(std::string)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setRealm(std::string)" << endl;
  }

  m_realm = a_realm;
}

std::string
FieldSolver::getRealm() const
{
  CH_TIME("FieldSolver::getRealm()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getRealm()" << endl;
  }

  return m_realm;
}

void
FieldSolver::parseVerbosity()
{
  CH_TIME("FieldSolver::parseVerbosity()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::parseVerbosity()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("verbosity", m_verbosity);
}

void
FieldSolver::parsePlotVariables()
{
  CH_TIME("FieldSolver::parsePlotVariables()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::parsePlotVariables()" << endl;
  }

  m_plotPotential          = false;
  m_plotRho                = false;
  m_plotElectricField      = false;
  m_plotResidue            = false;
  m_plotPermittivity       = false;
  m_plotSigma              = false;
  m_plotElectricFieldSolid = false;

  ParmParse pp(m_className.c_str());
  const int num = pp.countval("plt_vars");

  if (num > 0) {
    Vector<std::string> str(num);
    pp.getarr("plt_vars", str, 0, num);

    for (int i = 0; i < num; i++) {
      if (str[i] == "phi")
        m_plotPotential = true;
      else if (str[i] == "rho")
        m_plotRho = true;
      else if (str[i] == "resid")
        m_plotResidue = true;
      else if (str[i] == "E")
        m_plotElectricField = true;
      else if (str[i] == "Esol")
        m_plotElectricFieldSolid = true;
      else if (str[i] == "sigma")
        m_plotSigma = true;
      else if (str[i] == "perm")
        m_plotPermittivity = true;
    }
  }
}

void
FieldSolver::parseRegridSlopes()
{
  CH_TIME("FieldSolver::parseRegridSlopes()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::parseRegridSlopes()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("use_regrid_slopes", m_regridSlopes);
}

std::string
FieldSolver::makeBcString(const int a_dir, const Side::LoHiSide a_side) const
{
  CH_TIME("FieldSolver::makeBcString(int, Side::LoHiSide)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::makeBcString(int, Side::LoHiSide)" << endl;
  }

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

ElectrostaticDomainBc::BcType
FieldSolver::parseBcString(const std::string a_str) const
{
  CH_TIME("FieldSolver::parseBcString(std::string)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::parseBcString(std::string)" << endl;
  }

  ElectrostaticDomainBc::BcType ret;

  if (a_str == "dirichlet") {
    ret = ElectrostaticDomainBc::BcType::Dirichlet;
  }
  else if (a_str == "neumann") {
    ret = ElectrostaticDomainBc::BcType::Neumann;
  }
  else {
    MayDay::Error("ElectrostaticDomainBc::BcType - unknown BC type!");
  }

  return ret;
}

void
FieldSolver::setDefaultDomainBcFunctions()
{
  CH_TIME("FieldSolver::setDefaultDomainBcFunctions()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setDefaultDomainBcFunctions()" << endl;
  }

  // Default space/time dependency of domain BCs.
  auto defaultDomainBcFunction = [](const RealVect a_position, const Real a_time) -> Real {
    return 1.0;
  };

  m_domainBcFunctions.clear();
  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      m_domainBcFunctions.emplace(std::make_pair(dir, sit()), defaultDomainBcFunction);
    }
  }
}

void
FieldSolver::setDefaultEbBcFunctions()
{
  CH_TIME("FieldSolver::setDefaultEbBcFunctions()");
  CH_assert(!m_computationalGeometry.isNull());
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setDefaultEbBcFunctions()" << endl;
  }

  // TLDR: This sets the default potential function on the electrodes. We use a lambda for passing in
  //       some member variables (m_voltage, m_time) and get other parameters directly from the electrodes.
  //       We create one such voltage function for each electrode so that the voltage can be obtained
  //       anywhere on the electrode.

  const Vector<Electrode> electrodes = m_computationalGeometry->getElectrodes();

  for (int i = 0; i < electrodes.size(); i++) {
    const Electrode& elec = electrodes[i];

    const Real val  = (elec.isLive()) ? 1.0 : 0.0;
    const Real frac = elec.getFraction();

    ElectrostaticEbBc::BcFunction curFunc =
      [&, &time = this->m_time, &voltage = this->m_voltage, val, frac](const RealVect a_position, const Real a_time) {
        return voltage(time) * val * frac;
      };

    m_ebBc.addEbBc(elec, curFunc);
  }
}

void
FieldSolver::setDomainSideBcFunction(const int                                a_dir,
                                     const Side::LoHiSide                     a_side,
                                     const ElectrostaticDomainBc::BcFunction& a_function)
{
  CH_TIME("FieldSolver::setDomainSideBcFunction(int, Side::LoHiSide, ElectrostaticDomainBc::BcFunction)");
  if (m_verbosity > 4) {
    pout() << "FieldSolver::setDomainSideBcFunction(int, Side::LoHiSide, ElectrostaticDomainBc::BcFunction)" << endl;
  }

  const ElectrostaticDomainBc::DomainSide domainSide = std::make_pair(a_dir, a_side);

  m_domainBcFunctions.at(domainSide) = a_function;
}

void
FieldSolver::parseDomainBc()
{
  CH_TIME("FieldSolver::parseDomainBc()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::parseDomainBc()" << endl;
  }

  // TLDR: This routine might seem big and complicated. What we are doing is that we are creating one function object which returns some value
  //       anywhere in space and time on a domain edge (face). The FieldSolver class supports Dirichlet and Neumann, and the below code simply
  //       creates those functions and associates them with an edge.
  //
  //       For flexibility we want to be able to specify the potential directy without invoking m_voltage, while at the same time we want to offer
  //       the simplistic method of setting a domain side to be "grounded", "live", or otherwise given by some fraction of m_voltage.
  //       We thus make a distinction between "dirichlet" and "dirichlet_custom". The difference between these is that for "dirichlet_custom" the contents
  //       of m_domainBcFunctions are used as boundary conditions. For "dirichlet 0.5" the contents of m_domainBcFunctions are multiplied by 0.5*m_voltage.
  //       The same approach is used for Neumann boundary conditions.

  ParmParse pp(m_className.c_str());

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {

      const ElectrostaticDomainBc::DomainSide domainSide = std::make_pair(dir, sit());
      const std::string                       bcString   = this->makeBcString(dir, sit());
      const int                               num        = pp.countval(bcString.c_str());

      std::string str;

      ElectrostaticDomainBc::BcType      bcType;
      ElectrostaticDomainBc::BcFunction& bcFunc = m_domainBcFunctions.at(domainSide); // = F(x,t) in the comments below

      std::function<Real(const RealVect, const Real)> curFunc;

      if (
        num ==
        1) { // If we only had one argument the user has asked for "custom" boundary conditions, and we then pass in F(x,t) directly (evaluated at m_time)
        pp.get(bcString.c_str(), str, 0);

        curFunc = [&bcFunc, &time = this->m_time](const RealVect a_pos, const Real a_time) {
          return bcFunc(a_pos, time);
        };

        if (str == "dirichlet_custom") {
          bcType = ElectrostaticDomainBc::BcType::Dirichlet;
        }
        else if (str == "neumann_custom") {
          bcType = ElectrostaticDomainBc::BcType::Neumann;
        }
        else {
          MayDay::Error(
            "FieldSolver::parseDomainBc -- got only one argument but this argument was not dirichlet/neumann_custom. Maybe you have the wrong BC specification?");
        }
      }
      else if (
        num ==
        2) { // If we had two arguments the user has asked to run with less verbose specifications. E.g. "dirichlet 0.5" => V(x,t) = V(t) * 0.5 * F(x,t)
        Real val;

        pp.get(bcString.c_str(), str, 0);
        pp.get(bcString.c_str(), val, 1);

        bcType = this->parseBcString(str);

        // Build a function computing the value at the boundary.
        switch (bcType) {
        case ElectrostaticDomainBc::BcType::Dirichlet:
          curFunc = [&bcFunc, &voltage = this->m_voltage, &time = this->m_time, val](const RealVect a_pos,
                                                                                     const Real     a_time) {
            return bcFunc(a_pos, time) * voltage(time) * val;
          };
          break;
        case ElectrostaticDomainBc::BcType::Neumann:
          curFunc = [&bcFunc, &time = this->m_time, val](const RealVect a_pos, const Real a_time) {
            return bcFunc(a_pos, time) * val;
          };
          break;
        default:
          MayDay::Error("FieldSolver::parseDomainBc -- unsupported boundary condition requested!");
          break;
        }
      }
      else {
        const std::string errorString = "FieldSolver::parseDomainBc -- bad or no input parameter for " + bcString;
        MayDay::Error(errorString.c_str());
      }

      m_domainBc.setBc(domainSide, std::make_pair(bcType, curFunc));
    }
  }
}

void
FieldSolver::setPermittivities()
{
  CH_TIME("FieldSolver::setPermittivities()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::setPermittivities()" << endl;
  }

  const Real relEpsGas = m_computationalGeometry->getGasPermittivity();

  DataOps::setValue(m_permittivityCell, relEpsGas);
  DataOps::setValue(m_permittivityFace, relEpsGas);
  DataOps::setValue(m_permittivityEB, relEpsGas);

  const Vector<Dielectric>& dielectrics = m_computationalGeometry->getDielectrics();

  if (dielectrics.size() > 0 && m_multifluidIndexSpace->numPhases() > 1) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
      const Real               dx     = m_amr->getDx()[lvl];
      const RealVect           probLo = m_amr->getProbLo();

      LevelData<EBCellFAB>       cellPerm;
      LevelData<EBFluxFAB>       facePerm;
      LevelData<BaseIVFAB<Real>> ebPerm;

      MultifluidAlias::aliasMF(cellPerm, phase::solid, *m_permittivityCell[lvl]);
      MultifluidAlias::aliasMF(facePerm, phase::solid, *m_permittivityFace[lvl]);
      MultifluidAlias::aliasMF(ebPerm, phase::solid, *m_permittivityEB[lvl]);

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        EBCellFAB&       cellPermFAB = cellPerm[dit()];
        EBFluxFAB&       facePermFAB = facePerm[dit()];
        BaseIVFAB<Real>& ebPermFAB   = ebPerm[dit()];
        const Box        cellBox     = dbl[dit()];
        const EBISBox&   ebisbox     = cellPermFAB.getEBISBox();

        this->setCellPermittivities(cellPermFAB, cellBox, ebisbox, probLo, dx);
        this->setFacePermittivities(facePermFAB, cellBox, ebisbox, probLo, dx);
        this->setEbPermittivities(ebPermFAB, cellBox, ebisbox, probLo, dx);
      }
    }
  }
}

void
FieldSolver::setCellPermittivities(EBCellFAB&      a_relPerm,
                                   const Box&      a_cellBox,
                                   const EBISBox&  a_ebisbox,
                                   const RealVect& a_probLo,
                                   const Real&     a_dx)
{
  CH_TIME("FieldSolver::setCellPermittivities(EBCellFAB, Box, EBISBox, RealVect, Real)");
  if (m_verbosity > 10) {
    pout() << "FieldSolver::setCellPermittivities(EBCellFAB, Box, EBISBox, RealVect, Real)" << endl;
  }

  CH_assert(a_relPerm.nComp() == 1);

  const EBGraph&   ebgraph = a_ebisbox.getEBGraph();
  const IntVectSet irreg   = a_ebisbox.getIrregIVS(a_cellBox);

  BaseFab<Real>& relPermFAB = a_relPerm.getSingleValuedFAB();

  // Regular kernel
  auto regularKernel = [&](const IntVect& iv) -> void {
    if (a_ebisbox.isRegular(iv)) {
      const RealVect pos     = a_probLo + (0.5 * RealVect::Unit + RealVect(iv)) * a_dx;
      relPermFAB(iv, m_comp) = this->getDielectricPermittivity(pos);
    }
  };

  // Irregular kernel
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const RealVect pos     = Location::position(m_dataLocation, vof, a_ebisbox, a_dx);
    a_relPerm(vof, m_comp) = this->getDielectricPermittivity(pos);
  };

  // Kernel regions.
  VoFIterator irregRegion(irreg, ebgraph);

  // Launch kernels.
  BoxLoops::loop(a_cellBox, regularKernel);
  BoxLoops::loop(irregRegion, irregularKernel);
}

void
FieldSolver::setFacePermittivities(EBFluxFAB&      a_relPerm,
                                   const Box&      a_cellBox,
                                   const EBISBox&  a_ebisbox,
                                   const RealVect& a_probLo,
                                   const Real&     a_dx)
{
  CH_TIME("FieldSolver::setFacePermittivities(EBFluxFAB, Box, EBISBox, RealVect, Real)");
  if (m_verbosity > 10) {
    pout() << "FieldSolver::setFacePermittivities(EBFluxFAB, Box, EBISBox, RealVect, Real)" << endl;
  }

  CH_assert(a_relPerm.nComp() == 1);

  const EBGraph&   ebgraph = a_ebisbox.getEBGraph();
  const IntVectSet irreg   = a_ebisbox.getIrregIVS(a_cellBox);

  for (int dir = 0; dir < SpaceDim; dir++) {

    // Kernel regions.
    const Box    facebox = surroundingNodes(a_cellBox, dir);
    FaceIterator faceit  = FaceIterator(irreg, ebgraph, dir, FaceStop::SurroundingWithBoundary);

    // Single-valued data.
    BaseFab<Real>& relPermFAB = a_relPerm[dir].getSingleValuedFAB();

    // Regular kernel
    auto regularKernel = [&](const IntVect& iv) -> void {
      const RealVect pos     = a_probLo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit) - 0.5 * a_dx * BASISREALV(dir);
      relPermFAB(iv, m_comp) = this->getDielectricPermittivity(pos);
    };

    // Irregular kernel.
    auto irregularKernel = [&](const FaceIndex& face) -> void {
      const RealVect pos           = a_probLo + Location::position(m_faceLocation, face, a_ebisbox, a_dx);
      a_relPerm[dir](face, m_comp) = this->getDielectricPermittivity(pos);
    };

    // Launch kernels.
    BoxLoops::loop(facebox, regularKernel);
    BoxLoops::loop(faceit, irregularKernel);
  }
}

void
FieldSolver::setEbPermittivities(BaseIVFAB<Real>& a_relPerm,
                                 const Box&       a_cellBox,
                                 const EBISBox&   a_ebisbox,
                                 const RealVect&  a_probLo,
                                 const Real&      a_dx)
{
  CH_TIME("FieldSolver::setEbPermittivities(BaseIVFAB<Real>, Box, EBISBox, RealVect, Real)");
  if (m_verbosity > 10) {
    pout() << "FieldSolver::setEbPermittivities(BaseIVFAB<Real>, Box, EBISBox, RealVect, Real)" << endl;
  }

  CH_assert(a_relPerm.nComp() == 1);

  const IntVectSet& ivs     = a_relPerm.getIVS();
  const EBGraph&    ebgraph = a_relPerm.getEBGraph();

  VoFIterator vofit(ivs, ebgraph);

  auto kernel = [&](const VolIndex& vof) -> void {
    const RealVect pos = Location::position(Location::Cell::Boundary, vof, a_ebisbox, a_dx);

    a_relPerm(vof, m_comp) = this->getDielectricPermittivity(pos);
  };

  BoxLoops::loop(vofit, kernel);
}

void
FieldSolver::writePlotFile()
{
  CH_TIME("FieldSolver::writePlotFile()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::writePlotFile()" << endl;
  }

  // Number of output components
  const int numPlotVars = this->getNumberOfPlotVariables();

  // Get plot variable names
  const Vector<std::string> plotVarNames = getPlotVariableNames();

  // Allocate storage for output
  EBAMRCellData output;
  m_amr->allocate(output, m_realm, phase::gas, numPlotVars);

  // Copy internal data to be plotted over to 'output'
  int icomp = 0;
  this->writePlotData(output, icomp, true);

  // Filename
  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_className.c_str(), m_timeStep, SpaceDim);

  // Need to alias because EBHDF5 is not too smart.
  Vector<LevelData<EBCellFAB>*> outputPtr(1 + m_amr->getFinestLevel());
  m_amr->alias(outputPtr, output);

  Vector<Real> covered_values(numPlotVars, 0.0);
  string       fname(file_char);

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

#ifdef CH_USE_HDF5
void
FieldSolver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const
{
  CH_TIME("FieldSolver::writeCheckpointLevel(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::writeCheckpointLevel(HDF5Handle, int)" << endl;
  }

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Used for aliasing phases
  LevelData<EBCellFAB> potentialGas;
  LevelData<EBCellFAB> potentialSol;

  if (!ebisGas.isNull()) {
    MultifluidAlias::aliasMF(potentialGas, phase::gas, *m_potential[a_level]);
  }
  if (!ebisSol.isNull()) {
    MultifluidAlias::aliasMF(potentialSol, phase::solid, *m_potential[a_level]);
  }

  // Write data
  if (!ebisGas.isNull()) {
    write(a_handle, potentialGas, "FieldSolver::m_potential(gas)");
  }
  if (!ebisSol.isNull()) {
    write(a_handle, potentialSol, "FieldSolver::m_potential(solid)");
  }
}
#endif

#ifdef CH_USE_HDF5
void
FieldSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("FieldSolver::readCheckpointLevel(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::readCheckpointLevel(HDF5Handle, int)" << endl;
  }

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Used for aliasing phases
  LevelData<EBCellFAB> potentialGas;
  LevelData<EBCellFAB> potentialSol;

  if (!ebisGas.isNull()) {
    MultifluidAlias::aliasMF(potentialGas, phase::gas, *m_potential[a_level]);
  }
  if (!ebisSol.isNull()) {
    MultifluidAlias::aliasMF(potentialSol, phase::solid, *m_potential[a_level]);
  }

  // Read data
  if (!ebisGas.isNull()) {
    read<EBCellFAB>(a_handle,
                    potentialGas,
                    "FieldSolver::m_potential(gas)",
                    m_amr->getGrids(m_realm)[a_level],
                    Interval(0, 0),
                    false);
  }
  if (!ebisSol.isNull()) {
    read<EBCellFAB>(a_handle,
                    potentialSol,
                    "FieldSolver::m_potential(solid)",
                    m_amr->getGrids(m_realm)[a_level],
                    Interval(0, 0),
                    false);
  }
}
#endif

void
FieldSolver::postCheckpoint()
{
  CH_TIME("FieldSolver::postCheckpoint()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::postCheckpoint()" << endl;
  }
}

void
FieldSolver::writePlotData(EBAMRCellData& a_output, int& a_comp, const bool a_forceNoInterp)
{
  CH_TIME("FieldSolver::writePlotData(EBAMRCellData, int, bool)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::writePlotData(EBAMRCellData, int, bool)" << endl;
  }

  // This routine always outputs data on the centroid. If the data was defined on the center we move it to the centroid (but not forcefully)
  const bool doInterp = (m_dataLocation == Location::Cell::Center) && !a_forceNoInterp;

  // Add phi to output
  if (m_plotPotential) {
    this->writeMultifluidData(a_output, a_comp, m_potential, phase::gas, doInterp);
  }
  if (m_plotRho) {
    this->writeMultifluidData(a_output, a_comp, m_rho, phase::gas, false);
  }
  if (m_plotSigma) {
    this->writeSurfaceData(a_output, a_comp, m_sigma);
  }
  if (m_plotResidue) {
    this->writeMultifluidData(a_output, a_comp, m_residue, phase::gas, false);
  }
  if (m_plotPermittivity) {
    this->writeMultifluidData(a_output, a_comp, m_permittivityCell, phase::gas, false);
  }
  if (m_plotElectricField) {
    this->writeMultifluidData(a_output, a_comp, m_electricField, phase::gas, doInterp);
  }
  if (m_plotElectricFieldSolid) {
    this->writeMultifluidData(a_output, a_comp, m_electricField, phase::solid, doInterp);
  }
}

void
FieldSolver::writeMultifluidData(EBAMRCellData&           a_output,
                                 int&                     a_comp,
                                 const MFAMRCellData&     a_data,
                                 const phase::which_phase a_phase,
                                 const bool               a_interp) const
{
  CH_TIME("FieldSolver::writeMultifluidData(EBAMRCellData, int, MFAMRCellData, bool)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::writeMultifluidData(EBAMRCellData, int, MFAMRCellData, bool)" << endl;
  }

  // So the problem with the Chombo HDF5 I/O routines is that they are not really designed for multiphase. We happen to know that a_data can be multifluid data
  // which we want to put onto a single-phase data holder. There is ambiguity only in the cut-cells because the data there has multiple degrees of freedom. We also
  // happen to know that a_output is on the gas phase. Our issue is that we need to decide if the gas-side or solid-side data goes into the output data holder.

  const int numComp = a_data[0]->nComp();

  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  const bool reallyMultiPhase = !(ebisSol.isNull());

  // Aliases of the input multifluid data.
  EBAMRCellData aliasGas;
  EBAMRCellData aliasSolid;
  if (reallyMultiPhase) {
    aliasGas   = m_amr->alias(phase::gas, a_data);
    aliasSolid = m_amr->alias(phase::solid, a_data);
  }
  else {
    aliasGas = m_amr->alias(phase::gas, a_data);
  }

  // Allocate some scratch data that we can use. We happen to know that the
  // output data is always on the gas phase realm.
  EBAMRCellData scratchSolid;
  EBAMRCellData scratchGas;
  if (reallyMultiPhase) {
    m_amr->allocate(scratchGas, m_realm, phase::gas, numComp);
    m_amr->allocate(scratchSolid, m_realm, phase::solid, numComp);
  }
  else {
    m_amr->allocate(scratchGas, m_realm, phase::gas, numComp);
  }

  // Copy into data holders and interpolate to centroids if we need to.
  if (reallyMultiPhase) {
    DataOps::copy(scratchGas, aliasGas);
    DataOps::copy(scratchSolid, aliasSolid);
  }
  else {
    DataOps::copy(scratchGas, aliasGas);
  }

  // Coarsen data.
  if (reallyMultiPhase) {
    m_amr->conservativeAverage(scratchGas, m_realm, phase::gas);
    m_amr->conservativeAverage(scratchSolid, m_realm, phase::solid);
  }
  else {
    m_amr->conservativeAverage(scratchGas, m_realm, phase::gas);
  }

  // Interpolate ghost cells.
  if (reallyMultiPhase) {
    m_amr->interpGhost(scratchGas, m_realm, phase::gas);
    m_amr->interpGhost(scratchSolid, m_realm, phase::solid);
  }
  else {
    m_amr->interpGhost(scratchGas, m_realm, phase::gas);
  }

  // Put cell-centered data on centroid.
  if (a_interp) {
    if (reallyMultiPhase) {
      m_amr->interpToCentroids(scratchGas, m_realm, phase::gas);
      m_amr->interpToCentroids(scratchSolid, m_realm, phase::solid);
    }
    else {
      m_amr->interpToCentroids(scratchGas, m_realm, phase::gas);
    }
  }

  // Go through all levels and grid patches and replace the covered gas-side scratch data with the regular solid-side scratch data. On cut-cells
  // we determine the data based on the a_phase input flag.
  if (reallyMultiPhase) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl        = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebislGas   = m_amr->getEBISLayout(m_realm, phase::gas)[lvl];
      const EBISLayout&        ebislSolid = m_amr->getEBISLayout(m_realm, phase::solid)[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const EBISBox& ebisBoxGas   = ebislGas[dit()];
        const EBISBox& ebisBoxSolid = ebislSolid[dit()];

        const bool isGasRegular   = ebisBoxGas.isAllRegular();
        const bool isGasCovered   = ebisBoxGas.isAllCovered();
        const bool isGasIrregular = !isGasRegular && !isGasCovered;

        const bool isSolidRegular   = ebisBoxSolid.isAllRegular();
        const bool isSolidCovered   = ebisBoxSolid.isAllCovered();
        const bool isSolidIrregular = !isSolidRegular && !isSolidCovered;

        FArrayBox&       fabGas   = (*scratchGas[lvl])[dit()].getFArrayBox();
        const FArrayBox& fabSolid = (*scratchSolid[lvl])[dit()].getFArrayBox();

        // This is true if the current EBISBox covers a region that is either inside the gas phase
        // or inside the solid phase.
        const bool validData = !(isGasCovered && isSolidCovered);
        if (validData) {

          if (isSolidRegular) {
            // In this case we are purely inside the solid region -- take the data from the solid phase.
            fabGas.copy(fabSolid);
          }
          else if (isGasIrregular && isSolidIrregular) {
            // In this case we are looking at a grid patch that lies on the gas-solid boundary. We need to determine which cells
            // go into the output region. We happen to know that all gas-side data is already filled, so we only need to grok
            // the solid-side data.

            for (int comp = 0; comp < numComp; comp++) {

              auto kernel = [&](const IntVect& iv) -> void {
                const bool coveredGas   = ebisBoxGas.isCovered(iv);
                const bool irregGas     = ebisBoxGas.isIrregular(iv);
                const bool regularSolid = ebisBoxSolid.isRegular(iv);
                const bool irregSolid   = ebisBoxSolid.isIrregular(iv);

                if (regularSolid && coveredGas) {
                  fabGas(iv, comp) = fabSolid(iv, comp);
                }
                else if (irregGas && irregSolid) {
                  if (a_phase == phase::solid) {
                    fabGas(iv, comp) = fabSolid(iv, comp);
                  }
                }
              };

              BoxLoops::loop(dbl[dit()], kernel);
            }
          }
        }
      }
    }
  }

  // Copy the single-phase data to the output data holder.
  const Interval srcInterv(0, numComp - 1);
  const Interval dstInterv(a_comp, numComp - 1 + a_comp);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    if (m_realm == a_output.getRealm()) {
      scratchGas[lvl]->localCopyTo(srcInterv, *a_output[lvl], dstInterv);
    }
    else {
      scratchGas[lvl]->copyTo(srcInterv, *a_output[lvl], dstInterv);
    }
  }

  a_comp += numComp;
}

void
FieldSolver::writeSurfaceData(EBAMRCellData& a_output, int& a_comp, const EBAMRIVData& a_data)
{
  CH_TIME("FieldSolver::writeSurfaceData(EBAMRCellData, int, EBAMRIVData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::writeSurfaceData(EBAMRCellData, int, EBAMRIVData)" << endl;
  }

  // Put a_data in volume format.
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, phase::gas, 1);

  DataOps::setValue(scratch, 0.0);
  DataOps::incr(scratch, a_data, 1.0);

  // Copy to a_output
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    if (m_realm == a_output.getRealm()) {
      scratch[lvl]->localCopyTo(Interval(0, 0), *a_output[lvl], Interval(a_comp, a_comp));
    }
    else {
      scratch[lvl]->copyTo(Interval(0, 0), *a_output[lvl], Interval(a_comp, a_comp));
    }
  }

  a_comp++;
}

int
FieldSolver::getNumberOfPlotVariables() const
{
  CH_TIME("FieldSolver::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getNumberOfPlotVariables()" << endl;
  }

  int numPltVars = 0;

  if (m_plotPotential) {
    numPltVars = numPltVars + 1;
  }
  if (m_plotRho) {
    numPltVars = numPltVars + 1;
  }
  if (m_plotSigma) {
    numPltVars = numPltVars + 1;
  }
  if (m_plotResidue) {
    numPltVars = numPltVars + 1;
  }
  if (m_plotPermittivity) {
    numPltVars = numPltVars + 1;
  }
  if (m_plotElectricField) {
    numPltVars = numPltVars + SpaceDim;
  }
  if (m_plotElectricFieldSolid) {
    numPltVars = numPltVars + SpaceDim;
  }

  return numPltVars;
}

Vector<std::string>
FieldSolver::getPlotVariableNames() const
{
  CH_TIME("FieldSolver::getPlotVariableNames()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getPlotVariableNames()" << endl;
  }

  Vector<std::string> pltVarNames(0);

  if (m_plotPotential) {
    pltVarNames.push_back("Electrostatic potential");
  }
  if (m_plotRho) {
    pltVarNames.push_back("Space charge density");
  }
  if (m_plotSigma) {
    pltVarNames.push_back("Electrostatic sigma");
  }
  if (m_plotResidue) {
    pltVarNames.push_back("Electrostatic potential_residue");
  }
  if (m_plotPermittivity) {
    pltVarNames.push_back("Electrostatic permittivity");
  }

  if (m_plotElectricField) {
    pltVarNames.push_back("x-Electric field");
    pltVarNames.push_back("y-Electric field");
    if (SpaceDim == 3) {
      pltVarNames.push_back("z-Electric field");
    }
  }

  if (m_plotElectricFieldSolid) {
    pltVarNames.push_back("x-Electric field solid");
    pltVarNames.push_back("y-Electric field solid");
    if (SpaceDim == 3) {
      pltVarNames.push_back("z-Electric field solid");
    }
  }

  return pltVarNames;
}

Vector<long long>
FieldSolver::computeLoads(const DisjointBoxLayout& a_dbl, const int a_level)
{
  CH_TIME("FieldSolver::computeLoads(DisjointBoxLayout, int)");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::computeLoads(DisjointBoxLayout, int)" << endl;
  }

  // Compute number of cells
  Vector<int> numCells(a_dbl.size(), 0);
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit) {
    numCells[dit().intCode()] = a_dbl[dit()].numPts();
  }

#ifdef CH_MPI
  int         count = numCells.size();
  Vector<int> tmp(count);
  MPI_Allreduce(&(numCells[0]), &(tmp[0]), count, MPI_INT, MPI_SUM, Chombo_MPI::comm);
  numCells = tmp;
#endif

  Vector<long long> loads(numCells.size());
  for (int i = 0; i < loads.size(); i++) {
    loads[i] = (long long)numCells[i];
  }

  return loads;
}

const std::function<Real(const Real a_time)>&
FieldSolver::getVoltage() const
{
  CH_TIME("FieldSolver::getVoltage()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getVoltage()" << endl;
  }

  return m_voltage;
}

Real
FieldSolver::getCurrentVoltage() const
{
  CH_TIME("FieldSolver::getCurrentVoltage()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getCurrentVoltage()" << endl;
  }

  return m_voltage(m_time);
}

Real
FieldSolver::getTime() const
{
  CH_TIME("FieldSolver::getTime()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getTime()" << endl;
  }

  return m_time;
}

MFAMRCellData&
FieldSolver::getPotential()
{
  CH_TIME("FieldSolver::getPotential()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getPotential()" << endl;
  }

  return m_potential;
}

MFAMRCellData&
FieldSolver::getElectricField()
{
  CH_TIME("FieldSolver::getElectricField()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getElectricField()" << endl;
  }

  return m_electricField;
}

MFAMRCellData&
FieldSolver::getRho()
{
  CH_TIME("FieldSolver::getRho()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getRho()" << endl;
  }

  return m_rho;
}

EBAMRIVData&
FieldSolver::getSigma()
{
  CH_TIME("FieldSolver::getSigma()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getSigma()" << endl;
  }

  return m_sigma;
}

MFAMRCellData&
FieldSolver::getResidue()
{
  CH_TIME("FieldSolver::getResidue()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getResidue()" << endl;
  }

  return m_residue;
}

MFAMRCellData&
FieldSolver::getPermittivityCell()
{
  CH_TIME("FieldSolver::getPermittivityCell()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getPermittivityCell()" << endl;
  }

  return m_permittivityCell;
}

MFAMRFluxData&
FieldSolver::getPermittivityFace()
{
  CH_TIME("FieldSolver::getPermittivityFace()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getPermittivityFace()" << endl;
  }

  return m_permittivityFace;
}

MFAMRIVData&
FieldSolver::getPermittivityEB()
{
  CH_TIME("FieldSolver::getPermittivityEB()");
  if (m_verbosity > 5) {
    pout() << "FieldSolver::getPermittivityEB()" << endl;
  }

  return m_permittivityEB;
}

#include <CD_NamespaceFooter.H>
