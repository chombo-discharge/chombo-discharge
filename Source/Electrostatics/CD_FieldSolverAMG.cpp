/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_FieldSolverAMG.cpp
  @brief  Implementation of CD_FieldSolver.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_FieldSolverAMG.H>
#include <CD_MFHelmholtzPETSc.H>
#include <CD_NamespaceHeader.H>

FieldSolverAMG::FieldSolverAMG()
{
  CH_TIME("FieldSolverAMG::FieldSolverAMG");

  m_className              = "FieldSolverAMG";
  m_realm                  = Realm::Primal;
  m_isVoltageSet           = false;
  m_regridSlopes           = true;
  m_verbosity              = 10;
  m_plotPotential          = false;
  m_plotRho                = false;
  m_plotElectricField      = false;
  m_plotResidue            = false;
  m_plotPermittivity       = false;
  m_plotSigma              = false;
  m_plotElectricFieldSolid = false;
}

FieldSolverAMG::~FieldSolverAMG()
{
  CH_TIME("FieldSolverAMG::~FieldSolverAMG");
}

bool
FieldSolverAMG::solve(MFAMRCellData&       a_potential,
                      const MFAMRCellData& a_rho,
                      const EBAMRIVData&   a_sigma,
                      const bool           a_zeroPhi)
{
  CH_TIME("FieldSolverAMG::solve");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::solve" << endl;
  }

  bool converged = false;

  return converged;
}

void
FieldSolverAMG::parseOptions()
{
  CH_TIME("FieldSolverAMG::parseOptions");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::parseOptions" << endl;
  }

  //  this->parseVerbosity();
  //  this->parsePlotVariables();
}

void
FieldSolverAMG::parseRuntimeOptions()
{
  CH_TIME("FieldSolverAMG::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::parseRuntimeOptions" << endl;
  }
}

void
FieldSolverAMG::computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential) const
{
  CH_TIME("FieldSolverAMG::computeElectricField(MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::computeElectricField(MFAMRCellData)" << endl;
  }
}

void
FieldSolverAMG::computeElectricField(MFAMRFluxData& a_E, const MFAMRCellData& a_potential) const
{
  CH_TIME("FieldSolverAMG::computeElectricField(MFAMRFluxData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::computeElectricField(MFAMRFluxData)" << endl;
  }
}

void
FieldSolverAMG::computeElectricField(EBAMRCellData&           a_E,
                                     const phase::which_phase a_phase,
                                     const MFAMRCellData&     a_potential) const
{
  CH_TIME("FieldSolverAMG::computeElectricField(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::computeElectricField(EBAMRCellData)" << endl;
  }
}

void
FieldSolverAMG::computeElectricField(EBAMRFluxData&           a_E,
                                     const phase::which_phase a_phase,
                                     const MFAMRCellData&     a_potential) const
{
  CH_TIME("FieldSolverAMG::computeElectricField(EBAMRFluxData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::computeElectricField(EBAMRFluxData)" << endl;
  }
}

void
FieldSolverAMG::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("FieldSolverAMG::preRegrid");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::preRegrid" << endl;
  }

  FieldSolver::preRegrid(a_lbase, a_oldFinestLevel);
}

void
FieldSolverAMG::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("FieldSolverAMG::regrid");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::regrid" << endl;
  }
}

void
FieldSolverAMG::registerOperators()
{
  CH_TIME("FieldSolverAMG::registerOperators");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::registerOperators" << endl;
  }

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
FieldSolverAMG::setupSolver()
{
  CH_TIME("FieldSolverAMG::setupSolver");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::setupSolver" << endl;
  }
}

void
FieldSolverAMG::setSolverPermittivities(const MFAMRCellData& a_permittivityCell,
                                        const MFAMRFluxData& a_permittivityFace,
                                        const MFAMRIVData&   a_permittivityEB)
{
  CH_TIME("FieldSolverAMG::setSolverPermittivities");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::setSolverPermittivities" << endl;
  }
}

void
FieldSolverAMG::setPermittivities()
{
  CH_TIME("FieldSolverAMG::setPermittivities");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::setPermittivities" << endl;
  }
}

#include <CD_NamespaceFooter.H>
