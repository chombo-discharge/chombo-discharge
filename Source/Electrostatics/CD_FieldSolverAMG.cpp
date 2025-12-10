/* chombo-discharge
 * Copyright © 2026 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_FieldSolverAMG.cpp
  @brief  Implementation of CD_FieldSolver.H
  @author Robert Marskar
*/

#if CH_USE_PETSC

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_FieldSolverAMG.H>
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

  DataOps::setValue(a_potential, 0.0);

#if 1 // Check that PETSc works.
  const RefCountedPtr<PetscGrid>& petscGrid = m_amr->getPetscGrid(m_realm);

  Vec x;
  petscGrid->create(x);
  petscGrid->setValue(x, 1.0 * procID());
  petscGrid->putPetscInChombo(a_potential, x);
  petscGrid->destroy(x);
#endif

  return converged;
}

void
FieldSolverAMG::parseOptions()
{
  CH_TIME("FieldSolverAMG::parseOptions");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::parseOptions" << endl;
  }

  this->parseVerbosity();
  this->parsePlotVariables();
  this->parseRegridSlopes();
}

void
FieldSolverAMG::parseRuntimeOptions()
{
  CH_TIME("FieldSolverAMG::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::parseRuntimeOptions" << endl;
  }

  this->parseVerbosity();
  this->parsePlotVariables();
  this->parseRegridSlopes();
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
FieldSolverAMG::computeElectricField(MFAMRFluxData& a_electricField, const MFAMRCellData& a_potential) const
{
  CH_TIME("FieldSolverAMG::computeElectricField(MFAMRFluxData)");
  if (m_verbosity > 5) {
    pout() << "FieldSolverAMG::computeElectricField(MFAMRFluxData)" << endl;
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

  FieldSolver::regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
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

  if(!m_helmholtzPetsc->isDefined()){

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

#endif
