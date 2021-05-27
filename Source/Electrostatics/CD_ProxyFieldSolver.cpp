/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ProxyFieldSolver.cpp
  @brief  Implementation of CD_ProxyFieldSolver.H
  @author Robert Marskar
*/

// Chombo includes
#include <AMRMultiGrid.H>
#include <BiCGStabSolver.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <MFSimpleSolver.H>
#include <GMRESSolver.H>
#include <BaseBCFuncEval.H>

// Our includes
#include <CD_ProxyFieldSolver.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

bool ProxyFieldSolver::solve(MFAMRCellData&       a_potential,
			     const MFAMRCellData& a_rho,
			     const EBAMRIVData&   a_sigma,
			     const bool           a_zerophi) {
  DataOps::setValue(a_potential, 0.0);

  EBAMRCellData gasData = m_amr->alias(phase::gas, a_potential);

  this->solveOnePhase(gasData);
  
  return true;
}

void ProxyFieldSolver::registerOperators()  {

}

void ProxyFieldSolver::solveOnePhase(EBAMRCellData& a_phi){
  DataOps::setValue(a_phi, 1.23456789);
}

#include <CD_NamespaceFooter.H>
