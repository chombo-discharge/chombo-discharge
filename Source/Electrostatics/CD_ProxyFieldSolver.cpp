/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ProxyFieldSolver.cpp
  @brief  Implementation of CD_ProxyFieldSolver.H
  @author Robert Marskar
*/


// Our includes
#include <CD_ProxyFieldSolver.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

ProxyFieldSolver::ProxyFieldSolver(){

}
ProxyFieldSolver::~ProxyFieldSolver(){

}

bool ProxyFieldSolver::solve(MFAMRCellData&       a_potential,
			     const MFAMRCellData& a_rho,
			     const EBAMRIVData&   a_sigma,
			     const bool           a_zerophi = false) override final{
  DataOps::setValue(a_potential, 0.0);
  
  return true;
}

void ProxyFieldSolver::registerOperators() override final {

}
