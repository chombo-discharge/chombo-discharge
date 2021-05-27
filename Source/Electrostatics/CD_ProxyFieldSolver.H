/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ProxyFieldSolver.H
  @brief  Development proxy for new multigrid solver
  @author Robert Marskar
*/

#ifndef CD_ProxyFieldSolver_H
#define CD_ProxyFieldSolver_H

// Chombo includes
#include <AMRMultiGrid.H>
#include <BiCGStabSolver.H>
#include <BaseDomainBC.H>
#include <BaseEBBC.H>
#include <MFSimpleSolver.H>
#include <GMRESSolver.H>
#include <BaseBCFuncEval.H>

// Our includes
#include <CD_FieldSolver.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Proxy class for field solver. This class ignores the other phase and uses a single-phase solver for the "gas phase".
*/
class ProxyFieldSolver : public FieldSolver {
public:

  ProxyFieldSolver();
  ~ProxyFieldSolver();

  bool solve(MFAMRCellData&       a_potential,
	     const MFAMRCellData& a_rho,
	     const EBAMRIVData&   a_sigma,
	     const bool           a_zerophi = false) override final;

  void parseOptions() override final{}
  
  void parseRuntimeOptions() override final{ }
  void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override final{
    MayDay::Abort("ProxyFieldSolver::regrid - should never be called");
  }
  
  void registerOperators() override final;
  
protected:


  
};


#endif
