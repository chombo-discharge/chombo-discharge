/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ProxyFieldSolver.H
  @brief  Declaration of a debugging class for Chombo multigrid codes. 
  @author Robert Marskar
*/

#ifndef CD_ProxyFieldSolver_H
#define CD_ProxyFieldSolver_H

// Chombo includes
#include <NWOEBQuadCFInterp.H>
#include <EBMGAverage.H>

// Our includes
#include <CD_EBMultigridInterpolator.H>
#include <CD_FieldSolver.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Proxy class for field solver. This class ignores the other phase and uses a single-phase solver for the "gas phase".
  @note  This class is not usable as a field solver - it is only meant for debugging and development of multigrid code/stencils. 
*/
class ProxyFieldSolver : public FieldSolver {
public:

  ProxyFieldSolver(){m_className = "ProxyFieldSolver";}
  ~ProxyFieldSolver(){}

  bool solve(MFAMRCellData&       a_potential,
	     const MFAMRCellData& a_rho,
	     const EBAMRIVData&   a_sigma,
	     const bool           a_zerophi = false) override final;

  void parseOptions() override final{
    parsePlotVariables();
  }
  
  void parseRuntimeOptions() override final{ }
  
  void registerOperators() override final;
  
protected:

  void solveEBCond   (EBAMRCellData& a_solution, EBAMRCellData& a_residue);
  void solveHelmholtz(EBAMRCellData& a_solution, EBAMRCellData& a_residue);
  void setupMF();
  
  Vector<EBLevelGrid> getEBLevelGrids();
  Vector<RefCountedPtr<EBQuadCFInterp> > getQuadCFI();
  Vector<RefCountedPtr<EBMultigridInterpolator> > getMultigridInterpolators(const phase::which_phase a_phase);
};

#include <CD_NamespaceFooter.H>


#endif
