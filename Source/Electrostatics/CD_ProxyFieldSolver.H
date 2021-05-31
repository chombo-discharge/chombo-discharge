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
#include <NWOEBQuadCFInterp.H>
#include <EBMGAverage.H>

// Our includes
#include <CD_FieldSolver.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Proxy class for field solver. This class ignores the other phase and uses a single-phase solver for the "gas phase".
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

  void solveOnePhase(EBAMRCellData& a_solution, EBAMRCellData& a_residue);
  void coarsenMG(EBAMRCellData& a_solution);
  void coarsenConservative(EBAMRCellData& a_solution);
  
  Vector<EBLevelGrid> getEBLevelGrids();
  Vector<RefCountedPtr<NWOEBQuadCFInterp> > getInterpNWO();
  Vector<RefCountedPtr<EBQuadCFInterp> >    getInterpOld();
};

#include <CD_NamespaceFooter.H>


#endif
