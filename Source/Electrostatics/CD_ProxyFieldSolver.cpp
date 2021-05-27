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
#include <CD_NwoEbQuadCfInterp.H>
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
  if(m_amr.isNull()){
    MayDay::Abort("FieldSolverMultigrid::registerOperators - need to set AmrMesh!");
  }
  else{
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, phase::gas);
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, phase::solid);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, phase::gas);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, phase::solid);
    m_amr->registerOperator(s_eb_pwl_interp,   m_realm, phase::gas);
    m_amr->registerOperator(s_eb_pwl_interp,   m_realm, phase::solid);
    m_amr->registerOperator(s_eb_quad_cfi,     m_realm, phase::gas);
    m_amr->registerOperator(s_eb_quad_cfi,     m_realm, phase::solid);
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, phase::gas);
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, phase::solid);
    m_amr->registerOperator(s_eb_flux_reg,     m_realm, phase::gas);
    m_amr->registerOperator(s_eb_flux_reg,     m_realm, phase::solid);
  }
}

Vector<EBLevelGrid> ProxyFieldSolver::getEBLevelGrids(){
  Vector<EBLevelGrid> ret;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    ret.push_back(*(m_amr->getEBLevelGrid(m_realm, phase::gas)[lvl]));
  }

  return ret;
}

Vector<RefCountedPtr<NWOEBQuadCFInterp> > ProxyFieldSolver::getInterp(){
  const int finestLevel = m_amr->getFinestLevel();
  
  Vector<RefCountedPtr<NWOEBQuadCFInterp> > ret(1 + finestLevel);

  for (int lvl = 1; lvl <= finestLevel; lvl++){
    if(lvl > 0){
      ret[lvl] = RefCountedPtr<NWOEBQuadCFInterp> (new NWOEBQuadCFInterp(m_amr->getGrids(m_realm)[lvl],
									 m_amr->getGrids(m_realm)[lvl-1],
									 m_amr->getEBISLayout(m_realm, phase::gas)[lvl],
									 m_amr->getEBISLayout(m_realm, phase::gas)[lvl-1],
									 m_amr->getDomains()[lvl-1],
									 m_amr->getRefinementRatios()[lvl-1],
									 1,
									 m_amr->getDx()[lvl],
									 m_amr->getNumberOfGhostCells()*IntVect::Unit,
									 *(m_amr->getEBLevelGrid(m_realm, phase::gas)[lvl]->getCFIVS()),
									 m_amr->getEBISLayout(m_realm, phase::gas)[0].getEBIS()));
									 
    }
  }

  
  return ret;
}

void ProxyFieldSolver::solveOnePhase(EBAMRCellData& a_phi){
  DataOps::setValue(a_phi, 1.23456789);


  // Define coefficients
  EBAMRCellData aco;
  EBAMRFluxData bco;
  EBAMRIVData   bcoIrreg;

  m_amr->allocate(aco,      m_realm, phase::gas, 1);
  m_amr->allocate(bco,      m_realm, phase::gas, 1);
  m_amr->allocate(bcoIrreg, m_realm, phase::gas, 1);

  DataOps::setValue(aco, 0.0);
  DataOps::setValue(bco, 1.0);
  DataOps::setValue(bcoIrreg, 1.0);

  const Real alpha =  0.0;
  const Real beta  = -1.0;

  auto levelGrids    = this->getEBLevelGrids();
  auto interpolators = this->getInterp();
}

#include <CD_NamespaceFooter.H>
