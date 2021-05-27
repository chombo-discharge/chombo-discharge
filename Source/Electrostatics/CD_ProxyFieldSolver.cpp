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
#include <NWOEBConductivityOpFactory.H>
#include <EBConductivityOpFactory.H>
#include <DirichletConductivityDomainBC.H>
#include <DirichletConductivityEBBC.H>
#include <ParmParse.H>

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

Vector<RefCountedPtr<NWOEBQuadCFInterp> > ProxyFieldSolver::getInterpNWO(){
  const int finestLevel = m_amr->getFinestLevel();
  
  Vector<RefCountedPtr<NWOEBQuadCFInterp> > ret(1 + finestLevel);

  for (int lvl = 1; lvl <= finestLevel; lvl++){
    if(lvl > 0){
      ret[lvl] = RefCountedPtr<NWOEBQuadCFInterp> (new NwoEbQuadCfInterp(m_amr->getGrids(m_realm)[lvl],
									 m_amr->getGrids(m_realm)[lvl-1],
									 m_amr->getEBISLayout(m_realm, phase::gas)[lvl],
									 m_amr->getEBISLayout(m_realm, phase::gas)[lvl-1],
									 m_amr->getDomains()[lvl-1],
									 m_amr->getRefinementRatios()[lvl-1],
									 1,
									 m_amr->getDx()[lvl],
									 m_amr->getNumberOfGhostCells(),
									 *(m_amr->getEBLevelGrid(m_realm, phase::gas)[lvl]->getCFIVS()),
									 m_amr->getEBIndexSpace(phase::gas)));
									 
    }
  }

  
  return ret;
}

Vector<RefCountedPtr<EBQuadCFInterp> > ProxyFieldSolver::getInterpOld(){
  const int finestLevel = m_amr->getFinestLevel();
  
  Vector<RefCountedPtr<EBQuadCFInterp> > ret(1 + finestLevel);

  for (int lvl = 1; lvl <= finestLevel; lvl++){
    if(lvl > 0){
      ret[lvl] = RefCountedPtr<EBQuadCFInterp> (new EBQuadCFInterp(m_amr->getGrids(m_realm)[lvl],
								   m_amr->getGrids(m_realm)[lvl-1],
								   m_amr->getEBISLayout(m_realm, phase::gas)[lvl],
								   m_amr->getEBISLayout(m_realm, phase::gas)[lvl-1],
								   m_amr->getDomains()[lvl-1],
								   m_amr->getRefinementRatios()[lvl-1],
								   1,
								   *(m_amr->getEBLevelGrid(m_realm, phase::gas)[lvl]->getCFIVS()),
								   m_amr->getEBIndexSpace(phase::gas)));
									 
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
  EBAMRCellData rho;

  m_amr->allocate(aco,      m_realm, phase::gas, 1);
  m_amr->allocate(bco,      m_realm, phase::gas, 1);
  m_amr->allocate(bcoIrreg, m_realm, phase::gas, 1);
  m_amr->allocate(rho,      m_realm, phase::gas, 1);

  DataOps::setValue(aco, 1.0);
  DataOps::setValue(bco, 1.0);
  DataOps::setValue(bcoIrreg, 1.0);
  DataOps::setValue(rho, 0.0);

  const Real alpha =  0.0;
  const Real beta  = -1.0;

  auto levelGrids    = this->getEBLevelGrids();


  const IntVect ghostPhi = m_amr->getNumberOfGhostCells()*IntVect::Unit;
  const IntVect ghostRHS = m_amr->getNumberOfGhostCells()*IntVect::Unit;


  auto ebbcFactory = RefCountedPtr<DirichletConductivityEBBCFactory> (new DirichletConductivityEBBCFactory());
  ebbcFactory->setValue(1.0);
  ebbcFactory->setOrder(1);

  auto domainFactory = RefCountedPtr<DirichletConductivityDomainBCFactory> (new DirichletConductivityDomainBCFactory());
  domainFactory->setValue(0.0);


  //  auto interpolators = this->getInterpNWO();
  auto factoryNWO = RefCountedPtr<NWOEBConductivityOpFactory> (new NWOEBConductivityOpFactory(levelGrids,
											      this->getInterpNWO(),
											      alpha,
											      beta,
											      aco.getData(),
											      bco.getData(),
											      bcoIrreg.getData(),
											      m_amr->getDx()[0],
											      m_amr->getRefinementRatios(),
											      domainFactory,
											      ebbcFactory,
											      ghostPhi,
											      ghostRHS,
											      2));

  auto factoryOld = RefCountedPtr<EBConductivityOpFactory> (new EBConductivityOpFactory(levelGrids,
											this->getInterpOld(),
											alpha,
											beta,
											aco.getData(),
											bco.getData(),
											bcoIrreg.getData(),
											m_amr->getDx()[0],
											m_amr->getRefinementRatios(),
											domainFactory,
											ebbcFactory,
											ghostPhi,
											ghostRHS,
											2));

  
  BiCGStabSolver<LevelData<EBCellFAB> > bicgstab;
  AMRMultiGrid<LevelData<EBCellFAB> >   multigridSolver;

  ParmParse pp(m_className.c_str());
  bool useNWO;
  pp.get("use_nwo", useNWO);
  if(useNWO){
    multigridSolver.define(m_amr->getDomains()[0], *factoryNWO, &bicgstab, 1 + m_amr->getFinestLevel());
  }
  else{
    multigridSolver.define(m_amr->getDomains()[0], *factoryOld, &bicgstab, 1 + m_amr->getFinestLevel());
  }
  multigridSolver.setSolverParameters(16, 16, 16, 1, 32, 1E-30, 1E-30, 1E-60);


  // Solve
  Vector<LevelData<EBCellFAB>* > phi;
  Vector<LevelData<EBCellFAB>* > rhs;

  m_amr->alias(phi, a_phi);
  m_amr->alias(rhs, rho);

  multigridSolver.init( phi, rhs, m_amr->getFinestLevel(), 0);
  multigridSolver.solveNoInit(phi, rhs, m_amr->getFinestLevel(), 0);
  multigridSolver.m_verbosity = 10;

  m_amr->averageDown(a_phi, m_realm, phase::gas);
  m_amr->interpGhost(a_phi, m_realm, phase::gas);
  
}

#include <CD_NamespaceFooter.H>
