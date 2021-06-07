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
#include <EBConductivityOpFactory.H>
#include <NWOEBConductivityOpFactory.H>
#include <slowEBCOFactory.H>
#include <EBAMRPoissonOpFactory.H>
#include <DirichletConductivityDomainBC.H>
#include <DirichletConductivityEBBC.H>
#include <DirichletPoissonDomainBC.H>
#include <DirichletPoissonEBBC.H>
#include <NeumannConductivityDomainBC.H>
#include <NeumannConductivityEBBC.H>
#include <EBSimpleSolver.H>
#include <ParmParse.H>

// Our includes
#include <CD_ProxyFieldSolver.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_NwoEbQuadCfInterp.H>
#include <CD_DataOps.H>
#include <CD_EBHelmholtzOpFactory.H>
#include <CD_NamespaceHeader.H>


bool ProxyFieldSolver::solve(MFAMRCellData&       a_potential,
			     const MFAMRCellData& a_rho,
			     const EBAMRIVData&   a_sigma,
			     const bool           a_zerophi) {
  DataOps::setValue(a_potential, 0.0);
  DataOps::setValue(m_residue,   0.0);

  EBAMRCellData gasData = m_amr->alias(phase::gas, a_potential);
  EBAMRCellData gasResi = m_amr->alias(phase::gas, m_residue);

  this->setupHelmholtz(); // Just for testing, right now. 

  this->solveOnePhase(gasData, gasResi);
  
  return true;
}

void ProxyFieldSolver::registerOperators()  {
  if(m_amr.isNull()){
    MayDay::Abort("FieldSolverMultigrid::registerOperators - need to set AmrMesh!");
  }
  else{
    // For regridding
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, phase::gas);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, phase::solid);

    // For coarsening
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, phase::gas);
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, phase::solid);

    // For linearly filling ghost cells
    m_amr->registerOperator(s_eb_pwl_interp,   m_realm, phase::gas);
    m_amr->registerOperator(s_eb_pwl_interp,   m_realm, phase::solid);

    // For interpolating to cell centroids
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, phase::gas);
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, phase::solid);

    // For making helmholtz happy
    m_amr->registerOperator(s_eb_flux_reg, m_realm, phase::gas);
    m_amr->registerOperator(s_eb_flux_reg, m_realm, phase::solid);
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

  const int nGhosts = m_amr->getNumberOfGhostCells();
  
  for (int lvl = 1; lvl <= finestLevel; lvl++){
    if(lvl > 0){

      // Make the CFIVS
      const DisjointBoxLayout& fineGrids = m_amr->getGrids(m_realm)[lvl];
      LayoutData<IntVectSet> ghostCells(fineGrids);
      for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit){

	Box grownBox      = grow(fineGrids[dit()], nGhosts);
	grownBox         &= m_amr->getDomains()[lvl];
	ghostCells[dit()] = IntVectSet(grownBox);

	const Vector<LayoutIndex>& neighbors = (*m_amr->getNeighbors(m_realm, phase::gas)[lvl])[dit()];
	for (int i = 0; i < neighbors.size(); i++){
	  ghostCells[dit()] -= fineGrids[neighbors[i]];
	}
	ghostCells[dit()] -= fineGrids[dit()];
      }

      // Define interpolator. 
      ret[lvl] = RefCountedPtr<NWOEBQuadCFInterp> (new NwoEbQuadCfInterp(m_amr->getGrids(m_realm)[lvl],
									 m_amr->getGrids(m_realm)[lvl-1],
									 m_amr->getEBISLayout(m_realm, phase::gas)[lvl],
									 m_amr->getEBISLayout(m_realm, phase::gas)[lvl-1],
									 m_amr->getDomains()[lvl-1],
									 m_amr->getRefinementRatios()[lvl-1],
									 1,
									 m_amr->getDx()[lvl],
									 nGhosts*IntVect::Unit,
									 ghostCells,
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

void ProxyFieldSolver::solveOnePhase(EBAMRCellData& a_phi, EBAMRCellData& a_residue){
  ParmParse pp(m_className.c_str());

  // Define coefficients
  EBAMRCellData rho;
  EBAMRCellData aco;
  EBAMRCellData zero;
  EBAMRFluxData bco;
  EBAMRIVData   bcoIrreg;


  m_amr->allocate(aco,      m_realm, phase::gas, 1);
  m_amr->allocate(bco,      m_realm, phase::gas, 1);
  m_amr->allocate(bcoIrreg, m_realm, phase::gas, 1);
  m_amr->allocate(rho,      m_realm, phase::gas, 1);
  m_amr->allocate(zero,     m_realm, phase::gas, 1);  

  DataOps::setValue(zero,     0.0);
  DataOps::setValue(aco,      0.0);
  DataOps::setValue(bco,      1.0);
  DataOps::setValue(bcoIrreg, 1.0);
  DataOps::setValue(rho,      0.0);

  const Real alpha =  0.;
  const Real beta  =  1.;

  auto levelGrids    = this->getEBLevelGrids();
  auto interpNWO     = this->getInterpNWO();
  auto interpOld     = this->getInterpOld();

  const IntVect ghostPhi = m_amr->getNumberOfGhostCells()*IntVect::Unit;
  const IntVect ghostRHS = m_amr->getNumberOfGhostCells()*IntVect::Unit;

  int         relaxType;
  int         eb_order;
  Real        eb_value;
  Real        dom_value;
  Real        phi_init;
  std::string str;

  pp.get("relax",    relaxType);
  pp.get("eb_order", eb_order);
  pp.get("eb_val",   eb_value);
  pp.get("dom_val",  dom_value);
  pp.get("phi_ini",  phi_init);

  DataOps::setValue(a_phi, phi_init);

  // BC factories for conductivity ops
  auto ebbcFactory   = RefCountedPtr<DirichletConductivityEBBCFactory>     (new DirichletConductivityEBBCFactory());
  auto domainFactory = RefCountedPtr<DirichletConductivityDomainBCFactory> (new DirichletConductivityDomainBCFactory());
  ebbcFactory  ->setValue(eb_value);
  ebbcFactory  ->setOrder(eb_order);
  domainFactory->setValue(dom_value);
  

  // BC factories for EBAMRPoissonOp
  auto poissonEBFactory  = RefCountedPtr<DirichletPoissonEBBCFactory>     (new DirichletPoissonEBBCFactory());
  auto poissonDomFactory = RefCountedPtr<DirichletPoissonDomainBCFactory> (new DirichletPoissonDomainBCFactory());
  poissonEBFactory ->setValue(eb_value);
  poissonEBFactory ->setOrder(eb_order);
  poissonDomFactory->setValue(dom_value);


  auto factoryNWO = RefCountedPtr<NWOEBConductivityOpFactory> (new NWOEBConductivityOpFactory(levelGrids,
											      interpNWO,
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
											      relaxType));

  auto factoryOld = RefCountedPtr<EBConductivityOpFactory> (new EBConductivityOpFactory(levelGrids,
											interpOld,
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
											relaxType));

  
  auto factoryPoiss = RefCountedPtr<EBAMRPoissonOpFactory> (new EBAMRPoissonOpFactory(levelGrids,
										      m_amr->getRefinementRatios(),
										      interpOld,
										      m_amr->getDx()[0]*RealVect::Unit,
										      m_amr->getProbLo(),
										      40,
										      relaxType,
										      poissonDomFactory,
										      poissonEBFactory,
										      alpha,
										      beta,
										      0.0,
										      ghostPhi,
										      ghostRHS));


  auto factorySlow = RefCountedPtr<slowEBCOFactory> (new slowEBCOFactory(levelGrids,
									 interpOld,
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
									 relaxType));

										   

  BiCGStabSolver<LevelData<EBCellFAB> > bicgstab;
  AMRMultiGrid<LevelData<EBCellFAB> >   multigridSolver;

  bool useCond;
  bool useNWO;
  int  numSmooth;
  Real tolerance;
  std::string whichFactory;

  pp.get("mg_tol",   tolerance);
  pp.get("smooth",   numSmooth);
  pp.get("factory",  whichFactory);

  if(whichFactory == "nwo"){
    pout() << "using nwo ebconductivityop" << endl;
    multigridSolver.define(m_amr->getDomains()[0], *factoryNWO, &bicgstab, 1 + m_amr->getFinestLevel());
  }
  else if(whichFactory == "cond"){
    pout() << "using old ebconductivityop" << endl;
    multigridSolver.define(m_amr->getDomains()[0], *factoryOld, &bicgstab, 1 + m_amr->getFinestLevel());
  }
  else if(whichFactory == "poiss"){
    pout() << "using ebamrpoissonopfactory" << endl;
    multigridSolver.define(m_amr->getDomains()[0], *factoryPoiss, &bicgstab, 1 + m_amr->getFinestLevel());
  }
  else if(whichFactory == "slowco"){
    pout() << "using slowconductivity" << endl;
    Chombo_EBIS::alias(&(*m_amr->getEBIndexSpace(phase::gas)));
    multigridSolver.define(m_amr->getDomains()[0], *factorySlow, &bicgstab, 1 + m_amr->getFinestLevel());
  }
  else{
    MayDay::Abort("Bad argument to 'factory'");
  }

  multigridSolver.setSolverParameters(numSmooth, numSmooth, numSmooth, 1, 32, tolerance, 1E-60, 1E-60);


  // Solve
  Vector<LevelData<EBCellFAB>* > phi;
  Vector<LevelData<EBCellFAB>* > rhs;
  Vector<LevelData<EBCellFAB>* > res;
  Vector<LevelData<EBCellFAB>* > zer;

  m_amr->alias(phi, a_phi);
  m_amr->alias(rhs, rho);
  m_amr->alias(res, a_residue);
  m_amr->alias(zer, zero);

  const int finestLevel = m_amr->getFinestLevel();
  const int baseLevel   = 0;

  multigridSolver.m_verbosity = 10;
  multigridSolver.init(phi, rhs, finestLevel, baseLevel);
  multigridSolver.m_convergenceMetric = multigridSolver.computeAMRResidual(zer, rhs, finestLevel, baseLevel);

  Real zerResid = multigridSolver.computeAMRResidual(res, zer, rhs, finestLevel, baseLevel);
  Real phiResid = multigridSolver.computeAMRResidual(res, phi, rhs, finestLevel, baseLevel);

  int iter = 0;
  while(phiResid >= zerResid*tolerance && iter < 10 ){
    multigridSolver.m_convergenceMetric = multigridSolver.computeAMRResidual(zer, rhs, finestLevel, baseLevel);

    multigridSolver.solveNoInit(phi, rhs, finestLevel, baseLevel, false, false);
    coarsenConservative(a_phi);
    phiResid = multigridSolver.computeAMRResidual(res, phi, rhs, finestLevel, baseLevel);

    iter++;
  }
  
  this->computeElectricField();
  this->writePlotFile();
}

void ProxyFieldSolver::coarsenMG(EBAMRCellData& a_phi){
  for (int lvl = m_amr->getFinestLevel(); lvl > 0; lvl--){
    EBMGAverage aveOp(m_amr->getGrids(m_realm)[lvl],
		      m_amr->getGrids(m_realm)[lvl-1],
		      m_amr->getEBISLayout(m_realm, phase::gas)[lvl],
		      m_amr->getEBISLayout(m_realm, phase::gas)[lvl-1],
		      m_amr->getDomains()[lvl-1],
		      m_amr->getRefinementRatios()[lvl-1],
		      1,
		      m_amr->getEBIndexSpace(phase::gas),
		      m_amr->getNumberOfGhostCells()*IntVect::Unit,
		      true);
	  
    aveOp.average(*a_phi[lvl-1], *a_phi[lvl], Interval(0,0));
  }
}

void ProxyFieldSolver::coarsenConservative(EBAMRCellData& a_phi){
  m_amr->averageDown(a_phi, m_realm, phase::gas);
}

void ProxyFieldSolver::setupHelmholtz(){

  const int finestLevel = m_amr->getFinestLevel();
  
  Vector<RefCountedPtr<EBMultigridInterpolator> > interpolators(1+finestLevel);
  for (int lvl = 0; lvl <= finestLevel; lvl++){

    if(lvl > 0){
      const EBLevelGrid eblgFine = *m_amr->getEBLevelGrid(m_realm, phase::gas)[lvl];
      const EBLevelGrid eblgCoar = *m_amr->getEBLevelGrid(m_realm, phase::gas)[lvl-1];
      const LayoutData<IntVectSet>& cfivs = *eblgFine.getCFIVS();
      
      interpolators[lvl] = RefCountedPtr<EBMultigridInterpolator> (new EBMultigridInterpolator(eblgFine,
											       eblgCoar,
											       m_amr->getRefinementRatio(lvl, lvl-1),
											       1,
											       cfivs));
    }
  }

  const Real alpha = 1.0;
  const Real beta  = 1.0;

  EBAMRCellData Aco;
  EBAMRFluxData Bco;
  EBAMRIVData   BcoIrreg;



  m_amr->allocate(Aco,      m_realm, phase::gas, 1);
  m_amr->allocate(Bco,      m_realm, phase::gas, 1);
  m_amr->allocate(BcoIrreg, m_realm, phase::gas, 1);

  auto ebbcFactory   = RefCountedPtr<DirichletConductivityEBBCFactory>     (new DirichletConductivityEBBCFactory());
  auto domainFactory = RefCountedPtr<DirichletConductivityDomainBCFactory> (new DirichletConductivityDomainBCFactory());
  ebbcFactory  ->setValue(1);
  ebbcFactory  ->setOrder(1);
  domainFactory->setValue(-1);

  // Set the bottom domain. Don't go below 8x cells in any direction
  ProblemDomain bottomDomain = m_amr->getDomains()[0];
  while(bottomDomain.domainBox().longside() >= 16){
    bottomDomain.coarsen(2);
  }

  EBHelmholtzOpFactory fact(alpha,
			    beta,
			    m_amr->getEBLevelGrid(m_realm, phase::gas),
			    interpolators,
			    m_amr->getFluxRegister(m_realm, phase::gas),
			    m_amr->getCoarseAverage(m_realm, phase::gas),
			    m_amr->getRefinementRatios(),
			    m_amr->getDx(),
			    Aco.getData(),
			    Bco.getData(),
			    BcoIrreg.getData(),
			    domainFactory,
			    ebbcFactory,
			    m_amr->getNumberOfGhostCells()*IntVect::Unit,
			    m_amr->getNumberOfGhostCells()*IntVect::Unit,
			    EBHelmholtzOp::RelaxationMethod::Jacobi,
			    bottomDomain,
			    m_amr->getBlockingFactor());

  BiCGStabSolver<LevelData<EBCellFAB> > bicgstab;
  AMRMultiGrid<LevelData<EBCellFAB> > multigridSolver;
  
  multigridSolver.define(m_amr->getDomains()[0], fact, &bicgstab, 1 + m_amr->getFinestLevel());

}

#include <CD_NamespaceFooter.H>
