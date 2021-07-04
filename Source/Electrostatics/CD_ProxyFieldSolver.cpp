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
#include <MFSimpleSolver.H>
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
#include <CD_EBHelmholtzOpFactory.H>
#include <CD_EBHelmholtzNeumannEBBCFactory.H>
#include <CD_EBHelmholtzDirichletEBBCFactory.H>
#include <CD_EBHelmholtzRobinEBBCFactory.H>
#include <CD_EBHelmholtzDirichletDomainBCFactory.H>
#include <CD_EBHelmholtzNeumannDomainBCFactory.H>
#include <CD_EBHelmholtzRobinDomainBCFactory.H>

// Includes for MFHelmholtzOpFactory
#include <CD_MFHelmholtzDirichletEBBCFactory.H>
#include <CD_MFHelmholtzNeumannEBBCFactory.H>
#include <CD_MFHelmholtzRobinEBBCFactory.H>
#include <CD_MFHelmholtzDirichletDomainBCFactory.H>
#include <CD_MFHelmholtzNeumannDomainBCFactory.H>
#include <CD_MFHelmholtzRobinDomainBCFactory.H>
#include <CD_DataOps.H>
#include <CD_MFHelmholtzOpFactory.H>
#include <CD_MFQuadCFInterp.H>
#include <CD_MFLevelGrid.H>
#include <CD_MFFluxReg.H>
#include <CD_MFCoarAve.H>
#include <CD_Units.H>

// Old news
#include <CD_RobinConductivityDomainBcFactory.H>
#include <CD_RobinConductivityEbBcFactory.H>
#include <CD_NamespaceHeader.H>


bool ProxyFieldSolver::solve(MFAMRCellData&       a_potential,
			     const MFAMRCellData& a_rho,
			     const EBAMRIVData&   a_sigma,
			     const bool           a_zerophi) {
  if(m_verbosity > 3){
    pout() << "ProxyFieldSolver::solve(...)" << endl;
  }
  
  DataOps::setValue(a_potential, 0.0);
  DataOps::setValue(m_residue,   0.0);

  EBAMRCellData gasData = m_amr->alias(phase::gas, a_potential);
  EBAMRCellData gasResi = m_amr->alias(phase::gas, m_residue);

  ParmParse pp(m_className.c_str());
  std::string str;
  pp.get("solver", str);

  Real t1 = -MPI_Wtime();
  if(str == "ebcond"){
    this->solveEBCond(   gasData, gasResi);
  }
  else if(str == "helm"){
    this->solveHelmholtz(gasData, gasResi);
  }
  else if(str == "mfhelm"){
    this->solveMF(a_potential, a_rho, a_sigma, a_zerophi);
  }
  else{
    MayDay::Error("ProxyFieldSolver::solve - bad argument to ProxyFieldSolver.solver");
  }
  t1 += MPI_Wtime();
  if(procID() == 0) std::cout << "solve time = " << t1 << "\n";

  m_amr->averageDown(a_potential, m_realm);
  m_amr->interpGhost(a_potential, m_realm);
  
  this->computeElectricField();
  
  this->writePlotFile();
  
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

    // For multigrid interpolation
    m_amr->registerOperator(s_eb_quad_cfi,     m_realm, phase::gas);
    m_amr->registerOperator(s_eb_quad_cfi,     m_realm, phase::solid);

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

Vector<RefCountedPtr<EBQuadCFInterp> > ProxyFieldSolver::getQuadCFI(){
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

Vector<RefCountedPtr<EBMultigridInterpolator> > ProxyFieldSolver::getMultigridInterpolators(const phase::which_phase a_phase){
  const int finestLevel = m_amr->getFinestLevel();

  ParmParse pp("ProxyFieldSolver");

  int ghostCF;
  int weight;
  int order;

  pp.get("interp_weight", weight);
  pp.get("interp_num", ghostCF);
  pp.get("interp_order", order);

  Vector<RefCountedPtr<EBMultigridInterpolator> > interpolators(1+finestLevel);
  
  for (int lvl = 0; lvl <= finestLevel; lvl++){

    if(lvl > 0){
      const EBLevelGrid eblgFine = *m_amr->getEBLevelGrid(m_realm, a_phase)[lvl];
      const EBLevelGrid eblgCoar = *m_amr->getEBLevelGrid(m_realm, a_phase)[lvl-1];
      const LayoutData<IntVectSet>& cfivs = *eblgFine.getCFIVS();
      
      interpolators[lvl] = RefCountedPtr<EBMultigridInterpolator> (new EBMultigridInterpolator(eblgFine,
											       eblgCoar,
											       Location::Cell::Center,
											       m_amr->getNumberOfGhostCells()*IntVect::Unit,
											       m_amr->getRefinementRatios()[lvl-1],
											       1, // Variables
											       ghostCF, // # of ghost cells to fill
											       order, // Order
											       weight)); // Weight
    }
  }

  return interpolators;
}

void ProxyFieldSolver::solveEBCond(EBAMRCellData& a_phi, EBAMRCellData& a_residue){
  ParmParse pp(m_className.c_str());
  
  Real alpha;
  Real beta;

  Real Acoef;
  Real Bcoef;

  pp.get("alpha", alpha);
  pp.get("beta",  beta);
  pp.get("aco",   Acoef);
  pp.get("bco",   Bcoef);

  // Define coefficients
  EBAMRCellData zero;
  EBAMRCellData aco;
  EBAMRFluxData bco;
  EBAMRIVData   bcoIrreg;
  EBAMRCellData rho;

  m_amr->allocate(aco,      m_realm, phase::gas, 1);
  m_amr->allocate(bco,      m_realm, phase::gas, 1);
  m_amr->allocate(bcoIrreg, m_realm, phase::gas, 1);
  m_amr->allocate(rho,      m_realm, phase::gas, 1);
  m_amr->allocate(zero,     m_realm, phase::gas, 1);  

  DataOps::setValue(zero,     0.0);
  DataOps::setValue(aco,      Acoef);
  DataOps::setValue(bco,      Bcoef);
  DataOps::setValue(bcoIrreg, Bcoef);
  DataOps::setValue(rho,      0.0);

  const IntVect ghostPhi = m_amr->getNumberOfGhostCells()*IntVect::Unit;
  const IntVect ghostRHS = m_amr->getNumberOfGhostCells()*IntVect::Unit;

  int         eb_order;
  Real        eb_value;
  Real        dom_value;
  std::string str;

  // Set up boundary conditions. 
  RefCountedPtr<BaseEBBCFactory>     ebbcFactory;
  RefCountedPtr<BaseDomainBCFactory> domainBcFactory;

  // EBBC
  pp.get("eb_bc", str);
  if(str == "dirichlet"){
    pp.get("eb_order", eb_order);
    pp.get("eb_val",   eb_value);
    
    auto bcFactory = new DirichletConductivityEBBCFactory();
    bcFactory->setValue(eb_value);
    bcFactory->setOrder(eb_order);

    ebbcFactory = RefCountedPtr<BaseEBBCFactory> (bcFactory);

    if(eb_order == 2 && m_amr->getNumberOfGhostCells() < 3)  MayDay::Abort("ProxyFieldSolver::solveEBCond - not enough ghost cells!");
  }
  else if(str == "neumann"){
    pp.get("eb_val",   eb_value);
    
    auto bcFactory = new NeumannConductivityEBBCFactory();
    bcFactory->setValue(eb_value);

    ebbcFactory = RefCountedPtr<BaseEBBCFactory>(bcFactory);
  }
  else if (str == "robin"){
    pp.get("eb_val", eb_value);

    // Coefficients for radiative transfer with Robin. 
    const Real A =  1.5*eb_value;
    const Real B = -1.0*eb_value;
    const Real C =  0.0;

    auto bcFactory = new RobinConductivityEbBcFactory(m_amr->getProbLo());
    bcFactory->setCoefficients(A, B, C);
    
    ebbcFactory = RefCountedPtr<BaseEBBCFactory> (bcFactory);
  }
  else{
    MayDay::Error("ProxyFieldSolver::solveEBCond - uknown EBBC factory requested");
  }

  // Domain BC
  pp.get("domain_bc", str);
  pp.get("dom_val", dom_value);  
  if(str == "dirichlet"){
    auto bcFactory = new DirichletConductivityDomainBCFactory();
    bcFactory->setValue(dom_value);

    domainBcFactory = RefCountedPtr<BaseDomainBCFactory>(bcFactory);
  }
  else if(str == "neumann"){
    auto bcFactory = new NeumannConductivityDomainBCFactory();
    bcFactory->setValue(dom_value);

    domainBcFactory = RefCountedPtr<BaseDomainBCFactory>(bcFactory);
  }
  else if(str == "robin"){
    auto bcFactory = new RobinConductivityDomainBcFactory();

    // Coeffs for radiative transfer with Robin
    const Real A =  1.5*dom_value;
    const Real B = -1.0*dom_value;
    const Real C = 0.0;

    bcFactory->setCoefficients(A, B, C);

    domainBcFactory = RefCountedPtr<BaseDomainBCFactory>(bcFactory);
  }
  else{
    MayDay::Error("ProxyFieldSolver::solveEBCond - uknown domain factory requested");
  }

  int relaxType;
  pp.get("relax",relaxType);

  BiCGStabSolver<LevelData<EBCellFAB> > bicgstab;
  AMRMultiGrid<LevelData<EBCellFAB> >   multigridSolver;  
  EBConductivityOpFactory condFactory(this->getEBLevelGrids(),
				      this->getQuadCFI(),
				      alpha,
				      beta,
				      aco.getData(),
				      bco.getData(),
				      bcoIrreg.getData(),
				      m_amr->getDx()[0],
				      m_amr->getRefinementRatios(),
				      domainBcFactory,
				      ebbcFactory,
				      ghostPhi,
				      ghostRHS,
				      relaxType);
										   



  int  numSmooth;
  Real tolerance;

  pp.get("mg_tol",   tolerance);
  pp.get("smooth",   numSmooth);

  const int finestLevel = m_amr->getFinestLevel();
  const int baseLevel   = 0;

  // Define
  multigridSolver.define(m_amr->getDomains()[0], condFactory, &bicgstab, 1 + finestLevel);
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

  //  multigridSolver.m_verbosity = 10;
  multigridSolver.init(phi, rhs, finestLevel, baseLevel);
  multigridSolver.m_convergenceMetric = multigridSolver.computeAMRResidual(zer, rhs, finestLevel, baseLevel);
  multigridSolver.solveNoInit(phi, rhs, finestLevel, baseLevel, false);
}

void ProxyFieldSolver::solveHelmholtz(EBAMRCellData& a_phi, EBAMRCellData& a_residue){
  ParmParse pp(m_className.c_str());
  
  Real alpha;
  Real beta;

  Real Acoeff;
  Real Bcoeff;

  pp.get("alpha", alpha);
  pp.get("beta",  beta);
  pp.get("aco",   Acoeff);
  pp.get("bco",   Bcoeff);

  EBAMRCellData zero;
  EBAMRCellData Aco;
  EBAMRFluxData Bco;
  EBAMRIVData   BcoIrreg;
  EBAMRCellData rho;

  m_amr->allocate(zero,     m_realm, phase::gas, 1);
  m_amr->allocate(Aco,      m_realm, phase::gas, 1);
  m_amr->allocate(Bco,      m_realm, phase::gas, 1);
  m_amr->allocate(BcoIrreg, m_realm, phase::gas, 1);
  m_amr->allocate(rho,      m_realm, phase::gas, 1);

  DataOps::setValue(zero,     0.0);
  DataOps::setValue(Aco,      Acoeff);
  DataOps::setValue(Bco,      Bcoeff);
  DataOps::setValue(BcoIrreg, Bcoeff);
  DataOps::setValue(rho,      0.0);

  const IntVect ghostPhi = m_amr->getNumberOfGhostCells()*IntVect::Unit;
  const IntVect ghostRhs = m_amr->getNumberOfGhostCells()*IntVect::Unit;

  int         eb_order;
  int         eb_weight;
  Real        eb_value;
  Real        dom_value;
  std::string str;

  RefCountedPtr<EBHelmholtzEBBCFactory> ebbcFactory;
  RefCountedPtr<EBHelmholtzDomainBCFactory> domainBcFactory;

  // EBBC
  pp.get("eb_bc", str);
  if(str == "dirichlet"){
    pp.get("eb_order",  eb_order);
    pp.get("eb_val",    eb_value);
    pp.get("eb_weight", eb_weight);

    ebbcFactory = RefCountedPtr<EBHelmholtzEBBCFactory> (new EBHelmholtzDirichletEBBCFactory(eb_order, eb_weight, eb_value));

    // Make sure we have enough ghost cells.
    if(eb_order > m_amr->getNumberOfGhostCells()) MayDay::Abort("ProxyFieldSolver::solveHelm - not enough ghost cells!");
  }
  else if(str == "neumann"){
    pp.get("eb_val", eb_value);

    ebbcFactory = RefCountedPtr<EBHelmholtzEBBCFactory> (new EBHelmholtzNeumannEBBCFactory(eb_value));    
  }
  else if(str == "robin"){
    pp.get("eb_val", eb_value);

    // Coefficients for radiative transfer with Robin. 
    const Real A =  1.5*eb_value;
    const Real B = -1.0*eb_value;
    const Real C =  0.0;
    ebbcFactory = RefCountedPtr<EBHelmholtzEBBCFactory> (new EBHelmholtzRobinEBBCFactory(A, B, C));
  }
  else{
    MayDay::Error("ProxyFieldSolver::solveEBCond - uknown EBBC factory requested");
  }

  // Domain bc
  pp.get("domain_bc", str);
  pp.get("dom_val", dom_value);

  // BC function. Spatially dependent. 
  auto bcFunction = [dom_value](const RealVect& a_pos) -> Real {
    return dom_value;
  };
  
  if(str == "dirichlet"){
    domainBcFactory = RefCountedPtr<EBHelmholtzDomainBCFactory>(new EBHelmholtzDirichletDomainBCFactory(bcFunction));
  }
  else if(str == "neumann"){
    domainBcFactory = RefCountedPtr<EBHelmholtzDomainBCFactory>(new EBHelmholtzNeumannDomainBCFactory(bcFunction));
  }
  else if(str == "robin"){
    
    // Coeffs for radiative transfer with Robin
    const Real A =  1.5*dom_value;
    const Real B = -1.0*dom_value;
    const Real C =  0.0;
    
    domainBcFactory = RefCountedPtr<EBHelmholtzDomainBCFactory>(new EBHelmholtzRobinDomainBCFactory(A, B, C));
  }
  else{
    MayDay::Abort("ProxyFieldSolver::solveHelm - unknown domain bc requested");
  }


  // Set the bottom domain. Don't go below 8x cells in any direction
  int minCells;
  pp.get("min_cells", minCells);
  ProblemDomain bottomDomain = m_amr->getDomains()[0];
  while(bottomDomain.domainBox().shortside() >= 2*minCells){
    bottomDomain.coarsen(2);
  }

  // Relax stuff consistent with EBConductivityOp (with exception of GauSaiRedBlack...)
  EBHelmholtzOp::RelaxationMethod relaxType;
  int relax;
  pp.get("relax",relax);
  switch(relax){
  case 0:
    relaxType = EBHelmholtzOp::RelaxationMethod::PointJacobi;
    break;
  case 1:
    relaxType = EBHelmholtzOp::RelaxationMethod::GauSaiMultiColor;
    break;
  case 2:
    relaxType = EBHelmholtzOp::RelaxationMethod::GauSaiRedBlack;
    break;
  default:
    MayDay::Abort("ProxyFieldSolver::solveHelmholtz - unknown relaxation method requested");
    break;
  }


  EBHelmholtzOpFactory fact(alpha,
			    beta,
			    m_amr->getProbLo(),
			    m_amr->getEBLevelGrid(m_realm, phase::gas),
			    this->getMultigridInterpolators(phase::gas),
			    m_amr->getFluxRegister(m_realm, phase::gas),
			    m_amr->getCoarseAverage(m_realm, phase::gas),
			    m_amr->getRefinementRatios(),
			    m_amr->getDx(),
			    Aco.getData(),
			    Bco.getData(),
			    BcoIrreg.getData(),
			    domainBcFactory,
			    ebbcFactory,
			    ghostPhi,
			    ghostRhs,
			    relaxType,
			    bottomDomain,
			    m_amr->getMaxBoxSize());

  BiCGStabSolver<LevelData<EBCellFAB> > bicgstab;
  EBSimpleSolver simpleSolver;
  AMRMultiGrid<LevelData<EBCellFAB> > multigridSolver;

  int  numSmooth;
  int  botSmooth;
  Real tolerance;

  pp.get("mg_tol",   tolerance);
  pp.get("smooth",   numSmooth);
  pp.get("botsolver", str);
  
  pp.get("botsmooth", botSmooth);

  simpleSolver.setNumSmooths(botSmooth);
  
  const int baseLevel   = 0;
  const int finestLevel = m_amr->getFinestLevel();

  // Define
  //  bicgstab.m_verbosity=10;
  if(str == "bicgstab"){
    multigridSolver.define(m_amr->getDomains()[0], fact, &bicgstab, 1 + finestLevel);
  }
  else if(str == "simple"){
    multigridSolver.define(m_amr->getDomains()[0], fact, &simpleSolver, 1 + finestLevel);
  }
  multigridSolver.setSolverParameters(numSmooth, numSmooth, numSmooth, 1, 32, tolerance, 1E-60, 1E-60);
  
  // Solve
  Vector<LevelData<EBCellFAB>* > zer;
  Vector<LevelData<EBCellFAB>* > phi;
  Vector<LevelData<EBCellFAB>* > res;
  Vector<LevelData<EBCellFAB>* > rhs; 

  m_amr->alias(phi, a_phi);
  m_amr->alias(rhs, rho);
  m_amr->alias(res, a_residue);
  m_amr->alias(zer, zero);

  multigridSolver.m_verbosity=10;
  multigridSolver.init(phi, rhs, finestLevel, 0);
  multigridSolver.m_convergenceMetric = multigridSolver.computeAMRResidual(zer, rhs, finestLevel, baseLevel);
  Real t1 = -MPI_Wtime();
  multigridSolver.solveNoInit(phi, rhs, finestLevel, baseLevel, false);
  t1 += MPI_Wtime();
  if(procID() == 0) std::cout << "Multigrid solve time for Helm = " << t1 << std::endl;
}

void ProxyFieldSolver::solveMF(MFAMRCellData&       a_potential,
			       const MFAMRCellData& a_rho,
			       const EBAMRIVData&   a_sigma,
			       const bool           a_zerophi){
  const RefCountedPtr<EBIndexSpace>& ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);
  
  ParmParse pp("ProxyFieldSolver");
  ParmParse pp2("RodDielectric");
  
  Real alpha;
  Real beta;

  Real aco;
  Real bco;
  Real perm;

  MFAMRCellData Aco;
  MFAMRFluxData Bco;
  MFAMRIVData   BcoIrreg;

  m_amr->allocate(Aco,      m_realm, 1);
  m_amr->allocate(Bco,      m_realm, 1);
  m_amr->allocate(BcoIrreg, m_realm, 1);

  pp.get("alpha", alpha);
  pp.get("beta",  beta);
  
  pp.get("aco", aco);
  pp.get("bco",  bco);

  DataOps::setValue(Aco,      aco);
  DataOps::setValue(Bco,      bco);
  DataOps::setValue(BcoIrreg, bco);

  if(!ebisSol.isNull()){
    EBAMRFluxData B2    = m_amr->alias(phase::solid, Bco);
    EBAMRIVData   B2Irr = m_amr->alias(phase::solid, BcoIrreg);
    
    pp2.get("dielectric.permittivity", perm);

    DataOps::setValue(B2,    perm);
    DataOps::setValue(B2Irr, perm);
  }

  const int numPhases = m_multifluidIndexSpace->numPhases();


  const int finestLevel = m_amr->getFinestLevel();
  Vector<MFLevelGrid>             mflg(1 + finestLevel);
  Vector<MFMultigridInterpolator> mfInterp(1 + finestLevel);
  Vector<MFFluxReg>               mfFluxReg(1 + finestLevel);
  Vector<MFCoarAve>               mfCoarAve(1 + finestLevel);

  for (int lvl = 0; lvl <= finestLevel; lvl++){
    Vector<EBLevelGrid>                    eblgPhases(numPhases);
    Vector<RefCountedPtr<EBMultigridInterpolator> > quadPhases(numPhases);
    Vector<RefCountedPtr<EBFluxRegister> > fluxPhases(numPhases);
    Vector<RefCountedPtr<EbCoarAve> >      avePhases(numPhases);

    if(!ebisGas.isNull()) eblgPhases[phase::gas]   = *(m_amr->getEBLevelGrid(m_realm, phase::gas)  [lvl]);
    if(!ebisSol.isNull()) eblgPhases[phase::solid] = *(m_amr->getEBLevelGrid(m_realm, phase::solid)[lvl]);

    if(!ebisGas.isNull()) quadPhases[phase::gas]   = this->getMultigridInterpolators(phase::gas)  [lvl];
    if(!ebisSol.isNull()) quadPhases[phase::solid] = this->getMultigridInterpolators(phase::solid)[lvl];

    if(!ebisGas.isNull()) fluxPhases[phase::gas]   = (m_amr->getFluxRegister(m_realm, phase::gas)  [lvl]);
    if(!ebisSol.isNull()) fluxPhases[phase::solid] = (m_amr->getFluxRegister(m_realm, phase::solid)[lvl]);

    if(!ebisGas.isNull()) avePhases[phase::gas]   = (m_amr->getCoarseAverage(m_realm, phase::gas)  [lvl]);
    if(!ebisSol.isNull()) avePhases[phase::solid] = (m_amr->getCoarseAverage(m_realm, phase::solid)[lvl]);

    mflg[lvl].define(m_multifluidIndexSpace, eblgPhases);
    mfInterp[lvl].define(quadPhases);
    mfFluxReg[lvl].define(fluxPhases);
    mfCoarAve[lvl].define(avePhases);
  }

  // Set up BC factories
  int         eb_order;
  int         eb_weight;
  Real        eb_value;
  Real        dom_value;
  std::string str;

  RefCountedPtr<MFHelmholtzEBBCFactory> ebbcFactory;
  RefCountedPtr<MFHelmholtzDomainBCFactory> domainBcFactory;

  // EBBC
  pp.get("eb_bc", str);
  if(str == "dirichlet"){
    pp.get("eb_order",  eb_order);
    pp.get("eb_weight", eb_weight);
    pp.get("eb_val",    eb_value);

    ebbcFactory = RefCountedPtr<MFHelmholtzEBBCFactory> (new MFHelmholtzDirichletEBBCFactory(eb_order, eb_weight, eb_value));

    // Make sure we have enough ghost cells.
    if(eb_order > m_amr->getNumberOfGhostCells()) MayDay::Abort("ProxyFieldSolver::solveHelm - not enough ghost cells!");
  }
  else if(str == "neumann"){
    pp.get("eb_order",  eb_order);
    pp.get("eb_weight", eb_weight);
    pp.get("eb_val",    eb_value);

    ebbcFactory = RefCountedPtr<MFHelmholtzEBBCFactory> (new MFHelmholtzNeumannEBBCFactory(eb_order, eb_weight, eb_value));
  }
  else if(str == "robin"){
    pp.get("eb_order",  eb_order);
    pp.get("eb_weight", eb_weight);
    pp.get("eb_val",    eb_value);



    // Coefficients for radiative transfer with Robin. 
    const Real A =  1.5*eb_value;
    const Real B = -1.0*eb_value;
    const Real C =  0.0;

    ebbcFactory = RefCountedPtr<MFHelmholtzEBBCFactory> (new MFHelmholtzRobinEBBCFactory(eb_order, eb_weight, A, B, C));
  }
  else{
    MayDay::Error("ProxyFieldSolver::solveEBCond - uknown EBBC factory requested");
  }

  // Domain bc
  pp.get("domain_bc", str);
  pp.get("dom_val", dom_value);

  // BC function. Spatially dependent. 
  auto bcFunction = [dom_value](const RealVect& a_pos) -> Real {
    return dom_value;
  };
  
  if(str == "dirichlet"){
    domainBcFactory = RefCountedPtr<MFHelmholtzDomainBCFactory>(new MFHelmholtzDirichletDomainBCFactory(bcFunction));
  }
  else if(str == "neumann"){
    domainBcFactory = RefCountedPtr<MFHelmholtzDomainBCFactory>(new MFHelmholtzNeumannDomainBCFactory(bcFunction));
  }
  else if(str == "robin"){
    
    // Coeffs for radiative transfer with Robin
    const Real A =  1.5*dom_value;
    const Real B = -1.0*dom_value;
    const Real C =  0.0;
    
    domainBcFactory = RefCountedPtr<MFHelmholtzDomainBCFactory>(new MFHelmholtzRobinDomainBCFactory(A, B, C));
  }
  else{
    MayDay::Abort("ProxyFieldSolver::solveHelm - unknown domain bc requested");
  }

  int minCells;
  pp.get("min_cells", minCells);
  ProblemDomain bottomDomain = m_amr->getDomains()[0];
  while(bottomDomain.domainBox().shortside() >= 2*minCells){
    bottomDomain.coarsen(2);
  }


  const IntVect ghostPhi = m_amr->getNumberOfGhostCells()*IntVect::Unit;
  const IntVect ghostRhs = m_amr->getNumberOfGhostCells()*IntVect::Unit;

  MFHelmholtzOp::RelaxationMethod relaxType;
  int relax;
  pp.get("relax",relax);
  switch(relax){
  case 0:
    relaxType = MFHelmholtzOp::RelaxationMethod::PointJacobi;
    break;
  case 1:
    relaxType = MFHelmholtzOp::RelaxationMethod::GauSaiMultiColor;
    break;
  case 2:
    relaxType = MFHelmholtzOp::RelaxationMethod::GauSaiRedBlack;
    break;
  default:
    MayDay::Error("ProxyFieldSolver::solveMF - bad relaxation method");
  }

  int jumpOrder;
  int jumpWeight;

  pp.get("jump_order",  jumpOrder);
  pp.get("jump_weight", jumpWeight);

  MFHelmholtzOpFactory* fact = new MFHelmholtzOpFactory(m_multifluidIndexSpace,
							alpha,
							beta,
							m_amr->getProbLo(),
							mflg,
							mfInterp,
							mfFluxReg,
							mfCoarAve,
							m_amr->getRefinementRatios(),
							m_amr->getDx(),
							Aco.getData(),
							Bco.getData(),
							BcoIrreg.getData(),
							domainBcFactory,
							ebbcFactory,
							ghostPhi,
							ghostRhs,
							relaxType,
							bottomDomain,
							jumpOrder,
							jumpWeight,
							m_amr->getMaxBoxSize());

  fact->setJump(a_sigma, 1./Units::eps0);

  
  // Define the multigrid solver
  AMRMultiGrid<LevelData<MFCellFAB> >   multigridSolver;  


  int  numSmooth;
  Real tolerance;

  pp.get("mg_tol",   tolerance);
  pp.get("smooth",   numSmooth);

  const int baseLevel   = 0;

  MFSimpleSolver mfSolver;
  mfSolver.setNumSmooths(300);

  BiCGStabSolver<LevelData<MFCellFAB> > bicgstab;
  //  bicgstab.m_verbosity=4;

  // Define
  //  multigridSolver.define(m_amr->getDomains()[0], *fact, &mfSolver, 1 + finestLevel);
  multigridSolver.define(m_amr->getDomains()[0], *fact, &bicgstab, 1 + finestLevel);
  multigridSolver.setSolverParameters(numSmooth, numSmooth, numSmooth, 1, 32, tolerance, 1E-60, 1E-60);

  // Solve
  Vector<LevelData<MFCellFAB>* > phi;
  Vector<LevelData<MFCellFAB>* > rhs;
  Vector<LevelData<MFCellFAB>* > res;
  //  Vector<LevelData<MFCellFAB>* > zer;

  m_amr->alias(phi, a_potential);
  m_amr->alias(rhs, a_rho);
  m_amr->alias(res, m_residue);
  //  m_amr->alias(zer, zero);

  multigridSolver.m_verbosity = 10;
  multigridSolver.init(phi, rhs, finestLevel, baseLevel);
  Real t1 = -MPI_Wtime();
  multigridSolver.solveNoInit(phi, rhs, finestLevel, baseLevel, false);
  t1 += MPI_Wtime();
  if(procID() == 0) std::cout << "Multigrid solve time for MFHelm = " << t1 << std::endl;
}

#include <CD_NamespaceFooter.H>
