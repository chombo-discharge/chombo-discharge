/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*! 
  @file   CD_FieldSolver.cpp
  @brief  Implementation of CD_FieldSolver.H
  @author Robert Marskar
*/

// Std includes
#include <iostream>

// Chombo includes
#include <EBArith.H>
#include <ParmParse.H>
#include <MFAMRIO.H>
#include <EBAMRIO.H>

// Our includes
#include <CD_FieldSolver.H>
#include <CD_MultifluidAlias.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

Real FieldSolver::s_defaultDomainBcFunction(const RealVect a_position, const Real a_time){
  return 1.0;
}

Real FieldSolver::s_voltageOne(const Real a_time){
  return 1.0;
}

FieldSolver::FieldSolver(){
  m_className = "FieldSolver";
  m_realm     = Realm::Primal;

  this->setVerbosity(-1);
}

FieldSolver::~FieldSolver(){
  
}

bool FieldSolver::solve(const bool a_zerophi) {
  CH_TIME("FieldSolver::solve(bool)");
  if(m_verbosity > 5){
    pout() << "FieldSolver::solve(bool)" << endl;
  }

  const bool converged = this->solve(m_potential, m_rho, a_zerophi);
  
  return converged;
}

bool FieldSolver::solve(MFAMRCellData& a_potential, const bool a_zerophi){
  CH_TIME("FieldSolver::solve(MFAMRCellData, bool)");
  if(m_verbosity > 5){
    pout() << "FieldSolver::solve(MFAMRCellData, bool)" << endl;
  }

  const bool converged = this->solve(a_potential, m_rho, a_zerophi);

  return converged;
}

void FieldSolver::computeElectricField(){
  CH_TIME("FieldSolver::computeElectricField()");
  if(m_verbosity > 5){
    pout() << "FieldSolver::computeElectricField()" << endl;
  }

  this->computeElectricField(m_electricField, m_potential);
}

void FieldSolver::computeElectricField(MFAMRCellData& a_electricField, const MFAMRCellData& a_potential){
  CH_TIME("FieldSolver::computeElectricField(MFAMRCellData, MFAMRCellData)");
  if(m_verbosity > 5){
    pout() << "FieldSolver::computeElectricField(MFAMRCellData, MFAMRCellData)" << endl;
  }

  m_amr->computeGradient(a_electricField, a_potential, m_realm);
  DataOps::scale(a_electricField, -1.0);

  m_amr->averageDown(a_electricField, m_realm);
  m_amr->interpGhost(a_electricField, m_realm);
}

void FieldSolver::allocateInternals(){
  CH_TIME("FieldSolver::allocateInternals");
  if(m_verbosity > 5){
    pout() << "FieldSolver::allocateInternals" << endl;
  }

  const int ncomp = 1;

  m_amr->allocate(m_potential,     m_realm, ncomp);
  m_amr->allocate(m_rho,           m_realm, ncomp);
  m_amr->allocate(m_residue,       m_realm, ncomp);
  m_amr->allocate(m_sigma,         m_realm, phase::gas, ncomp);
  m_amr->allocate(m_electricField, m_realm, SpaceDim);

  DataOps::setValue(m_potential,     0.0);
  DataOps::setValue(m_rho,           0.0);
  DataOps::setValue(m_sigma,         0.0);
  DataOps::setValue(m_residue,       0.0);
  DataOps::setValue(m_electricField, 0.0);
}

void FieldSolver::preRegrid(const int a_lbase, const int a_oldFinestLevel){
  CH_TIME("FieldSolver::preRegrid");
  if(m_verbosity > 5){
    pout() << "FieldSolver::preRegrid" << endl;
  }

  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();
  
  m_amr->allocate(m_cache, m_realm, ncomp);
  
  for (int lvl = 0; lvl <= a_oldFinestLevel; lvl++){
    m_potential[lvl]->localCopyTo(*m_cache[lvl]);
  }
}

void FieldSolver::computeDisplacementField(MFAMRCellData& a_displacementField, const MFAMRCellData& a_electricField){
  CH_TIME("FieldSolver::computeDisplacementField");
  if(m_verbosity > 5){
    pout() << "FieldSolver::computeDisplacementField" << endl;
  }

  const Vector<Dielectric>& dielectrics = m_computationalGeometry->getDielectrics();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    LevelData<MFCellFAB>& D = *a_displacementField[lvl];
    LevelData<MFCellFAB>& E = *a_electricField[lvl];

    LevelData<EBCellFAB> D_gas, D_solid;
    LevelData<EBCellFAB> E_gas, E_solid;

    // This is all we need for the gas phase
    MultifluidAlias::aliasMF(D_gas,   phase::gas, D);
    MultifluidAlias::aliasMF(E_gas,   phase::gas, E);
    E_gas.localCopyTo(D_gas);
    DataOps::scale(D_gas,   Units::eps0);

    // For the solid phase, we multiply by epsilon. 
    if(m_multifluidIndexSpace->numPhases() > 1){
      MultifluidAlias::aliasMF(D_solid, phase::solid, D);
      MultifluidAlias::aliasMF(E_solid, phase::solid, E);
      E_solid.localCopyTo(D_solid);
      DataOps::scale(D_solid, Units::eps0);

      // Now scale by relative epsilon
      if(dielectrics.size() > 0){
	const RealVect dx            = m_amr->getDx()[lvl]*RealVect::Unit;
	const RealVect origin        = m_amr->getProbLo();
	const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	  EBCellFAB& dg = D_gas[dit()];

	  const Box box = dbl.get(dit());
	  for (VoFIterator vofit(IntVectSet(box), dg.getEBISBox().getEBGraph()); vofit.ok(); ++vofit){
	    const VolIndex& vof = vofit();
	    const RealVect& pos = EBArith::getVofLocation(vof, origin, dx);

	    Real dist = std::numeric_limits<Real>::infinity();
	    int closest;
	    
	    for (int i = 0; i < dielectrics.size(); i++){
	      const RefCountedPtr<BaseIF> func = dielectrics[i].getImplicitFunction();

	      const Real cur_dist = func->value(pos);
	
	      if(std::abs(cur_dist) <= std::abs(dist)){
		dist    = cur_dist;
		closest = i;
	      }
	    }
	    
	    const Real eps = dielectrics[closest].getPermittivity(pos);

	    for (int comp = 0; comp < dg.nComp(); comp++){
	      dg(vof, comp) *= eps;
	    }
	  }
	}
      }
    }
  }
}

Real FieldSolver::computeEnergyDensity(const MFAMRCellData& a_electricField){
  CH_TIME("FieldSolver::computeEnergyDensity");
  if(m_verbosity > 5){
    pout() << "FieldSolver::computeEnergyDensity" << endl;
  }

  MFAMRCellData D, EdotD;   
  m_amr->allocate(D, m_realm, SpaceDim);
  m_amr->allocate(EdotD, m_realm, 1);
  this->computeDisplacementField(D, a_electricField);
  DataOps::dotProduct(EdotD, D, a_electricField);

  Real U_g = 0.0;
  Real U_s = 0.0;

  // Energy in gas phase
  EBAMRCellData data_g;
  m_amr->allocatePointer(data_g);
  m_amr->alias(data_g, phase::gas, EdotD);
  m_amr->averageDown(data_g, m_realm, phase::gas);
  DataOps::norm(U_g, *data_g[0], m_amr->getDomains()[0], 1);

  if(m_multifluidIndexSpace->numPhases() > 1){
    EBAMRCellData data_s;
    m_amr->allocatePointer(data_s);
    m_amr->alias(data_s, phase::solid, EdotD);
    m_amr->averageDown(data_s, m_realm, phase::solid);
    DataOps::norm(U_s, *data_s[0], m_amr->getDomains()[0], 1);
  }

  return 0.5*(U_g + U_s);
}

Real FieldSolver::computeCapacitance(){
  CH_TIME("FieldSolver::computeCapacitance");
  if(m_verbosity > 5){
    pout() << "FieldSolver::computeCapacitance" << endl;
  }

  // TLDR; We MUST compute the energy density with the Laplace field, so no sources here...
  Real C;

  MFAMRCellData phi;
  MFAMRCellData source;
  EBAMRIVData   sigma;

  m_amr->allocate(phi,    m_realm, 1);
  m_amr->allocate(source, m_realm, 1);
  m_amr->allocate(sigma,  m_realm, phase::gas, 1);

  DataOps::setValue(phi,    0.0);
  DataOps::setValue(source, 0.0);
  DataOps::setValue(sigma,  0.0);

  this->solve(phi, source, sigma);

  // Solve and compute energy density
  MFAMRCellData E;
  m_amr->allocate(E, m_realm, SpaceDim);
  m_amr->computeGradient(E, phi, m_realm); // -E
  const Real U = this->computeEnergyDensity(E); // Energy density

  // U = 0.5*CV^2
  const Real pot = m_voltage(m_time);
  if(pot == 0.0){
    MayDay::Abort("FieldSolver::compute_capacitance - error, can't compute energy density with V = 0");
  }
  C = 2.0*U/(pot*pot);

  return C;
}

void FieldSolver::deallocateInternals(){
  CH_TIME("FieldSolver::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "FieldSolver::deallocateInternals" << endl;
  }
  
  m_amr->deallocate(m_potential);
  m_amr->deallocate(m_rho);
  m_amr->deallocate(m_residue);
  m_amr->deallocate(m_sigma);
}

void FieldSolver::regrid(const int a_lmin, const int a_old_finest, const int a_new_finest){
  CH_TIME("FieldSolver::regrid");
  if(m_verbosity > 5){
    pout() << "FieldSolver::regrid" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  this->allocateInternals();

  for (int i = 0; i < phase::numPhases; i++){
    phase::which_phase cur_phase;    
    if(i == 0){
      cur_phase = phase::gas;
    }
    else{
      cur_phase = phase::solid;
    }

    const RefCountedPtr<EBIndexSpace>& ebis = m_multifluidIndexSpace->getEBIndexSpace(cur_phase);

    if(!ebis.isNull()){
      EBAMRCellData potential_phase = m_amr->alias(cur_phase, m_potential);
      EBAMRCellData scratch_phase   = m_amr->alias(cur_phase, m_cache);

      Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->getPwlInterpolator(m_realm, cur_phase);

      // These levels have not changed
      for (int lvl = 0; lvl <= Max(0, a_lmin-1); lvl++){
	scratch_phase[lvl]->copyTo(*potential_phase[lvl]); // Base level should never change, but ownership might
      }
      for (int lvl = Max(1,a_lmin); lvl <= a_new_finest; lvl++){
	interpolator[lvl]->interpolate(*potential_phase[lvl], *potential_phase[lvl-1], interv);
	if(lvl <= a_old_finest){
	  scratch_phase[lvl]->copyTo(*potential_phase[lvl]);
	}

	potential_phase[lvl]->exchange();
      }
    }
  }

  m_amr->averageDown(m_potential, m_realm);
  m_amr->interpGhost(m_potential, m_realm);

  // Now recompute E
  this->computeElectricField();
}

void FieldSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry){
  CH_TIME("FieldSolver::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setComputationalGeometry" << endl;
  }

  m_computationalGeometry = a_computationalGeometry;
  m_multifluidIndexSpace     = m_computationalGeometry->getMfIndexSpace();

  this->setDefaultEbBcFunctions();
}

void FieldSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("FieldSolver::setAmr");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setAmr" << endl;
  }

  m_amr = a_amr;
}

void FieldSolver::setVoltage(std::function<Real(const Real a_time)> a_voltage){
  m_voltage = a_voltage;
}

void FieldSolver::setDomainBcWallFunction(const int a_dir, const Side::LoHiSide a_side, const ElectrostaticDomainBc::BcFunction& a_function){
  CH_TIME("FieldSolver::setDomainBcWallFunction");
  if(m_verbosity > 4){
    pout() << "FieldSolver::setDomainBcWallFunction" << endl;
  }

  const ElectrostaticDomainBc::Wall curWall = std::make_pair(a_dir, a_side);

  m_domainBcFunctions.at(curWall) = a_function;
}

void FieldSolver::setElectrodeDirichletFunction(const int a_electrode, const ElectrostaticEbBc::BcFunction& a_function){
  CH_TIME("FieldSolver::setElectrodeDirichletFunction");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setElectrodeDirichletFunction" << endl;
  }

  m_ebBc.setEbBc(a_electrode, a_function);
}

void FieldSolver::setVerbosity(const int a_verbosity){
  CH_TIME("FieldSolver::setVerbosity");
  m_verbosity = a_verbosity;

  if(m_verbosity > 4){
    pout() << "FieldSolver::setVerbosity" << endl;
  }
}

void FieldSolver::setTime(const int a_timeStep, const Real a_time, const Real a_dt) {
  CH_TIME("FieldSolver::setTime");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setTime" << endl;
  }
  
  m_timeStep = a_timeStep;
  m_time     = a_time;
  m_dt       = a_dt;
}

void FieldSolver::setRealm(const std::string a_realm){
  CH_TIME("FieldSolver::setRealm");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setRealm" << endl;
  }  
  m_realm = a_realm;
}

const std::string FieldSolver::getRealm() const{
  CH_TIME("FieldSolver::getRealm");
  if(m_verbosity > 5){
    pout() << "FieldSolver::getRealm" << endl;
  }
  
  return m_realm;
}

void FieldSolver::parsePlotVariables(){
  CH_TIME("FieldSolver::parsePlotVariables");
  if(m_verbosity > 5){
    pout() << "FieldSolver::parsePlotVariables" << endl;
  }

  m_plotPotential     = false;
  m_plotRho           = false;
  m_plotElectricField = false;
  m_plotResidue       = false;

  ParmParse pp(m_className.c_str());
  const int num = pp.countval("plt_vars");

  if(num > 0){
    Vector<std::string> str(num);
    pp.getarr("plt_vars", str, 0, num);

    for (int i = 0; i < num; i++){
      if(     str[i] == "phi")   m_plotPotential     = true;
      else if(str[i] == "rho")   m_plotRho           = true;
      else if(str[i] == "resid") m_plotResidue       = true;
      else if(str[i] == "E")     m_plotElectricField = true;
    }
  }
}

std::string FieldSolver::makeBcString(const int a_dir, const Side::LoHiSide a_side) const {
  CH_TIME("FieldSolver::makeBcString");
  if(m_verbosity > 5){
    pout() << "FieldSolver::makeBcString" << endl;
  }

  std::string strDir;
  std::string strSide;
  
  if(a_dir == 0){
    strDir = "x";
  }
  else if(a_dir == 1){
    strDir = "y";
  }
  else if(a_dir == 2){
    strDir = "z";
  }

  if(a_side == Side::Lo){
    strSide = "lo";
  }
  else if(a_side == Side::Hi){
    strSide = "hi";
  }

  const std::string ret = std::string("bc.") + strDir + std::string(".") + strSide;

  return ret;
}

ElectrostaticDomainBc::BcType FieldSolver::parseBcString(const std::string a_str) const {
  CH_TIME("FieldSolver::parseBcString");
  if(m_verbosity > 5){
    pout() << "FieldSolver::parseBcString" << endl;
  }

  ElectrostaticDomainBc::BcType ret;

  if(a_str == "dirichlet"){
    ret = ElectrostaticDomainBc::BcType::Dirichlet;
  }
  else if(a_str == "neumann"){
    ret = ElectrostaticDomainBc::BcType::Neumann;
  }
  else{
    MayDay::Abort("ElectrostaticDomainBc::BcType - unknown BC type!");
  }

  return ret;
}

void FieldSolver::setDefaultDomainBcFunctions(){
  CH_TIME("FieldSolver::setDefaultDomainBcFunctions");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setDefaultDomainBcFunctions" << endl;
  }

  m_domainBcFunctions.clear();
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      m_domainBcFunctions.emplace(std::make_pair(dir, sit()), FieldSolver::s_defaultDomainBcFunction);
    }
  }
}

void FieldSolver::setDefaultEbBcFunctions() {
  CH_TIME("FieldSolver::setDefaultEbBcFunctions");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setDefaultEbBcFunctions" << endl;
  }

  const Vector<Electrode> electrodes = m_computationalGeometry->getElectrodes();

  for (int i = 0; i < electrodes.size(); i++){
    const Electrode& elec = electrodes[i];

    Real val  = (elec.isLive()) ? 1.0 : 0.0;
    Real frac = elec.getFraction();
    
    ElectrostaticEbBc::BcFunction curFunc = [&, val, frac](const RealVect a_position, const Real a_time){
      return m_voltage(m_time)*val*frac;
    };

    m_ebBc.addEbBc(elec, curFunc);
  }
}


void FieldSolver::parseDomainBc(){
  CH_TIME("FieldSolver::parseDomainBc");
  if(m_verbosity > 5){
    pout() << "FieldSolver::parseDomainBc" << endl;
  }
  
  ParmParse pp(m_className.c_str());

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){

      const ElectrostaticDomainBc::Wall curWall = std::make_pair(dir, sit());
      const std::string bcString = this->makeBcString(dir, sit());
      const int num = pp.countval(bcString.c_str());

      std::string str;

      ElectrostaticDomainBc::BcType      bcType;
      ElectrostaticDomainBc::BcFunction& bcFunc = m_domainBcFunctions.at(curWall);

      std::function<Real(const RealVect, const Real)> curFunc;

      if(num == 1){
	pp.get(bcString.c_str(), str, 0);

	curFunc = [&] (const RealVect a_pos, const Real a_time) {
	  return bcFunc(a_pos, m_time);
	};

	if(str == "dirichlet_custom"){
	  bcType  = ElectrostaticDomainBc::BcType::Dirichlet;
	}
	else if(str == "neumann_custom"){
	  bcType  = ElectrostaticDomainBc::BcType::Neumann;
	}
	else{
	  MayDay::Abort("FieldSolver::parseDomainBc -- got only one argument but this argument was not dirichlet_custom. Maybe you have the wrong BC specification?");
	}
      }
      else if(num == 2){
	Real val;

	pp.get(bcString.c_str(), str, 0);
	pp.get(bcString.c_str(), val, 1);

	bcType = this->parseBcString(str);

	// Build a function computing the value at the boundary. 
	switch (bcType){
	case ElectrostaticDomainBc::BcType::Dirichlet:
	  curFunc = [&, val] (const RealVect a_pos, const Real a_time){
	    return bcFunc(a_pos, m_time)*m_voltage(m_time)*val;
	  };
	  break;
	case ElectrostaticDomainBc::BcType::Neumann:
	  curFunc = [&, val] (const RealVect a_pos, const Real a_time){
	    return bcFunc(a_pos, m_time)*val;
	  };
	  break;
	default:
	  MayDay::Abort("FieldSolver::parseDomainBc -- unsupported boundary condition requested!");
	  break;
	}
      }
      else{
	const std::string errorString = "FieldSolver::parseDomainBc -- bad or no input parameter for " + bcString;
	MayDay::Abort(errorString.c_str());
      }

      m_domainBc.setBc(curWall, std::make_pair(bcType, curFunc));
    }
  }
}

void FieldSolver::writePlotFile(){
  CH_TIME("FieldSolver::writePlotFile");
  if(m_verbosity > 5){
    pout() << "FieldSolver::writePlotFile" << endl;
  }

  // Number of output components and their names
  const int ncomps = this->getNumberOfPlotVariables();
  const Vector<std::string> names = getPlotVariableNames();

  // Allocate storage for output
  EBAMRCellData output;
  m_amr->allocate(output, m_realm, phase::gas, ncomps);

  // Copy internal data to be plotted over to 'output'
  int icomp = 0;
  writePlotData(output, icomp);

  // Filename
  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_className.c_str(), m_timeStep, SpaceDim);

  // Alias
  Vector<LevelData<EBCellFAB>* > output_ptr(1+m_amr->getFinestLevel());
  m_amr->alias(output_ptr, output);

  Vector<Real> covered_values(ncomps, 0.0);
  string fname(file_char);
  writeEBHDF5(fname,
	      m_amr->getGrids(m_realm),
	      output_ptr,
	      names,
	      m_amr->getDomains()[0].domainBox(),
	      m_amr->getDx()[0],
	      m_dt,
	      m_time,
	      m_amr->getRefinementRatios(),
	      m_amr->getFinestLevel() + 1,
	      false,
	      covered_values,
	      IntVect::Unit);
}

void FieldSolver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("FieldSolver::writeCheckpointLevel");
  if(m_verbosity > 5){
    pout() << "FieldSolver::writeCheckpointLevel" << endl;
  }

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Used for aliasing phases
  LevelData<EBCellFAB> potential_gas;
  LevelData<EBCellFAB> potential_sol;

  if(!ebis_gas.isNull()) MultifluidAlias::aliasMF(potential_gas,  phase::gas,   *m_potential[a_level]);
  if(!ebis_sol.isNull()) MultifluidAlias::aliasMF(potential_sol,  phase::solid, *m_potential[a_level]);
  
  // Write data
  if(!ebis_gas.isNull()) write(a_handle, potential_gas, "poisson_g");
  if(!ebis_sol.isNull()) write(a_handle, potential_sol, "poisson_s");
}

void FieldSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level){
  CH_TIME("FieldSolver::readCheckpointLevel");
  if(m_verbosity > 5){
    pout() << "FieldSolver::readCheckpointLevel" << endl;
  }

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Used for aliasing phases
  LevelData<EBCellFAB> potential_gas;
  LevelData<EBCellFAB> potential_sol;

  if(!ebis_gas.isNull()) MultifluidAlias::aliasMF(potential_gas,  phase::gas,   *m_potential[a_level]);
  if(!ebis_sol.isNull()) MultifluidAlias::aliasMF(potential_sol,  phase::solid, *m_potential[a_level]);
  
  // Read data
  if(!ebis_gas.isNull()) read<EBCellFAB>(a_handle, potential_gas, "poisson_g", m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);
  if(!ebis_sol.isNull()) read<EBCellFAB>(a_handle, potential_sol, "poisson_s", m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);
}

void FieldSolver::postCheckpoint(){
  CH_TIME("FieldSolver::postCheckpoint");
  if(m_verbosity > 5){
    pout() << "FieldSolver::postCheckpoint" << endl;
  }
}

void FieldSolver::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("FieldSolver::writePlotData");
  if(m_verbosity > 5){
    pout() << "FieldSolver::writePlotData" << endl;
  }

  // Add phi to output
  if(m_plotPotential) {
    this->writeMultifluidData(a_output, a_comp, m_potential,  true);
  }
  if(m_plotRho) {
    this->writeMultifluidData(a_output, a_comp, m_rho, false);
    //    DataOps::setCoveredValue(a_output, a_comp-1, 0.0); // Why was this here? If the dielectric side contains charge, it should not be here. 
  }
  if(m_plotResidue) {
    this->writeMultifluidData(a_output, a_comp, m_residue,  false);
  }
  if(m_plotElectricField) {
    //    this->computeElectricField();
    this->writeMultifluidData(a_output, a_comp, m_electricField, true);
  }
}

void FieldSolver::writeMultifluidData(EBAMRCellData& a_output, int& a_comp, const MFAMRCellData& a_data, const bool a_interp){
  CH_TIME("FieldSolver::writeMultifluidData");
  if(m_verbosity > 5){
    pout() << "FieldSolver::writeMultifluidData" << endl;
  }

  const int ncomp = a_data[0]->nComp();

  const RefCountedPtr<EBIndexSpace> ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Allocate some scratch data that we can use
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, phase::gas, ncomp);

  //
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    LevelData<EBCellFAB> data_gas;
    LevelData<EBCellFAB> data_sol;

    if(!ebis_gas.isNull()) MultifluidAlias::aliasMF(data_gas,  phase::gas,   *a_data[lvl]);
    if(!ebis_sol.isNull()) MultifluidAlias::aliasMF(data_sol,  phase::solid, *a_data[lvl]);

    data_gas.localCopyTo(*scratch[lvl]);

    // Copy all covered cells from the other phase
    if(!ebis_sol.isNull()){
      for (DataIterator dit = data_sol.dataIterator(); dit.ok(); ++dit){
	const Box box = data_sol.disjointBoxLayout().get(dit());
	const IntVectSet ivs(box);
	const EBISBox& ebisb_gas = data_gas[dit()].getEBISBox();
	const EBISBox& ebisb_sol = data_sol[dit()].getEBISBox();

	FArrayBox& scratch_gas    = (*scratch[lvl])[dit()].getFArrayBox();
	const FArrayBox& fab_gas = data_gas[dit()].getFArrayBox();
	const FArrayBox& fab_sol = data_sol[dit()].getFArrayBox();

	// TLDR: There are four cases
	// 1. Both are covered => inside electrode
	// 2. Gas is covered, solid is regular => inside solid phase only
	// 3. Gas is regular, solid is covered => inside gas phase only
	// 4. Gas is !(covered || regular) and solid is !(covered || regular) => on dielectric boundary
	if(!(ebisb_gas.isAllCovered() && ebisb_sol.isAllCovered())){ // Outside electrode, lets do stuff
	  if(ebisb_gas.isAllCovered() && ebisb_sol.isAllRegular()){ // Case 2, copy directly. 
	    scratch_gas += fab_sol;
	  }
	  else if(ebisb_gas.isAllRegular() && ebisb_sol.isAllCovered()) { // Case 3. Inside gas phase. Already did this. 
	  }
	  else{ // Case 4, needs special treatment. 
	    for (BoxIterator bit(box); bit.ok(); ++bit){ // Loop through all cells here
	      const IntVect iv = bit();
	      if(ebisb_gas.isCovered(iv) && !ebisb_sol.isCovered(iv)){   // Regular cells from phase 2
		for (int comp = 0; comp < ncomp; comp++){
		  scratch_gas(iv, comp) = fab_sol(iv, comp);
		}
	      }
	      else if(ebisb_sol.isIrregular(iv) && ebisb_gas.isIrregular(iv)){ // Irregular cells. Use gas side data
		for (int comp = 0; comp < ncomp; comp++){
		  scratch_gas(iv, comp) = fab_gas(iv,comp);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Average down shit and interpolate to centroids
  m_amr->averageDown(scratch, m_realm, phase::gas);
  m_amr->interpGhost(scratch, m_realm, phase::gas);
  if(a_interp){
    m_amr->interpToCentroids(scratch, m_realm, phase::gas);
  }

  const Interval src_interv(0, ncomp-1);
  const Interval dst_interv(a_comp, a_comp + ncomp - 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(m_realm == a_output.getRealm()){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else{
      scratch[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }

  a_comp += ncomp;
}

int FieldSolver::getNumberOfPlotVariables() const {
  CH_TIME("FieldSolver::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << "FieldSolver::getNumberOfPlotVariables" << endl;
  }

  int num_output = 0;

  if(m_plotPotential)     num_output = num_output + 1;
  if(m_plotRho)           num_output = num_output + 1;
  if(m_plotResidue)       num_output = num_output + 1;
  if(m_plotElectricField) num_output = num_output + SpaceDim;

  return num_output;
}

Vector<std::string> FieldSolver::getPlotVariableNames() const {
  CH_TIME("FieldSolver::getPlotVariableNames");
  if(m_verbosity > 5){
    pout() << "FieldSolver::getPlotVariableNames" << endl;
  }
  Vector<std::string> names(0);
  
  if(m_plotPotential) names.push_back("Electrostatic potential");
  if(m_plotRho)       names.push_back("Space charge density");
  if(m_plotResidue)   names.push_back("Electrostatic potential_residue");
  if(m_plotElectricField){
    names.push_back("x-Electric field"); 
    names.push_back("y-Electric field"); 
    if(SpaceDim == 3){
      names.push_back("z-Electric field");
    }
  }
  
  return names;
}

Vector<long long> FieldSolver::computeLoads(const DisjointBoxLayout& a_dbl, const int a_level){
  CH_TIME("FieldSolver::computeLoads");
  if(m_verbosity > 5){
    pout() << "FieldSolver::computeLoads" << endl;
  }

  // Compute number of cells
  Vector<int> numCells(a_dbl.size(), 0);
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit){
    numCells[dit().intCode()] = a_dbl[dit()].numPts();
  }

#ifdef CH_MPI
  int count = numCells.size();
  Vector<int> tmp(count);
  MPI_Allreduce(&(numCells[0]),&(tmp[0]), count, MPI_INT, MPI_SUM, Chombo_MPI::comm);
  numCells = tmp;
#endif

  Vector<long long> loads(numCells.size());
  for (int i = 0; i < loads.size(); i++){
    loads[i] = (long long) numCells[i];
  }

  return loads;
}

Real FieldSolver::getTime() const{
  return m_time;
}

MFAMRCellData& FieldSolver::getPotential(){
  return m_potential;
}

MFAMRCellData& FieldSolver::getElectricField(){
  return m_electricField;
}

MFAMRCellData& FieldSolver::getRho(){
  return m_rho;
}

MFAMRCellData& FieldSolver::getResidue(){
  return m_residue;
}

#include <CD_NamespaceFooter.H>
