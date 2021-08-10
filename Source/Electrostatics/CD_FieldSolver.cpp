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

constexpr int FieldSolver::m_comp;
constexpr int FieldSolver::m_nComp;

Real FieldSolver::s_defaultDomainBcFunction(const RealVect a_position, const Real a_time){
  return 1.0;
}

FieldSolver::FieldSolver(){

  // Default settings.
  m_className    = "FieldSolver";
  m_realm        = Realm::Primal;
  m_isVoltageSet = false;
  m_verbosity    = -1;

  this->setDataLocation(Location::Cell::Center);
  this->setDefaultDomainBcFunctions();
}

FieldSolver::~FieldSolver(){
  
}

void FieldSolver::setDataLocation(const Location::Cell a_dataLocation){
  CH_TIME("FieldSolver::setDataLocation(Location::Cell)");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setDataLocation(Location::Cell)" << endl;
  }
  
  switch(a_dataLocation){
  case Location::Cell::Center:{
    m_dataLocation = a_dataLocation;
    m_faceLocation = Location::Face::Center;
    break;
  }
  case Location::Cell::Centroid:{
    MayDay::Error("FieldSolver::setDataLocation - centroid discretization is not yet fully supported (due to the multigrid interpolator)");
    
    m_dataLocation = a_dataLocation;
    m_faceLocation = Location::Face::Centroid;
    break;
  }
  default:
    MayDay::Error("FieldSolver::setDataLocation - location must be either cell center or cell centroid");
    break;
  }
}

bool FieldSolver::solve(const bool a_zerophi) {
  CH_TIME("FieldSolver::solve(bool)");
  CH_assert(m_isVoltageSet);
  if(m_verbosity > 5){
    pout() << "FieldSolver::solve(bool)" << endl;
  }

  const bool converged = this->solve(m_potential, m_rho, a_zerophi);
  
  return converged;
}

bool FieldSolver::solve(MFAMRCellData& a_potential, const bool a_zerophi){
  CH_TIME("FieldSolver::solve(MFAMRCellData, bool)");
  CH_assert(m_isVoltageSet);  
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
  
#if 1 // Temporary warning showing that centroid-centered discretizations aren't ready (yet)
  if(m_dataLocation == Location::Cell::Centroid){
    pout() << "FieldSolver::computeElectricField - AmrMesh does not support centroid data for gradients (yet)" << endl;
  }
#endif

  // Update ghost cells. Use scratch storage for this. 
  MFAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_nComp);
  scratch.copy(a_potential);
  m_amr->averageDown(scratch, m_realm);
  m_amr->interpGhost(scratch, m_realm);

  // Compute the cell-centered gradient everywhere. 
  m_amr->computeGradient(a_electricField, scratch, m_realm);
  DataOps::scale(a_electricField, -1.0);

  // Coarsen solution and update ghost cells. 
  m_amr->averageDown(a_electricField, m_realm);
  m_amr->interpGhost(a_electricField, m_realm);
}

void FieldSolver::allocateInternals(){
  CH_TIME("FieldSolver::allocateInternals");
  if(m_verbosity > 5){
    pout() << "FieldSolver::allocateInternals" << endl;
  }

  m_amr->allocate(m_potential,        m_realm, m_nComp );
  m_amr->allocate(m_rho,              m_realm, m_nComp );
  m_amr->allocate(m_residue,          m_realm, m_nComp );
  m_amr->allocate(m_permittivityCell, m_realm, m_nComp);
  m_amr->allocate(m_permittivityFace, m_realm, m_nComp);
  m_amr->allocate(m_permittivityEB,   m_realm, m_nComp);
  m_amr->allocate(m_electricField,    m_realm, SpaceDim);

  m_amr->allocate(m_sigma,         m_realm, phase::gas, m_nComp);

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

  // TLDR: This version does a backup of m_potential on the old grid. This is later used for interpolating onto
  //       the new grid (in the regrid routine).
  
  m_amr->allocate(m_cache, m_realm, m_nComp);
  
  for (int lvl = 0; lvl <= a_oldFinestLevel; lvl++){
    m_potential[lvl]->localCopyTo(*m_cache[lvl]);
  }
}

void FieldSolver::computeDisplacementField(MFAMRCellData& a_displacementField, const MFAMRCellData& a_electricField){
  CH_TIME("FieldSolver::computeDisplacementField");
  if(m_verbosity > 5){
    pout() << "FieldSolver::computeDisplacementField" << endl;
  }

  // TLDR: This computes the displacement field D = eps*E on both phases. The data is either cell-centered or centroid-centered, and the
  //       permittivity can be spatially varying (in the dielectric). So, for the gas phase we only need to compute eps0*E while for the
  //       dielectric phase we have to iterate through each cell and find the corresponding permittivity.

  const Vector<Dielectric>& dielectrics = m_computationalGeometry->getDielectrics();

  // Make D = eps0*E
  DataOps::copy (a_displacementField, a_electricField);
  DataOps::scale(a_displacementField, Units::eps0);  
  //
  if(m_multifluidIndexSpace->numPhases() > 1 && dielectrics.size() > 0){

    // Lambda which returns the relative permittivity at some position
    auto permittivity = [&dielectrics](const RealVect& pos) -> Real{
      Real minDist = std::numeric_limits<Real>::infinity();
      int closest;
	    
      for (int i = 0; i < dielectrics.size(); i++){
	const RefCountedPtr<BaseIF> func = dielectrics[i].getImplicitFunction();

	const Real curDist = func->value(pos);
	
	if(std::abs(curDist) <= std::abs(minDist)){
	  minDist = curDist;
	  closest = i;
	}
      }
	    
      return dielectrics[closest].getPermittivity(pos);
    };

    // Iterate through data
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const Real dx                = m_amr->getDx()[lvl];
      const RealVect probLo        = m_amr->getProbLo();
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, phase::solid)[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit){
	const Box cellBox      = dbl[dit()];
	const EBISBox& ebisbox = ebisl[dit()];
	const EBGraph& ebgraph = ebisbox.getEBGraph();
	const IntVectSet ivs   = ebisbox.getIrregIVS(cellBox);

	// Get handle to data on the solid phase
	MFCellFAB& D    = (*a_displacementField[lvl])[dit()];
	EBCellFAB& Dsol = D.getPhase(phase::solid);
	FArrayBox& Dreg = Dsol.getFArrayBox();

	// Regular cells
	for (BoxIterator bit(cellBox); bit.ok(); ++bit){
	  const IntVect iv   = bit();
	  
	  if(ebisbox.isRegular(iv)){
	    const RealVect pos = probLo + dx*RealVect(iv);
	    const Real epsRel  = permittivity(pos);
	  
	    for (int comp = 0; comp < SpaceDim; comp++){
	      Dreg(iv, comp) *= epsRel;
	    }
	  }
	}

	// Irregular cells
	for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
	  const VolIndex& vof = vofit();
	  const RealVect pos  = probLo + Location::position(m_dataLocation, vof, ebisbox, dx);
	  const Real epsRel   = permittivity(pos);

	  for (int comp = 0; comp < SpaceDim; comp++){
	    Dsol(vof, comp) *= epsRel;
	  }
	}
      }
    }
  }
}

Real FieldSolver::computeEnergy(const MFAMRCellData& a_electricField){
  CH_TIME("FieldSolver::computeEnergy");
  if(m_verbosity > 5){
    pout() << "FieldSolver::computeEnergy" << endl;
  }

  // TLDR: This routine computes Int(E*D dV) over the entire domain. Since we use conservative averaging, we coarsen E*D onto
  //       the coarsest grid level and do the integratino there.

  const bool reallyMultiPhase = (m_multifluidIndexSpace->numPhases() > 1);

  // Allocate storage first, because we have no scratch storage for this. 
  MFAMRCellData D;
  MFAMRCellData EdotD;

  m_amr->allocate(D,     m_realm,  SpaceDim);
  m_amr->allocate(EdotD, m_realm,  1);  

  // Compute D = eps*E and then the dot product E*D. 
  this->computeDisplacementField(D, a_electricField);
  DataOps::dotProduct(EdotD, D, a_electricField);
  m_amr->averageDown(EdotD, m_realm);


  // This defines a lambda which computes the energy in a specified phase
  auto phaseEnergy = [&EdotD, &amr=this->m_amr] (const phase::which_phase a_phase) -> Real {
    const EBAMRCellData phaseEdotD = amr->alias(a_phase, EdotD);

    Real energy;
    
    DataOps::norm(energy, *phaseEdotD[0], amr->getDomains()[0], 1);

    return energy;
  };

  // Contributions from each phase. 
  const Real Ug = phaseEnergy(phase::gas);
  const Real Us = reallyMultiPhase ? phaseEnergy(phase::solid) : 0.0;

  return 0.5*(Ug + Us);
}

Real FieldSolver::computeCapacitance(){
  CH_TIME("FieldSolver::computeCapacitance");
  CH_assert(m_isVoltageSet);
  if(m_verbosity > 5){
    pout() << "FieldSolver::computeCapacitance" << endl;
  }

  // TLDR: The energy density U = 0.5*C*V^2, or U = Int(E*D dV). We set the potential to one and solve the Poisson equation. 
  //       We then compute U from E*D without sources and use C = 2*U/(V*V). The caveat to this approach is that the solver
  //       may have been set with a zero voltage, non-zero space charge and non-zero surface charge. Here, we do a "clean" solve
  //       with voltage set to 1, and without charges.

  MFAMRCellData phi;
  MFAMRCellData E;  
  MFAMRCellData source;
  EBAMRIVData   sigma;

  m_amr->allocate(phi,    m_realm, 1);
  m_amr->allocate(E,      m_realm, SpaceDim);
  m_amr->allocate(source, m_realm, 1);
  m_amr->allocate(sigma,  m_realm, phase::gas, 1);

  DataOps::setValue(phi,    0.0);
  DataOps::setValue(E,      0.0);
  DataOps::setValue(source, 0.0);
  DataOps::setValue(sigma,  0.0);

  // Do a backup of the voltage.
  auto voltageBackup = m_voltage;
  auto voltageOne    = [](const Real a_time) -> Real { return 100.0; };

  // Set the voltage to one and use that to compute the potential/E-field
  this->setVoltage(voltageOne);

  // Do a solve without a source term. 
  this->solve(phi, source, sigma, true);
  this->computeElectricField(E, phi);
  
  // The energy is U = 0.5*CV^2 so we can just compute the electrostatic energy and revert that expression.
  const Real U = this->computeEnergy(E); // Electrostatic energy
  const Real C = 2.0*U;                  // C = 2*U/(V*V) but voltage is 1.

  // Set the voltage back to what it was. 
  this->setVoltage(voltageBackup);

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

  const Interval interv(m_comp, m_comp);

  // Reallocate internals
  this->allocateInternals();

  for (int i = 0; i < phase::numPhases; i++){
    phase::which_phase curPhase;
    if(i == 0){
      curPhase = phase::gas;
    }
    else{
      curPhase = phase::solid;
    }

    const RefCountedPtr<EBIndexSpace>& ebis = m_multifluidIndexSpace->getEBIndexSpace(curPhase);

    if(!ebis.isNull()){
      EBAMRCellData potential_phase = m_amr->alias(curPhase, m_potential);
      EBAMRCellData scratch_phase   = m_amr->alias(curPhase, m_cache);

      Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->getPwlInterpolator(m_realm, curPhase);

      // These levels have not changed so we can do a direct copy. 
      for (int lvl = 0; lvl <= Max(0, a_lmin-1); lvl++){
	scratch_phase[lvl]->copyTo(*potential_phase[lvl]); // Base level should never change, but ownership might. So no localCopyTo here...
      }

      // These levels might have changed so we interpolate data here. 
      for (int lvl = Max(1,a_lmin); lvl <= a_new_finest; lvl++){
	interpolator[lvl]->interpolate(*potential_phase[lvl], *potential_phase[lvl-1], interv);
	if(lvl <= a_old_finest){
	  scratch_phase[lvl]->copyTo(*potential_phase[lvl]);
	}

	potential_phase[lvl]->exchange();
      }
    }
  }

  // Synchronize over levels. 
  m_amr->averageDown(m_potential, m_realm);
  m_amr->interpGhost(m_potential, m_realm);

  // Recompute E. 
  this->computeElectricField();

  // Set permittivities
  this->setPermittivities();
}

void FieldSolver::setRho(const Real a_rho){
  CH_TIME("FieldSolver::setRho(Real");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setRho(Real)" << endl;
  }
  
  DataOps::setValue(m_rho, a_rho);
}

void FieldSolver::setRho(const std::function<Real(const RealVect)>& a_rho){
  CH_TIME("FieldSolver::setRho(std::function");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setRho(std::function)" << endl;
  }
  
  DataOps::setValue(m_rho, a_rho, m_amr->getProbLo(), m_amr->getDx(), 0);
}

void FieldSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry){
  CH_TIME("FieldSolver::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setComputationalGeometry" << endl;
  }

  m_computationalGeometry = a_computationalGeometry;
  m_multifluidIndexSpace  = m_computationalGeometry->getMfIndexSpace();


  this->setDefaultEbBcFunctions();
}

void FieldSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("FieldSolver::setAmr(amr)");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setAmr(amr)" << endl;
  }

  m_amr = a_amr;
}

void FieldSolver::setVoltage(std::function<Real(const Real a_time)> a_voltage){
  CH_TIME("FieldSolver::setVoltage(std::function)");
  if(m_verbosity > 4){
    pout() << "FieldSolver::setVoltage(std::function)" << endl;
  }
  
  m_voltage      = a_voltage;
  m_isVoltageSet = true;
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

std::string FieldSolver::getRealm() const{
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
    MayDay::Error("ElectrostaticDomainBc::BcType - unknown BC type!");
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
  CH_assert(!m_computationalGeometry.isNull());
  if(m_verbosity > 5){
    pout() << "FieldSolver::setDefaultEbBcFunctions" << endl;
  }

  // TLDR: This sets the default potential function on the electrodes. We use a lambda for passing in
  //       some member variables (m_voltage, m_time) and get other parameters directly from the electrodes.
  //       We create one such voltage function for each electrode so that the voltage can be obtained
  //       anywhere on the electrode.

  const Vector<Electrode> electrodes = m_computationalGeometry->getElectrodes();

  for (int i = 0; i < electrodes.size(); i++){
    const Electrode& elec = electrodes[i];

    const Real val  = (elec.isLive()) ? 1.0 : 0.0;
    const Real frac = elec.getFraction();
    
    ElectrostaticEbBc::BcFunction curFunc = [&, &time=this->m_time, &voltage=this->m_voltage, val, frac](const RealVect a_position, const Real a_time){
      return voltage(time)*val*frac;
    };

    m_ebBc.addEbBc(elec, curFunc);
  }
}

void FieldSolver::setDomainSideBcFunction(const int a_dir, const Side::LoHiSide a_side, const ElectrostaticDomainBc::BcFunction& a_function){
  CH_TIME("FieldSolver::setDomainSideBcFunction");
  if(m_verbosity > 4){
    pout() << "FieldSolver::setDomainSideBcFunction" << endl;
  }

  const ElectrostaticDomainBc::DomainSide domainSide = std::make_pair(a_dir, a_side);

  m_domainBcFunctions.at(domainSide) = a_function;
}

void FieldSolver::parseDomainBc(){
  CH_TIME("FieldSolver::parseDomainBc");
  if(m_verbosity > 5){
    pout() << "FieldSolver::parseDomainBc" << endl;
  }

  // TLDR: This routine might seem big and complicated. What we are doing is that we are creating one function object which returns some value
  //       anywhere in space and time on a domain edge (face). The FieldSolver class supports Dirichlet and Neumann, and the below code simply
  //       creates those functions and associates them with an edge.
  //
  //       For flexibility we want to be able to specify the potential directy without invoking m_voltage, while at the same time we want to offer
  //       the simplistic method of setting a domain side to be "grounded", "live", or otherwise given by some fraction of m_voltage. 
  //       We thus make a distinction between "dirichlet" and "dirichlet_custom". The difference between these is that for "dirichlet_custom" the contents
  //       of m_domainBcFunctions are used as boundary conditions. For "dirichlet 0.5" the contents of m_domainBcFunctions are multiplied by 0.5*m_voltage. 
  //       The same approach is used for Neumann boundary conditions. 
  
  ParmParse pp(m_className.c_str());

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){

      const ElectrostaticDomainBc::DomainSide domainSide = std::make_pair(dir, sit());
      const std::string bcString = this->makeBcString(dir, sit());
      const int num = pp.countval(bcString.c_str());

      std::string str;

      ElectrostaticDomainBc::BcType      bcType;
      ElectrostaticDomainBc::BcFunction& bcFunc = m_domainBcFunctions.at(domainSide); // = F(x,t) in the comments below

      std::function<Real(const RealVect, const Real)> curFunc;

      if(num == 1){ // If we only had one argument the user has asked for "custom" boundary conditions, and we then pass in F(x,t) directly (evaluated at m_time)
	pp.get(bcString.c_str(), str, 0);

	curFunc = [&bcFunc, &time = this->m_time] (const RealVect a_pos, const Real a_time) {
	  return bcFunc(a_pos, time);
	};

	if(str == "dirichlet_custom"){
	  bcType  = ElectrostaticDomainBc::BcType::Dirichlet;
	}
	else if(str == "neumann_custom"){
	  bcType  = ElectrostaticDomainBc::BcType::Neumann;
	}
	else{
	  MayDay::Error("FieldSolver::parseDomainBc -- got only one argument but this argument was not dirichlet/neumann_custom. Maybe you have the wrong BC specification?");
	}
      }
      else if(num == 2){ // If we had two arguments the user has asked to run with less verbose specifications. E.g. "dirichlet 0.5" => V(x,t) = V(t) * 0.5 * F(x,t) 
	Real val;

	pp.get(bcString.c_str(), str, 0);
	pp.get(bcString.c_str(), val, 1);

	bcType = this->parseBcString(str);

	// Build a function computing the value at the boundary. 
	switch (bcType){
	case ElectrostaticDomainBc::BcType::Dirichlet:
	  curFunc = [&bcFunc, &voltage=this->m_voltage, &time=this->m_time, val] (const RealVect a_pos, const Real a_time){
	    return bcFunc(a_pos, time)*voltage(time)*val;
	  };
	  break;
	case ElectrostaticDomainBc::BcType::Neumann:
	  curFunc = [&bcFunc, &time=this->m_time, val] (const RealVect a_pos, const Real a_time){
	    return bcFunc(a_pos, time)*val;
	  };
	  break;
	default:
	  MayDay::Error("FieldSolver::parseDomainBc -- unsupported boundary condition requested!");
	  break;
	}
      }
      else{
	const std::string errorString = "FieldSolver::parseDomainBc -- bad or no input parameter for " + bcString;
	MayDay::Error(errorString.c_str());
      }

      m_domainBc.setBc(domainSide, std::make_pair(bcType, curFunc));
    }
  }
}

void FieldSolver::setPermittivities(){
  CH_TIME("FieldSolver::setPermittivities");
  if(m_verbosity > 5){
    pout() << "FieldSolver::setPermittivities" << endl;
  }

  const Real relEpsGas = m_computationalGeometry->getGasPermittivity();

  DataOps::setValue(m_permittivityCell, relEpsGas);
  DataOps::setValue(m_permittivityFace, relEpsGas); 
  DataOps::setValue(m_permittivityEB,   relEpsGas);

  const Vector<Dielectric>& dielectrics = m_computationalGeometry->getDielectrics();

  if(dielectrics.size() > 0 && m_multifluidIndexSpace->numPhases() > 1){
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const Real dx                = m_amr->getDx()[lvl];
      const RealVect probLo        = m_amr->getProbLo();

      LevelData<EBCellFAB>        cellPerm;
      LevelData<EBFluxFAB>        facePerm;
      LevelData<BaseIVFAB<Real> > ebPerm  ;

      MultifluidAlias::aliasMF(cellPerm, phase::solid, *m_permittivityCell[lvl]);
      MultifluidAlias::aliasMF(facePerm, phase::solid, *m_permittivityFace[lvl]);
      MultifluidAlias::aliasMF(ebPerm,   phase::solid, *m_permittivityEB   [lvl]);

      for (DataIterator dit(dbl); dit.ok(); ++dit){
	EBCellFAB& cellPermFAB     = cellPerm[dit()];
	EBFluxFAB& facePermFAB     = facePerm[dit()];
	BaseIVFAB<Real>& ebPermFAB = ebPerm  [dit()];

	
	const Box cellBox          = dbl.get(dit());
	const EBISBox& ebisbox     = cellPermFAB.getEBISBox();	

	this->setCellPermittivities(cellPermFAB,  cellBox, ebisbox, probLo, dx);
	this->setFacePermittivities(facePermFAB,  cellBox, ebisbox, probLo, dx);
	this->setEbPermittivities  (ebPermFAB,    cellBox, ebisbox, probLo, dx);
      }
    }
  }
}

void FieldSolver::setCellPermittivities(EBCellFAB&      a_relPerm,
					const Box&      a_cellBox,
					const EBISBox&  a_ebisbox,					
					const RealVect& a_probLo,
					const Real&     a_dx){
  CH_TIME("FieldSolver::setCellPermittivities");
  if(m_verbosity > 10){
    pout() << "FieldSolver::setCellPermittivities" << endl;
  }

  const EBGraph& ebgraph = a_ebisbox.getEBGraph();
  const IntVectSet irreg = a_ebisbox.getIrregIVS(a_cellBox);
  
  // Regular cells
  BaseFab<Real>& relPermFAB = a_relPerm.getSingleValuedFAB();
  for (BoxIterator bit(a_cellBox); bit.ok(); ++bit){
    const IntVect iv = bit();
    if(a_ebisbox.isRegular(iv)){
      const RealVect pos = a_probLo + a_dx*RealVect(iv) + 0.5*RealVect::Unit;
      relPermFAB(iv, m_comp) = this->getDielectricPermittivity(pos);
    }
  }

  // Irregular cells
  for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect pos  = Location::position(m_dataLocation, vof, a_ebisbox, a_dx);
    a_relPerm(vof, m_comp) = this->getDielectricPermittivity(pos);
  }
}

void FieldSolver::setFacePermittivities(EBFluxFAB&      a_relPerm,
					const Box&      a_cellBox,
					const EBISBox&  a_ebisbox,
					const RealVect& a_probLo,
					const Real&     a_dx){
  CH_TIME("FieldSolver::setFacePermittivities");
  if(m_verbosity > 10){
    pout() << "FieldSolver::setFacePermittivities" << endl;
  }

  const EBGraph& ebgraph         = a_ebisbox.getEBGraph();
  const IntVectSet irreg         = a_ebisbox.getIrregIVS(a_cellBox);
  FaceStop::WhichFaces stop_crit = FaceStop::SurroundingWithBoundary;
  
  for (int dir = 0; dir < SpaceDim; dir++){

    // Regular faces
    Box facebox = a_cellBox;
    facebox.surroundingNodes(dir);
    BaseFab<Real>& relPermFAB = a_relPerm[dir].getSingleValuedFAB();
    for (BoxIterator bit(facebox); bit.ok(); ++bit){
      const IntVect iv   = bit();
      const RealVect pos = a_probLo + a_dx*(RealVect(iv) + 0.5*RealVect::Unit) - 0.5*a_dx*BASISREALV(dir);
      
      relPermFAB(iv, m_comp) = this->getDielectricPermittivity(pos);
    }
    
    // Irregular faces
    for (FaceIterator faceit(irreg, ebgraph, dir, stop_crit); faceit.ok(); ++faceit){
      const FaceIndex& face  = faceit();
      const RealVect pos     = a_probLo + Location::position(m_faceLocation, face, a_ebisbox, a_dx);

      a_relPerm[dir](face, m_comp) = this->getDielectricPermittivity(pos);
    }
  }
}

void FieldSolver::setEbPermittivities(BaseIVFAB<Real>& a_relPerm,
				      const Box&       a_cellBox,
				      const EBISBox&   a_ebisbox,
				      const RealVect&  a_probLo,
				      const Real&      a_dx){
  CH_TIME("FieldSolver::setEbPermittivities");
  if(m_verbosity > 10){
    pout() << "FieldSolver::setEbPermittivities" << endl;
  }
  
  const int comp         = 0;
  const IntVectSet& ivs  = a_relPerm.getIVS();
  const EBGraph& ebgraph = a_relPerm.getEBGraph();

  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();
    const RealVect pos  = Location::position(Location::Cell::Boundary, vof, a_ebisbox, a_dx);

    a_relPerm(vof, m_comp) = this->getDielectricPermittivity(pos);
  }
}

void FieldSolver::writePlotFile(){
  CH_TIME("FieldSolver::writePlotFile");
  if(m_verbosity > 5){
    pout() << "FieldSolver::writePlotFile" << endl;
  }

  // Number of output components
  const int numPlotVars = this->getNumberOfPlotVariables();

  // Get plot variable names
  const Vector<std::string> plotVarNames = getPlotVariableNames();

  // Allocate storage for output
  EBAMRCellData output;
  m_amr->allocate(output, m_realm, phase::gas, numPlotVars);

  // Copy internal data to be plotted over to 'output'
  int icomp = 0; 
  this->writePlotData(output, icomp, true);

  // Filename
  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_className.c_str(), m_timeStep, SpaceDim);

  // Need to alias because EBHDF5 is not too smart. 
  Vector<LevelData<EBCellFAB>* > outputPtr(1+m_amr->getFinestLevel());
  m_amr->alias(outputPtr, output);

  Vector<Real> covered_values(numPlotVars, 0.0);
  string fname(file_char);
  writeEBHDF5(fname,
	      m_amr->getGrids(m_realm),
	      outputPtr,
	      plotVarNames,
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

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Used for aliasing phases
  LevelData<EBCellFAB> potentialGas;
  LevelData<EBCellFAB> potentialSol;

  if(!ebisGas.isNull()) MultifluidAlias::aliasMF(potentialGas,  phase::gas,   *m_potential[a_level]);
  if(!ebisSol.isNull()) MultifluidAlias::aliasMF(potentialSol,  phase::solid, *m_potential[a_level]);
  
  // Write data
  if(!ebisGas.isNull()) write(a_handle, potentialGas, "poisson_g");
  if(!ebisSol.isNull()) write(a_handle, potentialSol, "poisson_s");
}

void FieldSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level){
  CH_TIME("FieldSolver::readCheckpointLevel");
  if(m_verbosity > 5){
    pout() << "FieldSolver::readCheckpointLevel" << endl;
  }

  const RefCountedPtr<EBIndexSpace> ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Used for aliasing phases
  LevelData<EBCellFAB> potentialGas;
  LevelData<EBCellFAB> potentialSol;

  if(!ebisGas.isNull()) MultifluidAlias::aliasMF(potentialGas,  phase::gas,   *m_potential[a_level]);
  if(!ebisSol.isNull()) MultifluidAlias::aliasMF(potentialSol,  phase::solid, *m_potential[a_level]);
  
  // Read data
  if(!ebisGas.isNull()) read<EBCellFAB>(a_handle, potentialGas, "poisson_g", m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);
  if(!ebisSol.isNull()) read<EBCellFAB>(a_handle, potentialSol, "poisson_s", m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);
}

void FieldSolver::postCheckpoint(){
  CH_TIME("FieldSolver::postCheckpoint");
  if(m_verbosity > 5){
    pout() << "FieldSolver::postCheckpoint" << endl;
  }
}

void FieldSolver::writePlotData(EBAMRCellData& a_output, int& a_comp, const bool a_forceNoInterp){
  CH_TIME("FieldSolver::writePlotData");
  if(m_verbosity > 5){
    pout() << "FieldSolver::writePlotData" << endl;
  }

  // This routine always outputs data on the centroid. If the data was defined on the center we move it to the centroid (but not forcefully)
  const bool doInterp = (m_dataLocation == Location::Cell::Center) && !a_forceNoInterp;

  // Add phi to output
  if(m_plotPotential)     this->writeMultifluidData(a_output, a_comp, m_potential,     doInterp); // Possibly on cell center, so-recenter to centroid if needed
  if(m_plotRho)           this->writeMultifluidData(a_output, a_comp, m_rho,           false);    // Always centroid
  if(m_plotResidue)       this->writeMultifluidData(a_output, a_comp, m_residue,       false);    // Always centroid
  if(m_plotElectricField) this->writeMultifluidData(a_output, a_comp, m_electricField, doInterp); // Possibly on cell center, so-recenter to centroid if needed
}

void FieldSolver::writeMultifluidData(EBAMRCellData& a_output, int& a_comp, const MFAMRCellData& a_data, const bool a_interp){
  CH_TIME("FieldSolver::writeMultifluidData");
  if(m_verbosity > 5){
    pout() << "FieldSolver::writeMultifluidData" << endl;
  }

  // So the problem with the Chombo HDF5 I/O routines is that they are not really designed for multiphase. 

  const int numComp = a_data[0]->nComp();

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Allocate some scratch data that we can use
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, phase::gas, numComp);

  //
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    LevelData<EBCellFAB> gasData;
    LevelData<EBCellFAB> solData;

    if(!ebisGas.isNull()) MultifluidAlias::aliasMF(gasData,  phase::gas,   *a_data[lvl]);
    if(!ebisSol.isNull()) MultifluidAlias::aliasMF(solData,  phase::solid, *a_data[lvl]);

    gasData.localCopyTo(*scratch[lvl]);

    // Copy all covered cells from the other phase
    if(!ebisSol.isNull()){
      for (DataIterator dit = solData.dataIterator(); dit.ok(); ++dit){
	const Box box = solData.disjointBoxLayout().get(dit());
	const IntVectSet ivs(box);
	const EBISBox& ebisboxGas = gasData[dit()].getEBISBox();
	const EBISBox& ebisboxSol = solData[dit()].getEBISBox();

	FArrayBox& scratchGas    = (*scratch[lvl])[dit()].getFArrayBox();
	const FArrayBox& fabGas = gasData[dit()].getFArrayBox();
	const FArrayBox& fabSol = solData[dit()].getFArrayBox();

	// TLDR: There are four cases
	// 1. Both are covered => inside electrode
	// 2. Gas is covered, solid is regular => inside solid phase only
	// 3. Gas is regular, solid is covered => inside gas phase only
	// 4. Gas is !(covered || regular) and solid is !(covered || regular) => on dielectric boundary
	if(!(ebisboxGas.isAllCovered() && ebisboxSol.isAllCovered())){ // Outside electrode, lets do stuff
	  if(ebisboxGas.isAllCovered() && ebisboxSol.isAllRegular()){ // Case 2, copy directly. 
	    scratchGas += fabSol;
	  }
	  else if(ebisboxGas.isAllRegular() && ebisboxSol.isAllCovered()) { // Case 3. Inside gas phase. Already did this. 
	  }
	  else{ // Case 4, needs special treatment. 
	    for (BoxIterator bit(box); bit.ok(); ++bit){ // Loop through all cells here
	      const IntVect iv = bit();
	      if(ebisboxGas.isCovered(iv) && !ebisboxSol.isCovered(iv)){   // Regular cells from phase 2
		for (int comp = 0; comp < numComp; comp++){
		  scratchGas(iv, comp) = fabSol(iv, comp);
		}
	      }
	      else if(ebisboxSol.isIrregular(iv) && ebisboxGas.isIrregular(iv)){ // Irregular cells. Use gas side data
		for (int comp = 0; comp < numComp; comp++){
		  scratchGas(iv, comp) = fabGas(iv,comp);
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

  const Interval src_interv(0, numComp-1);
  const Interval dst_interv(a_comp, a_comp + numComp - 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    if(m_realm == a_output.getRealm()){
      scratch[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else{
      scratch[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }

  a_comp += numComp;
}

int FieldSolver::getNumberOfPlotVariables() const {
  CH_TIME("FieldSolver::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << "FieldSolver::getNumberOfPlotVariables" << endl;
  }

  int numPltVars = 0;

  if(m_plotPotential)     numPltVars = numPltVars + 1;
  if(m_plotRho)           numPltVars = numPltVars + 1;
  if(m_plotResidue)       numPltVars = numPltVars + 1;
  if(m_plotElectricField) numPltVars = numPltVars + SpaceDim;

  return numPltVars;
}

Vector<std::string> FieldSolver::getPlotVariableNames() const {
  CH_TIME("FieldSolver::getPlotVariableNames");
  if(m_verbosity > 5){
    pout() << "FieldSolver::getPlotVariableNames" << endl;
  }
  
  Vector<std::string> pltVarNames(0);
  
  if(m_plotPotential) pltVarNames.push_back("Electrostatic potential");
  if(m_plotRho)       pltVarNames.push_back("Space charge density");
  if(m_plotResidue)   pltVarNames.push_back("Electrostatic potential_residue");
  
  if(m_plotElectricField){
    pltVarNames.push_back("x-Electric field"); 
    pltVarNames.push_back("y-Electric field"); 
    if(SpaceDim == 3){
      pltVarNames.push_back("z-Electric field");
    }
  }
  
  return pltVarNames;
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

const std::function<Real(const Real a_time)>& FieldSolver::getVoltage() const{
  return m_voltage;
}

Real FieldSolver::getCurrentVoltage() const{
  return m_voltage(m_time);
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

MFAMRCellData& FieldSolver::getPermittivityCell(){
  return m_permittivityCell;
}

MFAMRFluxData& FieldSolver::getPermittivityFace(){
  return m_permittivityFace;
}

MFAMRIVData& FieldSolver::getPermittivityEB(){
  return m_permittivityEB;
}

#include <CD_NamespaceFooter.H>
