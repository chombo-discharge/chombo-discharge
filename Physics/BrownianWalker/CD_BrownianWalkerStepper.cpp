/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BrownianWalkerStepper.cpp
  @brief  Implementation of CD_BrownianWalkerStepper.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <PolyGeom.H>
#include <BinFab.H>

// Our includes
#include <CD_BrownianWalkerStepper.H>
#include <CD_BrownianWalkerSpecies.H>
#include <CD_PolyUtils.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::BrownianWalker;

BrownianWalkerStepper::BrownianWalkerStepper(){
  ParmParse pp("BrownianWalker");

  m_phase = phase::gas;

  pp.get("realm",          m_realm);
  pp.get("diffco",         m_faceCenteredDiffusionCoefficient);
  pp.get("omega",          m_omega);
  pp.get("verbosity",      m_verbosity);
  pp.get("ppc",            m_ppc);
  pp.get("max_cells_hop",  m_max_cells_hop);
  pp.get("load_balance",   m_LoadBalancing);
}

BrownianWalkerStepper::BrownianWalkerStepper(RefCountedPtr<ItoSolver>& a_solver) : BrownianWalkerStepper() {
  m_solver = a_solver;
}

BrownianWalkerStepper::~BrownianWalkerStepper(){

}

void BrownianWalkerStepper::parseRuntimeOptions() {

  ParmParse pp("BrownianWalker");
  
  pp.get("verbosity",      m_verbosity);
  pp.get("ppc",            m_ppc);
  pp.get("max_cells_hop",  m_max_cells_hop);
  pp.get("load_balance",   m_LoadBalancing);
  
  m_solver->parseRuntimeOptions();
}

void BrownianWalkerStepper::initialData(){
  CH_TIME("BrownianWalkerStepper::initialData");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::initialData" << endl;
  }
  
  m_solver->initialData();
  if(m_ppc > 0){
    m_solver->sortParticlesByCell(ItoSolver::WhichContainer::Bulk);
    m_solver->makeSuperparticles(ItoSolver::WhichContainer::Bulk,m_ppc);
    m_solver->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  }

  if(m_solver->isDiffusive()){
    m_solver->setDiffusionFunction(m_faceCenteredDiffusionCoefficient);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }
}

void BrownianWalkerStepper::postInitialize(){
  CH_TIME("BrownianWalkerStepper::postInitialize");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::postInitialize" << endl;
  }
}

void BrownianWalkerStepper::setVelocity(){
  CH_TIME("BrownianWalkerStepper::setVelocity");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::setVelocity" << endl;
  }
  m_solver->setParticleMobility(1.0);
  
  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    this->setVelocity(lvl);
  }

  EBAMRCellData& vel = m_solver->getVelocityFunction();
  m_amr->averageDown(vel, m_realm, m_phase);
  m_amr->interpGhost(vel, m_realm, m_phase);
}

bool BrownianWalkerStepper::loadBalanceThisRealm(const std::string a_realm) const {
  CH_TIME("BrownianWalkerStepper::loadBalanceThisRealm");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::loadBalanceThisRealm" << endl;
  }

  bool ret = false;

  if(m_LoadBalancing && a_realm == m_realm){
    ret = true;
  }

  return ret;
}

void BrownianWalkerStepper::loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
					     Vector<Vector<Box> >&            a_boxes,
					     const std::string                a_realm,
					     const Vector<DisjointBoxLayout>& a_grids,
					     const int                        a_lmin,
					     const int                        a_finestLevel){
  CH_TIME("BrownianWalkerStepper::loadBalanceBoxes");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::loadBalanceBoxes" << endl;
  }
  
  if(m_LoadBalancing && a_realm == m_realm){
    ParticleContainer<ItoParticle>& particles = m_solver->getParticles(ItoSolver::WhichContainer::Bulk);
  
    particles.regrid(a_grids, m_amr->getDomains(), m_amr->getDx(), m_amr->getRefinementRatios(), a_lmin, a_finestLevel);

    a_procs.resize(1 + a_finestLevel);
    a_boxes.resize(1 + a_finestLevel);
  
    // Compute loads on each level
    for (int lvl = 0; lvl < a_lmin; lvl++){
      a_procs[lvl] = a_grids[lvl].procIDs();
      a_boxes[lvl] = a_grids[lvl].boxArray();
    }

    for (int lvl = a_lmin; lvl <= a_finestLevel; lvl++){
      Vector<long int> loads;
      a_boxes[lvl] = a_grids[lvl].boxArray();
    
      m_solver->computeLoads(loads, a_grids[lvl], lvl);

#ifdef CH_MPI
      int count = loads.size();
      Vector<long int> tmp(count);
      MPI_Allreduce(&(loads[0]),&(tmp[0]), count, MPI_LONG, MPI_SUM, Chombo_MPI::comm);
      loads = tmp;
#endif

      LoadBalance(a_procs[lvl], loads, a_boxes[lvl]);
    }

    // Put particles back
    particles.preRegrid(a_lmin);
  }
  else{
    MayDay::Abort("BrownianWalkerStepper::loadBalanceBoxes - logic bust");
  }
}

void BrownianWalkerStepper::setVelocity(const int a_level){
  CH_TIME("BrownianWalkerStepper::setVelocity(level)");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::setVelocity(level)" << endl;
  }

  // TLDR: This code goes down to each cell on grid level a_level and sets the velocity to omega*r
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box = dbl.get(dit());

    EBCellFAB& vel = (*(m_solver->getVelocityFunction())[a_level])[dit()];
    BaseFab<Real>& vel_reg = vel.getSingleValuedFAB();

    vel.setVal(0.0);
      
    // Regular cells
    for (BoxIterator bit(box); bit.ok(); ++bit){
      const IntVect iv = bit();
      const RealVect pos = m_amr->getProbLo() + (RealVect(iv) + 0.5*RealVect::Unit)*m_amr->getDx()[a_level];

      const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
      const Real theta = atan2(pos[1],pos[0]);

      vel_reg(iv,0) = -r*m_omega*sin(theta);
      vel_reg(iv,1) =  r*m_omega*cos(theta);
    }

    // Irregular and multicells
    const EBISBox& ebisbox = vel.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){

      const VolIndex vof = vofit();
      const IntVect iv   = vof.gridIndex();
      const RealVect pos = m_amr->getProbLo() + (RealVect(iv) + 0.5*RealVect::Unit)*m_amr->getDx()[a_level];

      const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
      const Real theta = atan2(pos[1],pos[0]);

      vel(vof,0) = -r*m_omega*sin(theta);
      vel(vof,1) =  r*m_omega*cos(theta);
    }

    // Now set the mobility for all the particles
    List<ItoParticle>& particles = m_solver->getParticles(ItoSolver::WhichContainer::Bulk)[a_level][dit()].listItems();
    for (ListIterator<ItoParticle> lit(particles); lit.ok(); ++lit){
      lit().mobility() = 1.0;
    }
  }
}

#ifdef CH_USE_HDF5
void BrownianWalkerStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const{
  CH_TIME("BrownianWalkerStepper::writeCheckpointData");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::writeCheckpointData" << endl;
  }

  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
void BrownianWalkerStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl) {
  CH_TIME("BrownianWalkerStepper::readCheckpointData");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::readCheckpointData" << endl;
  }
  
  m_solver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

void BrownianWalkerStepper::postCheckpointSetup() {
  CH_TIME("BrownianWalkerStepper::postCheckpointSetup");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::postCheckpointSetup" << endl;
  }

  m_solver->remap();
  if(m_ppc > 0){
    m_solver->sortParticlesByCell(ItoSolver::WhichContainer::Bulk);
    m_solver->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_ppc);
    m_solver->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  }
  m_solver->depositParticles();

  if(m_solver->isDiffusive()){
    m_solver->setDiffusionFunction(m_faceCenteredDiffusionCoefficient);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }
}

int BrownianWalkerStepper::getNumberOfPlotVariables() const {
  CH_TIME("BrownianWalkerStepper::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::getNumberOfPlotVariables" << endl;
  }

  return m_solver->getNumberOfPlotVariables();
}

void BrownianWalkerStepper::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {
  CH_TIME("BrownianWalkerStepper::writePlotData");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::writePlotData" << endl;
  }
  a_plotVariableNames.append(m_solver->getPlotVariableNames());
  m_solver->writePlotData(a_output, a_icomp);
}

void BrownianWalkerStepper::computeDt(Real& a_dt, TimeCode& a_timeCode) {
  CH_TIME("BrownianWalkerStepper::computeDt");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::computeDt" << endl;
  }

  m_solver->setParticleMobility(1.0); 
  m_solver->interpolateVelocities();
  m_solver->interpolateDiffusion();

  a_dt = m_max_cells_hop*m_solver->computeDt();
}

void BrownianWalkerStepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("BrownianWalkerStepper::synchronizeSolverTimes");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::synchronizeSolverTimes" << endl;
  }
  
  m_solver->setTime(a_step, a_time, a_dt);

  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void BrownianWalkerStepper::printStepReport() {
  CH_TIME("BrownianWalkerStepper::printStepReport");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::printStepReport" << endl;
  }

  // Do nothing
  const size_t local_particles  = m_solver->getNumParticles(ItoSolver::WhichContainer::Bulk, true);
  const size_t global_particles = m_solver->getNumParticles(ItoSolver::WhichContainer::Bulk, false);

  pout() << "                                   #part = " << local_particles << " (" << global_particles << ")" << endl;
  
}

bool BrownianWalkerStepper::needToRegrid() {
  return false;
}

void BrownianWalkerStepper::preRegrid(const int a_lbase, const int a_oldFinestLevel){
  CH_TIME("BrownianWalkerStepper::preRegrid");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::preRegrid" << endl;
  }

  // TLDR: base is the finest level that DOES NOT CHANGE
  const int base = 0;
  const int finest_level = m_amr->getFinestLevel();
  
  m_solver->preRegrid(base, finest_level);
}

void BrownianWalkerStepper::setupSolvers() {
  CH_TIME("BrownianWalkerStepper::setupSolvers");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::setupSolvers" << endl;
  }

  m_species = RefCountedPtr<ItoSpecies> (new BrownianWalkerSpecies());

  m_solver->setVerbosity(m_verbosity);
  m_solver->parseOptions();
  m_solver->setAmr(m_amr);
  m_solver->setSpecies(m_species);
  m_solver->setPhase(m_phase);
  m_solver->setComputationalGeometry(m_computationalGeometry);
  m_solver->setRealm(m_realm);
}

void BrownianWalkerStepper::registerRealms() {
  CH_TIME("BrownianWalkerStepper::registerRealms");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_realm);
}

void BrownianWalkerStepper::registerOperators() {
  CH_TIME("BrownianWalkerStepper::registerOperators");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::registerOperators" << endl;
  }

  m_solver->registerOperators();
}

void BrownianWalkerStepper::allocate() {
  m_solver->allocateInternals(); // Allocate some internal storage
}

Real BrownianWalkerStepper::advance(const Real a_dt) {
  CH_TIME("BrownianWalkerStepper::advance");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::advance" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();
  const RealVect origin  = m_amr->getProbLo();
  
  //  m_solver->remap();

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const RealVect dx                     = m_amr->getDx()[lvl]*RealVect::Unit;
    const DisjointBoxLayout& dbl          = m_amr->getGrids(m_realm)[lvl];
    ParticleData<ItoParticle>& particles = m_solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

    const EBISLayout& ebisl = m_amr->getEBISLayout(m_realm, m_solver->getPhase())[lvl];

    if(m_solver->isMobile() || m_solver->isDiffusive()){
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){

	// Create a copy. 
	List<ItoParticle>& particleList = particles[dit()].listItems();
	List<ItoParticle>  particleCopy = List<ItoParticle>(particleList);

	// The list iterator is NOT an indexing iterator but iterates over the list given
	// in the constructor. So, we need one for velocities and one for the copy
	ListIterator<ItoParticle> lit(particleList);
	ListIterator<ItoParticle> litC(particleCopy);

	// Half Euler step and evaluate velocity at half step
	if(m_solver->isMobile()){
	  for (lit.rewind(); lit; ++lit){ 
	    ItoParticle& p = particleList[lit];
	    p.position() += 0.5*p.velocity()*a_dt;
	  }

	  m_solver->interpolateVelocities(lvl, dit());

	  // Final stage
	  for (lit.rewind(), litC.rewind(); lit, litC; ++lit, ++litC){
	    ItoParticle& p    = particleList[lit];
	    ItoParticle& oldP = particleCopy[litC];
	    p.position() = oldP.position() + p.velocity()*a_dt;
	  }
	}

	// Add a diffusion hop
	if(m_solver->isDiffusive()){
	  for (lit.rewind(); lit; ++lit){ 
	    ItoParticle& p = particleList[lit];
	    const RealVect ran = m_solver->randomGaussian();
	    const RealVect hop = ran*sqrt(2.0*p.diffusion()*a_dt);
	    p.position() += hop;
	  }
	}

	// Do particle bounceback on the EB
	const EBISBox& ebisbox = ebisl[dit()];
	if(!ebisbox.isAllRegular() || !ebisbox.isAllCovered()){ 
	  
	  // This is the implicit function
	  RefCountedPtr<BaseIF> func;
	  if(m_solver->getPhase() == phase::gas){
	    func = m_computationalGeometry->getGasImplicitFunction();
	  }
	  else {
	    func = m_computationalGeometry->getSolidImplicitFunction();
	  }

	  for (lit.rewind(), litC.rewind(); lit, litC; ++lit, ++litC){
	    ItoParticle& p    = particleList[lit];
	    ItoParticle& oldP = particleCopy[litC];

	    const RealVect& oldPos = oldP.position();
	    const RealVect& newPos = p.position();

	    const Real fOld = func->value(oldPos);
	    const Real fNew = func->value(newPos);

	    // If the particle crossed the boundary, 
	    if(fOld*fNew <= 0.0){
	      const RealVect xb = PolyUtils::brentRootFinder(func, oldPos, newPos);
	      const IntVect iv = locateBin(xb, dx, origin);
	      const VolIndex vof(iv, 0);

	      // Do a dummy check
	      //	      const Box b(iv-IntVect::Unit, iv+IntVect::Unit);
	      if(!ebisbox.isIrregular(iv)){
		MayDay::Abort("BrownianWalkerStepper::advance - logic bust in EB bounce-back code");
	      }
	      else{
		const RealVect n     = ebisbox.normal(vof);
		const RealVect bback = 2.0*n*PolyGeom::dot((newPos - xb), n);
		p.position() -= bback;
	      }
	    }
	  }
	}
      }
    }
  }

  // Remap and deposit particles
  m_solver->remap();
  m_solver->depositParticles();

  return a_dt;
}

void BrownianWalkerStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  CH_TIME("BrownianWalkerStepper::regrid");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::regrid" << endl;
  }

  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  if(m_solver->isDiffusive()){
    m_solver->setDiffusionFunction(m_faceCenteredDiffusionCoefficient);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }

  if(m_ppc > 0){
    m_solver->sortParticlesByCell(ItoSolver::WhichContainer::Bulk);
    m_solver->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_ppc);
    m_solver->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk);
    
    m_solver->setParticleMobility(1.0); // Superparticle algorithm only conserves mass, energy. Diffusion and mobility needs to be reset.
  }
}

void BrownianWalkerStepper::postRegrid(){
  CH_TIME("BrownianWalkerStepper::postRegrid");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::postRegrid" << endl;
  }
}

#include <CD_NamespaceFooter.H>
