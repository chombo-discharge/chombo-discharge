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
#include <CH_Timer.H>
#include <BinFab.H>

// Our includes
#include <CD_BrownianWalkerStepper.H>
#include <CD_BrownianWalkerSpecies.H>
#include <CD_PolyUtils.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::BrownianWalker;

BrownianWalkerStepper::BrownianWalkerStepper(){
  CH_TIME("BrownianWalkerStepper::BrownianWalkerStepper");
  
  ParmParse pp("BrownianWalker");

  m_phase = phase::gas;

  pp.get("realm",        m_realm);
  pp.get("diffco",       m_diffCo);
  pp.get("mobility",     m_mobility);  
  pp.get("omega",        m_omega);
  pp.get("verbosity",    m_verbosity);
  pp.get("ppc",          m_ppc);
  pp.get("cfl",          m_cfl);
  pp.get("load_balance", m_loadBalance);
}

BrownianWalkerStepper::BrownianWalkerStepper(RefCountedPtr<ItoSolver>& a_solver) : BrownianWalkerStepper() {
  CH_TIME("BrownianWalkerStepper::BrownianWalkerStepper(full)");

  CH_assert(!a_solver.isNull());
  
  m_solver = a_solver;
}

BrownianWalkerStepper::~BrownianWalkerStepper(){
  CH_TIME("BrownianWalkerStepper::~BrownianWalkerStepper");
}

void BrownianWalkerStepper::parseRuntimeOptions() {
  CH_TIME("BrownianWalkerStepper::parseRuntimeOptions");
  
  ParmParse pp("BrownianWalker");
  
  pp.get("verbosity",    m_verbosity);
  pp.get("ppc",          m_ppc);
  pp.get("cfl",          m_cfl);
  pp.get("load_balance", m_loadBalance);
  
  m_solver->parseRuntimeOptions();
}

void BrownianWalkerStepper::initialData(){
  CH_TIME("BrownianWalkerStepper::initialData");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::initialData" << endl;
  }

  // Fill initial particles and then make the desired number of superparticles. 
  m_solver->initialData();
  this->makeSuperParticles();
}

void BrownianWalkerStepper::postInitialize(){
  CH_TIME("BrownianWalkerStepper::postInitialize");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::postInitialize" << endl;
  }

  // Set advection and diffusion fields.
  this->setAdvectionDiffusion();

  // Set particle diffusion coefficient and mobility
  m_solver->setParticleDiffusion(m_diffCo  );
  m_solver->setParticleMobility (m_mobility);

  m_solver->interpolateVelocities();  
}

void BrownianWalkerStepper::setAdvectionDiffusion(){
  CH_TIME("BrownianWalkerStepper::setAdvectionDiffusion");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::setAdvectionDiffusion" << endl;
  }

  if(m_solver->isDiffusive()){
    this->setDiffusion();
  }

  if(m_solver->isMobile()){
    this->setVelocity();
  }
}

void BrownianWalkerStepper::setDiffusion(){
  CH_TIME("BrownianWalkerStepper::setDiffusion");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::setDiffusion" << endl;
  }

  CH_assert(m_solver->isDiffusive());

  // Set something crazy for the diffusion field. This should not matter because we set the particle diffusion coefficients directly. 
  m_solver->setDiffusionFunction(std::numeric_limits<Real>::max());  
}

void BrownianWalkerStepper::setVelocity(){
  CH_TIME("BrownianWalkerStepper::setVelocity");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::setVelocity" << endl;
  }

  CH_assert(m_solver->isMobile());

  // TLDR: This just sets the velocity field everywhere. 

  // Velocity field in solver
  EBAMRCellData& vel = m_solver->getVelocityFunction();  

  // Nifty lambda describing the advective field
  auto veloFunc = [omega=this->m_omega](const RealVect pos) -> RealVect {
    const Real r     = pos.vectorLength();
    const Real theta = atan2(pos[1], pos[0]);

    return RealVect(D_DECL(-r*omega*sin(theta), r*omega*cos(theta), 0.));
  };

  DataOps::setValue(vel, 0.0);
  DataOps::setValue(vel, veloFunc, m_amr->getProbLo(), m_amr->getDx());

  // Coarsen and update ghost cells.
  m_amr->averageDown  (vel, m_realm, m_phase);
  m_amr->interpGhostMG(vel, m_realm, m_phase);

  DataOps::setCoveredValue(vel, 0, 0.0);
}

bool BrownianWalkerStepper::loadBalanceThisRealm(const std::string a_realm) const {
  CH_TIME("BrownianWalkerStepper::loadBalanceThisRealm");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::loadBalanceThisRealm" << endl;
  }

  return m_loadBalance && (a_realm == m_realm);
}

void BrownianWalkerStepper::loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
					     Vector<Vector<Box> >&            a_boxes,
					     const std::string                a_realm,
					     const Vector<DisjointBoxLayout>& a_grids,
					     const int                        a_lmin,
					     const int                        a_finestLevel) {
  CH_TIME("BrownianWalkerStepper::loadBalanceBoxes");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::loadBalanceBoxes" << endl;
  }

  CH_assert(m_loadBalance && a_realm == m_realm);

  // TLDR: This routine is called AFTER AmrMesh::regridAMR which means that we have all EB-related information we need for building operators. We happen to
  //       know that ItoSolver computed the number of particles per cell in the preRegrid method and that these values are returned by a call to
  //       EBAMRCellData& ItoSolver::getScratch(). We take that data and regrid it onto the new grids. This requires us to manually build an operator which
  //       can regrid that data.

  constexpr int comp     = 0;
  constexpr int numComps = 1;

  // This is the number of particles per cell, but it is stored on the old grids. 
  const EBAMRCellData& oldParticlesPerCell = m_solver->getScratch();

  // Make some storage for the number of particles per cell on the new grids.
  EBAMRCellData newParticlesPerCell;
  m_amr->allocate(newParticlesPerCell, a_realm, m_phase, 1);

  // These are the EB layouts.
  const Vector<ProblemDomain>& domains = m_amr->getDomains();
  const Vector<EBISLayout>&    ebisl   = m_amr->getEBISLayout(a_realm, m_phase);
  const Vector<int>&           refRat  = m_amr->getRefinementRatios();

  for (int lvl = 0; lvl <= std::max(0,a_lmin-1); lvl++){
    oldParticlesPerCell[lvl]->copyTo(*newParticlesPerCell[lvl]);
  }

  // Regrid onto the new mesh
  for (int lvl = 0; lvl <= a_finestLevel; lvl++){

    const bool hasCoar = lvl > 0;

    if(hasCoar){
      EBPWLFineInterp fineInterp(a_grids[lvl],
				 a_grids[lvl-1],
				 ebisl[lvl],
				 ebisl[lvl-1],
				 domains[lvl-1],
				 refRat[lvl-1],
				 numComps,
				 ebisl[lvl].getEBIS());

      fineInterp.interpolate(*newParticlesPerCell[lvl], *newParticlesPerCell[lvl-1], Interval(0,0));

      if(lvl < std::min(newParticlesPerCell.size(), oldParticlesPerCell.size())){
	oldParticlesPerCell[lvl]->copyTo(*newParticlesPerCell[lvl]);
      }
    }
  }

  // At this point we need to replace the data UNDERNEATH the fine grids. This might seem weird but recall that we don't really have control
  // over what exists on the invalid regions on the coarse grids. Our simple way is just to call DataOps and have it set the invalid data to zero, and the same
  // with the covered data.
  DataOps::setCoveredData(nweParticlesPerCell, 0.0);

  // 
  
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

  // Update particle distribution
  m_solver->remap();
  this->makeSuperParticles();

  // Update advection-diffusion fields
  this->setAdvectionDiffusion();

  // Set particle diffusion coefficient and mobility
  m_solver->setParticleDiffusion(m_diffCo  );
  m_solver->setParticleMobility (m_mobility);

  // Interpolate particle velocities. 
  m_solver->interpolateVelocities();      
}

int BrownianWalkerStepper::getNumberOfPlotVariables() const {
  CH_TIME("BrownianWalkerStepper::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::getNumberOfPlotVariables" << endl;
  }

  const int numPlotVars = m_solver->getNumberOfPlotVariables();

  return numPlotVars;
}

void BrownianWalkerStepper::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {
  CH_TIME("BrownianWalkerStepper::writePlotData");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::writePlotData" << endl;
  }
  
  // We only write solver data -- it knows what to do. 
  a_plotVariableNames.append(m_solver->getPlotVariableNames());
  m_solver->writePlotData(a_output, a_icomp);
}

void BrownianWalkerStepper::computeDt(Real& a_dt, TimeCode& a_timeCode) {
  CH_TIME("BrownianWalkerStepper::computeDt");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::computeDt" << endl;
  }

  a_dt = m_cfl*m_solver->computeDt();
}

void BrownianWalkerStepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("BrownianWalkerStepper::synchronizeSolverTimes");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::synchronizeSolverTimes" << endl;
  }

  // Solver needs to synchronize. 
  m_solver->setTime(a_step, a_time, a_dt);

  // TimeStepper needs to synchronize. 
  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;
}

void BrownianWalkerStepper::printStepReport() {
  CH_TIME("BrownianWalkerStepper::printStepReport");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::printStepReport" << endl;
  }

  // Do nothing
  const size_t localParticles  = m_solver->getNumParticles(ItoSolver::WhichContainer::Bulk, true);
  const size_t globalParticles = m_solver->getNumParticles(ItoSolver::WhichContainer::Bulk, false);

  pout() << "                                   #part = " << localParticles << " (" << globalParticles << ")" << endl;
}

bool BrownianWalkerStepper::needToRegrid() {
  CH_TIME("BrownianWalkerStepper::needToRegrid");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::needToRegrid" << endl;
  }
  
  return false;
}

void BrownianWalkerStepper::preRegrid(const int a_lbase, const int a_oldFinestLevel){
  CH_TIME("BrownianWalkerStepper::preRegrid");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::preRegrid" << endl;
  }

  // Solver knows what to do. 
  m_solver->preRegrid(a_lbase, a_oldFinestLevel);
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
  CH_TIME("BrownianWalkerStepper::allocate");

  // Allocate solver storage -- it knows what to do.
  m_solver->allocateInternals();
}

Real BrownianWalkerStepper::advance(const Real a_dt) {
  CH_TIME("BrownianWalkerStepper::advance");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::advance" << endl;
  }

  CH_assert(m_solver->isMobile() || m_solver->isDiffusive());

  // TLDR: This function advances the particles using an Euler-Maruyama kernel. The steps are simply:
  //
  //          1. Compute Xnew = Xold + V*dt + sqrt(2*D*dt)*N0 where N0 is a random number
  //          2. Remap the particles, assigning them to new grid boxes.
  //          3. Remove particles that struck the EB.
  //          4. Make super-particles. 
  //          5. Update the particle velocities and diffusion coefficients.
  //          6. Deposit particles on mesh.
  //

  // 1. Euler-Maruayma kernel on each patch. 
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_realm)[lvl];
    ParticleData<ItoParticle>& particles = m_solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];
    
    for (DataIterator dit(dbl); dit.ok(); ++dit){

      // Particles that we iterate through. 
      List<ItoParticle>& particleList = particles[dit()].listItems();

      // Euler step.
      if(m_solver->isMobile()){
	for (ListIterator<ItoParticle> lit(particleList); lit; ++lit){ 
	  ItoParticle& p  = lit();
	  p.oldPosition() = p.position();
	  p.position()   += p.velocity()*a_dt;
	}
      }

      // Diffusion hop
      if(m_solver->isDiffusive()){
	for (ListIterator<ItoParticle> lit(particleList); lit; ++lit){ 
	  ItoParticle& p      = lit();	    
	  const RealVect ran  = m_solver->randomGaussian();
	  const RealVect hop  = ran*sqrt(2.0*p.diffusion()*a_dt);
	  p.position()       += hop;
	}
      }
    }
  }

  // 2. Remap particles and assign them to correct patches. This discards particles outside the simulation domain. 
  m_solver->remap();

  // 3. Particles that strike the EB are absorbed on it, and removed from the simulation. 
  m_solver->removeCoveredParticles(EbRepresentation::ImplicitFunction, 0.0);

  // 4. Make new super-particles. 
  this->makeSuperParticles();

  // 5. Update particle diffusion and velocities. 
  m_solver->setParticleMobility(m_mobility);
  m_solver->setParticleDiffusion(m_diffCo);
  m_solver->interpolateVelocities();  

  // Deposit onto mesh. 
  m_solver->depositParticles();

  return a_dt;
}

void BrownianWalkerStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  CH_TIME("BrownianWalkerStepper::regrid");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::regrid" << endl;
  }

  // Solver regrids.
  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  // Make superparticles
  this->makeSuperParticles();
}

void BrownianWalkerStepper::postRegrid(){
  CH_TIME("BrownianWalkerStepper::postRegrid");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::postRegrid" << endl;
  }

  // Update advection-diffusion fields
  this->setAdvectionDiffusion();

  // Set particle diffusion coefficient and mobility
  m_solver->setParticleDiffusion(m_diffCo  );
  m_solver->setParticleMobility (m_mobility);

  // Interpolate particle velocities. 
  m_solver->interpolateVelocities();    
}

void BrownianWalkerStepper::makeSuperParticles() {
  CH_TIME("BrownianWalkerStepper::makeSuperParticles");
  if(m_verbosity > 5){
    pout() << "BrownianWalkerStepper::makeSuperParticles" << endl;
  }

  // TLDR: ItoSolver requires the particles to be sorted by cell when making superparticles. So we explicitly
  //       need to call cell/patch sorting methods. 

  if(m_ppc > 0){
    m_solver->sortParticlesByCell (ItoSolver::WhichContainer::Bulk       );
    m_solver->makeSuperparticles  (ItoSolver::WhichContainer::Bulk, m_ppc);
    m_solver->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk       );
  }
}

#include <CD_NamespaceFooter.H>
