/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_StreamerInceptionStepperImplem.H
  @brief  Implementation of CD_StreamerInceptionStepper.H
  @author Robert Marskar
*/

#ifndef CD_StreamerInceptionStepperImplem_H
#define CD_StreamerInceptionStepperImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_StreamerInceptionStepper.H>
#include <CD_NamespaceHeader.H>

// Tasks:
//
// 1. Compute inception integral with particle tracing
// 2. Compute critical volume
// 3. Compute negative ion drift times to electrode(s). Compute as time/(V * mobility)
// 4. Evaluate field emission rates.
// 5. Background ionization rates
using namespace Physics::StreamerInception;

template <typename P, typename F>
StreamerInceptionStepper<P, F>::StreamerInceptionStepper()
{
  CH_TIME("StreamerInceptionStepper::StreamerInceptionStepper");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::StreamerInceptionStepper" << endl;
  }

  m_realm = Realm::Primal;
  m_phase = phase::gas;

  this->parseOptions();

  m_verbosity = -1;

  ParmParse pp("StreamerInceptionStepper");
  
  Real voltageLo = 0.0;
  Real voltageHi = 0.0;
  int voltageSteps = 0;
  std::string str;

  // Read in parameters
  pp.get("verbosity", m_verbosity);
  pp.get("voltage_lo", voltageLo);
  pp.get("voltage_hi", voltageHi);
  pp.get("voltage_steps", voltageSteps);
  pp.get("integration", str);
  if (str == "euler"){
    m_integrationAlgorithm = IntegrationAlgorithm::Euler;
  }
  else{
    MayDay::Error("StreamerInceptionStepper::StreamerInceptionStepper -- logic bust");
  }

  // Default alpha
  m_alpha = [](const Real E) { return tanh(E); };
  
}

template <typename P, typename F>
StreamerInceptionStepper<P, F>::~StreamerInceptionStepper()
{
  CH_TIME("StreamerInceptionStepper::~StreamerInceptionStepper");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::~StreamerInceptionStepper" << endl;
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setupSolvers()
{
  CH_TIME("StreamerInceptionStepper::setupSolvers");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setupSolvers" << endl;
  }

  // Define a voltage function to be used in the simulation.
  auto voltage =  [](const Real a_time) -> Real { return 1.0; };
  
  // Instantiate fieldSolver
  m_fieldSolver = RefCountedPtr<FieldSolver>(new F());
  m_fieldSolver->setVerbosity(m_verbosity);
  m_fieldSolver->parseOptions();
  m_fieldSolver->setAmr(m_amr);
  m_fieldSolver->setComputationalGeometry(m_computationalGeometry);
  m_fieldSolver->setVoltage(voltage);
  m_fieldSolver->setRealm(m_realm);
  m_fieldSolver->setTime(0, 0.0, 0.0);

  // Instantiate particleSolver
  m_particleSolver = RefCountedPtr<TracerParticleSolver<P>>(new TracerParticleSolver<P>());
  m_particleSolver->parseOptions();
  m_particleSolver->setAmr(m_amr);
  m_particleSolver->setComputationalGeometry(m_computationalGeometry);
  m_particleSolver->setRealm(m_realm);
  m_particleSolver->setPhase(m_phase);
  m_particleSolver->setTime(0, 0.0, 0.0);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::allocate()
{
  CH_TIME("StreamerInceptionStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::allocate" << endl;
  }

  m_fieldSolver->allocateInternals();
  m_particleSolver->allocate();

  m_amr->allocate(m_velocity, m_realm, m_phase, SpaceDim);
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::initialData()
{
  CH_TIME("StreamerInceptionStepper::initialData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::initialData" << endl;
  }
  
  // Set rho = sigma = 0 for field solver and then solve the Poisson equation. 
  m_fieldSolver->setRho(0.0);
  m_fieldSolver->setSigma(0.0);
  const bool converged = m_fieldSolver->solve(m_fieldSolver->getPotential(), m_fieldSolver->getRho(), m_fieldSolver->getSigma());
  if(!converged) {
    MayDay::Error("StreamerInceptionStepper::initialize -- could not solve Poisson equation. Aborting application.");
  }

  // Get the gas-side electric field and give it to the particle solvers. 
  DataOps::copy(m_velocity, m_amr->alias(phase::gas, m_fieldSolver->getElectricField()));
  
  m_particleSolver->setVelocity(m_velocity);
  
  this->initialParticles();
  this->postInitialize();
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::initialParticles()
{
  CH_TIME("StreamerInceptionStepper::initialParticles");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::initialParticles" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];
    const LevelData<BaseFab<bool>>& validCells = *m_amr->getValidCells(m_realm)[lvl];

    const Real dx = m_amr->getDx()[lvl];
    
    ParticleData<P>& levelParticles = amrParticles[lvl];
      
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box box = dbl[dit()];

      List<P>& particles = levelParticles[dit()].listItems();

      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();

	const RealVect pos = m_amr->getProbLo() + (0.5*RealVect::Unit + RealVect(iv))*dx;

	if(validCells[dit()](iv)){
	  particles.add(P(pos));
	}
      }

      // Define particles' masses (i.e. integration tracker)
      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
	P& p = lit();

	p.mass() = 0.0;
	p.template vect<0>() = p.position();     
	p.template scalar<2>() = 1; // contIntegral criterion. =1 if integration shall continue, =-1 if not (i.e. alpha<0, particle outside boundaries, or particle inside electrode 
      }      
    }
  }
    
  // Remapping
  m_particleSolver->remap(); 

  // // Remove particles inside the EB
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);
}
  

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::postInitialize()
{
  CH_TIME("StreamerInceptionStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::postInitialize" << endl;
  }

  m_particleSolver->interpolateVelocities();
  
  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  while (keepIntegrating(amrParticles)){
    // Check contIntegral criteria and move particles if true
    moveParticles(amrParticles);
  
    // Add integration parts to p.mass()
    addIntegration(amrParticles);
  }
    
  // Remove particles inside the EB
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);  
  
  // Move particles back to initial position
  moveParticlesInitialPosition(amrParticles);
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::moveParticles(ParticleContainer<P>& amrParticles){

  // Lambda for converting RealVect position to IntVect (lower-left corner). Input vector is displacement from probLo. 
  auto        getCellIndex = [](const RealVect a_rv, const Real dx) -> IntVect {
    return IntVect(D_DECL(std::floor(a_rv[0] / dx), std::floor(a_rv[1] / dx), std::floor(a_rv[2] / dx)));
  };
  
  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Set distance moved per step f_dx
  const Real f_dx = 0.01;  //m_amr->getDx()[m_amr->getFinestLevel()];
  
  // Move particles
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    
    const Real dx = m_amr->getDx()[lvl]; // used for getCellPosition, not particleDistanceMoved
    const EBISLayout& ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];
    ParticleData<P>& levelParticles = amrParticles[lvl];
    
    for (DataIterator dit(dbl); dit.ok(); ++dit){
            
      List<P>& particles = levelParticles[dit()].listItems();
      
      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
	P& p = lit();
	
	// Move particle and check if contIntegral if contIntegral not already false
	if (p.template scalar<2>() > 0){
	  
	  const RealVect& E = p.velocity();
	  const Real particleDistanceMoved = f_dx;
	  const RealVect normalizedVelocity = E/E.vectorLength();
	  const Real Emagnitude = E.vectorLength();	  
	  const Real startPointAlpha = m_alpha(Emagnitude);

	  // Check if alpha<0
	  if (startPointAlpha < 0) {
	    p.template scalar<2>() = -1;
	    std::cout << "alpha<0" << endl;
	    continue;
	  }
	  
	  p.template scalar<0>() = startPointAlpha;
	  p.template scalar<1>() = particleDistanceMoved;

	  const RealVect newPosition = p.position() + normalizedVelocity * particleDistanceMoved;

	  const IntVect particleCellIndex = getCellIndex(newPosition - probLo, dx);	

	  // Check if particle is outside of grid boundaries
	  if (particleOutsideGrid(newPosition, probLo, probHi)){
	    p.template scalar<2>() = -1;
	    continue;
	  }	    

	  // Check if particle is inside electrode
	  if (particleInsideElectrode(newPosition)){
	    p.template scalar<2>() = -1;
	    continue;
	  }
		
	  // Move particle if contIntegral
	  p.position() = newPosition;
	}		
      }
    }
  }
  m_particleSolver->remap();
  m_particleSolver->interpolateVelocities();
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::addIntegration(ParticleContainer<P>& amrParticles){
 for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];
    
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit){
	P& p = lit();

	// Add integration parts if contIntegral true
	if (p.template scalar<2>() > 0){	 	
	  const RealVect& E = lit().velocity();
	  const Real Emagnitude = E.vectorLength();
	  const Real endPointAlpha = m_alpha(Emagnitude);
	  const Real startPointAlpha = p.template scalar<0>();
	  const Real distanceMoved = p.template scalar<1>();
	  p.mass() += 0.5*(endPointAlpha + startPointAlpha) * distanceMoved;	
	}	 	
      }
    }
  }
}

template <typename P, typename F>
bool
StreamerInceptionStepper<P, F>::keepIntegrating(ParticleContainer<P>& amrParticles){
  CH_TIME("StreamerInceptionStepper::keepIntegrating(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::keepIntegrating(ParticleContainer)" << endl;
  }
  
  int myParticlesStopped = 1;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];
    
    ParticleData<P>& levelParticles = amrParticles[lvl];
    
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      List<P>& particles = levelParticles[dit()].listItems();
      
      for (ListIterator<P> lit(particles); lit.ok(); ++lit){
	if (lit().template scalar<2>() > 0.0){
	  myParticlesStopped = -1;
	}
      }
    }
  }

  return (ParallelOps::min(myParticlesStopped) > 0.0) ? false : true;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::moveParticlesInitialPosition(ParticleContainer<P>& amrParticles){
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];
    
    ParticleData<P>& levelParticles = amrParticles[lvl];     
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      
      List<P>& particles = levelParticles[dit()].listItems();
      
      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
	P& p = lit();
	p.position() = p.template vect<0>();
      }
    }
  }
  m_particleSolver->remap();
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::registerRealms()
{
  CH_TIME("StreamerInceptionStepper::registerRealms");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_realm);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::registerOperators()
{
  CH_TIME("StreamerInceptionStepper::registerOperators");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::registerOperators" << endl;
  }

  m_fieldSolver->registerOperators();
  m_particleSolver->registerOperators();
}

// Done in the constructor
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::parseOptions()
{
  CH_TIME("StreamerInceptionStepper::parseOptions");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseOptions" << endl;
  }
}

// Not used in this program
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::parseRuntimeOptions()
{
  CH_TIME("StreamerInceptionStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseRuntimeOptions" << endl;
  }
}

#ifdef CH_USE_HDF5
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("StreamerInceptionStepper::writeCheckpointData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writeCheckpointData" << endl;
  }
}
#endif

#ifdef CH_USE_HDF5
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("StreamerInceptionStepper::readCheckpointData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::readCheckpointData" << endl;
  }
}
#endif

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::writePlotData(EBAMRCellData&       a_output,
                                              Vector<std::string>& a_plotVariableNames,
                                              int&                 a_icomp) const
{
  CH_TIME("StreamerInceptionStepper::writePlotData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writePlotData" << endl;
  }

  // Append plot variable names.
  a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());
  a_plotVariableNames.append(m_particleSolver->getPlotVariableNames());

  // Write into plot data holder memory.
  m_fieldSolver->writePlotData(a_output, a_icomp);
  m_particleSolver->writePlotData(a_output, a_icomp);
}

template <typename P, typename F>
int
StreamerInceptionStepper<P, F>::getNumberOfPlotVariables() const
{
  CH_TIME("StreamerInceptionStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::getNumberOfPlotVariables" << endl;
  }

  int ncomp = 0;

  ncomp += m_fieldSolver->getNumberOfPlotVariables();
  ncomp += m_particleSolver->getNumberOfPlotVariables();
 
  return ncomp;
}

// TO DO: Do something with this one?
template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::computeDt()
{
  CH_TIME("StreamerInceptionStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeDt" << endl;
  }

  return std::numeric_limits<Real>::max();
}

// StreamerInceptionStepper is static. Any call to advance gives an error.
template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::advance(const Real a_dt)
{
  CH_TIME("StreamerInceptionStepper::advance");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::advance" << endl;
  }

  MayDay::Error("StreamerInceptionStepper::advance -- this is a static solver. Calling the advance method is an error.");

  return a_dt;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("StreamerInceptionStepper::synchronizeSolverTimes");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_fieldSolver->setTime(a_step, a_time, a_dt);
  m_particleSolver->setTime(a_step, a_time, a_dt);
  
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::printStepReport()
{
  CH_TIME("StreamerInceptionStepper::printStepReport");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::printStepReport" << endl;
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("StreamerInceptionStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::preRegrid" << endl;
  }

  m_fieldSolver->preRegrid(a_lmin, a_oldFinestLevel);
  m_particleSolver->preRegrid(a_lmin, a_oldFinestLevel);
}

// Combined implementation from FieldStepper and TracerParticleStepper
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("StreamerInceptionStepper::regrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::regrid" << endl;
  }

  // Regrid tracer particles and field
  m_fieldSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_particleSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  // Regrid velocity field
  m_amr->reallocate(m_velocity, m_phase, a_lmin);  
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::postRegrid()
{
  CH_TIME("StreamerInceptionStepper::postRegrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::postRegrid" << endl;
  }

  // This can be left empty because initialData is called right after the regrid method and calling
  // the advance method is an error. 
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setAlpha(const std::function<Real(const Real& E)>& a_alpha) noexcept
{
  CH_TIME("StreamerInceptionStepper::setAlpha");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setAlpha" << endl;
  }

  m_alpha = a_alpha;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setVoltage(const Real a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::setVoltage");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setVoltage" << endl;
  }

  CH_assert(!m_fieldSolver.isNull());

  // Give field solver a constant voltage.
  m_fieldSolver->setVoltage([a_voltage](const Real& a_time) -> Real { return a_voltage; });
}


template <typename P, typename F>
bool
StreamerInceptionStepper<P, F>::particleOutsideGrid(const RealVect pos, const RealVect probLo, const RealVect probHi){
  for (int dir = 0; dir < SpaceDim; dir++) {
    if(pos[dir] <= probLo[dir] || pos[dir] >= probHi[dir]) {
      return true;
    }
  }
  return false;  
}


template <typename P, typename F>
bool
StreamerInceptionStepper<P, F>::particleInsideElectrode(const RealVect pos){
   
  const RefCountedPtr<BaseIF>& implicitFunction = m_amr->getBaseImplicitFunction(m_phase);
  const Real value = implicitFunction->value(pos);
  if (value > 0.0){
    return true;
  }
  return false;
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::scaleElectrostatics(const Real& a_scale) noexcept
{
  CH_TIME("StreamerInceptionStepper::scaleElectrostatics");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::scaleElectrostatics" << endl;
  }

  CH_assert(!m_fieldSolver.isNull());

  MFAMRCellData& potential     = m_fieldSolver->getPotential();
  MFAMRCellData& electricField = m_fieldSolver->getElectricField();

  DataOps::scale(potential, a_scale);
  DataOps::scale(electricField, a_scale);
}


#include <CD_NamespaceFooter.H>

#endif
