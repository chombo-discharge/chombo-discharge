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
  // for (int lvl = ....
  // 	 for (DataIterator dit(...)
  // 		const Box bx = dbl[dit()];
  // 		EBCellFAB& velocity = (*m_velocity[lvl])[dit()];
  // 	      BaseFab<Real>& vFAB = velocity.getSingleValuedFAB();
	      
  // 	      for (BoxIterator bit(bx); bit.ok(); ++bit) {
  // 		const IntVect iv = bit();
  // 		  RealVect v = RealVect(D_DECL(vFAB(iv,0), vFab(iv,1), vFab(iv,2)));
  // 		  v = v/v.vectorLength();

  // 		  for (int dir = 0; dir < SpaceDim; dir++) {
  // 		    vFab(iv,dir) = v[dir];
  // 		  }
  // 	      }
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

      // Define particles' masses
      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
	P& p = lit();

	p.mass() = 0.0;
	p.template vect<0>() = p.position();     
	p.template scalar<2>() = 1; // criterion for continued integration (set to -1 if alpha<0, particle outside boundary, or particle inside electrode)
      }      
    }
  }
    
  // Remapping
  m_particleSolver->remap(); 

  // // Remove particles inside the EB
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);
}
  


// TO DO
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::postInitialize()
{
  CH_TIME("StreamerInceptionStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::postInitialize" << endl;
  }

  // Lambda for converting RealVect position to IntVect (lower-left corner). Input vector is displacement from probLo. 
  auto        getCellIndex = [](const RealVect a_rv, const Real dx) -> IntVect {
    return IntVect(D_DECL(std::floor(a_rv[0] / dx), std::floor(a_rv[1] / dx), std::floor(a_rv[2] / dx)));
  };

  m_particleSolver->interpolateVelocities();
  
  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();
  
  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Use finest resolution as dx for all particles
  const Real f_dx = -0.1;  //m_amr->getDx()[m_amr->getFinestLevel()];
  
  // Move particles
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    
    const Real dx = m_amr->getDx()[lvl]; // for getCellPosition, not particleDistanceMoved
    const EBISLayout& ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];
    ParticleData<P>& levelParticles = amrParticles[lvl];
    
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      
      
      List<P>& particles = levelParticles[dit()].listItems();
      
      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
	const RealVect& E = lit().velocity();
	const Real particleDistanceMoved = f_dx;
	const RealVect normalizedVelocity = E/E.vectorLength();

	P& p = lit();

	const Real Emagnitude = E.vectorLength();
	const Real startPointAlpha = m_alpha(Emagnitude);
	p.template scalar<0>() = startPointAlpha;
	p.template scalar<1>() = particleDistanceMoved;

	const RealVect newPosition = p.position() + normalizedVelocity * particleDistanceMoved;

	const IntVect particleCellIndex = getCellIndex(newPosition - probLo, dx);	

	// Check if particle is outside of grid boundaries
	for (int dir = 0; dir < SpaceDim; dir++) {
	  if(newPosition[dir] <= probLo[dir] || newPosition[dir] >= probHi[dir]) {
	    p.template scalar<2>() = -1;	    
	  }
	}

	// Move particle if inside grid boundaries
	if (p.template scalar<2>() > 0){
	  p.position() = newPosition;	     
	}		
      }
    }
  }
  m_particleSolver->remap();
  m_particleSolver->interpolateVelocities();

  
  // Iterate through particles and check if particles are covered or irregular (must be done after previous remap)
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const Real dx = m_amr->getDx()[lvl];
    const EBISLayout& ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];
    ParticleData<P>& levelParticles = amrParticles[lvl];
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const EBISBox& ebisbox = ebisl[dit()];
      List<P>& particles = levelParticles[dit()].listItems();	
      for (ListIterator<P> lit(particles); lit.ok(); ++lit){
	P& p = lit();

	const IntVect particleCellIndex = getCellIndex(p.position() - probLo, dx);	
	
	// Check if particle is covered
	if (ebisbox.isCovered(particleCellIndex)){
	  p.template scalar<2>() = -1;
	  //std::cout << "Particle integration stopped: Particle is in covered cell" << endl;
	}
	
	// Check if particle is irregular
	else if (ebisbox.isIrregular(particleCellIndex)){
	  VolIndex vof(particleCellIndex, 0);     
	  
	  const RealVect n = ebisbox.normal(vof);
	  const RealVect xc = ebisbox.bndryCentroid(vof) * dx * probLo;

	  if (n.dotProduct(p.position() - xc) <= 0.0) {
	    p.template scalar<2>() = -1;	    
	    std::cout << "Particle Integration stopped: Particle is in covered part of multi-cell" << endl;
	  }	  
	}      
      }
    }
  }

  // Add integration parts to p.mass()
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];
    
    for (DataIterator dit(dbl); dit.ok(); ++dit){
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit){
	P& p = lit();
	const RealVect& E = lit().velocity();
	const Real Emagnitude = E.vectorLength();
	const Real endPointAlpha = m_alpha(Emagnitude);
	const Real startPointAlpha = p.template scalar<0>();
	const Real distanceMoved = p.template scalar<1>();
		 
	p.mass() += 0.5*(endPointAlpha + startPointAlpha) * distanceMoved;
      }
    }
  }
  
  // Remove particles inside the EB
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);  
  
  // Move particles back to initial position
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
